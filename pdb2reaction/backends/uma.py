# pdb2reaction/backends/uma.py

"""
UMA (fairchem) backend for pdb2reaction.

Provides ``UMACalculator`` (PySisyphus-compatible) and ``UMAASECalculator``
(ASE-compatible, for DMF path optimization).
"""

from __future__ import annotations

import time
from typing import Any, Dict, Optional, Sequence

import click
from pdb2reaction.core.output import emit
import numpy as np
import torch
import torch.nn as nn
from ase import Atoms

from fairchem.core import pretrained_mlip, FAIRChemCalculator
from fairchem.core.datasets.atomic_data import AtomicData
from fairchem.core.datasets import data_list_collater

# Optional: only needed when workers > 1
try:
    from fairchem.core.units.mlip_unit.predict import ParallelMLIPPredictUnit
    from fairchem.core.units.mlip_unit.api.inference import guess_inference_settings
except Exception:
    ParallelMLIPPredictUnit = None
    guess_inference_settings = None

# fp64 base precision: switching OMol-trained UMA from default fp32 to
# fp64 can have non-trivial impact on TSopt + Hessian numerics. Available
# via InferenceSettings(base_precision_dtype="float64") in fairchem ≥ 2.0.
try:
    from fairchem.core.units.mlip_unit.api.inference import InferenceSettings as _UMAInferenceSettings
except Exception:
    _UMAInferenceSettings = None

from pysisyphus.constants import BOHR2ANG, ANG2BOHR, AU2EV

from .base import MLIPCalculator, BackendError, normalize_hessian_calc_mode
from pdb2reaction.core.defaults import DEFAULT_UMA_MODEL


# Strict deterministic-mode setup for UMA is opt-in (default no-op).
# It lives in `pdb2reaction.backends._determinism` and is enabled via the
# `--deterministic` CLI flag or `PDB2REACTION_STRICT_DETERMINISTIC=1`; the
# backend factory and the CLI option callback drive it.
# Precision controls numerical accuracy but does not guarantee determinism.
# Strict same-stack repeatability is opt-in and must be verified with the
# end-to-end two-run smoke gate on the installed software/hardware stack.


# ---- Unit conversion constants ----
EV2AU = 1.0 / AU2EV
F_EVAA_2_AU = EV2AU / ANG2BOHR
H_EVAA_2_AU = EV2AU / ANG2BOHR / ANG2BOHR


def _positive_worker_count(value: Any, name: str) -> int:
    """Return a positive integer worker count."""
    if value is None:
        return 1
    if isinstance(value, bool):
        raise BackendError(f"{name} must be a positive integer, got {value!r}.")
    try:
        count = int(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise BackendError(f"{name} must be a positive integer, got {value!r}.") from exc
    try:
        if float(value) != count:
            raise BackendError(f"{name} must be a positive integer, got {value!r}.")
    except (TypeError, ValueError, OverflowError) as exc:
        raise BackendError(f"{name} must be a positive integer, got {value!r}.") from exc
    if count < 1:
        raise BackendError(f"{name} must be a positive integer, got {value!r}.")
    return count


class UMAcore:
    """Thin wrapper around fairchem-UMA predict_unit.

    If ``workers > 1``, uses ``ParallelMLIPPredictUnit`` which does not
    expose ``predictor.model``.  Analytical Hessians are unavailable in
    that mode.
    """

    def __init__(
        self,
        elem: Sequence[str],
        *,
        charge: int = 0,
        spin: int = 1,
        model: str = DEFAULT_UMA_MODEL,
        task_name: str = "omol",
        device: str = "auto",
        workers: int = 1,
        workers_per_node: int = 1,
        max_neigh: Optional[int] = None,
        radius: Optional[float] = None,
        r_edges: bool = False,
        precision: str = "fp32",
    ):
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device_str = device
        self.device = torch.device(device)

        self.workers = _positive_worker_count(workers, "workers")
        self.workers_per_node = _positive_worker_count(
            workers_per_node, "workers_per_node"
        )
        self.parallel_predict = self.workers > 1

        # fp32 is the established baseline; fp64 enables full-precision
        # base inference and can change TSopt + Hessian numerics non-trivially.
        # Selected per-call via `precision="fp32"|"fp64"`.
        self.precision = str(precision or "fp32").lower()
        if self.precision not in ("fp32", "fp64"):
            raise BackendError(f"UMA precision must be 'fp32' or 'fp64', got {precision!r}")
        if self.precision == "fp64" and _UMAInferenceSettings is None:
            raise BackendError(
                "UMA precision='fp64' requires fairchem-core's InferenceSettings; "
                "upgrade fairchem-core (≥ 2.0) or pass precision='fp32'."
            )
        uma_inference_settings = (
            _UMAInferenceSettings(base_precision_dtype="float64")
            if self.precision == "fp64" and _UMAInferenceSettings is not None
            else None
        )

        self._AtomicData = AtomicData
        self._collater = data_list_collater

        if self.parallel_predict:
            if (ParallelMLIPPredictUnit is None) or (guess_inference_settings is None):
                raise ImportError(
                    "workers>1 requested, but ParallelMLIPPredictUnit/guess_inference_settings "
                    "could not be imported from fairchem. Please ensure your FAIR-Chem installation "
                    "includes `fairchem-core[extras]`."
                )
            ckpt_path = pretrained_mlip.pretrained_checkpoint_path_from_name(model)
            # fp64 path: pass our InferenceSettings; fp32 path: fairchem default
            inference_settings = uma_inference_settings or guess_inference_settings("default")
            atom_refs = pretrained_mlip.get_reference_energies(model, reference_type="atom_refs")
            form_elem_refs = pretrained_mlip.get_reference_energies(model, reference_type="form_elem_refs")
            self.predict = ParallelMLIPPredictUnit(
                inference_model_path=str(ckpt_path),
                device=self.device_str,
                inference_settings=inference_settings,
                atom_refs=atom_refs,
                form_elem_refs=form_elem_refs,
                num_workers=self.workers,
                num_workers_per_node=self.workers_per_node,
            )
        else:
            predict_kwargs = dict(device=self.device_str, workers=self.workers)
            if uma_inference_settings is not None:
                predict_kwargs["inference_settings"] = uma_inference_settings
            self.predict = pretrained_mlip.get_predict_unit(model, **predict_kwargs)

        self.has_torch_model = hasattr(self.predict, "model") and isinstance(
            getattr(self.predict, "model", None), nn.Module
        )
        if self.has_torch_model:
            self.predict.model.eval()
            for m in self.predict.model.modules():
                if isinstance(m, nn.Dropout):
                    m.p = 0.0

        self.elem = [e.capitalize() for e in elem]
        self.charge = charge
        self.spin = spin
        self.task_name = task_name
        self._max_neigh = max_neigh
        self._radius = radius
        self._r_edges = r_edges

    def _model_backbone(self):
        if not self.has_torch_model:
            return None
        mdl = self.predict.model
        mod = getattr(mdl, "module", mdl)
        return getattr(mod, "backbone", None)

    def _ase_to_batch(self, atoms: Atoms):
        backbone = self._model_backbone()
        default_max_neigh = getattr(backbone, "max_neighbors", None) if backbone is not None else None
        default_radius = getattr(backbone, "cutoff", None) if backbone is not None else None
        if default_radius is None:
            default_radius = 6.0
        max_neigh = self._max_neigh if self._max_neigh is not None else default_max_neigh
        radius = self._radius if self._radius is not None else default_radius
        r_edges = self._r_edges

        atoms.info.update({"charge": self.charge, "spin": self.spin})
        # When the inference path uses fp64, hand AtomicData the matching
        # target dtype so it does not down-cast positions to fp32 only to be
        # re-upcasted (and emit the fairchem `Upcasting atomic coordinates`
        # WARNING on every call). fp32 path keeps fairchem's default.
        target_dtype = torch.float64 if self.precision == "fp64" else torch.float32
        # `r_data_keys=["spin", "charge"]` is REQUIRED: fairchem's AtomicData.from_ase
        # only reads charge/spin from atoms.info when the key is listed here — otherwise
        # it silently defaults BOTH to 0 (see fairchem AtomicData.from_ase). Older
        # fairchem read atoms.info unconditionally; a fairchem version bump gated it on
        # r_data_keys, so omitting this makes the whole pysisyphus path run at
        # charge=0/spin=0 regardless of the user-requested charge. The ASE/DMF path
        # (FAIRChemCalculator) already passes r_data_keys, which is why only this
        # (pysisyphus) path was affected.
        data = self._AtomicData.from_ase(
            atoms,
            max_neigh=max_neigh,
            radius=radius,
            r_edges=r_edges,
            target_dtype=target_dtype,
            r_data_keys=["spin", "charge"],
        )
        data.dataset = self.task_name
        batch = self._collater([data], otf_graph=True)
        if not self.parallel_predict:
            batch = batch.to(self.device)
        return batch

    def compute(
        self,
        coord_ang: np.ndarray,
        *,
        forces: bool = False,
        hessian: bool = False,
    ) -> Dict[str, Any]:
        atoms = Atoms(self.elem, positions=coord_ang)
        batch = self._ase_to_batch(atoms)

        if self.parallel_predict or (not self.has_torch_model):
            if hessian:
                raise RuntimeError(
                    "Analytical Hessian is not available when predictor workers > 1 "
                    "or when predictor.model is not exposed. Use FiniteDifference Hessian."
                )
            res = self.predict.predict(batch)
            energy = float(res["energy"].squeeze().detach().item())
            forces_np = res["forces"].detach().cpu().numpy() if forces else None
            return {"energy": energy, "forces": forces_np, "hessian": None}

        batch.pos.requires_grad_(True)
        res = self.predict.predict(batch)
        energy = float(res["energy"].squeeze().detach().item())
        forces_np = res["forces"].detach().cpu().numpy() if (forces or hessian) else None

        if hessian:
            p_flags = [p.requires_grad for p in self.predict.model.parameters()]
            for p in self.predict.model.parameters():
                p.requires_grad_(False)
            self.predict.model.train()
            try:
                def e_fn(flat):
                    batch.pos = flat.view(-1, 3)
                    return self.predict.predict(batch)["energy"].squeeze()

                H = torch.autograd.functional.hessian(e_fn, batch.pos.view(-1), vectorize=False)
                H = H.view(len(atoms), 3, len(atoms), 3).detach()
            finally:
                self.predict.model.eval()
                for p, flag in zip(self.predict.model.parameters(), p_flags):
                    p.requires_grad_(flag)
                if self.device.type == "cuda":
                    torch.cuda.empty_cache()
        else:
            H = None

        return {"energy": energy, "forces": forces_np, "hessian": H}

    def forces_tensor(self, coord_ang: np.ndarray) -> torch.Tensor:
        """Return device-native forces (eV/Ang) for the FD Hessian assembler.

        Returns the force tensor on the model's execution device, skipping the
        ``.cpu().numpy()`` round-trip and the scalar energy ``.item()`` sync that
        ``compute()`` performs for each displacement. Only valid on the native
        torch-model path (``has_torch_model and not parallel_predict``); other
        providers must use ``compute()``. The returned values are the same tensor
        ``compute()`` would convert to NumPy, so the assembled central-difference
        Hessian is bit-identical to the NumPy round-trip path.
        """
        atoms = Atoms(self.elem, positions=coord_ang)
        batch = self._ase_to_batch(atoms)
        res = self.predict.predict(batch)
        return res["forces"].detach()


class UMACalculator(MLIPCalculator):
    """UMA (fairchem) backend.

    Extends ``MLIPCalculator`` with UMA-specific features:
    multi-worker inference and VRAM logging. Analytical Hessian
    (autograd) is also available on the other backends.
    """

    def __init__(
        self,
        *,
        model: str = DEFAULT_UMA_MODEL,
        task_name: str = "omol",
        workers: int = 1,
        workers_per_node: int = 1,
        max_neigh: Optional[int] = None,
        radius: Optional[float] = None,
        r_edges: bool = False,
        precision: str = "fp32",
        out_hess_torch: bool = True,
        print_vram: bool = True,
        # Base class parameters
        charge: int = 0,
        spin: int = 1,
        device: str = "auto",
        freeze_atoms: Optional[Sequence[int]] = None,
        hessian_calc_mode: str = "FiniteDifference",
        return_partial_hessian: bool = True,
        hessian_double: bool = True,
        print_timing: bool = True,
        **kwargs,
    ):
        mode = normalize_hessian_calc_mode(hessian_calc_mode)
        worker_count = _positive_worker_count(workers, "workers")
        worker_count_per_node = _positive_worker_count(
            workers_per_node, "workers_per_node"
        )
        if worker_count > 1 and mode == "Analytical":
            raise BackendError(
                "Analytical Hessian cannot be combined with UMA workers>1: "
                "the parallel predictor exposes no autograd model. Use workers=1 "
                "or select hessian_calc_mode='FiniteDifference'."
            )
        super().__init__(
            charge=charge,
            spin=spin,
            device=device,
            freeze_atoms=freeze_atoms,
            hessian_calc_mode=mode,
            return_partial_hessian=return_partial_hessian,
            hessian_double=hessian_double,
            print_timing=print_timing,
            **kwargs,
        )
        self._core: Optional[UMAcore] = None
        self._core_kw = dict(
            charge=charge,
            spin=spin,
            model=model,
            task_name=task_name,
            device=device,
            workers=worker_count,
            workers_per_node=worker_count_per_node,
            max_neigh=max_neigh,
            radius=radius,
            r_edges=r_edges,
            precision=precision,
        )
        self.out_hess_torch = out_hess_torch
        self.print_vram = bool(print_vram)

    def _ensure_core(self, elem: Sequence[str]):
        normalized = [symbol.capitalize() for symbol in elem]
        if self._core is None or self._core.elem != normalized:
            self._core = UMAcore(elem, **self._core_kw)

    def _supports_analytical_hessian(self) -> bool:
        return True  # Analytical path available only when workers == 1 and predictor.model is accessible

    # Subclass hooks are NOT used for UMA; we override get_* directly
    # to retain the GPU FD Hessian and torch-based analytical Hessian.

    def _compute_energy_forces_ev(self, elem, coord_ang):
        self._ensure_core(elem)
        res = self._core.compute(coord_ang, forces=True, hessian=False)
        return res["energy"], res["forces"]

    # Hessian helpers (UMA-specific)

    def _emit_hessian_logs(
        self,
        *,
        core: UMAcore,
        mode_label: str,
        mode_elapsed_s: float,
        total_elapsed_s: float,
        vram_base_alloc: Optional[float] = None,
        vram_base_reserved: Optional[float] = None,
        vram_total: Optional[float] = None,
    ) -> None:
        hessian_vram_summary = ""
        hessian_vram_detail = None
        if self.print_vram and core.device.type == "cuda":
            dev = core.device
            torch.cuda.synchronize(device=dev)
            base_alloc = float(vram_base_alloc or 0.0)
            base_reserved = float(vram_base_reserved or 0.0)
            peak_alloc = max(float(torch.cuda.max_memory_allocated(device=dev)) - base_alloc, 0.0) / 1e9
            peak_reserved_abs = float(torch.cuda.max_memory_reserved(device=dev))
            peak_reserved = max(peak_reserved_abs - base_reserved, 0.0) / 1e9
            total_vram = float(vram_total or torch.cuda.get_device_properties(dev).total_memory) / 1e9
            remaining_vram = max((total_vram * 1e9) - peak_reserved_abs, 0.0) / 1e9
            hessian_vram_summary = f", peak VRAM {peak_alloc:.2f} GB"
            hessian_vram_detail = (
                f"[HessianVRAM] total={total_vram:.3f} GB | "
                f"peak_allocated={peak_alloc:.3f} GB | "
                f"peak_reserved={peak_reserved:.3f} GB | "
                f"remaining={remaining_vram:.3f} GB"
            )
        from pdb2reaction.core.utils import verbose_level
        if self.print_timing:
            emit(
                f"[hessian] Completed {mode_label} Hessian: {total_elapsed_s:.2f} s"
                f"{hessian_vram_summary}",
                detail=True,
            )
            if verbose_level() >= 3:
                click.echo(f"[HessianTiming] mode: {mode_label} | elapsed: {mode_elapsed_s:.2f} s")
        if hessian_vram_detail is not None and verbose_level() >= 3:
            click.echo(hessian_vram_detail)

    def _au_hessian(self, H: torch.Tensor):
        """Convert Hessian from eV/Å² to Hartree/Bohr² (torch version)."""
        n = H.size(0)
        H = H.view(n * 3, n * 3)
        if self.hessian_double:
            H = H.to(dtype=torch.float64)
        # Bounded-peak symmetrization (helper writes both triangles).
        from pdb2reaction.core.utils import symmetrize_inplace
        symmetrize_inplace(H)
        H.mul_(H_EVAA_2_AU)
        if self.out_hess_torch:
            return H.detach()
        else:
            # C-NEED: numpy return at boundary; GPU H freed implicitly when caller drops the return value.
            return H.detach().cpu().numpy()

    def _apply_analytical_active_trim(self, H: torch.Tensor) -> torch.Tensor:
        """Extract the active-DOF sub-block (n_act×3×n_act×3) when
        ``return_partial_hessian=True``, or zero frozen rows/columns and
        return (n_atoms×3×n_atoms×3) when ``False``."""
        n_atoms = H.size(0)
        if len(self.freeze_atoms) == 0:
            return H
        _, active_dof_idx, frozen_dof_idx = self._active_and_frozen_dof_idx(n_atoms)
        H = H.view(n_atoms * 3, n_atoms * 3)
        if self.return_partial_hessian:
            idx = torch.tensor(active_dof_idx, device=H.device, dtype=torch.long)
            H_sub = H.index_select(0, idx).index_select(1, idx)
            n_act = len(active_dof_idx) // 3
            return H_sub.view(n_act, 3, n_act, 3)
        else:
            if frozen_dof_idx:
                cols = torch.tensor(frozen_dof_idx, device=H.device, dtype=torch.long)
                H.index_fill_(1, cols, 0.0)
                H.index_fill_(0, cols, 0.0)
            return H.view(n_atoms, 3, n_atoms, 3)

    def _build_fd_hessian_gpu(
        self,
        elem: Sequence[str],
        coord_ang: np.ndarray,
        *,
        eps_ang: float = 1.0e-3,
    ) -> Dict[str, Any]:
        """Assemble FD Hessian on GPU (UMA-specific)."""
        self._ensure_core(elem)
        core = self._core
        assert core is not None
        dev = core.device

        n_atoms = len(elem)
        dof = n_atoms * 3
        frozen_set = set(self.freeze_atoms)
        active_atoms = [i for i in range(n_atoms) if i not in frozen_set]
        active_dof_idx = [3 * i + j for i in active_atoms for j in range(3)]

        res0 = core.compute(coord_ang, forces=True, hessian=False)
        energy0_eV = res0["energy"]
        F0 = res0["forces"]

        force_dtype = torch.from_numpy(F0).dtype
        hessian_dtype = torch.float64 if self.hessian_double else force_dtype
        n_active = len(active_dof_idx)
        if self.return_partial_hessian and n_active < dof:
            # Compact-from-the-start: rows = active DOFs, cols = full dof.
            # Saves (dof - n_active) * dof * itemsize bytes of GPU VRAM
            # (e.g. ~7 GB at 10k atoms / 200 active). Partial Hessian is
            # preserved as-is (no upper-triangle-only tricks); row-write
            # populates both triangles of the active sub-block symmetrically
            # after the post-loop extract.
            H = torch.zeros((n_active, dof), device=dev, dtype=hessian_dtype)
            _partial_alloc = True
        else:
            H = torch.zeros((dof, dof), device=dev, dtype=hessian_dtype)
            _partial_alloc = False

        coord_plus = coord_ang.copy()
        coord_minus = coord_ang.copy()

        # On the native torch-model path, read forces directly as a
        # device tensor and assemble each column without a GPU->CPU->GPU
        # round-trip (and without the per-displacement scalar energy sync). The
        # NumPy round-trip is retained for non-native providers.
        native_tensor = getattr(core, "has_torch_model", False) and not getattr(
            core, "parallel_predict", False
        )

        for i_local, k in enumerate(active_dof_idx):
            a = k // 3
            c = k % 3
            coord_plus[a, c] = coord_ang[a, c] + eps_ang
            coord_minus[a, c] = coord_ang[a, c] - eps_ang
            if native_tensor:
                Fp_t = core.forces_tensor(coord_plus).reshape(-1).to(dev, dtype=hessian_dtype)
                Fm_t = core.forces_tensor(coord_minus).reshape(-1).to(dev, dtype=hessian_dtype)
            else:
                Fp = core.compute(coord_plus, forces=True, hessian=False)["forces"].reshape(-1)
                Fm = core.compute(coord_minus, forces=True, hessian=False)["forces"].reshape(-1)
                Fp_t = torch.from_numpy(Fp).to(dev, dtype=hessian_dtype)
                Fm_t = torch.from_numpy(Fm).to(dev, dtype=hessian_dtype)
            col = -(Fp_t - Fm_t) / (2.0 * eps_ang)
            if _partial_alloc:
                H[i_local, :] = col
            else:
                H[:, k] = col

            coord_plus[a, c] = coord_ang[a, c]
            coord_minus[a, c] = coord_ang[a, c]

        if self.return_partial_hessian:
            n_active_atoms = len(active_atoms)
            idx = torch.tensor(active_dof_idx, device=dev, dtype=torch.long)
            if _partial_alloc:
                # H is already (n_active, dof) — single index_select(1) extracts to (n_active, n_active).
                H = H.index_select(1, idx)
            elif frozen_set:
                # Bounded row-chunk active-square extraction avoids
                # the full (n_active, dof) row temporary of the chained form.
                from pysisyphus._array import active_square

                H = active_square(H, idx)
            H = H.view(n_active_atoms, 3, n_active_atoms, 3)
        else:
            # Mirror the analytical/base path: zero frozen rows/columns so the
            # non-partial FD Hessian does not keep half-magnitude frozen-atom
            # coupling (unreachable on shipped workflows, which force partial).
            if frozen_set:
                fidx = torch.tensor(
                    [3 * i + j for i in sorted(frozen_set) for j in range(3)],
                    device=dev,
                    dtype=torch.long,
                )
                H.index_fill_(0, fidx, 0.0)
                H.index_fill_(1, fidx, 0.0)
            H = H.view(n_atoms, 3, n_atoms, 3)

        return {"energy": energy0_eV, "forces": np.asarray(F0, dtype=np.float64).reshape(-1, 3), "hessian": H}

    # PySisyphus API (override base to use GPU path + torch Hessian)

    def get_energy(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        res = self._core.compute(coord_ang, forces=False, hessian=False)
        return {"energy": res["energy"] * EV2AU}

    def get_forces(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        res = self._core.compute(coord_ang, forces=True, hessian=False)
        F_ev = self._zero_frozen_forces_ev(np.asarray(res["forces"], dtype=np.float64).reshape(-1, 3))
        return {
            "energy": res["energy"] * EV2AU,
            "forces": (F_ev * F_EVAA_2_AU).reshape(-1),
        }

    def get_hessian(self, elem, coords):
        mode = normalize_hessian_calc_mode(self.hessian_calc_mode)
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG

        core = self._core
        assert core is not None
        force_fd = (core.parallel_predict or (not core.has_torch_model))

        # Parallel predictor (workers>1) has no autograd model, so an explicit
        # analytical request cannot be honoured. Fail instead of silently
        # changing the requested numerical method.
        if force_fd and mode == "Analytical":
            raise BackendError(
                "Analytical Hessian is unavailable because the UMA predictor "
                "does not expose an autograd model. Use workers=1 or select "
                "hessian_calc_mode='FiniteDifference'."
            )

        vram_base_alloc: Optional[float] = None
        vram_base_reserved: Optional[float] = None
        vram_total: Optional[float] = None
        if self.print_vram and core.device.type == "cuda":
            torch.cuda.synchronize(device=core.device)
            vram_base_alloc = float(torch.cuda.memory_allocated(device=core.device))
            vram_base_reserved = float(torch.cuda.memory_reserved(device=core.device))
            vram_total = float(torch.cuda.get_device_properties(core.device).total_memory)
            torch.cuda.reset_peak_memory_stats(device=core.device)

        hess_total_start = time.perf_counter()
        mode_elapsed_s = 0.0
        mode_label = "FiniteDifference"

        if (not force_fd) and mode == "Analytical":
            mode_label = "Analytical"
            t0 = time.perf_counter()
            try:
                res = self._core.compute(coord_ang, forces=True, hessian=True)
            except (torch.cuda.OutOfMemoryError, RuntimeError) as e:
                msg = str(e).lower()
                if "out of memory" in msg and "cuda" in msg:
                    raise RuntimeError(
                        "Analytical Hessian computation failed due to CUDA out-of-memory. "
                        "Please switch to `--hessian-calc-mode FiniteDifference`."
                    ) from e
                raise
            mode_elapsed_s = time.perf_counter() - t0

            res_forces_ev = self._zero_frozen_forces_ev(np.asarray(res["forces"], dtype=np.float64).reshape(-1, 3))
            H = self._apply_analytical_active_trim(res["hessian"])
            out = {
                "energy": res["energy"] * EV2AU,
                "forces": (res_forces_ev * F_EVAA_2_AU).reshape(-1),
                "hessian": self._au_hessian(H),
            }
        else:
            t0 = time.perf_counter()
            res = self._build_fd_hessian_gpu(elem, coord_ang)
            mode_elapsed_s = time.perf_counter() - t0

            res_forces_ev = self._zero_frozen_forces_ev(np.asarray(res["forces"], dtype=np.float64).reshape(-1, 3))
            out = {
                "energy": res["energy"] * EV2AU,
                "forces": (res_forces_ev * F_EVAA_2_AU).reshape(-1),
                "hessian": self._au_hessian(res["hessian"]),
            }

        total_elapsed_s = time.perf_counter() - hess_total_start
        self._emit_hessian_logs(
            core=core,
            mode_label=mode_label,
            mode_elapsed_s=mode_elapsed_s,
            total_elapsed_s=total_elapsed_s,
            vram_base_alloc=vram_base_alloc,
            vram_base_reserved=vram_base_reserved,
            vram_total=vram_total,
        )

        partial_meta = self._build_partial_hessian_meta(len(elem))
        if partial_meta is not None:
            out["within_partial_hessian"] = partial_meta

        return out


class UMAASECalculator(FAIRChemCalculator):
    """UMA ASE calculator (wraps FAIRChemCalculator) for use with ASE-based path optimizers (DMF/NEB/growing-string)."""

    def __init__(
        self,
        *,
        model: str = DEFAULT_UMA_MODEL,
        device: str = "auto",
        task_name: str = "omol",
        workers: int = 1,
        workers_per_node: int = 1,
        precision: str = "fp32",
    ):
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        num_workers = _positive_worker_count(workers, "workers")
        num_workers_per_node = _positive_worker_count(
            workers_per_node, "workers_per_node"
        )

        precision = str(precision or "fp32").lower()
        if precision not in ("fp32", "fp64"):
            raise BackendError(f"UMA precision must be 'fp32' or 'fp64', got {precision!r}")
        if precision == "fp64" and _UMAInferenceSettings is None:
            raise BackendError(
                "UMA precision='fp64' requires fairchem-core's InferenceSettings; "
                "upgrade fairchem-core (>= 2.0) or pass precision='fp32'."
            )
        uma_inference_settings = (
            _UMAInferenceSettings(base_precision_dtype="float64")
            if precision == "fp64" and _UMAInferenceSettings is not None
            else None
        )

        if num_workers > 1:
            if (ParallelMLIPPredictUnit is None) or (guess_inference_settings is None):
                raise ImportError(
                    "workers>1 requested, but ParallelMLIPPredictUnit/guess_inference_settings "
                    "could not be imported from fairchem."
                )
            ckpt_path = pretrained_mlip.pretrained_checkpoint_path_from_name(model)
            inference_settings = uma_inference_settings or guess_inference_settings("default")
            atom_refs = pretrained_mlip.get_reference_energies(model, reference_type="atom_refs")
            form_elem_refs = pretrained_mlip.get_reference_energies(model, reference_type="form_elem_refs")
            predictor = ParallelMLIPPredictUnit(
                inference_model_path=str(ckpt_path),
                device=device,
                inference_settings=inference_settings,
                atom_refs=atom_refs,
                form_elem_refs=form_elem_refs,
                num_workers=num_workers,
                num_workers_per_node=num_workers_per_node,
            )
        else:
            predict_kwargs = dict(device=device, workers=num_workers)
            if uma_inference_settings is not None:
                predict_kwargs["inference_settings"] = uma_inference_settings
            predictor = pretrained_mlip.get_predict_unit(model, **predict_kwargs)
        super().__init__(predictor, task_name=str(task_name))
