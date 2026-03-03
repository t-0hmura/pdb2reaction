# pdb2reaction/backends/uma.py

"""
UMA (fairchem) backend for pdb2reaction.

Provides ``UMACalculator`` (PySisyphus-compatible) and ``UMAASECalculator``
(ASE-compatible, for DMF path optimization).
"""

from __future__ import annotations

import time
from typing import Any, Dict, List, Optional, Sequence

import click
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

from pysisyphus.constants import BOHR2ANG, ANG2BOHR, AU2EV

from .base import MLIPCalculator, BackendError

# ---- Unit conversion constants ----
EV2AU = 1.0 / AU2EV
F_EVAA_2_AU = EV2AU / ANG2BOHR
H_EVAA_2_AU = EV2AU / ANG2BOHR / ANG2BOHR


# ===================================================================
#                         UMA core wrapper
# ===================================================================
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
        model: str = "uma-s-1p1",
        task_name: str = "omol",
        device: str = "auto",
        workers: int = 1,
        workers_per_node: int = 1,
        max_neigh: Optional[int] = None,
        radius: Optional[float] = None,
        r_edges: bool = False,
    ):
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device_str = device
        self.device = torch.device(device)

        self.workers = max(int(workers or 1), 1)
        self.workers_per_node = max(int(workers_per_node or 1), 1)
        self.parallel_predict = self.workers > 1

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
            inference_settings = guess_inference_settings("default")
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
            self.predict = pretrained_mlip.get_predict_unit(
                model, device=self.device_str, workers=self.workers
            )

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
        data = self._AtomicData.from_ase(atoms, max_neigh=max_neigh, radius=radius, r_edges=r_edges)
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


# ===================================================================
#                    PySisyphus calculator class
# ===================================================================
class UMACalculator(MLIPCalculator):
    """UMA (fairchem) backend — full feature set.

    Extends ``MLIPCalculator`` with UMA-specific features:
    GPU-accelerated FD Hessian, analytical Hessian (autograd),
    multi-worker inference, and VRAM logging.
    """

    def __init__(
        self,
        *,
        model: str = "uma-s-1p1",
        task_name: str = "omol",
        workers: int = 1,
        workers_per_node: int = 1,
        max_neigh: Optional[int] = None,
        radius: Optional[float] = None,
        r_edges: bool = False,
        out_hess_torch: bool = True,
        print_vram: bool = True,
        # Base class parameters
        charge: int = 0,
        spin: int = 1,
        device: str = "auto",
        freeze_atoms: Optional[Sequence[int]] = None,
        hessian_calc_mode: str = "FiniteDifference",
        return_partial_hessian: bool = False,
        hessian_double: bool = True,
        print_timing: bool = True,
        **kwargs,
    ):
        super().__init__(
            charge=charge,
            spin=spin,
            device=device,
            freeze_atoms=freeze_atoms,
            hessian_calc_mode=hessian_calc_mode,
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
            workers=workers,
            workers_per_node=workers_per_node,
            max_neigh=max_neigh,
            radius=radius,
            r_edges=r_edges,
        )
        self.out_hess_torch = out_hess_torch
        self.print_vram = bool(print_vram)

    def _ensure_core(self, elem: Sequence[str]):
        if self._core is None:
            self._core = UMAcore(elem, **self._core_kw)

    def _supports_analytical_hessian(self) -> bool:
        return True  # Conditional on workers == 1 at runtime

    # ------------------------------------------------------------------
    # Subclass hooks are NOT used for UMA; we override get_* directly
    # to retain the GPU FD Hessian and torch-based analytical Hessian.
    # ------------------------------------------------------------------

    def _compute_energy_forces_ev(self, elem, coord_ang):
        self._ensure_core(elem)
        res = self._core.compute(coord_ang, forces=True, hessian=False)
        return res["energy"], res["forces"]

    # ------------------------------------------------------------------
    # Hessian helpers (UMA-specific)
    # ------------------------------------------------------------------

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
        if self.print_timing:
            click.echo(f"[HessianTiming] mode: {mode_label} | elapsed: {mode_elapsed_s:.2f} s")
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
            click.echo(
                f"[HessianVRAM] total={total_vram:.3f} GB | "
                f"peak_allocated={peak_alloc:.3f} GB | "
                f"peak_reserved={peak_reserved:.3f} GB | "
                f"remaining={remaining_vram:.3f} GB\n"
            )

    def _au_hessian(self, H: torch.Tensor):
        """Convert Hessian from eV/Å² to Hartree/Bohr² (torch version)."""
        n = H.size(0)
        H = H.view(n * 3, n * 3)
        _t = H.T.clone()
        H.add_(_t).mul_(0.5)
        del _t
        H.mul_(H_EVAA_2_AU)
        if self.hessian_double:
            H = H.to(dtype=torch.float64)
        if self.out_hess_torch:
            return H.detach()
        else:
            return H.detach().cpu().numpy()

    def _apply_analytical_active_trim(self, H: torch.Tensor) -> torch.Tensor:
        """Apply Active-DOF trimming/column-zeroing to an Analytical Hessian."""
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
        H = torch.zeros((dof, dof), device=dev, dtype=hessian_dtype)

        coord_plus = coord_ang.copy()
        coord_minus = coord_ang.copy()

        for k in active_dof_idx:
            a = k // 3
            c = k % 3
            coord_plus[a, c] = coord_ang[a, c] + eps_ang
            res_p = core.compute(coord_plus, forces=True, hessian=False)
            Fp = res_p["forces"].reshape(-1)

            coord_minus[a, c] = coord_ang[a, c] - eps_ang
            res_m = core.compute(coord_minus, forces=True, hessian=False)
            Fm = res_m["forces"].reshape(-1)

            Fp_t = torch.from_numpy(Fp).to(dev, dtype=hessian_dtype)
            Fm_t = torch.from_numpy(Fm).to(dev, dtype=hessian_dtype)
            col = -(Fp_t - Fm_t) / (2.0 * eps_ang)
            H[:, k] = col

            coord_plus[a, c] = coord_ang[a, c]
            coord_minus[a, c] = coord_ang[a, c]

        if self.return_partial_hessian:
            idx = torch.tensor(active_dof_idx, device=dev, dtype=torch.long)
            H = H.index_select(0, idx).index_select(1, idx)
            n_active_atoms = len(active_atoms)
            H = H.view(n_active_atoms, 3, n_active_atoms, 3)
        else:
            H = H.view(n_atoms, 3, n_atoms, 3)

        return {"energy": energy0_eV, "forces": F0, "hessian": H}

    # ------------------------------------------------------------------
    # PySisyphus API (override base to use GPU path + torch Hessian)
    # ------------------------------------------------------------------

    def get_energy(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        res = self._core.compute(coord_ang, forces=False, hessian=False)
        return {"energy": res["energy"] * EV2AU}

    def get_forces(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        res = self._core.compute(coord_ang, forces=True, hessian=False)
        F_ev = self._zero_frozen_forces_ev(res["forces"])
        return {
            "energy": res["energy"] * EV2AU,
            "forces": (np.asarray(F_ev, dtype=np.float64) * F_EVAA_2_AU).reshape(-1),
        }

    def get_hessian(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG

        core = self._core
        assert core is not None
        force_fd = (core.parallel_predict or (not core.has_torch_model))

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

        mode = (self.hessian_calc_mode or "FiniteDifference").strip().lower()
        if (not force_fd) and (mode in ("analytical", "analytic")):
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

            res_forces_ev = self._zero_frozen_forces_ev(res["forces"])
            H = self._apply_analytical_active_trim(res["hessian"])
            out = {
                "energy": res["energy"] * EV2AU,
                "forces": (np.asarray(res_forces_ev, dtype=np.float64) * F_EVAA_2_AU).reshape(-1),
                "hessian": self._au_hessian(H),
            }
        else:
            t0 = time.perf_counter()
            res = self._build_fd_hessian_gpu(elem, coord_ang)
            mode_elapsed_s = time.perf_counter() - t0

            res_forces_ev = self._zero_frozen_forces_ev(res["forces"])
            out = {
                "energy": res["energy"] * EV2AU,
                "forces": (np.asarray(res_forces_ev, dtype=np.float64) * F_EVAA_2_AU).reshape(-1),
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
        return out


# ===================================================================
#                     ASE calculator class
# ===================================================================
class UMAASECalculator(FAIRChemCalculator):
    """UMA ASE calculator for DMF path optimization."""

    def __init__(
        self,
        *,
        model: str = "uma-s-1p1",
        device: str = "auto",
        task_name: str = "omol",
        workers: int = 1,
        workers_per_node: int = 1,
    ):
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        num_workers = max(int(workers or 1), 1)
        num_workers_per_node = max(int(workers_per_node or 1), 1)

        if num_workers > 1:
            if (ParallelMLIPPredictUnit is None) or (guess_inference_settings is None):
                raise ImportError(
                    "workers>1 requested, but ParallelMLIPPredictUnit/guess_inference_settings "
                    "could not be imported from fairchem."
                )
            ckpt_path = pretrained_mlip.pretrained_checkpoint_path_from_name(model)
            inference_settings = guess_inference_settings("default")
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
            predictor = pretrained_mlip.get_predict_unit(
                model, device=device, workers=num_workers,
            )
        super().__init__(predictor, task_name=str(task_name))
