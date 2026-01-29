"""UMA calculator wrapper for PySisyphus providing energy, forces, and Hessians.

Example:
    >>> from pdb2reaction.uma_pysis import uma_pysis
    >>> calc = uma_pysis(charge=0, spin=1, model="uma-s-1p1")
    >>> calc.get_energy(["C", "O"], coords)

For detailed documentation, see: docs/uma_pysis.md
"""

from __future__ import annotations
from typing import Sequence, Optional, Dict, Any, List

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

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import BOHR2ANG, ANG2BOHR, AU2EV
from pysisyphus import run

# ------------ unit conversion constants ----------------------------
EV2AU          = 1.0 / AU2EV                     # eV → Hartree
F_EVAA_2_AU    = EV2AU / ANG2BOHR                # eV Å⁻¹ → Hartree Bohr⁻¹
H_EVAA_2_AU    = EV2AU / ANG2BOHR / ANG2BOHR     # eV Å⁻² → Hartree Bohr⁻²

# Import defaults from centralized configuration
from .defaults import GEOM_KW_DEFAULT, UMA_CALC_KW

# Use centralized UMA calculator defaults (no local redefinition)
CALC_KW = UMA_CALC_KW

# ===================================================================
#                         UMA core wrapper
# ===================================================================
class UMAcore:
    """
    Thin wrapper around fairchem-UMA predict_unit.

    Notes on `workers`
    ------------------
    If `workers > 1`, this wrapper directly instantiates `ParallelMLIPPredictUnit`
    (with `num_workers=workers` and `num_workers_per_node=workers_per_node`) which
    typically does not expose `predictor.model`. In that situation, this wrapper will:
      - skip all `predictor.model` related operations (eval/dropout/etc.)
      - disallow Analytical Hessians (caller should use FD)
      - keep AtomicData/batch on CPU (predictor/worker handles device movement)
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
        # Select device ------------------------------------------------
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device_str = device
        self.device = torch.device(device)

        self.workers = int(workers) if workers is not None else 1
        if self.workers < 1:
            self.workers = 1

        self.workers_per_node = int(workers_per_node) if workers_per_node is not None else 1
        if self.workers_per_node < 1:
            self.workers_per_node = 1

        # If workers>1, we use ParallelMLIPPredictUnit directly.
        self.parallel_predict = self.workers > 1

        self._AtomicData = AtomicData
        self._collater = data_list_collater

        # Load predictor and place it on target device ----------------
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

        # Detect whether we can access an underlying torch model
        self.has_torch_model = hasattr(self.predict, "model") and isinstance(
            getattr(self.predict, "model", None), nn.Module
        )

        # If model is available, put it to eval and disable dropout
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

    # ----------------------------------------------------------------
    def _model_backbone(self):
        """
        Handle both plain and (D)DP-wrapped models.

        Returns
        -------
        backbone module or None
            None if predictor.model is not exposed (e.g. workers>1 predictor).
        """
        if not self.has_torch_model:
            return None
        mdl = self.predict.model
        mod = getattr(mdl, "module", mdl)
        return getattr(mod, "backbone", None)

    # ----------------------------------------------------------------
    def _ase_to_batch(self, atoms: Atoms):
        """
        Convert ASE Atoms → UMA AtomicData (batched).

        In parallel predictor mode (workers>1), keep AtomicData on CPU and let
        the predictor/workers handle device placement.
        """
        backbone = self._model_backbone()

        default_max_neigh = getattr(backbone, "max_neighbors", None) if backbone is not None else None
        default_radius = getattr(backbone, "cutoff", None) if backbone is not None else None

        # AtomicData.from_ase default radius is 6.0 Å; keep a sane fallback.
        if default_radius is None:
            default_radius = 6.0

        max_neigh = self._max_neigh if self._max_neigh is not None else default_max_neigh
        radius = self._radius if self._radius is not None else default_radius
        r_edges = self._r_edges

        atoms.info.update({"charge": self.charge, "spin": self.spin})
        data = self._AtomicData.from_ase(
            atoms,
            max_neigh=max_neigh,
            radius=radius,
            r_edges=r_edges,
        )
        data.dataset = self.task_name

        # collate
        batch = self._collater([data], otf_graph=True)

        # For non-parallel predictor, move batch to requested device here.
        if not self.parallel_predict:
            batch = batch.to(self.device)
        return batch

    # ----------------------------------------------------------------
    def compute(
        self,
        coord_ang: np.ndarray,
        *,
        forces: bool = False,
        hessian: bool = False,
    ) -> Dict[str, Any]:
        """
        Evaluate energy and (optionally) forces/Hessian at given coordinates.

        Parameters
        ----------
        coord_ang : (N, 3) array-like in Å
        forces : bool
            If True, return forces (and enable grads when using a local model).
        hessian : bool
            If True, keep autograd path so the analytical Hessian can be built.

        Returns
        -------
        dict with keys:
          - "energy"  : float, eV
          - "forces"  : np.ndarray, shape (N, 3), eV/Å or None
          - "hessian" : torch.Tensor, shape (N,3,N,3), eV/Å² or None

        Notes
        -----
        If the predictor does not expose `predictor.model` (e.g. workers>1),
        Analytical Hessian is not available. Use FD Hessian in the calculator.
        """
        atoms = Atoms(self.elem, positions=coord_ang)
        batch = self._ase_to_batch(atoms)

        # Parallel predictor mode: avoid touching predictor.model and avoid relying
        # on autograd graphs across process boundaries.
        if self.parallel_predict or (not self.has_torch_model):
            if hessian:
                raise RuntimeError(
                    "Analytical Hessian is not available when predictor workers > 1 "
                    "or when predictor.model is not exposed. Use FiniteDifference Hessian."
                )

            # Let the predictor decide how/when to enable grads for forces.
            res = self.predict.predict(batch)

            energy = float(res["energy"].squeeze().detach().item())
            forces_np = (
                res["forces"].detach().cpu().numpy()
                if forces
                else None
            )
            return {"energy": energy, "forces": forces_np, "hessian": None}

        batch.pos.requires_grad_(True)

        res = self.predict.predict(batch)

        energy = float(res["energy"].squeeze().detach().item())
        forces_np = (
            res["forces"].detach().cpu().numpy()
            if (forces or hessian)
            else None
        )

        if hessian:
            # Disable parameter grads during Hessian evaluation
            p_flags = [p.requires_grad for p in self.predict.model.parameters()]
            for p in self.predict.model.parameters():
                p.requires_grad_(False)

            self.predict.model.train()
            try:
                def e_fn(flat):
                    batch.pos = flat.view(-1, 3)
                    return self.predict.predict(batch)["energy"].squeeze()

                H = torch.autograd.functional.hessian(
                    e_fn, batch.pos.view(-1), vectorize=False
                )
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
class uma_pysis(Calculator):
    """
    PySisyphus-compatible UMA calculator.
    """

    implemented_properties = ["energy", "forces", "hessian"]

    def __init__(
        self,
        *,
        charge: int = CALC_KW["charge"],
        spin: int = CALC_KW["spin"],
        model: str = CALC_KW["model"],
        task_name: str = CALC_KW["task_name"],
        device: str = CALC_KW["device"],
        workers: int = CALC_KW["workers"],
        workers_per_node: int = CALC_KW["workers_per_node"],
        out_hess_torch: bool = CALC_KW["out_hess_torch"],
        max_neigh: Optional[int] = CALC_KW["max_neigh"],
        radius: Optional[float] = CALC_KW["radius"],
        r_edges: bool = CALC_KW["r_edges"],
        freeze_atoms: Optional[Sequence[int]] = CALC_KW["freeze_atoms"],
        hessian_calc_mode: str = CALC_KW["hessian_calc_mode"],
        return_partial_hessian: bool = CALC_KW["return_partial_hessian"],
        hessian_double: bool = CALC_KW["hessian_double"],
        # -------------------------------------------------------------------
        **kwargs,
    ):
        """
        Parameters
        ----------
        hessian_calc_mode : {"Analytical", "FiniteDifference"}, default "FiniteDifference"
            Select how to compute the Hessian in `get_hessian()`.
            - "Analytical": autograd-based analytical Hessian.
            - "FiniteDifference": central-difference Hessian; the matrix is
              assembled on the selected device.
            Note: if `workers > 1`, Analytical is disabled and FD is forced.
        workers : int, default 1
            Predictor worker count.
            If `workers > 1`, this wrapper instantiates `ParallelMLIPPredictUnit`
            directly (and `predictor.model` is not exposed); Hessian is forced
            to FiniteDifference.
        workers_per_node : int, default 1
            Passed as `num_workers_per_node` to `ParallelMLIPPredictUnit` when
            `workers > 1`.
        freeze_atoms : list[int], optional
            Atom indices (0-based). In both modes, DOFs of these atoms are
            treated as frozen.
        return_partial_hessian : bool, default False
            If True, return only the Active-DOF Hessian (submatrix for non-frozen atoms).
            If False, return a full (3N×3N) matrix where frozen-DOF columns are 0.
            Full Hessians avoid shape mismatches with pysisyphus optimizers.
        hessian_double : bool, default True
            If True, assemble/return the Hessian in float64 (double precision).
            Energy and forces are always returned as float64 regardless of this flag.
        """
        super().__init__(charge=charge, mult=spin, **kwargs)
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
        self.hessian_calc_mode = hessian_calc_mode
        if freeze_atoms is None:
            freeze_iter: List[int] = []
        else:
            try:
                freeze_iter = [int(i) for i in list(freeze_atoms)]
            except Exception:
                freeze_iter = []
        self.freeze_atoms = sorted(set(freeze_iter))
        self.return_partial_hessian = bool(return_partial_hessian)
        self.hessian_double = bool(hessian_double)

    # ---------- helpers ---------------------------------------------
    def _ensure_core(self, elem: Sequence[str]):
        if self._core is None:
            self._core = UMAcore(elem, **self._core_kw)

    def _au_hessian(self, H: torch.Tensor):
        """
        Convert Hessian from eV/Å² to Hartree/Bohr² and format the output.

        Parameters
        ----------
        H : torch.Tensor
            Shape (N,3,N,3), on device.

        Returns
        -------
        numpy.ndarray or torch.Tensor
            If `out_hess_torch` is False, returns a NumPy array on CPU.
            Otherwise returns a torch.Tensor on the current device.

        Notes
        -----
        The dtype of the returned Hessian is controlled solely by
        `self.hessian_double`:

        - If True, the result is float64 (double precision).
        - If False, the result stays in the model/native dtype.
        """
        n = H.size(0)
        H = H.view(n * 3, n * 3)
        H = 0.5 * (H + H.T)

        # Unit conversion
        H = H * H_EVAA_2_AU

        if self.hessian_double:
            H = H.to(dtype=torch.float64)

        if self.out_hess_torch:
            return H.detach()
        else:
            return H.detach().cpu().numpy()

    # ---- Common utilities for freeze/active DOF ----------------
    def _active_and_frozen_dof_idx(self, n_atoms: int):
        frozen_set = set(self.freeze_atoms)
        active_atoms = [i for i in range(n_atoms) if i not in frozen_set]
        active_dof_idx = [3 * i + j for i in active_atoms for j in range(3)]
        frozen_dof_idx = [3 * i + j for i in self.freeze_atoms for j in range(3)]
        return active_atoms, active_dof_idx, frozen_dof_idx

    def _zero_frozen_forces_ev(self, F: np.ndarray) -> np.ndarray:
        """Zero forces (eV/Å) on frozen atoms."""
        if (F is None) or (len(self.freeze_atoms) == 0):
            return F
        Fz = F.copy()
        Fz[np.asarray(self.freeze_atoms, dtype=int)] = 0.0
        return Fz

    def _apply_analytical_active_trim(self, H: torch.Tensor) -> torch.Tensor:
        """
        Apply Active‑DOF trimming/column-zeroing to an Analytical Hessian.

        If `return_partial_hessian=True`, reduce to Active‑DOF block.
        Else keep full size and zero columns corresponding to frozen DOF.
        """
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
                H.index_fill_(1, cols, 0.0)  # zero columns of frozen DOF
            return H.view(n_atoms, 3, n_atoms, 3)

    # ---------- Finite-Difference Hessian (GPU/CPU assembly) -------------
    def _build_fd_hessian_gpu(
        self,
        elem: Sequence[str],
        coord_ang: np.ndarray,
        *,
        eps_ang: float = 1.0e-3,  # central-difference step (Å)
    ) -> Dict[str, Any]:
        """
        Assemble the Hessian by central differences of forces, with the matrix
        stored and operated on the selected device:

            H_[:, k] = -(F(x + h e_k) - F(x - h e_k)) / (2 h)

        where F is the flattened force vector (eV/Å). Columns corresponding to
        `freeze_atoms` DOF are skipped. If `return_partial_hessian` is True, the
        returned matrix contains only the Active-DOF block (non-frozen atoms);
        otherwise the full (3N×3N) matrix is returned.

        Returns
        -------
        dict with keys:
          - "energy"  : float, eV
          - "forces"  : np.ndarray, shape (N, 3), eV/Å
          - "hessian" : torch.Tensor, shape (N_out,3,N_out,3), eV/Å² on device
        """
        self._ensure_core(elem)
        core = self._core
        assert core is not None

        dev = core.device

        n_atoms = len(elem)
        dof = n_atoms * 3

        # Active atoms / DOF
        frozen_set = set(self.freeze_atoms)
        active_atoms = [i for i in range(n_atoms) if i not in frozen_set]
        active_dof_idx = [3 * i + j for i in active_atoms for j in range(3)]

        # Base point evaluation
        res0 = core.compute(coord_ang, forces=True, hessian=False)
        energy0_eV = res0["energy"]
        F0 = res0["forces"]  # (N,3) eV/Å, native numpy dtype

        # Assemble Hessian in the same dtype as the model forces unless
        # hessian_double is requested.
        force_dtype = torch.from_numpy(F0).dtype
        hessian_dtype = torch.float64 if self.hessian_double else force_dtype

        # Device-side Hessian storage (2D for easy column insertion)
        H = torch.zeros((dof, dof), device=dev, dtype=hessian_dtype)

        # Host-side work arrays for coordinate perturbations
        coord_plus = coord_ang.copy()
        coord_minus = coord_ang.copy()

        # Compute columns only for active DOF
        for k in active_dof_idx:
            a = k // 3
            c = k % 3

            # x + h
            coord_plus[a, c] = coord_ang[a, c] + eps_ang
            res_p = core.compute(coord_plus, forces=True, hessian=False)
            Fp = res_p["forces"].reshape(-1)  # (3N,) eV/Å

            # x - h
            coord_minus[a, c] = coord_ang[a, c] - eps_ang
            res_m = core.compute(coord_minus, forces=True, hessian=False)
            Fm = res_m["forces"].reshape(-1)  # (3N,) eV/Å

            # Column on device
            Fp_t = torch.from_numpy(Fp).to(dev, dtype=hessian_dtype)
            Fm_t = torch.from_numpy(Fm).to(dev, dtype=hessian_dtype)
            col = -(Fp_t - Fm_t) / (2.0 * eps_ang)  # (3N,) eV/Å²

            H[:, k] = col

            # Restore
            coord_plus[a, c] = coord_ang[a, c]
            coord_minus[a, c] = coord_ang[a, c]

        # Reduce to Active-DOF block if requested
        if self.return_partial_hessian:
            idx = torch.tensor(active_dof_idx, device=dev, dtype=torch.long)
            H = H.index_select(0, idx).index_select(1, idx)
            n_active_atoms = len(active_atoms)
            H = H.view(n_active_atoms, 3, n_active_atoms, 3)
        else:
            H = H.view(n_atoms, 3, n_atoms, 3)

        return {"energy": energy0_eV, "forces": F0, "hessian": H}

    # ---------- PySisyphus API --------------------------------------
    def get_energy(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        res = self._core.compute(coord_ang, forces=False, hessian=False)
        return {"energy": res["energy"] * EV2AU}

    def get_forces(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        res = self._core.compute(coord_ang, forces=True, hessian=False)

        # Zero forces on frozen atoms (in eV/Å before conversion)
        F_ev = self._zero_frozen_forces_ev(res["forces"])

        return {
            "energy": res["energy"] * EV2AU,
            "forces": (np.asarray(F_ev, dtype=np.float64) * F_EVAA_2_AU).reshape(-1),
        }

    def get_hessian(self, elem, coords):
        """
        Compute the Hessian according to `hessian_calc_mode`.

        Modes
        -----
        - "Analytical": autograd-based analytical Hessian.
          If a CUDA out-of-memory error occurs, a clear message is raised
          instructing you to switch to finite differences.
          Active‑DOF trimming/column‑zeroing is applied to match FD semantics.
        - "FiniteDifference": central-difference Hessian assembled on the selected device.
          * Columns for `freeze_atoms` DOF are skipped.
          * If `return_partial_hessian` is True, return the Active-DOF submatrix;
            otherwise return a full-size matrix (frozen columns remain zero).

        If `hessian_calc_mode` is falsy or unrecognized, FiniteDifference is used.

        Note
        ---------
        If the underlying predictor does not expose `predictor.model` (typically
        when `workers > 1` was used), this method always uses FiniteDifference.
        """
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG

        core = self._core
        assert core is not None

        # Force FD Hessian when predictor.model is not accessible (e.g. workers>1)
        force_fd = (core.parallel_predict or (not core.has_torch_model))

        mode = (self.hessian_calc_mode or "FiniteDifference").strip().lower()
        if (not force_fd) and (mode in ("analytical", "analytic")):
            try:
                res = self._core.compute(
                    coord_ang,
                    forces=True,
                    hessian=True,
                )
            except (torch.cuda.OutOfMemoryError, RuntimeError) as e:
                msg = str(e).lower()
                if "out of memory" in msg and "cuda" in msg:
                    raise RuntimeError(
                        "Analytical Hessian computation failed due to CUDA out-of-memory. "
                        "Your GPU memory appears to be limited. Please switch to the finite-"
                        "difference Hessian by specifying `--hessian-calc-mode FiniteDifference` "
                        "in the external CLI, or `hessian_calc_mode=\"FiniteDifference\"` in this calculator."
                    ) from e
                raise

            # Zero forces for frozen atoms (eV/Å)
            res_forces_ev = self._zero_frozen_forces_ev(res["forces"])

            # Apply Active‑DOF trimming/column-zeroing for Analytical
            H = self._apply_analytical_active_trim(res["hessian"])

            return {
                "energy": res["energy"] * EV2AU,
                "forces": (np.asarray(res_forces_ev, dtype=np.float64) * F_EVAA_2_AU).reshape(-1),
                "hessian": self._au_hessian(H),
            }

        # FiniteDifference (default/fallback, also forced for workers>1 predictors)
        res = self._build_fd_hessian_gpu(elem, coord_ang)

        # Zero forces for frozen atoms (eV/Å)
        res_forces_ev = self._zero_frozen_forces_ev(res["forces"])

        return {
            "energy": res["energy"] * EV2AU,
            "forces": (np.asarray(res_forces_ev, dtype=np.float64) * F_EVAA_2_AU).reshape(-1),
            "hessian": self._au_hessian(res["hessian"]),
        }


class uma_ase(FAIRChemCalculator):
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
        predictor = pretrained_mlip.get_predict_unit(
            model,
            device=device,
            workers=int(workers),
            workers_per_node=int(workers_per_node),
        )
        super().__init__(predictor, task_name=str(task_name))