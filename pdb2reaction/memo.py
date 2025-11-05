# pdb2reaction/uma_pysis.py
"""
UMA Calculator Wrapper for PySisyphus.
Return Energy, Force, Analytic Hessian.

Notes (2025-11-04):
- Default compute precision is float32 (FP32) to reduce GPU memory and avoid CUDA OOM.
- You can now choose precision via `dtype` from `uma_pysis`, which is propagated to both
  torch and numpy so the whole pipeline uses a unified precision:
    dtype ∈ {"float32","float64","bfloat16","float16"} or the corresponding torch.dtype.
- Removed .double() casts; all unit scalings avoid implicit FP64 upcasts on GPU.
"""

from typing import Sequence, Optional, Union, Dict, Any

import numpy as np
import torch
import torch.nn as nn
from ase import Atoms

from fairchem.core import pretrained_mlip
from fairchem.core.datasets.atomic_data import AtomicData
from fairchem.core.datasets import data_list_collater

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import BOHR2ANG, ANG2BOHR, AU2EV
from pysisyphus import run

# ------------ unit conversion constants (scalar Python floats) ----------------
EV2AU          = 1.0 / AU2EV                     # eV → Hartree
F_EVAA_2_AU    = EV2AU / ANG2BOHR                # eV Å⁻¹ → Hartree Bohr⁻¹
H_EVAA_2_AU    = EV2AU / ANG2BOHR / ANG2BOHR     # eV Å⁻² → Hartree Bohr⁻²


# ===================================================================
#                         UMA core wrapper
# ===================================================================
class UMAcore:
    """Thin wrapper around fairchem‑UMA predict_unit with unified precision."""

    def __init__(
        self,
        elem: Sequence[str],
        *,
        charge: int = 0,
        spin: int = 1,
        model: str = "uma-s-1p1",
        task_name: str = "omol",
        device: str = "auto",
        max_neigh: Optional[int] = None,
        radius:    Optional[float] = None,
        r_edges:   bool = False,
        dtype: Union[str, torch.dtype] = "float32",
    ):
        # -------- device selection -----------------------------------
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device_str = device
        self.device = torch.device(device)

        # -------- dtype normalization --------------------------------
        self.dtype: torch.dtype = self._normalize_dtype(dtype)
        self.np_dtype: np.dtype = self._numpy_dtype(self.dtype)

        # Set torch default dtype only for (float32|float64)
        # (PyTorch does not allow setting default to fp16/bf16)
        self._prev_default_torch_dtype = torch.get_default_dtype()
        if self.dtype in (torch.float32, torch.float64) and self._prev_default_torch_dtype != self.dtype:
            torch.set_default_dtype(self.dtype)

        # -------- fairchem predictors & helpers ----------------------
        self._AtomicData = AtomicData
        self._collater   = data_list_collater

        # Load predictor and move model to target device/dtype
        self.predict = pretrained_mlip.get_predict_unit(model, device=self.device_str)
        self.predict.model.to(self.device, dtype=self.dtype)
        self.predict.model.eval()
        # Hard-disable any dropout
        for m in self.predict.model.modules():
            if isinstance(m, nn.Dropout):
                m.p = 0.0

        self.elem      = [e.capitalize() for e in elem]
        self.charge    = charge
        self.spin      = spin
        self.task_name = task_name

        self._max_neigh_user = max_neigh
        self._radius_user    = radius
        self._r_edges_user   = r_edges

    # ----------------------------------------------------------------
    @staticmethod
    def _normalize_dtype(dt: Union[str, torch.dtype]) -> torch.dtype:
        if isinstance(dt, torch.dtype):
            return dt
        if isinstance(dt, str):
            s = dt.lower()
            if s in ("float32", "fp32", "single"): return torch.float32
            if s in ("float64", "double", "fp64"): return torch.float64
            if s in ("bfloat16", "bf16"):          return torch.bfloat16
            if s in ("float16", "half", "fp16"):   return torch.float16
        raise ValueError(f"Unsupported dtype: {dt}")

    @staticmethod
    def _numpy_dtype(td: torch.dtype) -> np.dtype:
        """Map torch dtype → numpy dtype (best-effort for bf16)."""
        if td == torch.float32:
            return np.float32
        if td == torch.float64:
            return np.float64
        if td == torch.float16:
            return np.float16
        if td == torch.bfloat16:
            # numpy bfloat16 availability varies by version; fall back to float16 if unavailable
            try:
                return np.dtype("bfloat16")  # type: ignore[arg-type]
            except Exception:
                return np.float16
        # Shouldn't reach here due to normalization
        return np.float32

    # ----------------------------------------------------------------
    def _model_backbone(self):
        """Access backbone robustly for (D)DP-wrapped models."""
        mdl = self.predict.model
        mod = getattr(mdl, "module", mdl)
        return mod.backbone

    # ----------------------------------------------------------------
    def _ase_to_batch(self, atoms: Atoms):
        """Convert ASE Atoms → UMA AtomicData(Batch), respecting device/dtype."""
        backbone = self._model_backbone()
        default_max_neigh = getattr(backbone, "max_neighbors", None)
        default_radius    = getattr(backbone, "cutoff", None)

        max_neigh = self._max_neigh_user if self._max_neigh_user is not None else default_max_neigh
        radius    = self._radius_user    if self._radius_user    is not None else default_radius
        r_edges   = self._r_edges_user

        atoms.info.update({"charge": self.charge, "spin": self.spin})

        # Build AtomicData on target device
        data = self._AtomicData.from_ase(
            atoms,
            max_neigh=max_neigh,
            radius   =radius,
            r_edges  =r_edges,
        ).to(self.device)

        data.dataset = self.task_name

        # Collate -> Batch (otf_graph True replicates fairchem default behavior)
        batch = self._collater([data], otf_graph=True).to(self.device)

        # Ensure floating fields use the desired dtype (at least positions)
        if hasattr(batch, "pos") and isinstance(batch.pos, torch.Tensor) and batch.pos.dtype != self.dtype:
            batch.pos = batch.pos.to(self.dtype)

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
        coord_ang : (N,3) Å, numpy array (will be cast to the chosen np dtype)
        forces / hessian : toggles
        Returns dict with keys:
          - energy  (float, eV)
          - forces  (np.ndarray, eV/Å) or None
          - hessian (torch.Tensor, eV/Å^2) or None
        """
        # Make sure positions we pass to ASE are of unified numpy dtype.
        pos_ang = np.asarray(coord_ang, dtype=self.np_dtype)
        atoms = Atoms(self.elem, positions=pos_ang)

        batch = self._ase_to_batch(atoms)

        # Require gradients only when needed (saves memory)
        need_grad = bool(forces or hessian)
        pos = batch.pos.detach().clone()
        if pos.dtype != self.dtype:
            pos = pos.to(self.dtype)
        pos.requires_grad_(need_grad)
        batch.pos = pos

        if need_grad:
            res = self.predict.predict(batch)
        else:
            with torch.no_grad():
                res = self.predict.predict(batch)

        # Energy (tensor) → Python float
        energy_t = res["energy"].squeeze()
        energy   = float(energy_t.detach().item())

        # Forces (remain FP on device → CPU numpy of unified dtype)
        forces_np = None
        if forces or hessian:
            f_t = res["forces"]
            if f_t.dtype != self.dtype:
                f_t = f_t.to(self.dtype)
            forces_np = f_t.detach().cpu().numpy().astype(self.np_dtype, copy=False)

        # Hessian (computed as torch tensor in the chosen torch dtype)
        if hessian:
            # Temporarily freeze model params to avoid grad tracking
            p_flags = []
            for p in self.predict.model.parameters():
                p_flags.append(p.requires_grad)
                p.requires_grad_(False)

            self.predict.model.train()  # BN uses batch stats; dropout was disabled above

            def e_fn(flat):
                # flat: 3N vector (dtype == self.dtype)
                batch.pos = flat.view(-1, 3)
                out = self.predict.predict(batch)["energy"].squeeze()
                return out

            H_flat = torch.autograd.functional.hessian(
                e_fn,
                batch.pos.view(-1),
            )
            hess_t = H_flat.view(len(atoms), 3, len(atoms), 3).detach()

            self.predict.model.eval()
            # Restore flags
            for p, flag in zip(self.predict.model.parameters(), p_flags):
                p.requires_grad_(flag)

            if self.device.type == "cuda":
                torch.cuda.empty_cache()
        else:
            hess_t = None

        return {"energy": energy, "forces": forces_np, "hessian": hess_t}


# ===================================================================
#                    PySisyphus calculator class
# ===================================================================
class uma_pysis(Calculator):
    """PySisyphus‑compatible UMA calculator with unified torch/numpy precision."""

    implemented_properties = ["energy", "forces", "hessian"]

    def __init__(
        self,
        *,
        charge: int = 0,
        spin: int = 1,
        model: str = "uma-s-1p1",
        task_name: str = "omol",
        device: str = "auto",
        out_hess_torch: bool = False,
        max_neigh: Optional[int] = None,
        radius:    Optional[float] = None,
        r_edges:   bool = False,
        dtype: Union[str, torch.dtype] = "float32",
        **kwargs,
    ):
        super().__init__(charge=charge, mult=spin, **kwargs)

        # Normalize dtype once here so we can also unify numpy-side precision
        self._torch_dtype: torch.dtype = UMAcore._normalize_dtype(dtype)
        self._np_dtype:    np.dtype    = UMAcore._numpy_dtype(self._torch_dtype)

        # Prebuild numpy-typed unit constants to avoid implicit upcasts
        self._NP_EV2AU        = np.array(EV2AU, dtype=self._np_dtype)
        self._NP_F_EVAA_2_AU  = np.array(F_EVAA_2_AU, dtype=self._np_dtype)
        self._NP_BOHR2ANG     = np.array(BOHR2ANG, dtype=self._np_dtype)
        self._NP_ANG2BOHR     = np.array(ANG2BOHR, dtype=self._np_dtype)

        self._core: Optional[UMAcore] = None
        self._core_kw = dict(
            charge=charge,
            spin=spin,
            model=model,
            task_name=task_name,
            device=device,
            max_neigh=max_neigh,
            radius=radius,
            r_edges=r_edges,
            dtype=self._torch_dtype,  # ensure torch dtype object is used downstream
        )
        self.out_hess_torch = out_hess_torch

    # ---------- helpers ---------------------------------------------
    def _ensure_core(self, elem: Sequence[str]):
        if self._core is None:
            self._core = UMAcore(elem, **self._core_kw)

    def _au_energy(self, e_eV: float) -> float:
        # Return Python float (PySisyphus expects a scalar), but scale in unified numpy dtype
        return float(np.array(e_eV, dtype=self._np_dtype) * self._NP_EV2AU)

    def _au_forces(self, f_eV_A: np.ndarray) -> np.ndarray:
        # Flatten (3N,) array in unified numpy dtype
        f = np.asarray(f_eV_A, dtype=self._np_dtype)
        return (f * self._NP_F_EVAA_2_AU).reshape(-1).astype(self._np_dtype, copy=False)

    def _au_hessian(self, H_eV_AA: torch.Tensor):
        """
        H_eV_AA: (N,3,N,3) tensor on device with chosen torch dtype.
        - reshape to (3N,3N)
        - scale to atomic units using a same-dtype/device scalar to avoid FP64 promotion
        - symmetrize
        - return numpy (unified np dtype) or torch (unified torch dtype)
        """
        n = H_eV_AA.size(0)
        H = H_eV_AA.view(n * 3, n * 3)

        scale = H.new_tensor(H_EVAA_2_AU)  # same dtype & device as H
        H = H * scale
        H = 0.5 * (H + H.T)

        if not self.out_hess_torch:
            # to CPU → numpy with unified numpy dtype
            return (
                H.detach()
                 .to(dtype=self._torch_dtype)  # no-op if already same
                 .cpu()
                 .numpy()
                 .astype(self._np_dtype, copy=False)
            )
        else:
            # Keep as torch tensor in the unified torch dtype on the same device
            return H.detach().to(dtype=self._torch_dtype)

    # ---------- PySisyphus API --------------------------------------
    def get_energy(self, elem, coords):
        self._ensure_core(elem)
        # coords are Bohr in PySisyphus; convert to Å in unified numpy dtype
        coord_ang = np.asarray(coords, dtype=self._np_dtype).reshape(-1, 3) * self._NP_BOHR2ANG
        res = self._core.compute(coord_ang, forces=False, hessian=False)
        return {"energy": self._au_energy(res["energy"])}

    def get_forces(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=self._np_dtype).reshape(-1, 3) * self._NP_BOHR2ANG
        res = self._core.compute(coord_ang, forces=True, hessian=False)
        return {
            "energy": self._au_energy(res["energy"]),
            "forces": self._au_forces(res["forces"]),
        }

    def get_hessian(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=self._np_dtype).reshape(-1, 3) * self._NP_BOHR2ANG
        res = self._core.compute(coord_ang, forces=True, hessian=True)
        return {
            "energy": self._au_energy(res["energy"]),
            "forces": self._au_forces(res["forces"]),
            "hessian": self._au_hessian(res["hessian"]),
        }


# ---------- CLI ----------------------------------------
def run_pysis():
    """Enable `uma_pysis input.yaml`"""
    run.CALC_DICT["uma_pysis"] = uma_pysis
    run.run()


if __name__ == "__main__":
    run_pysis()

# pdb2reaction/uma_pysis.py
"""
UMA Calculator Wrapper for PySisyphus.
Return Energy, Force, Analytic Hessian.
"""
from typing import List, Sequence, Optional

import numpy as np
import torch
import torch.nn as nn
from ase import Atoms

from fairchem.core import pretrained_mlip
from fairchem.core.datasets.atomic_data import AtomicData
from fairchem.core.datasets import data_list_collater

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import BOHR2ANG, ANG2BOHR, AU2EV
from pysisyphus import run

# ------------ unit conversion constants ----------------------------
EV2AU          = 1.0 / AU2EV                     # eV → Hartree
F_EVAA_2_AU    = EV2AU / ANG2BOHR                # eV Å⁻¹ → Hartree Bohr⁻¹
H_EVAA_2_AU    = EV2AU / ANG2BOHR / ANG2BOHR     # eV Å⁻² → Hartree Bohr⁻²


# ===================================================================
#                         UMA core wrapper
# ===================================================================
class UMAcore:
    """Thin wrapper around fairchem‑UMA predict_unit."""

    def __init__(
        self,
        elem: Sequence[str],
        *,
        charge: int = 0,
        spin: int = 1,
        model: str = "uma-s-1p1",
        task_name: str = "omol",
        device: str = "auto",
        max_neigh: Optional[int] = None,
        radius:    Optional[float] = None,
        r_edges:   bool = False,
    ):
        # Select device ------------------------------------------------
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device_str = device
        self.device = torch.device(device)

        self._AtomicData = AtomicData
        self._collater   = data_list_collater

        # Predictor in double precision -------------------------------
        self.predict = pretrained_mlip.get_predict_unit(model, device=self.device_str)
        self.predict.model.eval()
        for m in self.predict.model.modules():
            if isinstance(m, nn.Dropout):
                m.p = 0.0

        self.elem      = [e.capitalize() for e in elem]
        self.charge    = charge
        self.spin      = spin
        self.task_name = task_name

        self._max_neigh_user = max_neigh
        self._radius_user    = radius
        self._r_edges_user   = r_edges

    # ----------------------------------------------------------------
    def _ase_to_batch(self, atoms: Atoms):
        """Convert ASE Atoms → UMA AtomicData(Batch)."""

        default_max_neigh = self.predict.model.module.backbone.max_neighbors
        default_radius    = self.predict.model.module.backbone.cutoff

        max_neigh = self._max_neigh_user if self._max_neigh_user is not None else default_max_neigh
        radius    = self._radius_user    if self._radius_user    is not None else default_radius
        r_edges   = self._r_edges_user

        atoms.info.update({"charge": self.charge, "spin": self.spin})
        data = self._AtomicData.from_ase(
            atoms,
            max_neigh=max_neigh,
            radius   =radius,
            r_edges  =r_edges,
        ).to(self.device)
        data.dataset = self.task_name
        return self._collater([data], otf_graph=True).to(self.device)

    # ----------------------------------------------------------------
    def compute(
        self,
        coord_ang: np.ndarray,
        *,
        forces: bool = False,
        hessian: bool = False,
    ):
        """
        coord_ang : (N,3) Å
        forces / hessian : return toggles
        Returns dict with keys energy (eV), forces (eV/Å), hessian (torch)
        """
        atoms = Atoms(self.elem, positions=coord_ang)
        batch = self._ase_to_batch(atoms)

        batch.pos = batch.pos.detach().clone().requires_grad_(True)

        res      = self.predict.predict(batch)
        energy   = float(res["energy"].double().squeeze().item())
        forces_np = res["forces"].double().cpu().numpy() if (forces or hessian) else None

        if hessian:
            self.predict.model.train()
            def e_fn(flat):
                batch.pos = flat.view(-1, 3)
                return self.predict.predict(batch)["energy"].double().squeeze()
            H_flat   = torch.autograd.functional.hessian(e_fn, batch.pos.view(-1))
            hess_t   = H_flat.view(len(atoms), 3, len(atoms), 3).detach()
            self.predict.model.eval()
        else:
            hess_t = None

        return {"energy": energy, "forces": forces_np, "hessian": hess_t}


# ===================================================================
#                    PySisyphus calculator class
# ===================================================================
class uma_pysis(Calculator):
    """PySisyphus‑compatible UMA calculator."""

    implemented_properties = ["energy", "forces", "hessian"]

    def __init__(
        self,
        *,
        charge: int = 0,
        spin: int = 1,
        model: str = "uma-s-1p1",
        task_name: str = "omol",
        device: str = "auto",
        out_hess_torch: bool = False,
        max_neigh: Optional[int] = None,
        radius:    Optional[float] = None,
        r_edges:   bool = False,
        **kwargs,
    ):
        super().__init__(charge=charge, mult=spin, **kwargs)
        self._core: Optional[UMAcore] = None
        self._core_kw = dict(
            charge=charge,
            spin=spin,
            model=model,
            task_name=task_name,
            device=device,
            max_neigh=max_neigh,
            radius=radius,
            r_edges=r_edges,
        )
        self.out_hess_torch = out_hess_torch

    # ---------- helpers ---------------------------------------------
    def _ensure_core(self, elem: Sequence[str]):
        if self._core is None:
            self._core = UMAcore(elem, **self._core_kw)

    @staticmethod
    def _au_energy(e_eV: float) -> float:
        return e_eV * EV2AU

    @staticmethod
    def _au_forces(f_eV_A: np.ndarray) -> np.ndarray:
        return (f_eV_A * F_EVAA_2_AU).reshape(-1)

    def _au_hessian(self, H_eV_AA: torch.Tensor):
        H = H_eV_AA.view(H_eV_AA.size(0)*3, H_eV_AA.size(2)*3) * H_EVAA_2_AU
        H = 0.5 * (H + H.T)
        return H.detach().cpu().numpy() if not self.out_hess_torch else H.detach()

    # ---------- PySisyphus API --------------------------------------
    def get_energy(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        res = self._core.compute(coord_ang)
        return {"energy": self._au_energy(res["energy"])}

    def get_forces(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        res = self._core.compute(coord_ang, forces=True)
        return {
            "energy": self._au_energy(res["energy"]),
            "forces": self._au_forces(res["forces"]),
        }

    def get_hessian(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        res = self._core.compute(coord_ang, forces=True, hessian=True)
        return {
            "energy": self._au_energy(res["energy"]),
            "forces": self._au_forces(res["forces"]),
            "hessian": self._au_hessian(res["hessian"]),
        }


# ---------- CLI ----------------------------------------
def run_pysis():
    """Enable `uma_pysis input.yaml`"""
    run.CALC_DICT["uma_pysis"] = uma_pysis
    run.run()


if __name__ == "__main__":
    run_pysis()


1番目はfloat32での計算をするように固定したスクリプト、2番目はfloat64でのもともとのスクリプトです。1番目では、精度を切り替えるオプションがありますが、可読性向上のために、float32固定にしたい。

変更は最低限にしつつ（float32での処理以外は変更なしで）、2番目のスクリプトを修正してください。

省略せずに完全なスクリプトを出力してください。