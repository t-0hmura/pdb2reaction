# pdb2reaction/uma_pysis.py
"""
UMA Calculator Wrapper for PySisyphus.
Return Energy, Force, Analytic Hessian.
"""

from __future__ import annotations
from typing import List, Sequence, Optional, Dict, Any

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

# Unified host/GPU dtypes for this implementation
_ndtype = np.float32
_tdtype = torch.float32


# ===================================================================
#                         UMA core wrapper
# ===================================================================
class UMAcore:
    """Thin wrapper around fairchem-UMA predict_unit."""

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

        # Load predictor and place it on *target device & dtype* ------
        self.predict = pretrained_mlip.get_predict_unit(model, device=self.device_str)
        self.predict.model.to(self.device, dtype=_tdtype)
        self.predict.model.eval()
        # Disable dropout hard
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
    def _model_backbone(self):
        """Handle both plain and (D)DP-wrapped models."""
        mdl = self.predict.model
        mod = getattr(mdl, "module", mdl)
        return mod.backbone

    # ----------------------------------------------------------------
    def _ase_to_batch(self, atoms: Atoms):
        """Convert ASE Atoms → UMA AtomicData(Batch)."""
        backbone = self._model_backbone()
        default_max_neigh = getattr(backbone, "max_neighbors", None)
        default_radius    = getattr(backbone, "cutoff", None)

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
        batch = self._collater([data], otf_graph=True).to(self.device)
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
        coord_ang : (N,3) Å
        forces / hessian : return toggles
        Returns dict with keys:
          - energy  (float, eV)
          - forces  (np.ndarray, eV/Å) or None
          - hessian (torch.Tensor, eV/Å^2) or None
        """
        atoms = Atoms(self.elem, positions=coord_ang)
        batch = self._ase_to_batch(atoms)

        need_grad = bool(forces or hessian)
        pos = batch.pos.detach().clone().to(self.device, dtype=_tdtype)
        pos.requires_grad_(need_grad)
        batch.pos = pos

        if need_grad:
            res = self.predict.predict(batch)
        else:
            with torch.no_grad():
                res = self.predict.predict(batch)

        energy = float(res["energy"].squeeze().detach().item())
        forces_np = (
            res["forces"].detach().cpu().numpy().astype(_ndtype, copy=False)
            if (forces or hessian)
            else None
        )

        if hessian:
            p_flags = [p.requires_grad for p in self.predict.model.parameters()]
            for p in self.predict.model.parameters():
                p.requires_grad_(False)

            self.predict.model.train()
            try:
                def e_fn(flat):
                    batch.pos = flat.view(-1, 3)
                    out = self.predict.predict(batch)["energy"].squeeze()
                    return out

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
class uma_pysis(Calculator):
    """PySisyphus-compatible UMA calculator."""

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
    def _au_energy(E: float) -> float:
        return E * EV2AU

    @staticmethod
    def _au_forces(F: np.ndarray) -> np.ndarray:
        return (F * F_EVAA_2_AU).reshape(-1)

    def _au_hessian(self, H: torch.Tensor):
        """
        H: (N,3,N,3) tensor on device.
        - reshape to (3N,3N)
        - symmetrize
        - scale with same-device/dtype scalar to avoid FP64 promotion
        - return numpy or torch depending on out_hess_torch
        """
        n = H.size(0)
        H = H.view(n * 3, n * 3)
        H = H * H_EVAA_2_AU
        return H.detach().cpu().numpy() if not self.out_hess_torch else H.detach()

    # ---------- PySisyphus API --------------------------------------
    def get_energy(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=_ndtype).reshape(-1, 3) * BOHR2ANG
        res = self._core.compute(coord_ang, forces=False, hessian=False)
        return {"energy": self._au_energy(res["energy"])}

    def get_forces(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=_ndtype).reshape(-1, 3) * BOHR2ANG
        res = self._core.compute(coord_ang, forces=True, hessian=False)
        return {
            "energy": self._au_energy(res["energy"]),
            "forces": self._au_forces(res["forces"]),
        }

    def get_hessian(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=_ndtype).reshape(-1, 3) * BOHR2ANG
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
