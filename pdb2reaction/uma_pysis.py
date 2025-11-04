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
