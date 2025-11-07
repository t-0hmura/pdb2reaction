# pdb2reaction/bond_changes.py

from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable, Tuple, Set, List, Optional

import warnings
import torch
import numpy as np

from pysisyphus.elem_data import (
    COVALENT_RADII as CR,
)

Pair = Tuple[int, int]


@dataclass
class BondChangeResult:
    formed_covalent: Set[Pair]
    broken_covalent: Set[Pair]
    distances_1: Optional[np.ndarray] = None
    distances_2: Optional[np.ndarray] = None


def _upper_pairs_from_mask(mask: torch.Tensor) -> Set[Pair]:
    idx = torch.nonzero(mask, as_tuple=False).detach().cpu().numpy()
    return set(map(tuple, idx))


def _element_arrays(atoms: Iterable[str]) -> Tuple[List[str], np.ndarray]:
    elems = [a.capitalize() for a in atoms]
    cov = np.array([CR[a.lower()] for a in elems], dtype=float)
    return elems, cov


def _resolve_device(device: str) -> torch.device:
    dev_str = (device or "cpu").lower()
    if dev_str.startswith("cuda"):
        if torch.cuda.is_available():
            try:
                _ = torch.device(dev_str)
                return torch.device(dev_str)
            except Exception:
                warnings.warn(
                    f"Requested device '{device}' is not available. Falling back to CPU.",
                    RuntimeWarning,
                )
                return torch.device("cpu")
        else:
            warnings.warn(
                "CUDA is not available. Falling back to CPU.",
                RuntimeWarning,
            )
            return torch.device("cpu")
    try:
        return torch.device(dev_str)
    except Exception:
        warnings.warn(
            f"Requested device '{device}' is not recognized. Falling back to CPU.",
            RuntimeWarning,
        )
        return torch.device("cpu")


@torch.no_grad()
def compare_structures(
    geom1,
    geom2,
    device: str = "cuda",
    bond_factor: float = 1.20,
    margin_fraction: float = 0.05,
    delta_fraction: float = 0.05,
) -> BondChangeResult:

    assert geom1.atoms == geom2.atoms, "原子の順序・種類が一致している必要があります。"
    N = len(geom1.atoms)

    elems, cov_np = _element_arrays(geom1.atoms)
    dev = _resolve_device(device)

    dtype = torch.float64
    R1 = torch.as_tensor(geom1.coords3d, dtype=dtype, device=dev)
    R2 = torch.as_tensor(geom2.coords3d, dtype=dtype, device=dev)
    cov = torch.as_tensor(cov_np, dtype=dtype, device=dev)

    T_cov = bond_factor * (cov[:, None] + cov[None, :])
    eps_cov = margin_fraction * T_cov

    D1 = torch.cdist(R1, R1)
    D2 = torch.cdist(R2, R2)

    up = torch.triu(torch.ones((N, N), dtype=torch.bool, device=dev), diagonal=1)

    A1 = (D1 <= (T_cov - eps_cov)) & up
    A2 = (D2 <= (T_cov - eps_cov)) & up

    dD = D2 - D1
    need_change = (dD.abs() >= (delta_fraction * T_cov)) & up

    formed_cov_mask = (~A1) & A2 & need_change
    broken_cov_mask = A1 & (~A2) & need_change

    formed_covalent = _upper_pairs_from_mask(formed_cov_mask)
    broken_covalent = _upper_pairs_from_mask(broken_cov_mask)

    return BondChangeResult(
        formed_covalent=formed_covalent,
        broken_covalent=broken_covalent,
        distances_1=D1.detach().cpu().numpy(),
        distances_2=D2.detach().cpu().numpy(),
    )


def _bond_str(i: int, j: int, elems: List[str], one_based: bool = True) -> str:
    ii = i + 1 if one_based else i
    jj = j + 1 if one_based else j
    return f"{elems[i]}{ii}-{elems[j]}{jj}"


def summarize_changes(geom, result: BondChangeResult, one_based: bool = True) -> str:
    elems = [a.capitalize() for a in geom.atoms]
    lines: List[str] = []

    def pairs_to_lines(title: str, pairs: Set[Pair]):
        if not pairs:
            lines.append(f"{title}: なし")
            return
        lines.append(f"{title} ({len(pairs)}):")
        for i, j in sorted(pairs):
            lines.append(f"  - {_bond_str(i, j, elems, one_based)}")

    pairs_to_lines("Bond formed", result.formed_covalent)
    pairs_to_lines("Bond broken", result.broken_covalent)

    return "\n".join(lines)


# if __name__ == "__main__":
#     from pysisyphus.helpers import geom_loader

#     geom1 = geom_loader("./reac.pdb", coord_type="cart")
#     geom2 = geom_loader("./prod.pdb", coord_type="cart")

#     res = compare_structures(
#         geom1, geom2,
#         device="cuda",
#         bond_factor=1.20,
#         margin_fraction=0.05,
#         delta_fraction=0.05,
#     )

#     print(summarize_changes(geom1, res, one_based=True))
#     pass
