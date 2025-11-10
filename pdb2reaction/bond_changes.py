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
from pysisyphus.constants import BOHR2ANG  # ← 追加：長さをÅ表示するため

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

    assert geom1.atoms == geom2.atoms, "Atom types and ordering must be identical."
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
    """
    結合の形成/解離を列挙。さらに、対応する結合長の変化を
    'x.xxx Å --> y.yyy Å' 形式で出力する（Å表示）。
    """
    elems = [a.capitalize() for a in geom.atoms]
    lines: List[str] = []

    # 距離行列（Bohr）→ Å に変換して使う（利用可能な場合のみ）
    D1 = result.distances_1
    D2 = result.distances_2
    have_lengths = (
        isinstance(D1, np.ndarray)
        and isinstance(D2, np.ndarray)
        and D1.shape == D2.shape
    )

    def _len_str(i: int, j: int) -> str:
        if not have_lengths:
            return ""
        # coords3d は Bohr 前提。Å に変換。
        d1 = float(D1[i, j]) * BOHR2ANG
        d2 = float(D2[i, j]) * BOHR2ANG
        return f" : {d1:.3f} Å --> {d2:.3f} Å"

    def pairs_to_lines(title: str, pairs: Set[Pair]):
        if not pairs:
            lines.append(f"{title}: None")
            return
        lines.append(f"{title} ({len(pairs)}):")
        for i, j in sorted(pairs):
            lines.append(f"  - {_bond_str(i, j, elems, one_based)}{_len_str(i, j)}")

    pairs_to_lines("Bond formed", result.formed_covalent)
    pairs_to_lines("Bond broken", result.broken_covalent)

    return "\n".join(lines)
