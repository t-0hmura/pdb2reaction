# pdb2reaction/utils.py
import sys
import math
from pathlib import Path

from ase.io import read, write

# ---------------------------------------------------------------------
# convert_xyz_to_pdb                                                         
# ---------------------------------------------------------------------
def convert_xyz_to_pdb(xyz_path: Path, ref_pdb_path: Path, out_pdb_path: Path) -> None:
    """Overlay coordinates from *xyz_path* onto *ref_pdb_path* topology and write *out_pdb_path*.

    *xyz_path* may contain one or many frames.  For multi‑frame trajectories a
    MODEL/ENDMDL block is appended for each subsequent frame.
    """
    ref_atoms = read(ref_pdb_path)  # Reference topology (single frame)
    traj = read(xyz_path, index=":", format="xyz")  # All XYZ frames
    if not traj:
        raise ValueError(f"No frames found in {xyz_path}.")

    for step, frame in enumerate(traj):
        atoms = ref_atoms.copy()
        atoms.set_positions(frame.get_positions())
        if step == 0:
            write(out_pdb_path, atoms)  # overwrite / create
        else:
            write(out_pdb_path, atoms, append=True)

# ---------------------------------------------------------------------
# freeze_links                                                            
# ---------------------------------------------------------------------
def parse_pdb_coords(pdb_path):
    with open(pdb_path, "r") as f:
        lines = f.readlines()

    others = []
    lkhs = []
    for line in lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        name    = line[12:16].strip()
        resname = line[17:20].strip()
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue

        if resname == "LKH" and name == "HL":
            lkhs.append((x, y, z, line))
        else:
            others.append((x, y, z, line))
    return others, lkhs

def nearest_index(point, pool):
    x, y, z = point
    best_i = -1
    best_d2 = float("inf")
    for i, (a, b, c, _) in enumerate(pool):
        d2 = (a - x) ** 2 + (b - y) ** 2 + (c - z) ** 2
        if d2 < best_d2:
            best_d2 = d2
            best_i = i
    return best_i, math.sqrt(best_d2)

def freeze_links(pdb_path):
    """
    Listup indices of link-parent atoms.
    """
    others, lkhs = parse_pdb_coords(pdb_path)

    if not lkhs:
        return []

    indices = []
    for (x, y, z, line) in lkhs:
        idx, dist = nearest_index((x, y, z), others)
        indices.append(idx)
    return indices

