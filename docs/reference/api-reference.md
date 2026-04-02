# Python API Reference

This page documents the most useful public Python APIs for programmatic use. These functions can be imported and used without the CLI.

---

## Bond Change Detection

```python
from pdb2reaction.bond_changes import compare_structures, has_bond_change, summarize_changes
```

| Function | Description |
|----------|-------------|
| `compare_structures(geom1, geom2, device="cuda", bond_factor=1.20)` | Detect covalent bonds formed or broken between two pysisyphus geometries. Returns a `BondChangeResult` with `formed_covalent` and `broken_covalent` (sets of 0-based index pairs). |
| `has_bond_change(geom_start, geom_end, bond_cfg)` | Convenience wrapper: returns `(has_changes: bool, summary_text: str)`. `bond_cfg` accepts keys: `device`, `bond_factor`, `margin_fraction`, `delta_fraction`. |
| `summarize_changes(geom, result, one_based=True)` | Build a human-readable report of bond changes with distances in Angstrom. |

---

## Active Site Model Extraction

```python
from pdb2reaction.extract import extract_api, compute_charge_summary, load_structure
```

| Function | Description |
|----------|-------------|
| `extract_api(complex_pdb, center, output=None, radius=2.6, ...)` | Extract active site model from PDB. Accepts PDB path(s), substrate center spec, and extraction options. Returns `{"outputs": [...], "counts": [...], "charge_summary": {...}}`. |
| `compute_charge_summary(structure, selected_ids, substrate_ids, ligand_charge=None)` | Compute total model charge from amino acid, ion, and ligand contributions. Returns dict with `total_charge`, `protein_charge`, etc. |
| `load_structure(path, name)` | Load a PDB file into a Biopython `Structure` object. |

---

## MLIP Calculator Factory

```python
from pdb2reaction.backends import create_calculator, create_ase_calculator
```

| Function | Description |
|----------|-------------|
| `create_calculator(backend="uma", **kwargs)` | Create a PySisyphus-compatible MLIP calculator. Backends: `uma`, `orb`, `mace`, `aimnet2`. Supports `solvent`, `solvent_model` kwargs. |
| `create_ase_calculator(backend="uma", **kwargs)` | Create an ASE-compatible MLIP calculator (used for DMF workflows). |

---

## Energy Diagram & Trajectory Utilities

```python
from pdb2reaction.utils import build_energy_diagram, read_xyz_energies, convert_xyz_to_pdb
```

| Function | Description |
|----------|-------------|
| `build_energy_diagram(energies, labels, ylabel="ΔE")` | Create a Plotly energy diagram figure with horizontal state segments and dotted connectors. |
| `read_xyz_energies(path)` | Extract energies (hartree) from the comment line of each frame in an XYZ trajectory. |
| `convert_xyz_to_pdb(xyz_path, ref_pdb_path, out_pdb_path)` | Overlay XYZ coordinates onto a reference PDB topology template. Supports multi-frame trajectories. |
| `read_xyz_first_last(trj_path)` | Read first and last frames from an XYZ trajectory. Returns `(elements, first_coords_A, last_coords_A)`. |

---

## Structure Alignment

```python
from pdb2reaction.align_freeze_atoms import kabsch_R_t, align_and_refine_pair_inplace
```

| Function | Description |
|----------|-------------|
| `kabsch_R_t(P, Q)` | Kabsch algorithm: find rotation `R` and translation `t` minimizing RMSD between (N,3) point sets. |
| `align_and_refine_pair_inplace(g_ref, g_mob, ...)` | Rigid Kabsch alignment + iterative scan-and-relax for a pair of geometries. |

---

## PDB Utilities

```python
from pdb2reaction.add_elem_info import guess_element, assign_elements
from pdb2reaction.fix_altloc import has_altloc, fix_altloc_file
```

| Function | Description |
|----------|-------------|
| `guess_element(atom_name, resname, is_het)` | Infer element symbol from PDB atom name + residue name. |
| `assign_elements(in_pdb, out_pdb)` | Populate element columns (77-78) for all ATOM/HETATM records. |
| `has_altloc(pdb_path)` | Check if a PDB file contains alternate location characters. |
| `fix_altloc_file(in_path, out_path)` | Resolve alternate locations by keeping highest-occupancy conformer. |

---

## Charge & Spin Resolution

```python
from pdb2reaction.utils import resolve_charge_spin, detect_freeze_links
```

| Function | Description |
|----------|-------------|
| `resolve_charge_spin(prepared_inputs, charge, spin, ...)` | Resolve charge/spin from CLI args, ligand metadata, and GJF templates. Returns `(charge, spin)`. |
| `detect_freeze_links(pdb_path)` | Identify 0-based atom indices of link-parent atoms for LKH/HL link hydrogens. |

---

## Configuration Defaults

```python
from pdb2reaction.defaults import CALC_KW_DEFAULT, LBFGS_KW, RFO_KW, IRC_KW, FREQ_KW
```

All default configuration dictionaries can be imported and overridden for programmatic use. See [YAML Reference](../yaml-reference.md) for the full list of keys.

---

## See Also

- [YAML Reference](../yaml-reference.md) — Complete YAML configuration options
- [CLI Command Reference](commands/index.md) — CLI subcommand help
- [MLIP Calculator](../uma-pysis.md) — Backend configuration details
