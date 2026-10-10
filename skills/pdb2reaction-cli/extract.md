# `pdb2reaction extract`

## When to use

Cut an active-site cluster model out of a protein–ligand PDB/mmCIF around the
substrate and catalytic residues, before running `opt`, `path-opt`,
`path-search`, `scan`, or `tsopt` on the cluster. Cut bonds at CA and CB are
capped with hydrogens (cap H), and the residue charges are summed (`-l` gives
the ligand charges). `all` with `-c` runs this step for you.

## Minimal run

```bash
pdb2reaction extract -i 1abc.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' -r 4.0 --out-json -o cluster.pdb
```

`-c` lists the centers; every match starts the radius expansion. It takes
residue names (`'GPP,SAM'`), chain:number IDs (`'A:44,B:321'`), chain-qualified
names (`'B:SAM'`), chain:name:number (`'B:SAM:321'`), or a PDB/mmCIF file of
the substrate. In a PDB with chains, prefer the chain:name:number form.

Without `-o`, the model is `model.pdb` (one input) or `model_<input name>.pdb`
(several inputs). For mmCIF input, or a PDB too large for the PDB columns,
`extract` also writes a file with the same stem as `.cif` that keeps the
original chain IDs and residue numbers.

Other inputs:

```bash
# Reactant and product cut with the same boundary
pdb2reaction extract -i 1.R.pdb -i 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -o 1.R_clu.pdb -o 3.P_clu.pdb

# Substrate given as a structure file
pdb2reaction extract -i complex.pdb -c substrate.pdb -r 3.5 -o cluster.pdb

# Repeated ligand name in a large mmCIF: writes cluster.pdb and cluster.cif
pdb2reaction extract -i complex.cif -c 'LONG_CHAIN:SAM:10001' -l 'SAM:1' -o cluster.pdb
```

## Judge success

The run exits with 0 and prints `[extract] Atoms after truncation: N` and
`[extract] Link-H to add: M`; the model has N + M atoms, M of them cap H. With
several inputs, the line is `[extract:multi] Atoms after truncation (model k): N`
for each model. With
`--out-json`, `result.json` (and its copy `summary.json`) is written next to
the first output file:

```python
import json
d = json.load(open("result.json"))
print(d["execution_status"], d["scientific_status"], d["total_charge"])
print(d["n_atoms_raw"], "->", d["n_atoms_extracted"], "+", d["n_link_hydrogens"], "cap H")
print(d["files"])              # {basename: full_path} for each written PDB/CIF
print(d["protein_charge"], d["ligand_total_charge"], d["ion_total_charge"])
```

`n_atoms_extracted` does not count the cap H. `extract` freezes nothing
itself; the downstream commands freeze the cap parents (below) and report the
count as `n_freeze_atoms`.

## Pitfalls and recovery

- When `-c` is a structure file, its atom names must match the complex exactly
  (case-sensitive).
- A blank element column stops the run with `Element symbols are missing in
  '…'`; run [`add-elem-info`](utilities.md#add-elem-info) first.
- `extract` keeps one coherent altLoc per residue on its own; run
  [`fix-altloc`](utilities.md#fix-altloc) only when a cleaned PDB file is
  itself needed.
- Several inputs must have the same atoms in the same order, or the run stops
  with `[multi] Atom count mismatch` or `[multi] Atom order mismatch`.
- Unknown ligand charges come only from `-l`; the internal tables cover
  standard and recognized modified amino acids and ions, not cofactors in
  general. A mapping name that matches no unknown selected residue gives a
  warning and is ignored. Ligand and ion details:
  [formats.md](../pdb2reaction-model-setup/formats.md#ligand-ion-and-metal-charges).
- A residue with non-standard amino-acid names (for example from MCPB.py) is
  neither cut nor capped, and a WARNING says that backbone truncation was not
  applied; register it with `--modified-residue 'HD1,HE1'` (a bare name counts
  as charge 0; `NAME:charge` sets it).
- `extract` takes no `-q`, `-m`, `--freeze-links`, `--config`,
  `--convert-files`, `--show-config`, or `--dry-run`; give them to the
  downstream command.
- No radius is safe for every system. How to choose and check the radius and
  the boundary: [model-setup](../pdb2reaction-model-setup/SKILL.md#check-the-boundary-and-the-charge).

## Freeze atoms at the cluster boundary

A cluster cut from a protein needs its boundary held in place, or the
optimizer pulls the dangling fragment into an unphysical geometry. `extract`
places each cap H 1.09 Å from the cut carbon along the old bond and writes it
as `HETATM` residue `LKH`, atom `HL`. The carbon that carries the cap (its
parent) is the atom to freeze. Other cut bonds get no cap. The record format
is in [formats.md](../pdb2reaction-model-setup/formats.md#cap-hydrogens).

### Three sources of frozen atoms

- `--freeze-links` (on by default) freezes the cap parents. It reads `LKH`
  from PDB/mmCIF input, or from XYZ/GJF input given with `--ref-pdb`.
- `--freeze-atoms 'i,j,k'` takes 1-based atom numbers in any input format.
- YAML `geom.freeze_atoms` (with `--config`) takes 1-based numbers in any
  input format.

A run freezes the union of the three; no source overrides another.

### Recipes

PDB from `extract` (`--freeze-links` is on by default):

```bash
pdb2reaction extract -i complex.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' -o model.pdb
pdb2reaction opt -i model.pdb -l 'SAM:1,GPP:-3' -m 1
pdb2reaction tsopt -i ts_model.pdb -l 'SAM:1,GPP:-3' -m 1 -o result_tsopt
pdb2reaction freq -i result_tsopt/final_geometry.pdb -l 'SAM:1,GPP:-3' -m 1 -o result_freq
```

XYZ or GJF input has no `LKH` residue, so give the atoms or a reference PDB:

```bash
pdb2reaction tsopt -i ts_candidate.xyz -q 0 -m 1 --freeze-atoms '12,15,28,29,42'
pdb2reaction tsopt -i ts_candidate.xyz -q 0 -m 1 --ref-pdb model.pdb
```

Keep a long list in YAML and add one-off atoms with `--freeze-atoms`:

```yaml
# tsopt.yaml
geom:
  freeze_atoms: [12, 15, 28, 29, 42, 88, 91, 92]
```

```bash
pdb2reaction tsopt -i ts.xyz -q 0 -m 1 --config tsopt.yaml
```

### Effect on the calculation

- `opt`, `tsopt`, `scan`, `freq`, `irc`, and GSM `path-opt`/`path-search` set
  the forces on frozen atoms to zero, so they keep their initial Cartesian
  coordinates.
- DMF (`--mep-mode dmf`) holds them with a stiff harmonic restraint
  (k = 300 eV/Å²) instead, so they can move slightly; inspect them.
- Frozen atoms are left out of the Hessian, and `freq` runs a partial Hessian
  vibrational analysis (PHVA) on the movable atoms.

### PHVA treatment

With frozen atoms, `freq`, `irc`, the TS check and Dimer direction in
`tsopt`, and `--flatten` use the Hessian of the movable atoms, and only rigid
motions that keep every frozen atom in place are removed. For a cluster with
three or more frozen atoms off one line, that removes nothing. This is
separate from `tsopt --ref-mode`. See the
[JSON record](../../docs/json-output.md#rigid-projection-provenance) and
[details](../../docs/freq.md#rigid-modes-with-frozen-boundaries).

### Which subcommands use it

`opt`, `tsopt`, `freq`, `irc`, `path-opt`, `path-search`, `scan`, `scan2d`,
`scan3d`, and `all` take all three sources. `sp` takes only `--freeze-atoms`
and YAML; with `--hess` it writes the Hessian of the movable atoms. `extract`
only writes the cap records.

### Pitfalls

- `LKH`/`HL` records deleted by hand leave `--freeze-links` nothing to find.
  Rerun `extract` or give `--freeze-atoms`.
- `--freeze-atoms` follows the input atom order. After a new extraction with a
  different selection, rebuild the list.
- XYZ/GJF input without `--ref-pdb` makes `--freeze-links` do nothing.
- `--no-freeze-links` is for diagnostic runs that let the boundary relax on
  purpose; keep `--freeze-links` on in production runs.
- If every atom is frozen, `freq` stops with an error; keep at least one atom
  movable.

## Next step

- Build, trim, or enlarge the model, and check its boundary and charge:
  [model-setup](../pdb2reaction-model-setup/SKILL.md).
- Residue selectors and PDB columns:
  [model-setup](../pdb2reaction-model-setup/SKILL.md#selecting-residues-and-atoms).
- Pre-clean a raw PDB: [`add-elem-info`](utilities.md#add-elem-info),
  [`fix-altloc`](utilities.md#fix-altloc).
- Full guide to freezing: [docs](../../docs/model-setup.md#freeze-atoms-and-restrain-distances).
- Flags and defaults: `--help-advanced` and [SKILL.md](SKILL.md#where-flags-and-defaults-live).
