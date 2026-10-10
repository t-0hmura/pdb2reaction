# Small utilities

Six commands that compute, clean, or plot one thing without running a
workflow. Each section gives a minimal run and what to check.

## Which utility

| Command | When to use | Output |
|---|---|---|
| `sp` | One energy and force evaluation, optionally a Hessian, to check a backend, a precision, or a stationary point without an optimizer | Energy on stdout, `forces.npy`; `hessian.npy` with `--hess` |
| `fix-altloc` | A PDB has alternate locations and a cleaned PDB file is itself needed | PDB with one altLoc per residue |
| `add-elem-info` | The PDB element column (77–78) is blank or wrong, as in many tleap PDBs | PDB with the element column filled |
| `bond-summary` | Check formed and broken bonds between R and P, or along a series of frames | Text on stdout; JSON with `--json` |
| `trj2fig` | Plot the energy profile of an XYZ trajectory from IRC, MEP, or scan | PNG, SVG, PDF, JPEG, HTML, or CSV |
| `energy-diagram` | Draw a diagram from energy values you already have, such as R, TS, and P from separate runs | `energy_diagram.png` |

## sp

### Minimal run

```bash
# Energy and forces
pdb2reaction sp -i geom.pdb -l 'SAM:1,GPP:-3' -b orb -o result_sp

# Analytical ORB Hessian at its default fp64 precision
pdb2reaction sp -i geom.xyz -q -2 -m 1 -b orb \
    --hess --hessian-calc-mode Analytical --out-json -o result_sp_hess

# Finite-difference Hessian as a cross-check
pdb2reaction sp -i geom.xyz -q -2 -m 1 -b mace \
    --hess --hessian-calc-mode FiniteDifference -o result_sp_fd
```

`-b dft` uses the PySCF/GPU4PySCF calculator instead of an MLIP. The Hessian
is finite-difference by default for every backend; UMA, ORB, MACE, and
AIMNet2 also take an explicit `--hessian-calc-mode Analytical`. Without
`--precision`, UMA and AIMNet2 run in fp32 and ORB and MACE in fp64.

### Check the result

The console prints `[sp] energy = … a.u.  |force|_max = … a.u./bohr`.
`forces.npy` holds the `(N, 3)` forces in Hartree/Bohr. With `--hess`,
`hessian.npy` holds the Cartesian Hessian without mass weighting
(Hartree/Bohr²): `(3N, 3N)`, or the block of the movable atoms when atoms are
frozen with `--freeze-atoms` or YAML `geom.freeze_atoms`. `sp` has no
`--freeze-links`, so it does not freeze the cap parents. With `--out-json`,
`result.json` and `summary.json` hold `execution_status`,
`scientific_status`, `energy_au`, `forces_path`, and `hessian_path` (null
without `--hess`); `sp` writes no `summary.log`.

- Compare Hessians only at the same geometry, atom order, frozen set, backend
  model, precision, charge, and multiplicity.
- AIMNet2 rejects fp64 instead of changing the precision silently.
- Without `--hess`, the [UMA-worker limit on Analytical Hessians](SKILL.md#cross-cutting-pitfalls)
  does not apply.
- `--show-config` prints the settings and continues; use `--dry-run` to check
  without computing.

## fix-altloc

`fix-altloc` blanks the altLoc column (column 17) without shifting any other
field and keeps one altLoc label per residue: the one with the highest mean
occupancy over that residue's labeled atoms, and on a tie the label that
appears first. A is not preferred. Atoms with a blank altLoc are kept, and
atoms that belong only to an unselected conformer are dropped, not merged into
a hybrid residue. `extract` and the geometry commands apply the same selection
on their own, so run `fix-altloc` only when a cleaned PDB file is itself
needed, to clean a directory, or to inspect the chosen conformer first.

### Minimal run

```bash
pdb2reaction fix-altloc -i raw.pdb -o cleaned.pdb
pdb2reaction fix-altloc -i raw_pdbs/ -o cleaned_pdbs/
pdb2reaction fix-altloc -i raw.pdb --inplace
```

`-i` is a file or a directory, and `-o` is then a file or a directory.
Without `-o`, the output is `<input>_clean.pdb` or `<input>_clean/`.
`--recursive` also cleans subdirectories, and `--inplace` overwrites each
input after saving `<name>.pdb.bak`.

### Check the result

A file prints `[fix-altloc] Fixed altLoc → …`; a directory prints
`[fix-altloc] Processed N file(s) → …` and `Skipped N file(s)` for files
without altLoc. Directory output keeps the input's relative paths. Column 17
is blank for every kept atom.

### Pitfalls

- The choice is an occupancy heuristic. When the active-site conformer must
  be chosen from ligand contacts or the mechanism, select it in a structure
  editor and check it.
- A file without altLoc is skipped and nothing is written; `--force`
  processes it anyway.
- An existing output stops the run with `Output exists: … (use --overwrite
  …)`. `--inplace` cannot be combined with `-o`, and an existing `.bak` is not
  replaced.
- Atom serial numbers are not renumbered, `CONECT` records are not updated,
  and each `MODEL` block is processed on its own.
- Run `fix-altloc` only when alternate locations are present, and
  `add-elem-info` only when element fields are missing or wrong; neither is a
  required step for every RCSB file.

## add-elem-info

`add-elem-info` fills or corrects the element column (77–78) of ATOM and
HETATM records, which tleap often leaves blank; `extract` stops without it.
`all` runs this step itself only when element fields are missing.

### Minimal run

```bash
pdb2reaction add-elem-info -i raw.pdb -o cleaned.pdb
```

Without `-o`, the output is `<input>_add_elem.pdb`; `--overwrite` without
`-o` replaces the input, and an `-o` equal to the input also needs
`--overwrite`. A valid existing symbol (or `EP` for a water virtual site) is
kept; `--overwrite-elem` infers every field again. Remove virtual sites before
`extract` ([Common edits](../pdb2reaction-model-setup/formats.md#common-edits)).

The element comes from the atom and residue names:

- Ion residues take the element from the residue name; `NH4` and `H3O+` are
  assigned per atom (H, N, O).
- Protein, nucleic acid, and water atoms use the usual PDB elements (H, C, N,
  O, S, P, Se).
- Other ligands use the column alignment of the atom name: ` NA ` is N and
  `NA  ` is Na. LEaP ` CL1` and ` BR1` are halogens, and `HG11` is H.

### Check the result

The console prints `[OK] Wrote: …` with the counts `total atoms`,
`assigned/updated`, and `kept existing`. Atoms that cannot be assigned are
left unchanged and listed after `[WARN] Could not confidently assign N atoms;
left unchanged.`, up to 50 of them; type their symbols into columns 77–78 by
hand. With no `[WARN]` line, every atom has an element.

### Pitfalls

- Only columns 77–78 of the repaired records change; diff the files to check
  the changes.
- Atom names that do not follow the usual conventions, as in some ligands, can
  be assigned the wrong element; spot-check them.

## bond-summary

`bond-summary` reports the bonds formed and broken between each consecutive
pair of structures, with the same criterion that `path-search` uses to split a
path and that `irc` and `all` use for `bond_changes`. Use it to check R
against P, or to see where a recursive `path-search` split the path.

### Minimal run

```bash
pdb2reaction bond-summary -i 1.R.pdb -i 3.P.pdb

# Several frames; files may also follow one -i or be given without -i
pdb2reaction bond-summary -i frame_01.xyz -i frame_05.xyz -i frame_10.xyz
```

`--bond-factor` (default 1.20) scales the covalent radii, `--device` defaults
to `cpu`, and atom labels are 1-based unless `--zero-based`.

### Check the result

Text output on stdout, one block per pair:

```text
============================================================
  1.R.pdb  →  3.P.pdb
============================================================
Bond formed (2):
  - C320-C321 : 3.170 Å --> 1.680 Å
  - O186-H321 : 2.230 Å --> 0.980 Å
Bond broken (2):
  - S320-C320 : 1.800 Å --> 3.430 Å
  - C321-H321 : 1.100 Å --> 2.320 Å
```

With `--json`, stdout holds one JSON object instead, and no file is written:

```json
{
  "execution_status": "completed",
  "scientific_status": "success",
  "comparisons": [
    {"structure_a": "1.R.pdb", "structure_b": "3.P.pdb",
     "bonds_formed": 2, "bonds_broken": 2}
  ]
}
```

A pair that cannot be compared prints `ERROR: …` to stderr and is left out of
`comparisons`; the run then exits with 1, `execution_status` is `failed`, and
`scientific_status` is `partial` when other pairs were compared or `failed`
when none were.

### Pitfalls

- All inputs need the same atoms in the same order, or the pair fails with
  `ERROR: Atom types and ordering must be identical.` `extract` checks this
  for several inputs but does not map atoms; reorder the full structures
  first, then extract them together.
- The cutoff is geometric only and does not tell covalent, ionic, and
  hydrogen bonds apart; metal–ligand contacts can sit near it. A larger
  `--bond-factor` (for example 1.30 for metal coordination) counts longer
  contacts but can add false bonds; report the value and inspect the pairs.
- The JSON gives counts only; the atom pairs are in the text output and in
  `bond_changes` of `summary.json`
  ([outputs.md](../pdb2reaction-overview/outputs.md#bond-changes)).

## trj2fig

`trj2fig` plots the energy along an XYZ trajectory, ΔE in kcal/mol relative to
the first frame by default.

### Minimal run

```bash
pdb2reaction trj2fig -i finished_irc_trj.xyz -o irc_profile.png

# Interactive HTML; the extension selects the format
pdb2reaction trj2fig -i mep_trj.xyz -o mep.html

# Several formats and result.json
pdb2reaction trj2fig -i mep_trj.xyz -o mep.png -o mep.html -o mep.csv --out-json
```

The formats are `.png`, `.jpg`, `.jpeg`, `.svg`, `.pdf`, `.html`, and `.csv`;
without `-o` the output is `energy.png`. `--unit hartree` changes the unit,
`-r` sets the reference (`init`, the left-end frame; a 0-based index; or
`none` for absolute energies), and `--reverse-x` flips the x-axis.

### Check the result

The console prints `[trj2fig] Saved figure -> …` for each figure. The CSV has
`frame`, `energy_hartree`, and the plotted value (`delta_kcal`,
`delta_hartree`, or with `-r none` `energy_kcal`). With `--out-json`,
`result.json` in the directory of the first output holds `execution_status`,
`scientific_status`, `n_frames`, the minimum and maximum energies in hartree,
`energy_source`, and `files`. `energy_source` is `trajectory_comment` when the
energies come from the comment lines, and then the calculator fields are null;
it is `mlip_recomputed` with `-q` or `-m`.

### Pitfalls

- Each comment line must give the energy. Trajectories written by
  pdb2reaction are read as they are. Otherwise write `E=<value>` (or
  `Energy:`) with an optional unit (`Ha`, `Eh`, `hartree`, `eV`, `kcal/mol`;
  hartree without one); a comment that is only one decimal number, such as
  `-12345.67`, is read as hartree. A comment with no number, only an integer,
  or several bare numbers stops the run with an error that names the frame.
- `-q` or `-m` makes the MLIP chosen with `-b` recompute every frame instead
  of reading the comments; with only one of them, the other is taken as charge 0 or
  multiplicity 1, so give both for a charged or open-shell system.
- For a labeled diagram of R, TS, IM, and P, use
  [`energy-diagram`](#energy-diagram).

## energy-diagram

`energy-diagram` draws an energy diagram from numbers, for example values
collected from several runs into one figure.

### Minimal run

```bash
# Two-step profile, one value per -i
pdb2reaction energy-diagram -i 0.0 -i 21.5 -i -0.7 -i 2.2 -i -18.2 -o diagram.png

# The same values as one list literal
pdb2reaction energy-diagram -i "[0.0, 21.5, -0.7, 2.2, -18.2]" -o diagram.png

# State labels, a y-axis label, and result.json
pdb2reaction energy-diagram -i "[0.0, 21.5, -18.2]" \
    --label-x R --label-x TS --label-x P \
    --label-y 'ΔG (kcal/mol)' --out-json -o diagram.png
```

The states are plotted in input order and named `S1`, `S2`, … without
`--label-x`. The y-axis label defaults to `ΔE (kcal/mol)`.

### Check the result

The console prints `[energy-diagram] Saved -> …`. With `--out-json`,
`result.json` and `summary.json` next to the image hold `execution_status`,
`scientific_status`, `n_points`, and `files`; they do not record the values or
the labels.

### Pitfalls

- The values are plotted as given, with no unit conversion; put them all in
  one unit (usually kcal/mol) first.
- `-i 0 12.5 4.3` and `--label-x R TS P` are rejected, because each option
  takes one value; repeat the option or quote one list. The number of labels
  must equal the number of values, and at least two values are needed.
- For a profile along a trajectory, use [`trj2fig`](#trj2fig).

## Next step

- Frequencies and thermochemistry from a Hessian: [freq.md](freq.md);
  saddle-point checks: [tsopt.md](tsopt.md).
- Backend installation and precision:
  [pdb2reaction-install](../pdb2reaction-install/SKILL.md).
- PDB columns and other cleanup:
  [formats.md](../pdb2reaction-model-setup/formats.md#pdb); cutting the
  cluster: [extract.md](extract.md).
- Trajectories to plot come from [irc.md](irc.md), [path.md](path.md), and
  [scan.md](scan.md).
- Per-segment energies in `summary.json` to feed `energy-diagram`:
  [outputs.md](../pdb2reaction-overview/outputs.md#energy-diagrams).
