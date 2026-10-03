# Building the cluster model

`pdb2reaction` computes on a cluster model: the residues around the substrate, cut out of the protein, with each cut bond capped by a hydrogen (a cap hydrogen).

## Quick guide

| Goal | What to do | Section |
| --- | --- | --- |
| Build the default model | Give the substrate, cofactors, and metals to `-c` | [Build the default model](#build-the-default-model) |
| See the atom count and charge before a long run | Run `extract` on its own, or `all --dry-run` | [Check the model](#check-the-model) |
| Make the calculation lighter | `--exclude-backbone`, `--no-include-h2o`, `-r 0 --selected-resn`, or trim the model in PyMOL | [Make the model smaller](#make-the-model-smaller) |
| A residue of the reaction is missing, or only a fragment of it is in | Raise `-r`, or add the residue to `-c` | [Make the model larger](#make-the-model-larger) |
| Use a model you built yourself | Omit `-c`, and freeze the boundary atoms with `--freeze-atoms` | [Use a model you built yourself](#use-a-model-you-built-yourself) |
| Decide or check the boundary by hand | Follow the checklist | [Building or auditing a cluster model manually](#building-or-auditing-a-cluster-model-manually) |
| Freeze atoms or restrain a distance | `--freeze-links` (on by default), `--freeze-atoms`, `--distance-restraint` | [Freeze atoms and restrain distances](#freeze-atoms-and-restrain-distances) |

## Build the default model

Give the substrate to `-c`, and put the cofactors and metals in the same list. Each residue in the list starts its own distance search.

```bash
pdb2reaction extract -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -o model.pdb --out-json
```

This is the bundled example in [`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples). The model is written to `model.pdb`, and `result.json` next to it records the atom count without the cap hydrogens (`n_atoms_extracted`), the number of cap hydrogens (`n_link_hydrogens`), and the charge (`total_charge`). The model has `n_atoms_extracted + n_link_hydrogens` atoms.

The default model is built as follows:

- **Radius**: a residue joins the model when any of its atoms lies within `-r` (default 2.6 Å) of an atom in `-c`; when that atom is a main-chain atom, the residues on both sides join too. Waters are included (`--include-h2o`).
- **Cuts**: a run of consecutive amino acids keeps its main chain and is cut at CA at both ends; a residue whose neighbors are not in the model is cut at CB and keeps only its side chain. Amino acids in `-c` keep all their atoms.
- **Cap hydrogens**: each bond cut at CA or CB gets a hydrogen 1.09 Å from the carbon (atom `HL` in residue `LKH`).
- **Several structures**: give the reactant and product in one run, and the union of the residues selected in each structure is used for all of them, so every model has the same boundary. The inputs must list the same atoms in the same order. In the bundled example, `1.R.pdb` alone gives 632 atoms with charge +1, and `1.R.pdb` with `3.P.pdb` gives 669 atoms with charge 0. In `3.P.pdb` a main-chain hydrogen of Cys167 lies within the radius, so Cys167, Val166, and Asp168 (−1) join every model of the run.

## Check the model

Check the atom count and the charge before a long run. Run `extract` on its own as above, or add `--dry-run` to `all`:

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --dry-run
```

With `-c`, `--dry-run` extracts the model in a temporary directory and checks the charge and whether the electron count is even or odd; it runs no calculation and writes no output. Read these lines on the console:

- **Atom count**: `Atoms after truncation: N` counts the atoms before the cap hydrogens. Add the number of cap hydrogens from `[extract] Link-H to add: M` (one input) or `[extract:multi] link-H targets common across models: M` (several inputs). For `1.R.pdb` alone, N + M = 608 + 24 = 632.
- **Charge**: `Total active site model charge`. With `--dry-run`, `[all] --dry-run extract: model total_charge = …` also appears.
- **Warnings**: `extract` warns when a bond other than C–C crosses the boundary, and when a residue looks like an amino acid but has an {ref}`unknown name <extract-modified-residue>`. Look at the boundary, the caps, and the charge.

After a real `all` run with `-c`, the model is kept as `_work/models/model_<input name>.pdb` in the output directory. Open it in a viewer and check that the residues of the reaction are in the model.

## Make the model smaller

DFT optimization is practical up to roughly 300 atoms. For `1.R.pdb` alone with the `-c` and `-l` above, the options below give these models; the atom counts include the cap hydrogens.

| How to build the model | Options | Atoms | Charge |
| --- | --- | --- | --- |
| Default extraction | — | 632 | +1 |
| Drop main-chain atoms | `--exclude-backbone` | 397 | +1 |
| …and waters | `--exclude-backbone --no-include-h2o` | 367 | +1 |
| Pick residues yourself | `-r 0 --selected-resn '44,63,186'` | 129 | −1 |
| Edit by hand | Open an extracted PDB in PyMOL, delete what you do not need, and save | — | — |

- **`--exclude-backbone`** removes the main-chain atoms (N, CA, C, O, and their hydrogens) from amino acids outside `-c`, so each side chain is cut at CB; prolines keep their ring. A residue whose only atoms near the center are main-chain atoms is not selected.
- **`--no-include-h2o`** leaves out the waters.
- **`-r 0 --selected-resn`** adds no residues by distance: the model holds the `-c` residues and the residues you list. In the bundled example, residues 44, 63, and 186 are the three closest to the methyl carbon of SAM (CS1); for your own system, pick the residues that take part in the reaction.
- **By hand**: delete residues in PyMOL or another viewer, save the file, and use it as in [Use a model you built yourself](#use-a-model-you-built-yourself).

`all` takes the same options, so the MLIP search can run on the small model from the start; see example 1 of [Refine an MLIP TS with DFT](dft-backend.md#examples). Removing residues changes the charge, so check the log each time.

## Make the model larger

Make the model larger when a residue of the reaction is missing, or when only a fragment of it is in the model.

- **Raise `-r`** (default 2.6 Å), for example to 3.5 Å or 4.0 Å.
- **`--radius-het2het`** (default 0, off) adds a second cutoff measured only between atoms other than C and H, on both the center side and the neighbor side. It picks up close N and O partners without enlarging the whole radius.
- **Add the residue to `-c`** to keep it whole. Amino acids in `-c` keep all their atoms, main chain included, and start their own distance search. A residue added with `--selected-resn` is not a center, so without a neighbor in the model it keeps only its side chain.
- **Residues of a partner chain** that lie outside the radius: add them to `-c` with their chain, such as `-c 'A:SAM,A:GPP,A:MG,B:MET:38'`.

A larger model costs more, and accuracy does not always improve with size. Check for your system that the barrier does not change when the model grows. To include the whole protein as the environment, use the ML/MM toolkit [mlmm-toolkit](https://github.com/t-0hmura/mlmm_toolkit).

## Use a model you built yourself

Omit `-c`, and `all` skips the extraction and uses the input as it is. Give the total charge with `-q`, or the charge of each residue name with `-l` (PDB/mmCIF input only).

```bash
pdb2reaction all -i cluster_R.pdb cluster_P.pdb -q 0 --tsopt --thermo
```

When you trim an extracted PDB, the cap hydrogens you keep still freeze their parents. Delete the cap hydrogens of the residues you removed; a cap hydrogen without its parent stops the run with `isolated LKH/HL`. Cap each new cut as in the checklist below. A model built without `extract` has no `LKH`/`HL` cap hydrogens unless you add them, so `--freeze-links` freezes nothing; freeze its boundary atoms with `--freeze-atoms`.

### Building or auditing a cluster model manually

Use this list when you decide the boundary yourself, or when you check a model that `extract` built.

- Choose the bonds to cut by their chemistry, not only by distance.
- Where you keep a main-chain fragment, stop both ends at CA and fill the valence with cap hydrogens.
- Cut side chains, ligands, and cofactors at nonpolar C–C single bonds where possible (usually CA–CB). Do not cut peptide C–N bonds, polar C–O or C–N bonds, aromatic or conjugated bonds, S–S bonds, or bonds to a metal, even when they reach past the radius; include the partner or move the boundary.
- Place one cap hydrogen on the kept atom of each cut bond, 1.09 Å along the old bond. Name it atom `HL` in residue `LKH` (chain `L`), and its parent atom is frozen automatically.
- Check that each boundary atom has the intended valence and that each cut bond has exactly one cap hydrogen.
- In the reactant, intermediates, and product, list the same atoms in the same order. Decide the selection and the caps once and apply them to every state; do not cut each state separately.
- After cutting, count the charge and the multiplicity again, and inspect the boundary before any minimum-energy path (MEP) or Hessian calculation.
- Covalently bound cofactors, modified residues, metal sites, and boundaries that do not fit these rules need hand edits after the automatic extraction.

(freeze-atoms-and-restraints)=
## Freeze atoms and restrain distances

Freezing keeps an atom in place by setting its force to zero. A restraint pulls the distance between two atoms toward a target with a harmonic spring. In an extracted model, freezing the parent atoms of the cap hydrogens keeps the boundary from deforming during optimizations, path searches, and IRC (intrinsic reaction coordinate).

### Three ways to freeze atoms

- **`--freeze-links/--no-freeze-links`** (default on) freezes the parent atoms of the cap hydrogens that `extract` added. It works on PDB and mmCIF input. XYZ and GJF input has no `LKH` residue, so pass a PDB with the same atoms with `--ref-pdb` to use it there; the coordinates still come from `-i`. `sp` does not freeze these atoms automatically.
- **`--freeze-atoms 'i,j,k'`** takes 1-based atom numbers and works with any input format. In `all` with `-c`, number the atoms of the original input; they are mapped onto the extracted model. Other commands use the numbering of the structure you pass.
- **YAML `geom.freeze_atoms`** (passed with `--config`) suits a long list, or a list you keep with the rest of the settings.

```bash
pdb2reaction extract -i complex.pdb -c 'A:GPP:301,A:SAM:302' -l 'GPP:-3,SAM:1' -o model.pdb
pdb2reaction opt -i model.pdb -l 'GPP:-3,SAM:1' -m 1   # the parents of the cap hydrogens are frozen
pdb2reaction tsopt -i ts_candidate.xyz -q 0 -m 1 --freeze-atoms '12,15,28,29,42'
```

```yaml
geom:
  freeze_atoms: [12, 15, 28, 29, 42]   # 1-based
```

A run freezes the union of the three lists. All three work in `all`, `opt`, `tsopt`, `freq`, `irc`, `path-opt`, `path-search`, `scan`, `scan2d`, and `scan3d`; `sp` takes only `--freeze-atoms` and YAML.

### What freezing does

- Frozen atoms get zero force in `opt`, `tsopt`, `scan`, `freq`, `irc`, and the GSM (growing string method) path search.
- Frozen atoms are left out of the Hessian, so `freq` runs a partial Hessian vibrational analysis (PHVA) on the movable atoms. For the rigid motions that are removed, see [Rigid modes with frozen boundaries](freq.md#rigid-modes-with-frozen-boundaries).
- The DMF (direct max flux) path search (`--mep-mode dmf`) holds frozen atoms with a harmonic restraint instead (k = 300 eV/Å²), so they can move slightly; see the [`path-opt` notes](path-opt.md#notes).

### Restrain a distance

`--distance-restraint` adds a harmonic restraint between two atoms. Give `(i, j, target)` with the target in Å, or `(i, j)` to keep the starting distance.

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 \
    --distance-restraint '[(1,5,2.0)]' --restraint-k 20.0 --out-dir ./result_opt_rest
```

- Atom numbers are 1-based; `--zero-based` counts from 0. `--restraint-k` sets the force constant (default 300 eV/Å²).
- In `scan`, `scan2d`, and `scan3d`, `--restraint-k` is the spring of each step (eV/Å² for distances, eV/rad² for angles) and takes priority over YAML `bias.k`. Each stage restrains only its own pairs. In `all`, pass the same value with `--scan-restraint-k`. For the format of the scan lists, see {ref}`Scan-list spec <scan-list-spec>`.
- Reported energies leave out the energy of the restraint.

## Notes

- **300 atoms** is a rule of thumb for DFT, not a limit in the code. GPU memory does not follow from the atom count alone; measure {ref}`GPU memory <troubleshooting-gpu-memory>` on a short trial with the same settings.
- **Chain column**: the bundled PDB has an empty chain column; in such a PDB, give residues by name or by number. In a PDB with chains, write each residue as [chain:name:number](cli-conventions.md#residue-selectors), such as `-c 'A:TYR:44'`.
- **Caps only at CA and CB**: `extract` adds cap hydrogens only where it cuts an amino acid at CA or CB. Other cut bonds get no cap.
- **CA–N warnings in the default model**: the default cut also prints the warning about bonds other than C–C at each main-chain CA–N cut. The cap hydrogen there goes on CA; these cuts are the main-chain ends of the checklist and need no action, so look at the other bonds that the warning lists.
- **Same atoms in every input**: inputs with different atom counts or order stop with `[multi] Atom count mismatch` or `[multi] Atom order mismatch`.
- **`--no-freeze-links`** is for diagnostic runs that let the boundary relax on purpose. Keep `--freeze-links` on for production runs.

## See also

- [`extract`](extract.md) — extraction options, cap hydrogens, and non-standard residue names
- [`all`](all.md) — the full workflow; extracts the model with `-c`
- [`opt`](opt.md) — optimization with distance restraints
- [`scan`](scan.md) — staged scans with restraints
- [`freq`](freq.md) — PHVA and rigid modes with frozen atoms
- [Tips for studying reaction mechanisms](mechanism-tips.md) — when to enlarge the model
- [Refine an MLIP TS with DFT](dft-backend.md) — keeping the model small enough for DFT
- [Common options and selectors](cli-conventions.md) — residue and atom selectors
- [YAML Reference](yaml-reference.md) — `geom.freeze_atoms` and the other settings
- [Troubleshooting](troubleshooting.md) — extraction errors
- [Glossary](glossary.md) — active-site model, cluster model, cap hydrogen
