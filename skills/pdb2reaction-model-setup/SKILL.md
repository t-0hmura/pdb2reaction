---
name: pdb2reaction-model-setup
description: "Input structures and the active-site cluster for pdb2reaction: PDB, mmCIF (including very large inputs), XYZ, and GJF layout and which subcommand reads which, residue and atom selectors, the charge and multiplicity decision (-q, -l, -m), and building and sizing the cluster with extract (-c centers, -r radius, --radius-het2het, waters, backbone, --selected-resn): boundary cuts and caps with a manual checklist, trimming (--exclude-backbone, --no-include-h2o, -r 0 with --selected-resn, a hand-edited cluster without -c), enlarging with a radius sensitivity check, and atom-set rules for R/IM/P paths and WT/mutant comparisons. TRIGGER on editing or inspecting a structure file, writing a selector, deciding -q, -l, or -m, choosing centers or a radius, a cluster that is too large, slow, or out of memory, a missing residue or water, a hand-edited cluster, or comparing variants. SKIP for subcommand flags and frozen-atom mechanics (pdb2reaction-cli), outputs and TS strategy (pdb2reaction-overview), install, or HPC."
---

# Structures and the cluster model

PDB or mmCIF keeps residue names, so `-l 'RES:Q'` derives the charge; XYZ needs `-q` (and `--ref-pdb` for residue context); GJF takes charge and spin from its header; set `-m` for open-shell systems.

All four formats use Å and standard element symbols. Column layouts, cap
hydrogens, and ligand, ion, and metal tables are in [formats.md](formats.md).

Start from `-c` (substrate and catalytic residues) at the default `-r 2.6`, check the boundary, caps, and charge, then trim (`--exclude-backbone`, `--no-include-h2o`, `-r 0 --selected-resn`) or enlarge (`-r`, `--radius-het2het`, more residues in `-c` or `--selected-resn`), and rerun the key step at another size when the boundary could change the mechanism.

## Which format

| Format | Carries | Use it for | How to set `-q` / `-m` |
|---|---|---|---|
| PDB | Atom and residue names, chain, occupancy, B-factor, element | A fresh extraction of a standard-size model | `-l 'RES:Q,...'` for unknown/non-standard ligand charges; standard amino acids and recognized ions use internal tables |
| mmCIF (`.cif`, `.mmcif`) | PDB fields without the one-character chain, four-digit residue, and five-digit atom-serial limits | PDB Bank mmCIF, multi-character chains, or 10,000 or more residues | Same residue rules as PDB; prefer mmCIF over a PDB with widened columns |
| XYZ | Element and Cartesian coordinates only | Optimized geometries, trajectories, IRC endpoints | `-q` and `-m` explicitly, or `--ref-pdb` to the original PDB so `-l` works |
| GJF | Element, coordinates, charge, and multiplicity | Gaussian-style inputs | Read from the header; `-q` / `-m` override it |

## Which subcommand reads which format

`path-search`, `path-opt`, `opt`, `tsopt`, `freq`, `irc`, `sp`, `dft`, `scan`,
`scan2d`, `scan3d`, `bond-summary`, and `all` read all four formats. `extract`
reads PDB and mmCIF only, so `all` cuts a cluster with `-c` only from PDB or
mmCIF; with XYZ or GJF, omit `-c`. `fix-altloc` and `add-elem-info` read PDB
only, `trj2fig` reads a trajectory XYZ, and `energy-diagram` takes no
structure.

mmCIF input, and a PDB too large for its fixed columns, is converted to a
`.pdb` with temporary chain IDs and residue numbers for the calculation. With
`--convert-files` enabled (the default), each coordinate output also gets a
`.cif` with the original chain IDs, residue numbers, insertion codes, residue
names, and atom names; `extract` writes that `.cif` automatically. Report
residues from the `.cif`, not from the `.pdb`.

An XYZ passed to a command that needs residue context, such as
`-l 'GLU:-1'`, needs `--ref-pdb` pointing to a PDB or mmCIF with the same atoms.
`sp` has no `--ref-pdb`: give it the PDB or mmCIF directly, or an XYZ with
explicit `-q` and `-m`.

## Selecting residues and atoms

`-c` on `extract` and `all`, and `--selected-resn`, take these residue forms:

| Form | Example | Selects |
|---|---|---|
| Chain + name + number (recommended) | `'A:SAM:321'`, `'A:TYR:44,A:SAM:321'` | Exactly one residue per entry |
| Chain + name | `'A:SAM'` | Every SAM in chain A; warns when more than one matches |
| Chain + number | `'A:44'`, `'A:123A'` | One residue; a trailing letter is the insertion code |
| Name | `'SAM,GPP,MG'` | Every residue with that name in any chain |
| Number | `'44,63,186'` | That number in every chain |
| Structure file | `substrate.pdb` | Residues whose coordinates match a separate PDB or mmCIF |

```bash
pdb2reaction extract -i complex.pdb -c 'A:SAM:321' -o cluster.pdb
pdb2reaction extract -i complex.pdb -c 'A:SAM' -o cluster.pdb
pdb2reaction extract -i complex.pdb -c 'SAM,GPP,MG' -o cluster.pdb
```

Use the chain-qualified form with a number when one residue is intended,
especially when a ligand name repeats across chains. A PDB with an empty chain
column, such as the bundled examples, takes only the name or number forms.
Use one kind per list: names only, chain + name with or without a number, or
numbers with or without a chain. `A:TYR:44,A:SAM` works; `A:SAM,SAM` stops with
an error. `TYR:44` is read as chain `TYR`, residue 44.

Atom selectors, used by scans and distance restraints, take three fields
(residue name, residue number, atom name) in any order, such as `'SAM 320 CS1'`
or `'SER:11:HG'`, or the four-field form `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` in
that order, such as `'A:SAM:12B:C1'`. Use three fields on a PDB with an empty
chain column; `_` does not stand for an empty chain, so `'_:SER:11:HG'` matches
nothing. Chain IDs are case-sensitive, and mmCIF chain IDs may be longer than
one character.

## Charge and multiplicity

Every run needs a total charge and a multiplicity. A wrong value gives
wrong-chemistry results, not an error.

For PDB or mmCIF, give `-l 'RES:Q'` for unknown/non-standard residues only and
let pdb2reaction sum the total. Standard amino acids and recognized ions come
from internal tables, and waters and cap hydrogens are neutral, so the total
follows the residues the model keeps. It is still only chemically correct when
residue names, protonation and oxidation states, ligand mappings, and the
cluster boundary are correct; read the charge breakdown printed as
`Total active site model charge` whenever the model changes. With `-c`, the
mapping feeds the extraction's charge summary; without `-c`, it is applied to
the whole input structure. A single number, as in `-l -3`, is the total ligand
charge.

Use `-q INTEGER` when there are no residues to sum, as for an XYZ without
`--ref-pdb`, to replace a GJF header on purpose, or to override the `-l`
derivation. In ordinary geometry commands `-q` wins over `-l` when both are
given. In `all -c`, extraction still derives and reports its value. Explicit
`-q` sets the total, and a mismatch produces a warning. `-l` alone on an XYZ
or GJF without `--ref-pdb` stops with an error.

The order is: explicit `-q`; then the `-l` total; then `calc.charge` from
`--config`; then, without `-l`, the standard residues and ions of the
extracted model in `all -c`, or the `.gjf` header; otherwise an error.

Ions follow the internal `ION` table. Listing an ion in `-l` with the table
value is accepted, and a different value is ignored with a warning; a mapping
does not override a standard amino acid or a recognized ion. For a non-default
protonation or oxidation state, use a distinct residue name when you build the
model, or give the verified total with `-q`. Dump both tables with:

```bash
python -c "from pdb2reaction.workflows.extract import AMINO_ACIDS, ION; print(dict(AMINO_ACIDS)); print(dict(ION))"
```

`-m 1` is the default but is valid only for a closed-shell state. Use `-m 2` for
a verified doublet and the verified integer (3 or more) for a high-spin state;
`-m` takes an integer, not `3+`. Derive charge and multiplicity for the whole
cluster from the mechanism and primary literature, especially for metals,
radicals, and antiferromagnetically coupled centers.

| `-m` | Typical use |
|---|---|
| 1 | Closed shell, only when the electron count and state are known to be closed-shell |
| 2 | Radicals, radical SAM enzymes, low-spin Fe(III) |
| 3 | O₂, some carbenes, high-spin Ni(II) |
| 4 | High-spin Co(II), Cr(III), V(II) |
| 5 | Mn(III), high-spin Fe(II) |
| 6 | High-spin Mn(II), high-spin Fe(III) |

These are examples, not a spin-state calculator; the metal table is in
[formats.md](formats.md#ligand-ion-and-metal-charges).

## Unknown substrate charge

Look up a ligand's formal charge in the mechanism's primary paper or deposited
structure documentation first, then PubChem or ChEBI (`Formal Charge`) and the
RCSB ligand page, for example `https://www.rcsb.org/ligand/SAM`. From SMILES,
use `sum(a.GetFormalCharge() for a in Chem.MolFromSmiles(smi).GetAtoms())`.
Record the source and the modeled protonation and oxidation state. If the
sources do not settle one state, stop and ask the user; never default a metal
or radical cluster to `-q 0 -m 1` without checking.

| Group at pH 7 | Typical charge |
|---|---|
| Carboxylate, Asp, Glu side chains | −1 each |
| Phosphate monoester (`-OPO₃²⁻`) | −2 |
| Phosphate diester (`-OPO₂⁻`) | −1 |
| Triphosphate, such as ATP | −4 |
| Sulfonium, such as the SAM cofactor | +1 |
| Lys and Arg side chains | +1 |
| His | 0 or +1 |

Mechanisms sometimes need an unusual protonation state; check the literature
for the system you model.

## Build the cluster

```bash
pdb2reaction extract -i complex.pdb -c 'A:SAM:321,A:GPP:322,A:MG:323' -l 'SAM:1,GPP:-3' -o model.pdb --out-json
```

Give `-c` the substrate, cofactors, metals, and catalytic residues as chain:name:number. A residue joins when any atom lies within `-r` of a `-c` atom; waters and main chains stay by default, and amino acids in `-c` stay whole when `-r` is above 0. `--radius-het2het` adds a second cutoff between atoms other than C and H; `--selected-resn` adds residues without a distance search. Extract R, IM, and P in one run so all models share the same residues and caps. `all` without `-c` uses the input as the cluster.

## Check the boundary and the charge

Success: the console prints `[extract] Atoms after truncation: N`, `[extract] Link-H to add: M`, and `Total active site model charge`; the model has N + M atoms (`n_atoms_extracted` + `n_link_hydrogens` in `result.json`). `all --dry-run` prints the same. Confirm in a viewer that hydrogen-bonding, catalytic, metal-coordination, and charge-compensating residues are in. For a boundary set or audited by hand:

- End each main-chain fragment at `CA` on both sides.
- Cut side chains, ligands, and cofactors at an aliphatic C–C single bond (`CA–CB` or farther out). Do not cut peptide C–N, polar C–N/C–O, aromatic, S–S, or metal-coordination bonds; include the partner or move the boundary.
- One cap H per cut bond, with the intended valence; cap parents are [frozen by default](../pdb2reaction-cli/extract.md#freeze-atoms-at-the-cluster-boundary).
- Recount charge and multiplicity; a wrong electron count makes the model invalid.

`extract` caps only CA and CB cuts; another cut bond between nonmetal atoms makes `extract` warn and `all` stop. `--no-freeze-links` is for diagnostics only.

Before a long job, use `--dry-run` on a calculation command. It exits before
the MLIP or DFT stages; `all -c ... --dry-run` runs the extraction in a
temporary directory to print and check the model's charge and electron parity,
then removes it. `extract` has no dry run but prints the per-residue charge
breakdown. `--show-config` prints the configuration and then runs the full job,
so it is not a preview.

## Trim to lower cost

DFT optimization is practical up to roughly 300 atoms (N + M). Remove main chains except between peptide-bonded `-c` residues (`--exclude-backbone`) or waters (`--no-include-h2o`); keep only chosen residues with `-r 0 --selected-resn` (disulfide partners and a proline's N-side neighbor still join); or trim an extracted PDB by hand and run `all` without `-c`, with `-q` or `-l`. Freezing distant atoms with `--freeze-atoms` makes `freq` a PHVA on the movable block, and the UMA finite-difference Hessian displaces only movable atoms.

Each removal changes the charge. Delete the cap H of removed residues, or the run stops with `isolated LKH/HL`; a model built without `extract` has no cap H, so freeze its boundary with `--freeze-atoms`, and regenerate those indices after re-extraction (`all` with `-c` numbers the original input). If `freq` runs out of CUDA memory, keep the default `--hessian-calc-mode FiniteDifference` or shrink the movable region. DFT cost does not follow atom count; time one structure before a batch.

## Enlarge when the model is too small

Raise `-r`, use `--radius-het2het` for hydrogen-bond and metal partners, and add catalytic or charge-compensating residues to `-c` (whole when `-r` is above 0) or `--selected-resn` (side chain unless a neighbor is in). Keep the main chain. No radius is universally safe: 2.6 Å is a starting value, not a chemically validated cutoff, and a larger model costs more without always being more accurate. When the boundary could change the mechanism, rerun the key step at another size and compare barriers.

## Same atoms across states and variants

Within one R/IM/P reaction path, all structures need the same atom identities and order, cluster boundary, and cap topology. Prefer one multi-input `extract`; mismatched inputs stop with `[multi] Atom count mismatch` or `[multi] Atom order mismatch`. States prepared separately must be harmonized first.

WT/mutant or other cross-variant models may differ in atoms. Keep residue positions, boundary policy, protonation, and cap rules the same except for the mutation, record the differences, and compare barriers, not total energies ([ΔΔG‡](../pdb2reaction-overview/ts-strategy.md#controlled-comparisons)). Separate automatic extractions can pick different residues near the cutoff; harmonize the positions. For one chemical composition, use one atom set and order.

## Editing approach

1. Read the file first: residues, atom count, and any charge and multiplicity.
2. Check that the change keeps the format rules, such as PDB column widths and
   the XYZ atom-count line.
3. For an unknown charge or multiplicity, confirm with the user or the
   literature before you choose ([Unknown substrate charge](#unknown-substrate-charge)).
4. Make the smallest edit, such as one residue rename or one charge change,
   not a full rewrite.

## Next step

- Formats, cap hydrogens, and ligand and ion tables: [formats.md](formats.md).
- `extract` options, outputs, and frozen atoms: [extract](../pdb2reaction-cli/extract.md).
- Pick an `all` mode: [pdb2reaction-overview](../pdb2reaction-overview/SKILL.md#pick-an-all-mode).
- Conventions shared by all commands: [pdb2reaction-cli](../pdb2reaction-cli/SKILL.md).
- What comes out of a run: [outputs](../pdb2reaction-overview/outputs.md).
- Docs: [Building the cluster model](../../docs/model-setup.md) and [Common options and selectors](../../docs/cli-conventions.md).
