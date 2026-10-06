---
name: pdb2reaction-model-setup
description: "Building and sizing the pdb2reaction active-site cluster: which residues extract keeps (-c centers, -r radius, --radius-het2het, waters, backbone, --selected-resn), how the boundary is cut and capped, how to trim the cluster to lower cost (--exclude-backbone, --no-include-h2o, -r 0 with --selected-resn, a hand-edited cluster run without -c), how to enlarge it (larger -r or --radius-het2het, added residues, a radius sensitivity check), the manual boundary checklist, and the atom-set rules for R/IM/P paths and WT/mutant comparisons. TRIGGER on choosing extraction centers or a radius, a cluster that is too large, slow, or out of memory, a missing catalytic residue or water, editing a cluster by hand, or comparing variants. SKIP for extract flag syntax and frozen-atom mechanics (pdb2reaction-cli), file formats or charge and multiplicity (pdb2reaction-structure-io), and TS strategy (pdb2reaction-overview)."
---

# Building the cluster model

Start from `-c` (substrate and catalytic residues) at the default `-r 2.6`, check the boundary, caps, and charge, then trim (`--exclude-backbone`, `--no-include-h2o`, `-r 0 --selected-resn`) or enlarge (`-r`, `--radius-het2het`, more residues in `-c` or `--selected-resn`), and rerun the key step at another size when the boundary could change the mechanism.

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

## Trim to lower cost

DFT optimization is practical up to roughly 300 atoms (N + M). Remove main chains except between peptide-bonded `-c` residues (`--exclude-backbone`) or waters (`--no-include-h2o`); keep only chosen residues with `-r 0 --selected-resn` (disulfide partners and a proline's N-side neighbor still join); or trim an extracted PDB by hand and run `all` without `-c`, with `-q` or `-l`. Freezing distant atoms with `--freeze-atoms` makes `freq` a PHVA on the movable block, and the UMA finite-difference Hessian displaces only movable atoms.

Each removal changes the charge. Delete the cap H of removed residues, or the run stops with `isolated LKH/HL`; a model built without `extract` has no cap H, so freeze its boundary with `--freeze-atoms`, and regenerate those indices after re-extraction (`all` with `-c` numbers the original input). If `freq` runs out of CUDA memory, keep the default `--hessian-calc-mode FiniteDifference` or shrink the movable region. DFT cost does not follow atom count; time one structure before a batch.

## Enlarge when the model is too small

Raise `-r`, use `--radius-het2het` for hydrogen-bond and metal partners, and add catalytic or charge-compensating residues to `-c` (whole when `-r` is above 0) or `--selected-resn` (side chain unless a neighbor is in). Keep the main chain. No radius is universally safe: 2.6 Å is a starting value, not a chemically validated cutoff, and a larger model costs more without always being more accurate. When the boundary could change the mechanism, rerun the key step at another size and compare barriers.

## Same atoms across states and variants

Within one R/IM/P reaction path, all structures need the same atom identities and order, cluster boundary, and cap topology. Prefer one multi-input `extract`; mismatched inputs stop with `[multi] Atom count mismatch` or `[multi] Atom order mismatch`. States prepared separately must be harmonized first.

WT/mutant or other cross-variant models may differ in atoms. Keep residue positions, boundary policy, protonation, and cap rules the same except for the mutation, record the differences, and compare barriers, not total energies ([ΔΔG‡](../pdb2reaction-overview/ts-strategy.md#controlled-comparisons)). Separate automatic extractions can pick different residues near the cutoff; harmonize the positions. For one chemical composition, use one atom set and order.

## Next step

- [`pdb2reaction-cli/extract.md`](../pdb2reaction-cli/extract.md): flags, outputs, frozen atoms.
- [`pdb2reaction-overview/SKILL.md`](../pdb2reaction-overview/SKILL.md#pick-an-all-mode): pick an `all` mode.
- [`pdb2reaction-structure-io/SKILL.md`](../pdb2reaction-structure-io/SKILL.md#charge-and-multiplicity): charge and multiplicity.
- [Building the cluster model](../../docs/model-setup.md): the full guide.
