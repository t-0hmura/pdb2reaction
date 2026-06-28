# Frozen Atoms

## Overview

Cluster models need a small set of atoms held in place at the truncation boundary so the optimizer cannot pull the dangling fragment into something unphysical. `pdb2reaction` handles this through **link hydrogens** (added at severed bonds by `extract`) and three layers of `freeze_atoms` specification.

When a residue is sliced out of a larger protein using the `extract` sub-command, the bond at the boundary is capped with a **link hydrogen** (residue `LKH`, atom `HL`, 1.09 Å along the original bond vector). If the parent atom of that link hydrogen is left free, gradient descent will quietly relax the cap+parent pair into a different geometry, deforming the boundary. Freezing the relevant atoms keeps the boundary stationary throughout optimization, MEP search, IRC, and vibrational analysis.

## Three ways to specify frozen atoms

### 1. `--freeze-links/--no-freeze-links` (default `True`)

Auto-freezes the parent atoms of link hydrogens added by `extract`. Active by default in every downstream subcommand. Recipe:

```bash
pdb2reaction extract -i complex.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' -o model.pdb
pdb2reaction opt -i model.pdb -q 0 -m 1   # --freeze-links is True; LKH parents auto-frozen
```

For XYZ/GJF inputs, no `LKH` atoms are present, so `--freeze-links` has no effect; use the next two methods instead. `--ref-pdb FILE` lets XYZ/GJF runs inherit a PDB topology and resurrect link-hydrogen detection.

### 2. `--freeze-atoms 'i,j,k,...'` (CLI explicit list)

Comma-separated **1-based** atom indices, applicable to any input format. Complements `--freeze-links` (the union is frozen at run time).

```bash
pdb2reaction tsopt -i ts_candidate.xyz -q 0 -m 1 \
  --freeze-atoms '12,15,28,29,42'
```

### 3. YAML `geom.freeze_atoms` (via `--config`)

```yaml
geom:
  freeze_atoms: [12, 15, 28, 29, 42]   # 1-based indices
```

Useful when the list is long or you want to ship it with the rest of the run configuration. CLI and YAML lists are merged, not replaced.

```bash
pdb2reaction tsopt -i ts.xyz -q 0 -m 1 --config tsopt.yaml
```

## How the three sources combine

The frozen-atom set used at run time is the **union** of:

- YAML `geom.freeze_atoms` (`--config FILE`), plus
- CLI `--freeze-atoms`, plus
- atoms detected as `LKH` parents when `--freeze-links` is on.

There is no mode that substitutes one for another; every entry that appears in any source is frozen.

## Effect on the calculation

- **Forces:** zeroed for every frozen DOF in `opt` / `tsopt` / `scan` / `freq` / `irc` and in `path-opt` / `path-search` `--mep-mode gsm` (hard freeze; the optimizer cannot move them).
- **Hessian:** rows and columns of frozen DOFs are either removed (`calc.return_partial_hessian: true`, the global calculator default; explicitly forced again by `opt` / `tsopt` / `scan` / `freq` / `irc`) or zeroed in the full matrix.
- **Vibrational analysis:** when frozen atoms are present, `freq` automatically performs partial Hessian vibrational analysis (PHVA) on the active block.
- **`path-opt --mep-mode dmf` and `path-search --mep-mode dmf` (soft restraint):** instead of zeroing forces, these stages add a `HarmonicFixAtoms` calculator (default `k_fix = 300 eV/Å²`, ASE units) per image so frozen atoms relax with a harmonic restraint, not a hard constraint. Coordinates may drift slightly from the input geometry.
- **MEP / IRC:** `path-opt` / `path-search` `--mep-mode gsm` and `irc` apply the hard freeze along the resolved path / IRC trajectory; `--mep-mode dmf` (path-opt or path-search) uses the soft restraint above.

## Subcommand coverage

| Subcommand | `--freeze-links` (PDB) | `--freeze-atoms` (any input) | YAML `geom.freeze_atoms` |
|---|:---:|:---:|:---:|
| [`extract`](extract.md) | (inserts `LKH/HL`; flag is for downstream stages) | n/a | n/a |
| [`opt`](opt.md) | yes | yes | yes |
| [`tsopt`](tsopt.md) | yes | yes | yes |
| [`freq`](freq.md) | yes (triggers PHVA) | yes | yes |
| [`irc`](irc.md) | yes | yes | yes |
| [`path-opt`](path-opt.md) | yes | yes | yes |
| [`path-search`](path-search.md) | yes | yes | yes |
| [`scan`](scan.md) / [`scan2d`](scan2d.md) / [`scan3d`](scan3d.md) | yes | yes | yes |
| [`all`](all.md) | yes | yes | yes |

## Common pitfalls

- **Manually deleted `LKH/HL` records.** `--freeze-links` finds nothing to freeze. Use `--freeze-atoms` to specify the boundary explicitly, or rerun `extract`.
- **Re-numbered atoms.** `--freeze-atoms` and `geom.freeze_atoms` are 1-based and tied to input atom order; if you re-extract and the order changes, regenerate the index list.
- **XYZ/GJF without a topology.** No `LKH` records exist, so `--freeze-links` is a no-op. Provide `--ref-pdb FILE` or an explicit `--freeze-atoms` list.
- **`--no-freeze-links`.** Disables the auto-freeze. Useful only for diagnostic runs that intentionally let the boundary relax; production cluster-model runs should leave `--freeze-links` on.

## See Also

- [`extract`](extract.md) — Where link hydrogens are inserted (residue `LKH`, atom `HL`, 1.09 Å along severed bond vector). Full algorithmic detail at {ref}`Link hydrogen and frozen atoms <link-hydrogen-and-frozen-atoms>`.
- [YAML Reference](yaml-reference.md) — `geom.freeze_atoms` schema and merge order.
- [CLI Conventions](cli-conventions.md) — flag conventions shared across subcommands.
- [Glossary](glossary.md) — Active Site Model, Cluster Model, Link Hydrogen.
