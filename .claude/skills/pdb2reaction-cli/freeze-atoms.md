# Frozen Atoms (cluster boundary)

Cluster models pulled from a larger protein need a small set of atoms held
in place at the truncation boundary, otherwise the optimizer pulls the
dangling fragment into something unphysical. `pdb2reaction` handles this
with **link hydrogens** + three layers of `freeze_atoms`.

## Background: link hydrogens (`LKH/HL`)

When `extract` cuts a covalent bond between an in-cluster atom (`A`) and
an out-of-cluster atom (`B`), it places a hydrogen along `A→B` at 1.09 Å.
The cap is written as a `HETATM` with residue name `LKH`, atom name `HL`.
The cluster-side parent atom `A` of each cap is what needs to be frozen.

## Three sources of frozen atoms (merged as a union at run time)

| Source | Where | Default | Applies to |
|---|---|---|---|
| `--freeze-links` / `--no-freeze-links` | CLI flag | `True` | PDB inputs (auto-detects `LKH` parents) |
| `--freeze-atoms 'i,j,k,...'` | CLI flag | empty | any input format; 1-based indices |
| `geom.freeze_atoms: [...]` | YAML via `--config` | empty | any input format; 1-based indices |

The three are unioned. There is no override mode — every atom mentioned
in any source is frozen.

## Recipes

### PDB extracted by pdb2reaction (typical case)

```bash
pdb2reaction extract -i complex.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' -o model.pdb
pdb2reaction opt -i model.pdb -q 0 -m 1            # --freeze-links is True by default
pdb2reaction tsopt -i model.pdb -q 0 -m 1
pdb2reaction freq -i tsopt_final.pdb -q 0 -m 1     # auto PHVA on frozen-atom set
```

### XYZ / GJF input (no LKH residue → must specify explicitly)

```bash
pdb2reaction tsopt -i ts_candidate.xyz -q 0 -m 1 \
  --freeze-atoms '12,15,28,29,42'
```

Or inherit topology from a reference PDB:

```bash
pdb2reaction tsopt -i ts_candidate.xyz -q 0 -m 1 \
  --ref-pdb model.pdb                              # --freeze-links resurrected
```

### Long list — keep it in YAML

```yaml
# tsopt.yaml
geom:
  freeze_atoms: [12, 15, 28, 29, 42, 88, 91, 92]
```

```bash
pdb2reaction tsopt -i ts.xyz -q 0 -m 1 --config tsopt.yaml
```

CLI `--freeze-atoms` and YAML `geom.freeze_atoms` are merged (union), so
you can put a long stable list in the YAML and add a few one-off indices
on the CLI.

## Effect on the calculation

- **Forces:** zeroed for every frozen DOF; optimizer cannot move them.
- **Hessian:** rows and columns of frozen DOFs are either removed
  (`calc.return_partial_hessian: true`, the default for `freq`, forced
  for `irc`) or zeroed in the full matrix.
- **`freq`:** when the frozen-atom set is non-empty, automatically runs
  Partial Hessian Vibrational Analysis (PHVA) on the active block.
  Imaginary-mode thresholds (5 cm⁻¹ for non-TS, 100 cm⁻¹ for TS) apply
  to the resulting eigenvalues.
- **`irc` / `path-opt` / `path-search`:** frozen atoms keep their
  initial Cartesian coordinates at every image and every step.

## Subcommand coverage

All subcommands that touch geometry honor all three sources:

| Subcommand | `--freeze-links` (PDB only) | `--freeze-atoms` | YAML `geom.freeze_atoms` |
|---|:---:|:---:|:---:|
| `extract` | inserts `LKH/HL`; flag is for downstream | n/a | n/a |
| `opt`, `tsopt`, `freq`, `irc` | ✓ | ✓ | ✓ |
| `path-opt`, `path-search` | ✓ | ✓ | ✓ |
| `scan`, `scan2d`, `scan3d` | ✓ | ✓ | ✓ |
| `all` | ✓ | ✓ | ✓ |

## Common pitfalls

- **`LKH/HL` records manually deleted** → `--freeze-links` finds nothing.
  Either rerun `extract` or specify `--freeze-atoms` explicitly.
- **Re-numbered atoms** → `--freeze-atoms` is tied to input atom order
  (1-based). If you re-extract with different residue selection,
  regenerate the index list.
- **XYZ / GJF without `--ref-pdb`** → `--freeze-links` is a no-op
  (no `LKH` residue to detect). Provide `--ref-pdb FILE` or pass an
  explicit `--freeze-atoms` list.
- **`--no-freeze-links`** is a diagnostic flag (let the boundary relax
  on purpose to inspect the result). Production cluster-model runs
  should leave `--freeze-links` on.

## See also

- `pdb2reaction-cli/extract.md` — link-hydrogen insertion details.
- `pdb2reaction-structure-io/pdb.md` — the `LKH/HL` record format.
- User-facing docs: `docs/freeze-atoms.md` (full guide with cross-refs).