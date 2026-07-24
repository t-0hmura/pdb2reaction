# Frozen Atoms (cluster boundary)

Cluster models pulled from a larger protein need a small set of atoms held
in place at the truncation boundary, otherwise the optimizer pulls the
dangling fragment into an unphysical geometry. `pdb2reaction` handles this
with **cap hydrogens** + three layers of `freeze_atoms`.

## Background: cap hydrogens (`LKH/HL`)

When `extract` recognizes one of its supported carbon-boundary truncations
between an in-cluster carbon (`A`) and an out-of-cluster atom (`B`), it places a
hydrogen along `A→B` at 1.09 Å. It does not generically cap arbitrary bond
types or non-carbon boundaries. The cap is written as a
`HETATM` with residue name `LKH`, atom name `HL`. The cluster-side parent
atom `A` of each cap is what needs to be frozen.

## Three sources of frozen atoms (merged as a union at run time)

| Source | Where | Default | Applies to |
|---|---|---|---|
| `--freeze-links` / `--no-freeze-links` | CLI flag | `True` | PDB/mmCIF topology (auto-detects `LKH` parents after bridging) |
| `--freeze-atoms 'i,j,k,...'` | CLI flag | empty | any input format; 1-based indices |
| `geom.freeze_atoms: [...]` | YAML via `--config` | empty | any input format; 1-based indices |

The three are unioned. There is no override mode — every atom mentioned
in any source is frozen.

## Recipes

### PDB extracted by pdb2reaction (typical case)

```bash
pdb2reaction extract -i complex.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' -o model.pdb
pdb2reaction opt -i model.pdb -l 'SAM:1,GPP:-3' -m 1  # --freeze-links defaults on
pdb2reaction tsopt -i model.pdb -l 'SAM:1,GPP:-3' -m 1 -o result_tsopt
pdb2reaction freq -i result_tsopt/final_geometry.pdb \
  -l 'SAM:1,GPP:-3' -m 1 -o result_freq  # auto PHVA
```

### XYZ / GJF input (no LKH residue → must specify explicitly)

```bash
pdb2reaction tsopt -i ts_candidate.xyz -q 0 -m 1 \
  --freeze-atoms '12,15,28,29,42'
```

Or inherit topology from a reference PDB/mmCIF:

```bash
pdb2reaction tsopt -i ts_candidate.xyz -q 0 -m 1 \
  --ref-pdb model.pdb                              # --freeze-links re-enabled
```

### Long list — keep it in YAML

```yaml
# tsopt.yaml
geom:
  freeze_atoms: [12, 15, 28, 29, 42, 88, 91, 92]
  tr_projection: constrained
```

```bash
pdb2reaction tsopt -i ts.xyz -q 0 -m 1 --config tsopt.yaml
```

CLI `--freeze-atoms` and YAML `geom.freeze_atoms` are merged (union), so
you can put a long stable list in the YAML and add a few one-off indices
on the CLI.

## Effect on the calculation

- **Forces:**
  - Hard freeze (`opt` / `tsopt` / `scan` / `freq` / `irc`, `path-opt` / `path-search` `--mep-mode gsm`): forces on frozen DOFs are zeroed; optimizer cannot move them.
  - Soft restraint (`path-opt` / `path-search` `--mep-mode dmf`): `HarmonicFixAtoms` adds a stiff harmonic penalty around the initial Cartesians; frozen atoms can drift slightly.
- **Hessian:** rows and columns of frozen DOFs are either removed
  (`calc.return_partial_hessian: true`, the global calculator default
  and explicitly re-asserted by `opt` / `tsopt` / `scan` / `freq` / `irc`)
  or zeroed in the full matrix.
- **`freq`:** when the frozen-atom set is non-empty, automatically runs
  Partial Hessian Vibrational Analysis (PHVA) on the active block.
- **`irc` and GSM `path-opt` / `path-search`:** hard-frozen atoms keep their
  initial Cartesian coordinates. DMF uses the soft restraint described above,
  so small displacement is possible and must be inspected.

## Rigid-mode treatment for PHVA

Use `geom.tr_projection` or `--tr-projection` in `freq`, `irc`, `opt`, `tsopt`, and `all`:

| Value | Behavior |
|---|---|
| `constrained` (default) | Remove only full-system rigid motions that leave every frozen anchor fixed. Generic effective ranks are 6 / 3 / 1 / 0 for 0 / 1 / 2 / at least 3 non-collinear anchors. A normal multi-anchor cluster boundary therefore usually has rank 0. |
| `legacy-active` | Deprecated: treat the active fragment as isolated for comparison only. Never use it for pass/HOSP transition-state certification. The current common kernel handles degeneracies; near-linear or degenerate cases do not guarantee bitwise replay of older results. |

This flag controls Cartesian PHVA-related eigensolvers: `freq`, the initial IRC mode, Dimer orientation, exact TS validation, and opt/TS flattening. It is unrelated to `tsopt --ref-mode`, which supplies an MEP reaction direction. An all-frozen structure is invalid because it has no active vibrational DOF.

With JSON output, inspect `result.json["rigid_projection"]` for `treatment`, `effective_rank`, and Hessian source/shape. `freq --dump` writes the same provenance to `thermoanalysis.yaml`.

## Subcommand coverage

The geometry workflows below honor the listed sources. `sp` has no freeze
CLI flags, but it honors 1-based YAML `geom.freeze_atoms`; with `--hess`, the
default partial-Hessian setting writes the active block. `extract` only creates
the cap records for downstream use.

| Subcommand | `--freeze-links` (PDB, or XYZ/GJF + `--ref-pdb`) | `--freeze-atoms` | YAML `geom.freeze_atoms` |
|---|:---:|:---:|:---:|
| `extract` | inserts `LKH/HL`; flag is for downstream | n/a | n/a |
| `opt`, `tsopt`, `freq`, `irc` | ✓ | ✓ | ✓ |
| `path-opt`, `path-search` | ✓ | ✓ | ✓ |
| `scan`, `scan2d`, `scan3d` | ✓ | ✓ | ✓ |
| `all` | ✓ | ✗ (use YAML `geom.freeze_atoms`) | ✓ |
| `sp` | ✗ | ✗ | ✓ |

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
- **All atoms frozen** → PHVA and IRC stop with an explicit error. Keep at least one active atom.

## See also

- `pdb2reaction-cli/extract.md` — cap-hydrogen insertion details.
- `pdb2reaction-structure-io/pdb.md` — the `LKH/HL` record format.
- User-facing docs: `docs/freeze-atoms.md` (full guide with cross-refs).
