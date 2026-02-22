# Quickstart: `pdb2reaction scan` with `--spec`

## Goal

Drive one or more bond distances from a single structure using a YAML scan specification.

## Prerequisites

- Input structure: `input.pdb`
- Charge and multiplicity (`-q`, `-m`) for the target state

## 1. Prepare `scan.yaml`

```yaml
one_based: true
stages:
 - [["TYR,285,CA", "MMT,309,C10", 1.35]]
 - [["TYR,285,CA", "MMT,309,C10", 2.20], ["TYR,285,CB", "MMT,309,C11", 1.80]]
```

## 2. Run scan

```bash
pdb2reaction scan -i input.pdb -q 0 -m 1 --spec scan.yaml --print-parsed --out-dir ./result_scan
```

## What to check

- `result_scan/stage_01/result.pdb`
- `result_scan/stage_02/result.pdb`
- Optional trajectories when `--dump` is enabled (`scan.trj`, `scan.pdb`)

## Notes

- `--spec` is the recommended format; `--scan-lists` remains available for advanced/compatibility use.
- Use `pdb2reaction scan --help-advanced` to inspect all scan controls.

## Next step

- Feed scan results to path refinement with [all](all.md) or [path-search](path_search.md).
