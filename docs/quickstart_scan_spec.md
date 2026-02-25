# Quickstart: `pdb2reaction scan` — staged scan

## Goal

Drive one or more bond distances from a single structure using staged scans.
Two input formats are available: **`--spec`** (YAML file, recommended) and **`--scan-lists`** (inline Python literals).

## Prerequisites

- Input structure: `input.pdb`
- Charge and multiplicity (`-q`, `-m`) for the target state

---

## Option A: `--spec` (YAML file — recommended)

### 1. Prepare `scan.yaml`

```yaml
one_based: true
stages:
 - [["TYR,285,CA", "MMT,309,C10", 1.35]]
 - [["TYR,285,CA", "MMT,309,C10", 2.20], ["TYR,285,CB", "MMT,309,C11", 1.80]]
```

### 2. Run scan

```bash
pdb2reaction scan -i input.pdb -q 0 -m 1 --spec scan.yaml --print-parsed --out-dir ./result_scan
```

---

## Option B: `--scan-lists` (inline syntax)

The same two-stage scan can be expressed directly on the command line without a YAML file:

```bash
pdb2reaction scan -i input.pdb -q 0 -m 1 \
 --scan-lists '[("TYR,285,CA","MMT,309,C10",1.35)]' \
              '[("TYR,285,CA","MMT,309,C10",2.20),("TYR,285,CB","MMT,309,C11",1.80)]' \
 --print-parsed --out-dir ./result_scan
```

Each Python-literal string after `--scan-lists` defines one stage.
**Do not repeat the flag** — list all stages as consecutive arguments.

---

## What to check

- `result_scan/stage_01/result.pdb`
- `result_scan/stage_02/result.pdb`
- Optional trajectories when `--dump` is enabled (`scan_trj.xyz`, `scan.pdb`)

## How to choose

| | `--spec` (YAML) | `--scan-lists` (inline) |
|---|---|---|
| **Best for** | Multi-stage, complex scans | Quick single-stage runs, scripting |
| **Readability** | Easy to read and version-control | Shell quoting can be tricky |
| **Reproducibility** | YAML file is self-documenting | Must copy the full command line |

## Notes

- Use `pdb2reaction scan --help-advanced` to inspect all scan controls.
- For full format details, see [scan — --spec format](scan.md#--spec-format-recommended) and [scan — --scan-lists format](scan.md#--scan-lists-format).

## Next step

- Feed scan results to path refinement with [all](all.md) or [path-search](path_search.md).
