# Quickstart: `pdb2reaction scan` — single-structure staged scan

## Goal

Run a staged scan from a single structure by driving one or more bond distances.

## Prerequisites

- Input structure: `input.pdb`
- Charge and multiplicity (`-q`, `-m`) for the target state

---

## Method A: YAML spec file (recommended for complex scans)

### 1. Prepare `scan.yaml`

Define each stage in order:

```yaml
one_based: true
stages:
 - [["TYR,285,CA", "SAM,309,C10", 1.35]]
 - [["TYR,285,CA", "SAM,309,C10", 2.20], ["TYR,285,CB", "SAM,309,C11", 1.80]]
```

### 2. Run scan

```bash
pdb2reaction scan -i input.pdb -q 0 -m 1 -s scan.yaml -o ./result_scan
```

---

## Method B: `-s/--scan-lists` inline literal (quick one-liners)

`-s/--scan-lists` also accepts Python-literal strings directly on the command line.

### Basic syntax

Each literal is a list of `(atom1, atom2, target_Å)` triples. One literal = one stage.

```bash
# Single stage, integer atom indices (1-based by default)
pdb2reaction scan -i input.pdb -q 0 -s '[(1, 5, 1.35)]' -o ./result_scan

# Single stage, PDB selector strings
pdb2reaction scan -i input.pdb -q 0 -s '[("TYR,285,CA", "SAM,309,C10", 1.35)]' -o ./result_scan
```

### PDB selectors

Atoms can be specified by residue name, residue number, and atom name. Token separators are flexible:

```bash
"TYR,285,CA"     # comma-separated
"TYR 285 CA"     # space-separated
"TYR/285/CA"     # slash-separated
"285,TYR,CA"     # order is flexible
```

### Multiple stages

Pass multiple literals — each becomes one sequential stage:

```bash
# Stage 1: drive one bond to 1.35 Å
# Stage 2: drive two bonds simultaneously
pdb2reaction scan -i input.pdb -q 0 -s \
  '[("TYR,285,CA","SAM,309,C10",1.35)]' \
  '[("TYR,285,CA","SAM,309,C10",2.20),("TYR,285,CB","SAM,309,C11",1.80)]' \
  -o ./result_scan
```

Stages run sequentially; each starts from the previous stage's relaxed result.

### Quoting rules

```bash
# Correct: single-quote the outer list, double-quote selector strings inside
-s '[("TYR,285,CA","SAM,309,C10",1.35)]'

# Correct: integer indices need no inner quotes
-s '[(1, 5, 2.0)]'

# Avoid: double-quoting the outer literal requires escaping
-s "[(\"TYR,285,CA\",\"SAM,309,C10\",1.35)]"
```

> **Tip:** Use `--print-parsed` to verify that your scan targets were parsed correctly before a full run.

---

## What to check

- `result_scan/stage_01/result.pdb`
- `result_scan/stage_02/result.pdb` (if multiple stages)
- Concatenated scan trajectories (`scan_trj.xyz`, `scan.pdb`) are always written

## Notes

- `-s/--scan-lists` accepts either a YAML/JSON file path or inline Python literals (not both at once).
- Use `pdb2reaction scan --help-advanced` to inspect all scan controls.
- For full input-format details, see [scan](scan.md).

## Next step

- Feed scan results to path refinement with [all](all.md) or [path-search](path_search.md).
