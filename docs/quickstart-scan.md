# Quickstart: single-structure scan workflow (`--scan-lists`)

## Goal

Run the full `pdb2reaction all` workflow from a single structure by driving one or more bond distances via `-s/--scan-lists`. This automatically chains: staged scan → MEP refinement → (optional) TS optimization + IRC.

## Prerequisites

- Input structure: `.pdb`
- Charge (`-q/--charge` or `-l/--ligand-charge`) and multiplicity (`-m`) for the target state

---

## Method A: `-s/--scan-lists` inline literal (quick one-liners)

`-s/--scan-lists` accepts Python-literal strings directly on the command line.

### Basic syntax

Each literal is a list of 3-tuples `(atom1, atom2, target_distance_Å)`. Exactly three elements per tuple are required; the third is always the target distance in **ångströms**. One literal = one stage.

```bash
# Single stage, integer atom indices (1-based by default)
pdb2reaction -i input.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -m 1 \
 -s '[(1, 5, 1.35)]' -o ./result_scan

# Single stage, PDB selector strings
pdb2reaction -i input.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -m 1 \
 -s '[("TYR,285,CA", "SAM,309,C10", 1.35)]' -o ./result_scan
```

The `-c/--center` cluster selector is required when running on a protein–ligand PDB; omit `-c` (and pass `-q` directly instead of `-l`) for small-molecule `.pdb` / `.xyz` inputs. `-m/--multiplicity` defaults to `1` (singlet) but is shown explicitly here for clarity.

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
pdb2reaction -i input.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -m 1 -s \
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

> **Tip:** To verify that your scan targets were parsed correctly before a full run, invoke the standalone `scan`/`scan2d`/`scan3d` subcommand with `--print-parsed --dry-run` (these flags are not exposed on `pdb2reaction all`).

---

## Method B: YAML spec file (for complex scans)

### 1. Prepare `scan.yaml`

Define each stage in order:

```yaml
one_based: true
stages:
 - [["TYR,285,CA", "SAM,309,C10", 1.35]]
 - [["TYR,285,CA", "SAM,309,C10", 2.20], ["TYR,285,CB", "SAM,309,C11", 1.80]]
```

### 2. Run

```bash
pdb2reaction -i input.pdb -q 0 -m 1 -s scan.yaml -o ./result_scan
```

---

## Expected output

```text
result_scan/
├── summary.log
├── summary.json
├── scan/
│   ├── preopt/                    # Pre-optimized structure
│   ├── stage_01/                  # Scan stage 1 results
│   │   ├── scan_trj.xyz          # Scan trajectory
│   │   └── scan.pdb
│   └── stage_02/                  # Scan stage 2 (if multi-stage)
└── path_search/                   # MEP search (default, recursive); path_opt/ with --refine-path False
    ├── mep.pdb
    └── energy_diagram_UMA_all.png
```

**What to check:**

1. `scan/stage_01/scan_trj.xyz` — open in PyMOL to verify bond distances change as expected
2. `path_search/mep.pdb` — the optimized MEP trajectory
3. `summary.log` — barrier heights and bond change summary

**Tip:** Use `--print-parsed --dry-run` to verify scan targets before a full run:

```bash
pdb2reaction scan -i input.pdb -q 0 -s '[(1, 5, 1.35)]' --print-parsed --dry-run
```

## Notes

- `-s/--scan-lists` accepts either a YAML/JSON file path or inline Python literals (not both at once).
- Use `pdb2reaction all --help-advanced` to inspect all options including scan controls.
- For the standalone `scan` subcommand (without MEP refinement), see [scan](scan.md).

## Next step

- Full option reference: [all](all.md)
- TS optimization and validation: [Quickstart: `pdb2reaction tsopt`](quickstart-tsopt-freq.md)
