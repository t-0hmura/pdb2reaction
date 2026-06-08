# Quickstart: `pdb2reaction all --scan-lists` (Scan mode)

## Goal

Run the full `pdb2reaction all` workflow from a single structure by driving one or more bond distances via `--scan-lists/-s`. This automatically chains: staged scan → MEP refinement → (optional) TS optimization + IRC.

## Prerequisites

- Input structure: `.pdb`
- Charge (`-q/--charge` or `--ligand-charge/-l`) and multiplicity (`-m`) for the target state

## Method A: `--scan-lists/-s` inline literal (default)

`--scan-lists/-s` accepts Python-literal strings directly on the command line. For atom selector syntax (residue/atom tokens, separators, ordering) and outer/inner quoting rules, see {ref}`CLI Conventions: Scan-list spec <scan-list-spec>`.

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

## Expected output

```text
result_scan/
├── summary.log
├── summary.json
├── mep.pdb                        # Full MEP path (promoted to the root)
├── energy_diagram_MEP.png         # MEP energy plot (promoted to the root)
├── segments/
│   └── seg_01/                    # Per-reactive-segment deliverables
└── _work/                         # Pipeline scratch (safe to delete)
    ├── scan/
    │   ├── preopt/                # Pre-optimized structure
    │   ├── stage_01/              # Scan stage 1 results
    │   │   ├── scan_trj.xyz       # Scan trajectory
    │   │   └── scan.pdb
    │   └── stage_02/              # Scan stage 2 (if multi-stage)
    └── path_search/               # MEP search (path_opt/ with --refine-path False)
```

**What to check:**

1. `_work/scan/stage_01/scan_trj.xyz` — open in PyMOL to verify bond distances change as expected
2. `mep.pdb` — the optimized MEP trajectory (promoted to the output root)
3. `summary.log` — barrier heights and bond change summary

**Tip:** Use `--print-parsed` (and abort with Ctrl-C) to verify scan targets before letting the full run proceed:

```bash
pdb2reaction scan -i input.pdb -q 0 -s '[(1, 5, 1.35)]' --print-parsed
```

## Notes

- `-s/--scan-lists` accepts inline Python literals when used with `all`. The standalone `scan` subcommand additionally accepts a YAML/JSON spec file path (see [scan](scan.md)).
- Default scan stepping: `--scan-max-step-size 0.20 Å` (per-distance step) and `--scan-bias-k 300 eV/Å²` (harmonic bias strength) when invoked via `pdb2reaction all`. The standalone `pdb2reaction scan` command exposes the same defaults under the un-prefixed `--max-step-size` / `--bias-k` names; override via either form or via the YAML `bias` block. See [scan](scan.md) and [yaml-reference](yaml-reference.md#bias) for the per-stage controls.
- Each scan stage ends with a bond-change check (`has_bond_change`) on the final relaxed geometry; the per-stage result is recorded under each `stage_XX/result.json` (when `--out-json` is set) and in the scan log. The recursive MEP refinement (`path-search`) consumes the scan endpoints unconditionally — it is gated by `--refine-path`, not by the scan-stage bond-change flag.
- Use `pdb2reaction all --help-advanced` to inspect all options including scan controls.
- For the standalone `scan` subcommand (without MEP refinement), see [scan](scan.md).

## Next step

- Full option reference: [all](all.md)
- TS candidate validation: [Quickstart: TS-only mode](quickstart-tsopt-freq.md)
