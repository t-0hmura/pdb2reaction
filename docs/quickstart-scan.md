# Quickstart: scan-seeded `pdb2reaction all` workflow

## Goal

Run `pdb2reaction all` from a single structure by defining one or more distance coordinates with `--scan-lists/-s`. The workflow performs the restrained scan, optimizes an MEP, and optionally continues to transition-state (TS) optimization and intrinsic reaction coordinate (IRC) analysis.

## Prerequisites

- Input structure: PDB/mmCIF; XYZ/GJF is also accepted when extraction is omitted and atom indices are used
- Charge (`-q/--charge` or `--ligand-charge/-l`) and multiplicity (`-m`) for the target state

## Selection among scan workflows

| Objective | Command | Coordinate treatment |
| --- | --- | --- |
| Generate restrained endpoint structures and trajectories only | `pdb2reaction scan` | Concerted multi-distance and multistage scans are supported |
| Continue from a scan-defined path to MEP and optional TS/IRC analysis | `pdb2reaction all -s ...` | The scan endpoints seed `path-opt` or, with `--refine-path`, `path-search` |
| Explore a two- or three-dimensional energy landscape and map a PES | `pdb2reaction scan2d` / `scan3d` | Cartesian product of independent distance grids |

For `scan` and scan-seeded `all`, one literal defines one stage. Multiple tuples
within that literal are advanced concertedly; multiple literals define
sequential stages.

## Scan targets via `--scan-lists/-s` inline literal (default)

`--scan-lists/-s` accepts Python-literal strings directly on the command line. For atom selector syntax (residue/atom tokens, separators, ordering) and outer/inner quoting rules, see {ref}`CLI Conventions: Scan-list spec <scan-list-spec>`.

### Basic syntax

Each literal is a list of 3-tuples `(atom1, atom2, target_distance_Å)`. Exactly three elements per tuple are required; the third is always the target distance in **ångströms**. One literal = one stage.

```bash
# Single stage, integer atom indices (1-based by default)
pdb2reaction all -i input.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -m 1 \
 -s '[(1, 5, 1.35)]' -o ./result_scan

# Single stage, PDB selector strings
pdb2reaction all -i input.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -m 1 \
 -s '[("TYR,285,CA", "SAM,309,C10", 1.35)]' -o ./result_scan
```

Use `-c/--center` when an active-site model must be extracted from a full-system
PDB/mmCIF. Omit it for an already extracted cluster model or for XYZ/GJF input;
in that case the supplied structure is analyzed as given. `-m/--multiplicity`
defaults to `1` (singlet) but is shown explicitly here for clarity.

### Multiple stages

Pass multiple literals — each becomes one sequential stage:

```bash
# Stage 1: drive one bond to 1.35 Å
# Stage 2: drive two bonds simultaneously
pdb2reaction all -i input.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -m 1 -s \
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
└── _work/                         # Pipeline scratch (safe to delete)
    ├── scan/
    │   ├── preopt/                # Pre-optimized structure
    │   ├── stage_01/              # Scan stage 1 results
    │   │   ├── result.{xyz,pdb}   # Final restrained endpoint (unbiased only with --scan-endopt)
    │   │   ├── scan_trj.xyz       # Scan trajectory
    │   │   └── scan.pdb
    │   └── stage_02/              # Scan stage 2 (if multi-stage)
    └── path_opt/                  # MEP search (path_search/ with --refine-path)
        └── hei_seg_01.{xyz,pdb}   # Highest-energy MEP image
```

This minimal command stops after the MEP stage and does **not** create
`segments/`. Add `--tsopt` to create canonical per-segment R/TS/P and IRC
outputs after successful validation; add `--thermo` for `freq/` outputs.

### Output validation

1. `_work/scan/stage_01/scan_trj.xyz` — open in PyMOL to verify bond distances change as expected
2. `mep.pdb` and `_work/path_opt/hei_seg_01.pdb` — inspect the optimized MEP and its highest-energy image
3. `summary.log` — barrier heights and bond change summary

Validate the complete `all` input, extraction, charge, and scan mapping before
execution with:

```bash
pdb2reaction all -i input.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
  -s '[(1, 5, 1.35)]' --dry-run
```

Standalone `scan --print-parsed` validates only the scan specification; it does
not validate `all`-specific extraction or atom-index remapping.

## Notes

- `-s/--scan-lists` accepts inline Python literals when used with `all`. The standalone `scan` subcommand additionally accepts a YAML/JSON spec file path (see [scan](scan.md)).
- The scan engine is shared, but the parent option names and preoptimization defaults differ:

  | Command | Step / bias | Relaxation limit | Preopt / endopt |
  | --- | --- | --- | --- |
  | `pdb2reaction all` | `--scan-max-step-size 0.20 Å`; `--scan-bias-k 300 eV/Å²` | `--scan-relax-max-cycles 100000` | scan preopt inherits parent `--preopt` (on by default); endopt is off (`--no-scan-endopt`) |
  | `pdb2reaction scan` | `--max-step-size 0.20 Å`; `--bias-k 300 eV/Å²` | `--relax-max-cycles 100000` | `--no-preopt`; `--no-endopt` |

  Override the step width with the corresponding CLI flag. Override the bias
  strength with that CLI flag or YAML `bias.k`. See [scan](scan.md) and
  [yaml-reference](yaml-reference.md#bias) for the per-stage controls.
- Each scan stage ends with a bond-change check (`has_bond_change`) on the final relaxed geometry. The per-stage result is recorded in the scan log, and — when `--out-json` is set — also in the aggregate `result.json` (under its `stages` array) written to the scan out-dir.
- The recursive MEP refinement (`path-search`) consumes the scan endpoints unconditionally. Whether it runs is gated by `--refine-path`, **not** by the scan-stage bond-change flag (`has_bond_change`).
- Successful scan completion establishes only that the restrained endpoint was generated. It does not establish an unconstrained minimum, an elementary reaction, or a transition state. Validate minima by unbiased optimization and validate TS candidates by saddle-order and IRC analysis.
- Use `pdb2reaction all --help-advanced` to inspect all options including scan controls.
- For the standalone `scan` subcommand (without MEP refinement), see [scan](scan.md).

## Next step

- Full option reference: [all](all.md)
- TS candidate validation: [Quickstart: TS-only mode](quickstart-tsopt-freq.md)
- Term definitions (MEP, TS, IRC, …): [Glossary](glossary.md)
