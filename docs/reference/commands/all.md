# pdb2reaction all

```
pdb2reaction ver. 0.2.1.dev116+gff86b942a

Usage: cli all [OPTIONS]

  Run pocket extraction → (optional single-structure staged scan) → MEP search →
  merge to full PDBs in one shot. If exactly one input is provided: (a) with
  --scan-lists, run staged scan on the pocket (or full structure when extraction
  is skipped) and use stage results as inputs for path_search; (b) with --tsopt
  True and no --scan-lists, run TSOPT-only mode.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Two or more **full structures** (PDB/XYZ/GJF)
                                  in reaction order (reactant [intermediates
                                  ...] product), or a single **full structure**
                                  (with --scan-lists or with --tsopt True).
                                  Extraction (-c/--center) requires PDB inputs.
                                  When using --scan-lists without extraction,
                                  the input may be PDB/XYZ/GJF (integer indices
                                  only for non-PDB inputs). Repeat -i/--input
                                  for each file.  [required]
  -c, --center TEXT               Substrate specification for the extractor: a
                                  PDB path, a residue-ID list like '123,124' or
                                  'A:123,B:456' (insertion codes OK: '123A' /
                                  'A:123A'), or a residue-name list like
                                  'GPP,SAM'. When omitted, extraction is skipped
                                  and the **full input structure(s)** are used
                                  directly as pockets.
  -o, --out-dir DIRECTORY         Top-level output directory for the pipeline.
                                  [default: result_all]
  -l, --ligand-charge TEXT        Total charge (number) or per-resname mapping
                                  like 'GPP:-3,SAM:1'. Used for extractor charge
                                  summaries; when extraction is skipped, PDB
                                  inputs derive the total charge and numeric
                                  values act as a total-charge fallback.
  -q, --charge INTEGER            Force the total system charge (overrides
                                  extractor/GJF/--ligand-charge-derived values;
                                  emits a warning when used).
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running any stage.  [default: no-dry-
                                  run]
  --tsopt BOOLEAN                 TS optimization + IRC per reactive segment (or
                                  TSOPT-only mode for single-structure), and
                                  build energy diagrams.  [default: False]
  --thermo BOOLEAN                Run freq on (R, TS, P) per reactive segment
                                  (or TSOPT-only mode) and build Gibbs free-
                                  energy diagram (MLIP).  [default: False]
  --dft BOOLEAN                   Run DFT single-point on (R, TS, P) and build
                                  DFT energy diagram. With --thermo True, also
                                  generate a DFT//MLIP Gibbs diagram.  [default:
                                  False]
  -s, --scan-lists TEXT           Scan targets: inline Python literal or a
                                  YAML/JSON spec file path. Multiple inline
                                  literals define sequential stages, e.g.
                                  '[(12,45,1.35)]'
                                  '[(10,55,2.20),(23,34,1.80)]'. Indices refer
                                  to the original full input PDB (1-based). When
                                  extraction is used, they are auto-mapped to
                                  the pocket after extraction.
  -h, --help                      Show this message and exit.
```
