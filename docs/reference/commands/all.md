# `pdb2reaction all`

```text

Usage: pdb2reaction all [OPTIONS]

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
  -r, --radius FLOAT              Inclusion cutoff (Å) around substrate atoms.
                                  [default: 2.6]
  --radius-het2het FLOAT          Independent hetero–hetero cutoff (Å) for
                                  non‑C/H pairs.  [default: 0.0]
  --include-h2o BOOLEAN           Include waters (HOH/WAT/TIP3/SOL) in the
                                  pocket.  [default: True]
  --exclude-backbone BOOLEAN      Remove backbone atoms on non‑substrate amino
                                  acids (with PRO/HYP safeguards).  [default:
                                  False]
  --add-linkh BOOLEAN             Add link hydrogens for severed bonds (carbon-
                                  only) in pockets.  [default: True]
  --selected-resn TEXT            Force-include residues (comma/space separated;
                                  chain/insertion codes allowed).  [default: ""]
  -l, --ligand-charge TEXT        Total charge (number) or per-resname mapping
                                  like 'GPP:-3,SAM:1'. Used for extractor charge
                                  summaries; when extraction is skipped, PDB
                                  inputs derive the total charge and numeric
                                  values act as a total-charge fallback.
  -q, --charge INTEGER            Force the total system charge (overrides
                                  extractor/GJF/--ligand-charge-derived values;
                                  emits a warning when used).
  --workers INTEGER               MLIP predictor workers; >1 spawns a parallel
                                  predictor (Hessian computation not supported
                                  with workers>1).  [default: 1]
  --workers-per-node INTEGER      Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  --solvent-model [alpb|cpcmx]    xTB solvent model.  [default: alpb]
  --verbose BOOLEAN               Enable INFO-level logging inside extractor.
                                  [default: True]
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1).  [default: 1]
  --freeze-links BOOLEAN          Freeze parent atoms of link hydrogens (PDB
                                  input or XYZ/GJF with --ref-pdb).  [default:
                                  True]
  --mep-mode [gsm|dmf]            MEP optimizer: Growing String Method (gsm) or
                                  Direct Max Flux (dmf).  [default: gsm]
  --max-nodes INTEGER             Max internal nodes for **segment** GSM (String
                                  has max_nodes+2 images including endpoints).
                                  [default: 20]
  --max-cycles INTEGER            Maximum GSM optimization cycles.  [default:
                                  300]
  --climb BOOLEAN                 Enable transition-state climbing after growth
                                  for the **first** segment in each pair.
                                  [default: True]
  --opt-mode [grad|hess]          Optimizer mode forwarded to scan/tsopt and
                                  used for single optimizations: grad
                                  (=LBFGS/Dimer) or hess (=RFO/RSIRFO).
                                  [default: grad]
  --opt-mode-post [grad|hess]     Optimizer mode override for TSOPT/post-IRC
                                  endpoint optimizations. If unset, uses --opt-
                                  mode when explicitly provided; otherwise falls
                                  back to the default ('hess').  [default: hess]
  --dump BOOLEAN                  Dump GSM/MEP trajectories. Always forwarded to
                                  path_search/path-opt; scan/tsopt receive it
                                  only when explicitly set here. The freq stage
                                  uses dump=True by default; set --dump False
                                  explicitly to disable it.  [default: False]
  --convert-files BOOLEAN         Convert XYZ/TRJ outputs into PDB/GJF
                                  companions based on the input format.
                                  [default: True]
  --refine-path BOOLEAN           If True, run recursive path_search on the full
                                  ordered series; if False, run a single-pass
                                  path-opt GSM between each adjacent pair and
                                  concatenate the segments (no path_search).
                                  [default: True]
  --thresh TEXT                   Convergence preset (gau_loose|gau|gau_tight|ga
                                  u_vtight|baker|never). Defaults to 'gau' when
                                  not provided.
  --thresh-post TEXT              Convergence preset for post-IRC endpoint
                                  optimizations (gau_loose|gau|gau_tight|gau_vti
                                  ght|baker|never).  [default: baker]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running any stage.  [default: no-dry-
                                  run]
  --resume / --no-resume          Resume a previous run from --out-dir.
                                  Completed stages whose output files already
                                  exist are skipped. Useful when a long pipeline
                                  was interrupted.  [default: no-resume]
  --preopt BOOLEAN                If True, run initial single-structure
                                  optimizations of the pocket inputs.  [default:
                                  True]
  --hessian-calc-mode [finitedifference|analytical]
                                  Common MLIP Hessian calculation mode forwarded
                                  to tsopt and freq. Defaults to
                                  'FiniteDifference'.
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
  --tsopt-max-cycles INTEGER      Override tsopt --max-cycles value. Defaults to
                                  10000 when not provided.
  --tsopt-out-dir DIRECTORY       Override tsopt output subdirectory (relative
                                  paths are resolved against the default).
  --flatten / --no-flatten        Enable the extra-imaginary-mode flattening
                                  loop in tsopt (grad: dimer loop, hess: post-
                                  RSIRFO); --no-flatten forces
                                  flatten_max_iter=0.  [default: no-flatten]
  --freq-out-dir DIRECTORY        Override freq output base directory (relative
                                  paths resolved against the default).
  --freq-max-write INTEGER        Override freq --max-write value. Defaults to
                                  10.
  --freq-amplitude-ang FLOAT      Override freq --amplitude-ang (Å). Defaults to
                                  0.8.
  --freq-n-frames INTEGER         Override freq --n-frames value. Defaults to
                                  20.
  --freq-sort [value|abs]         Override freq mode sorting. Defaults to
                                  'value'.
  --freq-temperature FLOAT        Override freq thermochemistry temperature (K).
                                  Defaults to 298.15 K.
  --freq-pressure FLOAT           Override freq thermochemistry pressure (atm).
                                  Defaults to 1.0 atm.
  --dft-out-dir DIRECTORY         Override dft output base directory (relative
                                  paths resolved against the default).
  --dft-func-basis TEXT           Override dft --func-basis value. Defaults to
                                  'wb97m-v/def2-tzvpd'.
  --dft-max-cycle INTEGER         Override dft --max-cycle value. Defaults to
                                  100.
  --dft-conv-tol FLOAT            Override dft --conv-tol value. Defaults to
                                  1e-9.
  --dft-grid-level INTEGER        Override dft --grid-level value. Defaults to
                                  3.
  --dft-engine [gpu|cpu|auto]     Preferred DFT backend: GPU (default), CPU, or
                                  auto (try GPU then CPU).  [default: gpu]
  -s, --scan-lists TEXT           Scan targets: inline Python literal or a
                                  YAML/JSON spec file path. Multiple inline
                                  literals define sequential stages, e.g.
                                  '[(12,45,1.35)]'
                                  '[(10,55,2.20),(23,34,1.80)]'. Indices refer
                                  to the original full input PDB (1-based). When
                                  extraction is used, they are auto-mapped to
                                  the pocket after extraction.
  --scan-out-dir DIRECTORY        Override the scan output directory (default:
                                  <out-dir>/scan/). Relative paths are resolved
                                  against the default parent.
  --scan-one-based BOOLEAN        Override the scan subcommand indexing
                                  interpretation (True = 1-based, False =
                                  0-based). Defaults to 1-based.
  --scan-max-step-size FLOAT      Override scan --max-step-size (Å). Defaults to
                                  0.20 Å.
  --scan-bias-k FLOAT             Override scan harmonic bias strength k
                                  (eV/Å^2). Defaults to 300.
  --scan-relax-max-cycles INTEGER
                                  Override scan relaxation max cycles per step.
                                  Defaults to 10000.
  --scan-preopt BOOLEAN           Override scan --preopt flag. When omitted,
                                  this follows --preopt (default False).
  --scan-endopt BOOLEAN           Override scan --endopt flag. Defaults to
                                  False.
  --ref-pdb FILE                  Reference PDB for topology when -i provides
                                  XYZ inputs. Enables PDB output conversion in
                                  TSOPT-only, scan, and path_search pipelines.
  -h, --help                      Show this message and exit.
```
