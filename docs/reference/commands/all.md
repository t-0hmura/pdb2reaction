# `pdb2reaction all`

```text
Usage: pdb2reaction all [OPTIONS]

  Run active site model extraction → (optional single-structure staged scan) →
  MEP search → merge to full PDBs in one shot. If exactly one input is provided:
  (a) with --scan-lists, run staged scan on the active site model (or full
  structure when extraction is skipped) and use stage results as inputs for
  path-opt (path_search with --refine-path True); (b) with --tsopt True and no
  --scan-lists, run TSOPT-only mode.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
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
                                  directly as active site models.
  -o, --out-dir DIRECTORY         Top-level output directory for the pipeline.
                                  [default: result_all]
  -r, --radius FLOAT              Inclusion cutoff (Å) around substrate atoms.
                                  [default: 2.6]
  --radius-het2het FLOAT          Independent hetero–hetero cutoff (Å) for
                                  non‑C/H pairs.  [default: 0.0]
  --include-h2o BOOLEAN           Include waters (HOH/WAT/TIP3/SOL) in the
                                  active site model.  [default: True]
  --exclude-backbone BOOLEAN      Remove backbone atoms on non‑substrate amino
                                  acids (with PRO/HYP safeguards).  [default:
                                  False]
  --add-linkh BOOLEAN             Add cap hydrogens for severed bonds (carbon
                                  boundaries only) in active site models.
                                  [default: True]
  --selected-resn TEXT            Force-include residues by residue ID (not
                                  name; e.g. '123', 'A:123A', 'B:456');
                                  comma/space separated. Use '-c/--center
                                  GPP,SAM' for residue-name selection.
                                  [default: ""]
  --modified-residue TEXT         Comma-separated residue names (with optional
                                  charge) to treat as amino acids for backbone
                                  truncation and charge assignment. Examples:
                                  'HD1,HD2,HD3' or 'HD1:0,SEP:-2'.  [default:
                                  ""]
  -l, --ligand-charge TEXT        Total charge (number) or per-resname mapping
                                  like 'GPP:-3,SAM:1'. Used for extractor charge
                                  summaries; when extraction is skipped, PDB
                                  inputs derive the total charge and numeric
                                  values act as a total-charge fallback.
  -q, --charge INTEGER            Force the total system charge (overrides
                                  extractor/GJF/--ligand-charge-derived values;
                                  emits a warning when used).
  --workers INTEGER               MLIP predictor workers; >1 spawns a parallel
                                  predictor. NOTE: analytical Hessian raises a
                                  RuntimeError when workers>1; pass --hessian-
                                  calc-mode FiniteDifference explicitly.
                                  [default: 1]
  --workers-per-node INTEGER      Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  --solvent-model [alpb|cpcmx]    xTB solvent model.  [default: alpb]
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1).  [default: 1]
  --freeze-links BOOLEAN          Freeze parent atoms of cap hydrogens (PDB
                                  input or XYZ/GJF with --ref-pdb).  [default:
                                  True]
  --mep-mode [gsm|dmf]            MEP optimizer: Growing String Method (gsm) or
                                  Direct Max Flux (dmf).  [default: gsm]
  --dmf-backend [cpu|gpu]         DMF compute backend (--mep-mode dmf only): gpu
                                  (dmf.torch / CUDA) or cpu (dmf / NumPy). On a
                                  GPU out-of-memory, retry with cpu.  [default:
                                  gpu]
  --max-nodes INTEGER             Max internal nodes for **segment** GSM (String
                                  has max_nodes+2 images including endpoints).
                                  [default: 20]
  --max-cycles INTEGER            Maximum GSM optimization cycles.  [default:
                                  300]
  --climb BOOLEAN                 Enable climbing image for standard GSM
                                  segments (bridge segments always disable
                                  climbing).  [default: True]
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
  --refine-path BOOLEAN           Run a single-pass path-opt GSM between each
                                  adjacent pair and concatenate the segments
                                  (default; no path_search). Use --refine-path
                                  True to run recursive path_search on the full
                                  ordered series for automatic multistep
                                  discovery.  [default: False]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset (gau_loose|gau|gau_tight|ga
                                  u_vtight|baker|never). Defaults to 'gau' when
                                  not provided.
  --thresh-post [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for post-IRC endpoint
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
  --preopt BOOLEAN                If True, run initial single-structure
                                  optimizations of the active site model inputs.
                                  [default: True]
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
  --dft-engine [gpu|cpu]          DFT backend: gpu (GPU4PySCF, raises error if
                                  unavailable) or cpu (PySCF).  [default: gpu]
  -s, --scan-lists TEXT           Scan targets: inline Python literal or a
                                  YAML/JSON spec file path. Multiple inline
                                  literals define sequential stages, e.g.
                                  '[(12,45,1.35)]'
                                  '[(10,55,2.20),(23,34,1.80)]'. Indices refer
                                  to the original full input PDB (1-based). When
                                  extraction is used, they are auto-mapped to
                                  the active site model after extraction.
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
  --scan-preopt BOOLEAN           Override scan --preopt flag. Inherits from
                                  --preopt when omitted.
  --scan-endopt BOOLEAN           Override scan --endopt flag. Defaults to
                                  False.
  --ref-pdb FILE                  Reference PDB for topology when -i provides
                                  XYZ inputs. Enables PDB output conversion in
                                  TSOPT-only, scan, and path_search pipelines.
  --coord-type [cart|dlc]         Optimization coordinate system (cart|dlc).
                                  cart is the robust default used in published
                                  numbers; dlc speeds up torsion-rich opts.
  --precision [fp32|fp64]         MLIP backend precision: fp32 (default) or
                                  fp64. Routed to backend-specific kwargs (UMA
                                  precision / ORB precision / MACE
                                  default_dtype). aimnet2: fp32 no-op; fp64
                                  rejected.
  --backend-model TEXT            Model variant for the selected --backend (e.g.
                                  uma-s-1p1 / uma-m-1p1 for uma,
                                  orb_v3_conservative_omol for orb, MACE-OMOL-0
                                  / MACE-OFF23_small for mace). Default: the
                                  backend's built-in model.
  --deterministic / --no-deterministic
                                  Strict bit-reproducible GPU runs
                                  (deterministic algorithms + index_reduce_
                                  shim). Slower; raises if unsupported. Default
                                  off.
  --allow-charge-mult-mismatch    Skip the cluster charge/multiplicity electron-
                                  parity check (logs that it was skipped). For
                                  an intentional open-shell or modified-residue
                                  cluster.
  -h, --help                      Show this message and exit.
```
