# `pdb2reaction all`

```text
Usage: pdb2reaction all [OPTIONS]

  Run active site model extraction → optional staged scan → MEP search → full-
  structure merge in one run. If exactly one input is provided: (a) with --scan-
  lists, run staged scan on the active site model (or full structure when
  extraction is skipped) and use stage results as inputs for path-opt
  (path_search with --refine-path); (b) with --tsopt and no --scan-lists, run
  TSOPT-only mode.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Two or more **full structures**
                                  (PDB/mmCIF/XYZ/GJF) in reaction order
                                  (reactant [intermediates ...] product), or a
                                  single **full structure** (with --scan-lists
                                  or with --tsopt). Extraction (-c/--center)
                                  accepts PDB/mmCIF. mmCIF is processed through
                                  an internal PDB bridge and emitted again as
                                  mmCIF; oversized PDBs use the same safe
                                  bridge. Pass multiple files after one
                                  -i/--input, or repeat -i/--input for each
                                  file.  [required]
  -c, --center TEXT               Substrate specification for the extractor: a
                                  PDB/mmCIF path, a residue-ID list like
                                  '123,124' or 'A:123,B:456' (insertion codes
                                  OK: '123A' / 'A:123A'), a residue-name list
                                  like 'GPP,SAM', or chain-qualified 'A:SAM' /
                                  'A:SAM:123'. When omitted, extraction is
                                  skipped and the **full input structure(s)**
                                  are used directly as active site models.
  -o, --out-dir DIRECTORY         Top-level output directory for the pipeline.
                                  [default: result_all]
  -r, --radius FLOAT RANGE        Inclusion cutoff (Å) around substrate atoms.
                                  Zero is accepted and evaluated internally as
                                  0.001 Å (effectively off for ordinary radius-
                                  based neighbors).  [default: 2.6; x>=0.0]
  --radius-het2het FLOAT RANGE    Independent hetero–hetero cutoff (Å) for
                                  non‑C/H pairs.  [default: 0.0; x>=0.0]
  --include-h2o BOOLEAN           Include waters (HOH/WAT/TIP3/SOL) in the
                                  active site model.  [default: True]
  --exclude-backbone BOOLEAN      Remove backbone atoms on non‑substrate amino
                                  acids (with PRO/HYP safeguards).  [default:
                                  False]
  --add-linkh BOOLEAN             Add cap hydrogens for severed bonds (carbon
                                  boundaries only) in active site models.
                                  [default: True]
  --selected-resn TEXT            Force-include residues using the same
                                  selectors as -c/--center: IDs ('123',
                                  'A:123A'), names ('SAM'), or chain-qualified
                                  names ('A:SAM', 'A:SAM:123'); comma/space
                                  separated.  [default: ""]
  --modified-residue TEXT         Comma-separated residue names to treat as
                                  amino acids for backbone truncation and charge
                                  assignment. NAME:charge adds or overrides the
                                  nominal charge for this extraction; bare NAME
                                  defaults to 0. Examples: 'HD1,HD2,HD3' or
                                  'HD1:0,SEP:-2'.  [default: ""]
  -l, --ligand-charge TEXT        Total charge (number) or per-resname mapping
                                  like 'GPP:-3,SAM:1'. The per-resname mapping
                                  is applied whether or not extraction runs:
                                  with -c/--center it feeds the extractor charge
                                  summary; with -c omitted (extraction skipped)
                                  the same mapping is applied to the full input
                                  PDB/mmCIF to derive the total system charge. A
                                  bare number sets the total directly. PDB/mmCIF
                                  inputs only. An explicit -q/--charge has
                                  highest priority over either derived value and
                                  emits a warning.
  -q, --charge INTEGER            Total system charge. This explicit value has
                                  highest priority over extractor/workflow-
                                  derived charge; a mismatch emits a warning.
                                  Omit it to use automatic charge derivation.
  --workers INTEGER               MLIP predictor workers; >1 spawns a parallel
                                  predictor. NOTE: with UMA, workers>1 plus an
                                  explicit Analytical Hessian request is an
                                  error; use workers=1 or FiniteDifference.
                                  [default: 1]
  --workers-per-node INTEGER      Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Experimental, computationally expensive xTB
                                  solvent delta correction. Examples: water,
                                  methanol, acetonitrile, dmso, thf, toluene.
                                  'none' disables it.  [default: none]
  --solvent-model [alpb|cpcmx]    xTB solvent model.  [default: alpb]
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1).  [default: 1]
  --freeze-links BOOLEAN          Freeze parent atoms of cap hydrogens
                                  (PDB/mmCIF input or XYZ/GJF with --ref-pdb).
                                  [default: True]
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  in every stage (e.g., '1,3,5'). With
                                  extraction, indices refer to the original full
                                  input and are mapped to the active site model;
                                  otherwise they refer to the current input.
                                  Merged with --freeze-links and YAML
                                  geom.freeze_atoms.
  --mep-mode [gsm|dmf]            MEP optimizer: Growing String Method (gsm) or
                                  Direct Max Flux (dmf).  [default: gsm]
  --dmf-backend [cpu|gpu]         DMF compute backend (--mep-mode dmf only): gpu
                                  (dmf.torch / CUDA) or cpu (dmf / NumPy). On a
                                  GPU out-of-memory error, retry with cpu.
                                  [default: gpu]
  --max-nodes INTEGER             Movable internal images per GSM/DMF segment;
                                  the complete segment has max_nodes+2 images
                                  including endpoints.  [default: 20]
  --gsm-param [equi|energy]       GSM node parameterization after string growth.
                                  The energy scheme concentrates nodes in high-
                                  energy regions and may be tried when an
                                  equidistant path skips the reaction-coordinate
                                  region near the HEI.  [default: (equi)]
  --max-cycles-gsm INTEGER RANGE  Maximum GSM string-optimizer cycles for the
                                  MEP stage.  [default: (300); x>=1]
  --max-cycles-dmf INTEGER RANGE  Maximum IPOPT iterations for the DMF MEP
                                  stage. This is a solver iteration count, not a
                                  string-optimizer cycle count.  [default:
                                  (300); x>=1]
  --climb BOOLEAN                 Enable climbing image for standard GSM
                                  segments (bridge segments always disable
                                  climbing).  [default: True]
  --opt-mode [grad|hess]          Optimizer mode forwarded to scan/tsopt and
                                  used for single optimizations: grad
                                  (=LBFGS/Dimer) or hess (=RFO for scan/opt; RS-
                                  P-RFO for tsopt).  [default: grad]
  --opt-mode-post [grad|hess]     Optimizer mode override for TSOPT/post-IRC
                                  endpoint optimizations. If unset, uses --opt-
                                  mode when explicitly provided; otherwise falls
                                  back to the default ('hess' = RS-P-RFO).
                                  [default: hess]
  --dump BOOLEAN                  Dump GSM/MEP trajectories. An explicit parent
                                  toggle is forwarded to path-search/path-opt
                                  and scan/tsopt; when omitted, child
                                  YAML/defaults apply. When --thermo is enabled,
                                  freq always retains thermoanalysis.yaml
                                  because the composite workflow consumes that
                                  file; --no-dump does not suppress it.
                                  [default: False]
  --convert-files BOOLEAN         Convert XYZ/TRJ outputs into PDB/CIF/GJF
                                  companions based on the input format.
                                  [default: True]
  --refine-path BOOLEAN           Run a single-pass path-opt GSM between each
                                  adjacent pair and concatenate the segments
                                  (default; no path_search). Use --refine-path
                                  to run recursive path_search on the full
                                  ordered series for automatic multistep
                                  discovery.  [default: False]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for single-structure
                                  optimizations and scan relaxations (gau_loose|
                                  gau|gau_tight|gau_vtight|baker|never). The MEP
                                  stage keeps its own --thresh-gsm / --thresh-
                                  dmf.  [default: (gau)]
  --thresh-post [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for post-IRC endpoint
                                  optimizations (gau_loose|gau|gau_tight|gau_vti
                                  ght|baker|never).  [default: baker]
  --thresh-gsm [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for the GSM string
                                  optimizer of the MEP stage (gau_loose|gau|gau_
                                  tight|gau_vtight|baker|never).  [default:
                                  (gau_loose)]
  --thresh-dmf TEXT               IPOPT dual-infeasibility tolerance for the DMF
                                  MEP stage: tight (0.04) | middle (0.10) |
                                  loose (0.20) or a positive float. This is not
                                  a Gaussian preset.  [default: (tight)]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan.
                                  With --center, runs extraction in a temporary
                                  directory to validate derived charge and
                                  electron parity; no computational stage or
                                  persistent output is produced.  [default: no-
                                  dry-run]
  --preopt BOOLEAN                If True, run initial single-structure
                                  optimizations of the active site model inputs.
                                  [default: True]
  --hessian-calc-mode [finitedifference|analytical]
                                  Common MLIP Hessian calculation mode forwarded
                                  to tsopt and freq.  [default:
                                  (FiniteDifference)]
  --tsopt BOOLEAN                 TS optimization + IRC per reactive segment (or
                                  TSOPT-only mode for single-structure), and
                                  build energy diagrams.  [default: False]
  --tsopt-from-mep-tan / --no-tsopt-from-mep-tan
                                  Guide Hessian-based TS root identity from MEP
                                  tangent candidate(s) at the highest-energy
                                  image. The CPU/file cache is not created or
                                  used when disabled. Dimer does not consume
                                  this Hessian reference mode.  [default: tsopt-
                                  from-mep-tan]
  --thermo BOOLEAN                Run freq on (R, TS, P) per reactive segment
                                  (or TSOPT-only mode) and build Gibbs free-
                                  energy diagram (MLIP).  [default: False]
  --dft BOOLEAN                   Run DFT single-point on (R, TS, P) and build
                                  DFT energy diagram. With --thermo, also
                                  generate a DFT//MLIP Gibbs diagram.  [default:
                                  False]
  --tsopt-max-cycles INTEGER RANGE
                                  Override tsopt --max-cycles.  [default:
                                  (100000); x>=1]
  --tsopt-out-dir DIRECTORY       Override tsopt output subdirectory (relative
                                  paths are resolved against the default).
                                  [default: (<segment>/ts)]
  --flatten / --no-flatten        Enable the extra-imaginary-mode flattening
                                  loop in tsopt (grad: dimer loop, hess: post-
                                  RS-P-RFO); --no-flatten forces
                                  flatten_max_iter=0.  [default: no-flatten]
  --reject-uphill / --no-reject-uphill
                                  Opt in to rejecting uphill RFO trials during
                                  post-IRC endpoint re-optimization only
                                  (tolerance: 1e-4 Hartree). Does not affect TS
                                  optimization or path search.  [default: no-
                                  reject-uphill]
  --stop-plateau / --no-stop-plateau
                                  Stop when the energy stops changing while the
                                  convergence criteria are still unmet, and
                                  report the run as stalled. It never signals
                                  convergence; each stage's own cycle limit
                                  remains the real bound.  [default: no-stop-
                                  plateau]
  --stop-plateau-thresh FLOAT     Energy range (hartree) below which --stop-
                                  plateau treats the window as flat.  [default:
                                  (1e-4)]
  --stop-plateau-window INTEGER   Number of consecutive cycles --stop-plateau
                                  inspects.  [default: (50)]
  --irc-step-size FLOAT           Override IRC --step-size (Bohr). If an IRC
                                  stops after only a few frames, retry with a
                                  smaller value such as 0.05.  [default: (0.10)]
  --irc-max-cycles INTEGER RANGE  Cycle cap for each post-TS IRC.  [default:
                                  (125); x>=1]
  --irc-never-stop / --no-irc-never-stop
                                  Ignore IRC RMS-gradient, hard-gradient,
                                  energy-rise, and energy-change stops and trace
                                  until the IRC max-cycle limit.
                                  Numerical/integration failures and external
                                  interruption still stop the run.  [default:
                                  (no-irc-never-stop)]
  --freq-out-dir DIRECTORY        Override freq output base directory (relative
                                  paths resolved against the default).
                                  [default: (<tsopt dir>/freq)]
  --freq-max-write INTEGER        Override freq --max-write value.  [default:
                                  (10)]
  --freq-amplitude-ang FLOAT      Override freq --amplitude-ang (Å).  [default:
                                  (0.8)]
  --freq-n-frames INTEGER         Override freq --n-frames value.  [default:
                                  (20)]
  --freq-sort [value|abs]         Override freq mode sorting.  [default:
                                  (value)]
  --freq-temperature FLOAT        Override freq thermochemistry temperature (K).
                                  [default: (298.15)]
  --freq-pressure FLOAT           Override freq thermochemistry pressure (atm).
                                  [default: (1.0)]
  --dft-out-dir DIRECTORY         Override dft output base directory (relative
                                  paths resolved against the default).
                                  [default: (<tsopt dir>/dft)]
  --dft-func-basis TEXT           Override dft --func-basis value.  [default:
                                  (wb97m-v/def2-tzvpd)]
  --dft-max-cycle INTEGER RANGE   Override dft --max-cycle value.  [default:
                                  (100); x>=1]
  --dft-conv-tol FLOAT            Override dft --conv-tol value.  [default:
                                  (1e-9)]
  --dft-grid-level INTEGER        Override dft --grid-level value.  [default:
                                  (3)]
  --dft-engine [gpu|cpu]          Override the DFT backend (gpu or cpu); omitted
                                  values inherit YAML/defaults.  [default:
                                  (gpu)]
  -s, --scan-lists TEXT           Scan targets: inline Python literal. Multiple
                                  inline literals define sequential stages, e.g.
                                  '[(12,45,1.35)]'
                                  '[(10,55,2.20),(23,34,1.80)]'. Indices refer
                                  to the original full input ordering (1-based);
                                  atom strings may use
                                  CHAIN:RESNAME:RESSEQ[ICODE]:ATOM. When
                                  extraction is used, selections are auto-mapped
                                  to the active site model after extraction.
  --scan-out-dir DIRECTORY        Override the scan output directory. Relative
                                  paths are resolved against the default parent.
                                  [default: (<out-dir>/_work/scan)]
  --scan-one-based BOOLEAN        Override the scan subcommand indexing
                                  interpretation (True = 1-based, False =
                                  0-based).  [default: (True)]
  --scan-max-step-size FLOAT      Override scan --max-step-size (Å).  [default:
                                  (0.20)]
  --scan-bias-k FLOAT             Override scan harmonic bias strength k
                                  (eV/Å^2).  [default: (300)]
  --scan-relax-max-cycles INTEGER RANGE
                                  Override scan relaxation max cycles per step.
                                  [default: (100000); x>=1]
  --scan-preopt BOOLEAN           Override scan --preopt flag. Inherits from
                                  --preopt when omitted.  [default: (inherits
                                  --preopt)]
  --scan-endopt BOOLEAN           Override scan --endopt flag.  [default:
                                  (False)]
  --ref-pdb FILE                  Reference PDB/mmCIF for topology when -i
                                  provides XYZ inputs. Enables topology-aware
                                  PDB/mmCIF output conversion in TSOPT-only,
                                  scan, and path_search pipelines.
  --coord-type [cart|dlc]         Optimization coordinate system (cart|dlc).
                                  [default: (cart)]
  --print-every INTEGER RANGE     Print optimizer status every N cycles.
                                  [default: (100); x>=1]
  --precision [fp32|fp64]         MLIP backend precision: fp32 or fp64. Unset
                                  defaults per backend (uma: fp32; orb, mace:
                                  fp64). Routed to backend-specific kwargs (UMA
                                  precision / ORB precision / MACE
                                  default_dtype). aimnet2: fp32 no-op; fp64
                                  rejected.  [default: (per backend: uma fp32;
                                  orb, mace fp64)]
  --backend-model TEXT            Model variant for the selected --backend (e.g.
                                  uma-s-1p2 / uma-m-1p1 for uma,
                                  orb_v3_conservative_omol for orb, MACE-OMOL-0
                                  / off:small for mace).  [default: (the
                                  selected backend's own model)]
  --calc-file FILE                Python file exposing get_calculator(...) -> an
                                  ASE Calculator, used as the energy/gradient
                                  backend (overrides --backend). Couples GFN-xTB
                                  / DFTB+ / any ASE engine. See --calc-file-
                                  func-name.
  --calc-file-func-name TEXT      Name of the callable in --calc-file that
                                  returns an ASE Calculator (or a module-level
                                  Calculator instance). CLI overrides config
                                  YAML; otherwise defaults to get_calculator.
                                  [default: (get_calculator)]
  --deterministic / --no-deterministic
                                  Request strict same-stack PyTorch determinism
                                  (deterministic algorithms + index_reduce_
                                  shim). Slower; raises for detected unsupported
                                  ops; custom calculators are outside its scope.
                                  [default: no-deterministic]
  --allow-charge-mult-mismatch    Skip the cluster charge/multiplicity electron-
                                  parity check (logs that it was skipped). Open-
                                  shell clusters need a matching multiplicity
                                  instead; use this only for an intentionally
                                  nonstandard electron count.
  -h, --help                      Show this message and exit.
```
