# `pdb2reaction all`

```text
Usage: pdb2reaction all [OPTIONS]

  Run active site model extraction → optional staged scan → MEP search in one
  command. If exactly one input is provided: (a) with --scan-lists, run staged
  scan on the active site model (or full structure when extraction is skipped)
  and use stage results as inputs for the MEP search; (b) with --tsopt and no
  --scan-lists, run TSOPT-only mode.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, output paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Two or more full structures
                                  (PDB/mmCIF/XYZ/GJF) in reaction order
                                  (reactant [intermediates ...] product), or a
                                  single full structure (with --scan-lists or
                                  with --tsopt). Extraction (-c/--center)
                                  accepts PDB/mmCIF. mmCIF and oversized PDB
                                  inputs give mmCIF outputs. Pass multiple files
                                  after one -i/--input, or repeat -i/--input for
                                  each file.  [required]
  -c, --center TEXT               Extraction centers (normally substrate +
                                  catalytic residues): a PDB/mmCIF path, a
                                  residue-ID list like '123,124' or
                                  'A:123,B:456' (insertion codes OK: '123A' /
                                  'A:123A'), a residue-name list like 'GPP,SAM',
                                  or chain-qualified 'A:SAM' / 'A:SAM:123'. Each
                                  match starts radius expansion. When omitted,
                                  extraction is skipped and the full input
                                  structure(s) are used directly as active site
                                  models.
  -o, --out-dir DIRECTORY         Top-level output directory.  [default:
                                  result_all]
  -r, --radius FLOAT RANGE        Inclusion cutoff (Å) around center atoms. 0
                                  disables radius-based expansion.  [default:
                                  2.6; x>=0.0]
  --radius-het2het FLOAT RANGE    Independent hetero-hetero cutoff (Å) for
                                  non-C/H pairs.  [default: 0.0; x>=0.0]
  --include-h2o / --no-include-h2o
                                  Include waters (HOH/WAT/H2O/DOD/TIP/TIP3/SOL)
                                  in the active site model.  [default:
                                  include-h2o]
  --exclude-backbone / --no-exclude-backbone
                                  Delete main-chain atoms from amino acids; only
                                  the main chain between peptide-bonded
                                  extraction centers is kept.  [default: no-
                                  exclude-backbone]
  --add-linkh / --no-add-linkh    Add cap hydrogens for severed bonds (carbon
                                  boundaries only) in active site models.
                                  [default: add-linkh]
  --selected-resn TEXT            Force-include residues without radius
                                  expansion using the same selectors as
                                  -c/--center: IDs ('123', 'A:123A'), names
                                  ('SAM'), or chain-qualified names ('A:SAM',
                                  'A:SAM:123'); comma/space separated.
                                  [default: ""]
  --modified-residue TEXT         Comma-separated residue names to treat as
                                  amino acids for backbone truncation and charge
                                  assignment. NAME:charge adds or overrides the
                                  nominal charge for this extraction; bare NAME
                                  defaults to 0. Examples: 'HD1,HD2,HD3' or
                                  'HD1:0,SEP:-2'.  [default: ""]
  -l, --ligand-charge TEXT        Total ligand charge (number) or per-resname
                                  mapping like 'GPP:-3,SAM:1'. The per-resname
                                  mapping is applied whether or not extraction
                                  runs: with -c/--center it feeds the extractor
                                  charge summary; with -c omitted (extraction
                                  skipped) the same mapping is applied to the
                                  full input PDB/mmCIF to derive the total
                                  system charge. A bare number specifies the
                                  total ligand charge; use -q/--charge to set
                                  the total system charge. PDB/mmCIF inputs
                                  only. An explicit -q/--charge has highest
                                  priority over either derived value and emits a
                                  warning.
  -q, --charge INTEGER            Total system charge. Overrides the derived
                                  charge; a mismatch emits a warning. Omit it to
                                  use automatic charge derivation.
  --uma-workers, --workers INTEGER
                                  MLIP predictor workers; >1 spawns a parallel
                                  predictor (needs fairchem-core[ray]). With
                                  UMA, workers>1 plus an explicit Analytical
                                  Hessian request is an error; use workers=1 or
                                  FiniteDifference.  [default: 1]
  --uma-workers-per-node, --workers-per-node INTEGER
                                  Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  -b, --backend [uma|orb|mace|aimnet2|dft]
                                  Energy/force calculator backend.  [default:
                                  uma]
  --solvent TEXT                  Computationally expensive xTB solvent delta
                                  correction for MLIP backends; dft uses native
                                  PySCF PCM/SMD. Examples: water, methanol,
                                  acetonitrile, dmso, thf, toluene. 'none'
                                  disables it.  [default: none]
  --solvent-model [alpb|cpcmx|pcm|smd]
                                  Solvent model: ALPB/CPCMx for MLIP backends;
                                  PCM/SMD for dft.  [default: alpb]
  --solvent-xtb-cmd TEXT          Command for the xTB solvent correction,
                                  including optional xTB arguments. If xTB SCC
                                  convergence is poor, increasing --etemp may
                                  help (for example: 'xtb --etemp 1000').
  -m, --multiplicity INTEGER RANGE
                                  Spin multiplicity (2S+1).  [default: 1; x>=1]
  --freeze-links / --no-freeze-links
                                  Freeze parent atoms of cap hydrogens
                                  (PDB/mmCIF input or XYZ/GJF with --ref-pdb).
                                  [default: freeze-links]
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
  --max-depth INTEGER RANGE       Maximum recursive subdivision levels; requires
                                  --refine-path. 0 disables subdivision. With a
                                  positive limit, an interval still undivided at
                                  the limit is named seg_NNN_maxdepth and may
                                  contain several steps.  [default: (10); x>=0]
  --gsm-param [equi|energy]       GSM node parameterization after string growth.
                                  The energy scheme concentrates nodes in high-
                                  energy regions and may be tried when an
                                  equidistant path skips the reaction-coordinate
                                  region near the HEI.  [default: (equi)]
  --max-cycles-gsm INTEGER RANGE  Maximum GSM string-optimizer cycles for the
                                  MEP stage.  [default: (300); x>=1]
  --dmf-max-iterations, --max-cycles-dmf INTEGER RANGE
                                  Maximum IPOPT iterations for the DMF MEP
                                  stage.  [default: (3000); x>=1]
  --climb / --no-climb            Enable climbing image for standard GSM
                                  segments (bridge segments always disable
                                  climbing).  [default: climb]
  --opt-mode [grad|hess]          Optimizer mode forwarded to scan/tsopt and
                                  used for single optimizations: grad
                                  (=LBFGS/Dimer) or hess (=RFO for scan/opt; RS-
                                  P-RFO for tsopt).  [default: grad]
  --opt-mode-post [grad|hess]     Optimizer mode override for TSOPT/post-IRC
                                  endpoint optimizations. If unset, uses --opt-
                                  mode when explicitly provided; otherwise falls
                                  back to the default ('hess' = RS-P-RFO).
                                  [default: hess]
  --dump / --no-dump              Dump trajectories from the MEP, scan and tsopt
                                  stages; when omitted, each stage's
                                  YAML/defaults apply. With --thermo, freq
                                  always keeps thermoanalysis.yaml.  [default:
                                  no-dump]
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB/CIF/GJF
                                  companions based on the input format.
                                  [default: convert-files]
  --refine-path / --no-refine-path
                                  By default, run a single-pass path-opt with
                                  the selected GSM/DMF method between each
                                  adjacent pair and concatenate the segments.
                                  Use --refine-path to run recursive path-search
                                  on the full ordered series: it proposes
                                  multistep reaction paths and also refines a
                                  single-step MEP, which can improve a poor HEI
                                  or TS estimate.  [default: no-refine-path]
  --write-ref-merge / --no-write-ref-merge
                                  Write the MEP and HEI merged into the full
                                  input (mep_w_ref/hei_w_ref) in --refine-path
                                  mode. Requires -c/--center and PDB/mmCIF
                                  input.  [default: no-write-ref-merge]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for single-structure
                                  optimizations and scan relaxations. The MEP
                                  stage keeps its own --thresh-gsm / --dmf-tol.
                                  [default: (gau)]
  --thresh-post [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for TS and post-IRC
                                  endpoint optimizations.  [default: baker]
  --thresh-gsm [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for the GSM string
                                  optimizer of the MEP stage.  [default:
                                  (gau_loose)]
  --dmf-tol, --thresh-dmf TEXT    IPOPT dual-infeasibility tolerance for the DMF
                                  MEP stage: tight (0.04) | middle (0.10) |
                                  loose (0.20) or a positive float.  [default:
                                  (tight)]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --resume / --no-resume          Resume from the first unfinished segment using
                                  the verified MEP in --out-dir.  [default: no-
                                  resume]
  --resume-segment INTEGER RANGE  Reuse the verified MEP in --out-dir and rerun
                                  post-processing from segment N. Repeat the
                                  original path/extraction/calculator options;
                                  post-processing options other than
                                  --precision-post may be changed.  [default:
                                  (disabled); x>=1]
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan.
                                  With --center, runs extraction in a temporary
                                  directory to validate derived charge and
                                  electron parity; no computational stage or
                                  persistent output is produced.  [default: no-
                                  dry-run]
  --preopt / --no-preopt          Run initial single-structure optimizations of
                                  the active site model inputs.  [default:
                                  preopt]
  --hessian-calc-mode [finitedifference|analytical]
                                  How the backend computes the Hessian in tsopt,
                                  irc and freq.  [default: (FiniteDifference)]
  --tsopt / --no-tsopt            Run TS optimization and IRC per reactive
                                  segment (or TSOPT-only mode for a single
                                  structure) and build energy diagrams.
                                  [default: no-tsopt]
  --tsopt-from-mep-tan / --no-tsopt-from-mep-tan
                                  Use the MEP tangent at the highest-energy
                                  image to choose which Hessian mode RS-P-RFO
                                  (hess) follows as the reaction mode. Not used
                                  with Dimer (grad).  [default: tsopt-from-mep-
                                  tan]
  --thermo / --no-thermo          Run freq on (R, TS, P) per reactive segment
                                  (or TSOPT-only mode) and build Gibbs free-
                                  energy diagram (MLIP).  [default: no-thermo]
  --dft / --no-dft                Run DFT single-point on (R, TS, P) and build
                                  DFT energy diagram. With --thermo, also
                                  generate a DFT//MLIP Gibbs diagram.  [default:
                                  no-dft]
  --tsopt-max-cycles INTEGER RANGE
                                  Override tsopt --max-cycles.  [default:
                                  (100000); x>=1]
  --tsopt-out-dir DIRECTORY       Override tsopt output subdirectory (relative
                                  paths are resolved against the default).
                                  [default: (<segment>/ts)]
  --flatten / --no-flatten        Enable the extra-imaginary-mode flattening
                                  loop in tsopt; --no-flatten forces
                                  flatten_max_iter=0.  [default: no-flatten]
  --skip-final-freq / --no-skip-final-freq
                                  Skip the final PHVA/frequency analysis in
                                  tsopt. The TS structure is retained with
                                  unverified saddle order, and IRC and the later
                                  stages of that segment are skipped; other
                                  segments continue.  [default: no-skip-final-
                                  freq]
  --reject-uphill / --no-reject-uphill
                                  Reject uphill RFO steps during post-IRC
                                  endpoint re-optimization only (tolerance: 1e-4
                                  Hartree). Does not affect TS optimization or
                                  path search.  [default: no-reject-uphill]
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
                                  [default: (<segment dir>/freq)]
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
                                  [default: (<segment dir>/dft)]
  --func-basis, --dft-func-basis TEXT
                                  Override dft --func-basis value.  [default:
                                  (wb97m-v/def2-svp)]
  --scf-max-cycles, --dft-max-cycle INTEGER RANGE
                                  Override the DFT SCF iteration limit.
                                  [default: (100); x>=1]
  --scf-tol, --dft-conv-tol FLOAT
                                  Override the DFT SCF convergence tolerance.
                                  [default: (1e-9)]
  --dft-grid-level INTEGER        Override dft --grid-level value.  [default:
                                  (3)]
  --dft-engine, --engine [gpu|cpu]
                                  Override the DFT backend (gpu or cpu); omitted
                                  values inherit YAML/defaults.  [default:
                                  (gpu)]
  --dft-solvent TEXT              Implicit solvent for post-processing DFT
                                  single points.  [default: (none)]
  --dft-solvent-model [pcm|smd]   Native PySCF solvent model for post-processing
                                  DFT single points.  [default: (smd)]
  -s, --scan-lists TEXT           Scan targets: distance (i,j,target), angle
                                  (i,j,k,target), or dihedral (i,j,k,l,target).
                                  Multiple inline literals define sequential
                                  stages. Distances use Å; angles and dihedrals
                                  use degrees. Indices refer to the original
                                  full input ordering (1-based); atom strings
                                  may use CHAIN:RESNAME:RESSEQ[ICODE]:ATOM. When
                                  extraction is used, selections are auto-mapped
                                  to the active site model after extraction.
  --scan-out-dir DIRECTORY        Override the scan output directory. Relative
                                  paths are resolved against the default parent.
                                  [default: (<out-dir>/_work/scan)]
  --scan-one-based / --scan-zero-based
                                  Override the scan subcommand indexing
                                  interpretation (True = 1-based, False =
                                  0-based).  [default: (True)]
  --scan-max-step-size FLOAT      Override scan --max-step-size (Å).  [default:
                                  (0.20)]
  --scan-restraint-k, --scan-bias-k FLOAT
                                  Override scan harmonic bias strength k
                                  (eV/Å^2).  [default: (300)]
  --scan-relax-max-cycles INTEGER RANGE
                                  Override scan relaxation max cycles per step.
                                  [default: (100000); x>=1]
  --scan-preopt / --no-scan-preopt
                                  Override scan --preopt flag.  [default:
                                  (inherits --preopt)]
  --scan-endopt / --no-scan-endopt
                                  Override scan --endopt flag.  [default:
                                  (False)]
  --ref-pdb FILE                  Reference PDB/mmCIF for topology when -i
                                  provides XYZ/GJF inputs. Enables PDB/mmCIF
                                  output conversion.
  --coord-type [cart|dlc]         Optimization coordinate system.  [default:
                                  (cart)]
  --print-every INTEGER RANGE     Print optimizer status every N cycles.
                                  [default: (100); x>=1]
  --precision [tf32|fp32|fp64]    MLIP backend precision: tf32, fp32 (IEEE) or
                                  fp64. tf32 uses TF32 matmul for speed (uma,
                                  orb only); avoid it for Hessians. aimnet2:
                                  fp32 only. --precision-post overrides it for
                                  TSOPT, IRC, endpoint optimization and freq.
                                  [default: (per backend: uma fp32; orb, mace
                                  fp64)]
  --precision-post [tf32|fp32|fp64]
                                  MLIP precision for TSOPT, IRC, the post-IRC
                                  endpoint optimizations and freq. Unset uses
                                  --precision.  [default: (inherits
                                  --precision)]
  --compile / --no-compile        torch.compile the UMA or ORB model: slower
                                  first call, slightly faster afterwards; worth
                                  it only for long runs.  [default: (off)]
  --backend-model TEXT            Model variant for the selected --backend (e.g.
                                  uma-s-1p2p1 / uma-s-1p2p1-ft-d4rev / uma-m-1p1
                                  for uma, orbmol-v2 / orb_v3_conservative_omol
                                  for orb, MACE-OMOL-0 / MACE-POLAR-1-M for
                                  mace, aimnet2 / aimnet2-nse for aimnet2).
                                  [default: (the selected backend's own model)]
  --calc-file FILE                Python file exposing get_calculator(...) -> an
                                  ASE Calculator, used as the energy/gradient
                                  backend (overrides --backend). Works with GFN-
                                  xTB, DFTB+ or any ASE engine. See --calc-file-
                                  func-name.
  --calc-factory, --calc-file-func-name TEXT
                                  Name of the callable in --calc-file that
                                  returns an ASE Calculator (or a module-level
                                  Calculator instance).  [default:
                                  (get_calculator)]
  --deterministic / --no-deterministic
                                  Make repeated runs on the same software
                                  versions and hardware give identical results
                                  (PyTorch deterministic algorithms). Slower;
                                  raises an error on unsupported ops; custom
                                  calculators are outside its scope.  [default:
                                  no-deterministic]
  --allow-charge-mult-mismatch    Skip the cluster charge/multiplicity electron-
                                  parity check (logs that it was skipped). Open-
                                  shell clusters need a matching multiplicity
                                  instead; use this only for an intentionally
                                  nonstandard electron count.
  --save-scf-checkpoint / --no-save-scf-checkpoint
                                  Save a structure-bound PySCF checkpoint.
                                  [default: (disabled)]
  --scf-checkpoint FILE           Load/save the optional structure-bound PySCF
                                  checkpoint.
  --scf-stepwise-grid / --no-scf-stepwise-grid
                                  Converge the first SCF on a coarse grid, then
                                  on the final grid (later SCFs reuse the
                                  previous density as usual).  [default:
                                  (enabled)]
  --dft-low-memory, --lowmem / --no-dft-low-memory, --no-lowmem
                                  Use GPU4PySCF rks_lowmem for closed-shell GPU
                                  DFT; open-shell GPU and CPU use standard
                                  direct JK. --no-dft-low-memory enables density
                                  fitting.  [default: (lowmem)]
  --dft-nprocs INTEGER RANGE      PySCF/OpenMP CPU threads; GPU count is
                                  unaffected.  [default: (auto); x>=1]
  --dft-memory, --dft-mem TEXT    PySCF host RAM limit (for example 64GB or
                                  120000MB).  [default: (auto)]
  -h, --help                      Show this message and exit.
```
