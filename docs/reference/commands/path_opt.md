# `pdb2reaction path-opt`

```text
Usage: pdb2reaction path-opt [OPTIONS]

  MEP optimization via GSM or DMF.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE...             Two endpoint structures (reactant and
                                  product); accepts .pdb or .xyz.  [required]
  --mep-mode [gsm|dmf]            MEP optimizer: Growing String Method (gsm) or
                                  Direct Max Flux (dmf).  [default: gsm]
  --dmf-backend [cpu|gpu]         DMF compute backend (--mep-mode dmf only): gpu
                                  (dmf.torch / CUDA) or cpu (dmf / NumPy). On a
                                  GPU out-of-memory error, retry with cpu.
                                  [default: gpu]
  -q, --charge INTEGER            Total charge. Required unless a .gjf template
                                  provides charge metadata or --ligand-charge is
                                  supplied for PDB inputs.
  --workers INTEGER               MLIP predictor workers; >1 spawns a parallel
                                  predictor. NOTE: when workers>1 the analytical
                                  Hessian is unavailable and auto-downgrades to
                                  finite differences with a warning; run with
                                  --workers 1 to force analytical.  [default: 1]
  --workers-per-node INTEGER      Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1).
  --freeze-links / --no-freeze-links
                                  Freeze parent atoms of cap hydrogens (PDB
                                  input or XYZ/GJF with --ref-pdb).  [default:
                                  freeze-links]
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., '1,3,5').
  --max-nodes INTEGER             Maximum number of internal nodes (string has
                                  up to max_nodes+2 images including endpoints).
                                  [default: 20]
  --max-cycles INTEGER            Maximum string optimizer cycles (GSM/DMF path
                                  optimization).  [default: 300]
  --climb / --no-climb            Search for a transition state (climbing image)
                                  after path growth.  [default: climb]
  --opt-mode [grad|hess]          Single-structure optimizer for endpoint
                                  preoptimization: grad (=LBFGS) or hess (=RFO).
                                  [default: grad]
  --dump / --no-dump              Write optimizer trajectory and restarts during
                                  the run.  [default: no-dump]
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB/GJF
                                  companions based on the input format.
                                  [default: convert-files]
  --ref-pdb FILE                  Reference PDB topology to use when the input
                                  is XYZ/GJF (keeps XYZ coordinates).
  -o, --out-dir TEXT              Output directory.  [default:
                                  ./result_path_opt/]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for endpoint
                                  preoptimization only (gau_loose|gau|gau_tight|
                                  gau_vtight|baker|never). Defaults to 'gau'
                                  when not provided.
  --thresh-stopt TEXT             Convergence preset for the string optimizer
                                  (stopt) (gau_loose|gau|gau_tight|gau_vtight|ba
                                  ker|never). Defaults to 'gau_loose' when not
                                  provided.
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running path optimization.  [default:
                                  no-dry-run]
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
  --preopt / --no-preopt          Preoptimize each endpoint via the selected
                                  single-structure optimizer (LBFGS/RFO) before
                                  alignment and path optimization.  [default:
                                  preopt]
  --preopt-max-cycles INTEGER     Maximum cycles for each endpoint
                                  preoptimization pass (LBFGS or RFO; only used
                                  when --preopt is enabled).  [default: 10000]
  --fix-ends / --no-fix-ends      Fix structures of input endpoints during path
                                  optimization (GSM/DMF).  [default: fix-ends]
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  --solvent-model [alpb|cpcmx]    xTB solvent model.  [default: alpb]
  --coord-type [cart|dlc]         Optimization coordinate system (cart|dlc).
                                  cart is the reliable default used for the
                                  published results; dlc speeds up torsion-rich
                                  optimizations.
  --precision [fp32|fp64]         MLIP backend precision: fp32 or fp64. Unset
                                  defaults per backend (uma: fp32; orb, mace:
                                  fp64). Routed to backend-specific kwargs (UMA
                                  precision / ORB precision / MACE
                                  default_dtype). aimnet2: fp32 no-op; fp64
                                  rejected.
  --backend-model TEXT            Model variant for the selected --backend (e.g.
                                  uma-s-1p2 / uma-m-1p1 for uma,
                                  orb_v3_conservative_omol for orb, MACE-OMOL-0
                                  / MACE-OFF23_small for mace). Default: the
                                  backend's built-in model.
  --calc-file FILE                Python file exposing get_calculator(...) -> an
                                  ASE Calculator, used as the energy/gradient
                                  backend (overrides --backend). Couples GFN-xTB
                                  / DFTB+ / any ASE engine. See --calc-factory.
  --calc-factory TEXT             Name of the callable in --calc-file that
                                  returns an ASE Calculator (or a module-level
                                  Calculator instance).  [default:
                                  get_calculator]
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
