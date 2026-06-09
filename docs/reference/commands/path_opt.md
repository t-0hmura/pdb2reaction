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
  -q, --charge INTEGER            Total charge. Required unless a .gjf template
                                  provides charge metadata or --ligand-charge is
                                  supplied for PDB inputs.
  --workers INTEGER               MLIP predictor workers; >1 spawns a parallel
                                  predictor. NOTE: the analytical Hessian raises
                                  a RuntimeError when workers>1; run with
                                  --workers 1 for Hessian-based modes.
                                  [default: 1]
  --workers-per-node INTEGER      Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1).
  --freeze-links / --no-freeze-links
                                  Freeze parent atoms of link hydrogens (PDB
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
                                  cart is the robust default used in published
                                  numbers; dlc speeds up torsion-rich opts.
  --precision [fp32|fp64]         MLIP backend precision: fp32 (default) or
                                  fp64. Routed to backend-specific kwargs (UMA
                                  precision / ORB precision / MACE
                                  default_dtype). aimnet2: fp32 no-op; fp64
                                  rejected.
  --deterministic / --no-deterministic
                                  Strict bit-reproducible GPU runs
                                  (deterministic algorithms + index_reduce_
                                  shim). Slower; raises if unsupported. Default
                                  off.
  -h, --help                      Show this message and exit.
```
