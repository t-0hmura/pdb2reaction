# `pdb2reaction opt`

```text
Usage: pdb2reaction opt [OPTIONS]

  Single-structure geometry optimization using LBFGS or RFO.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure file (.pdb, .xyz, _trj.xyz,
                                  ...).  [required]
  --workers INTEGER               MLIP predictor workers; >1 spawns a parallel
                                  predictor. NOTE: the analytical Hessian (used
                                  for RFO seeding) raises a RuntimeError when
                                  workers>1; run with --workers 1 for those
                                  modes.  [default: 1]
  --workers-per-node INTEGER      Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  --dist-freeze TEXT              Distance restraints: inline Python literal
                                  (e.g. '[(1,5,1.4)]') or a YAML/JSON spec file
                                  path. Same format as --scan-lists:
                                  (i,j,target_A) triples. Target may be omitted
                                  to freeze at the current distance: (i,j).
  --one-based / --zero-based      Interpret --dist-freeze / --scan-lists indices
                                  as 1-based (default) or 0-based.  [default:
                                  one-based]
  --bias-k FLOAT                  Harmonic restraint strength k [eV/Å^2] for
                                  --dist-freeze.  [default: 300]
  --freeze-links / --no-freeze-links
                                  Freeze parent atoms of cap hydrogens (PDB
                                  input or XYZ/GJF with --ref-pdb).  [default:
                                  freeze-links]
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., '1,3,5').
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB/GJF
                                  companions based on the input format.
                                  [default: convert-files]
  --ref-pdb FILE                  Reference PDB topology to use when the input
                                  is XYZ/GJF (keeps XYZ coordinates).
  --max-cycles INTEGER            Maximum number of optimization cycles.
                                  [default: 10000]
  --opt-mode [grad|hess|lbfgs|rfo]
                                  Optimization mode: grad (lbfgs) or hess (rfo).
                                  Aliases lbfgs/rfo are accepted.  [default:
                                  grad]
  --flatten / --no-flatten        Enable/disable imaginary-mode flatten loop
                                  after optimization.  [default: no-flatten]
  --dump / --no-dump              Write optimization trajectory to
                                  'optimization_trj.xyz'.  [default: no-dump]
  -o, --out-dir TEXT              Output directory.  [default: ./result_opt/]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset (gau_loose|gau|gau_tight|ga
                                  u_vtight|baker|never). Defaults to 'gau' when
                                  not provided.
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running optimization.  [default: no-
                                  dry-run]
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  --solvent-model [alpb|cpcmx]    xTB solvent model.  [default: alpb]
  --coord-type [cart|redund|dlc|tric]
                                  Optimization coordinate system
                                  (cart|redund|dlc|tric). cart is the reliable
                                  default used for the published results; dlc
                                  speeds up torsion-rich optimizations.
  -q, --charge INTEGER            Total charge. Required for non-.gjf inputs
                                  unless --ligand-charge is provided (.gjf
                                  templates inherit the charge automatically).
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1).
  --print-every INTEGER RANGE     Print optimizer status every N cycles (debug
                                  knob).  [x>=1]
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
