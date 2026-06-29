# `pdb2reaction scan3d`

```text
Usage: pdb2reaction scan3d [OPTIONS]

  3D distance scan with harmonic restraints.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure file (.pdb, .xyz, _trj.xyz,
                                  ...). Required unless --csv is provided.
  -s, --scan-lists TEXT           Scan targets: inline Python literal or a
                                  YAML/JSON spec file path. scan3d expects
                                  EXACTLY 3 quadruples (i, j, low, high) — one
                                  per scanned bond axis — e.g. '[(12,45,1.30,3.1
                                  0),(10,55,1.20,3.20),(15,60,1.10,3.00)]'. Atom
                                  indices may also be PDB-style strings like 'CE
                                  SAM   216'. Step count per axis is set via
                                  --max-step-size, NOT inside the tuple (scan3d
                                  does not accept a 5th element).
  --csv FILE                      If provided, skip the 3D scan and read a
                                  precomputed surface.csv from this path. When
                                  used, -i/--input and --scan-lists are
                                  optional.
  -q, --charge INTEGER            Total charge. Required for non-.gjf inputs
                                  unless --ligand-charge is provided (PDB inputs
                                  or XYZ/GJF with --ref-pdb).
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
  --one-based / --zero-based      Interpret (i,j) indices in --scan-lists as
                                  1-based (default) or 0-based.  [default: one-
                                  based]
  --max-step-size FLOAT           Maximum step size per scanned distance [Å].
                                  [default: 0.2]
  --bias-k FLOAT                  Harmonic well strength k [eV/Å^2]. Defaults to
                                  YAML bias.k (BIAS_KW['k']=300 in defaults.py)
                                  when omitted; explicit CLI value overrides
                                  YAML.
  --relax-max-cycles INTEGER      Maximum optimizer cycles per grid relaxation.
                                  When explicitly provided, used unless YAML
                                  sets opt.max_cycles.  [default: 10000]
  --opt-mode [grad|hess]          Relaxation mode: grad (=LBFGS) or hess (=RFO).
                                  [default: grad]
  --freeze-links / --no-freeze-links
                                  Freeze parent atoms of cap hydrogens (PDB
                                  input or XYZ/GJF with --ref-pdb).  [default:
                                  freeze-links]
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., '1,3,5').
  --dump / --no-dump              Write inner d3 scan trajectories per (d1,d2)
                                  as TRJ under result_scan3d/grid/.  [default:
                                  no-dump]
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB/GJF
                                  companions based on the input format.
                                  [default: convert-files]
  --ref-pdb FILE                  Reference PDB topology to use when the input
                                  is XYZ/GJF (keeps XYZ coordinates).
  -o, --out-dir TEXT              Base output directory.  [default:
                                  ./result_scan3d/]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset (gau_loose|gau|gau_tight|ga
                                  u_vtight|baker|never).  Defaults to 'baker'.
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --preopt / --no-preopt          Pre-optimize the initial structure without
                                  bias before the scan.  [default: no-preopt]
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  --solvent-model [alpb|cpcmx]    xTB solvent model.  [default: alpb]
  --baseline [min|first]          Reference for relative energy (kcal/mol):
                                  'min' or 'first' (i=0,j=0,k=0).  [default:
                                  min]
  --zmin FLOAT                    Lower bound of color scale for plots
                                  (kcal/mol).
  --zmax FLOAT                    Upper bound of color scale for plots
                                  (kcal/mol).
  --print-parsed / --no-print-parsed
                                  Print parsed scan targets after resolving
                                  --scan-lists.  [default: no-print-parsed]
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
  --dry-run / --no-dry-run        Resolve and validate options (input,
                                  charge/spin parity, --scan-lists parse) and
                                  print the planned scan, then exit without
                                  running any optimization.  [default: no-dry-
                                  run]
  --coord-type [cart|redund|dlc|tric]
                                  Optimization coordinate system
                                  (cart|redund|dlc|tric). cart is the reliable
                                  default used for the published results; dlc
                                  speeds up torsion-rich optimizations.
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
