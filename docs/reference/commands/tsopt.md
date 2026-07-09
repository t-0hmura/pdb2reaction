# `pdb2reaction tsopt`

```text
Usage: pdb2reaction tsopt [OPTIONS]

  Transition state optimization (Dimer or RS-P-RFO).

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
                                  predictor. NOTE: when workers>1 the analytical
                                  Hessian is unavailable (the parallel predictor
                                  exposes no autograd model); it auto-downgrades
                                  to finite differences with a warning.
                                  [default: 1]
  --workers-per-node INTEGER      Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
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
  --flatten / --no-flatten        Enable the extra-imaginary-mode flattening
                                  loop (grad: dimer loop, hess: post-RS-P-RFO).
                                  [default: no-flatten]
  --opt-mode [grad|hess|dimer|rsirfo|trim|rsprfo]
                                  TS optimizer: 'grad'/'dimer' → Hessian Guided
                                  Dimer; 'hess'/'rsprfo' → RS-P-RFO (Banerjee,
                                  default); 'trim' → TRIM (Helgaker); 'rsirfo' →
                                  RS-I-RFO.  [default: hess]
  --dump / --no-dump              Write the per-cycle optimization trajectory
                                  ('optimization_trj.xyz' for rsirfo/hess,
                                  'optimization_all_trj.xyz' for grad/dimer) in
                                  the output directory.  [default: no-dump]
  -o, --out-dir TEXT              Output directory.  [default: ./result_tsopt/]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for the active optimizer (g
                                  au_loose|gau|gau_tight|gau_vtight|baker|never)
                                  . Defaults to 'baker' when not provided.
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running TS optimization.  [default:
                                  no-dry-run]
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
  --hessian-calc-mode [finitedifference|analytical]
                                  Choose MLIP Hessian evaluation mode (used
                                  unless YAML sets calc.hessian_calc_mode).
                                  Defaults to 'FiniteDifference'.
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  --solvent-model [alpb|cpcmx]    xTB solvent model.  [default: alpb]
  -q, --charge INTEGER            Total charge. Required for non-.gjf inputs
                                  unless --ligand-charge is provided (.gjf
                                  templates inherit the charge automatically).
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1).
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
  --coord-type [cart|redund|dlc|tric]
                                  Optimization coordinate system
                                  (cart|redund|dlc|tric). cart is the reliable
                                  default used for the published results; dlc
                                  speeds up torsion-rich optimizations.
  --print-every INTEGER RANGE     Print optimizer status every N cycles (debug
                                  knob).  [x>=1]
  -h, --help                      Show this message and exit.
```
