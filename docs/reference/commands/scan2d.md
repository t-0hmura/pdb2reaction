# `pdb2reaction scan2d`

```text
Usage: pdb2reaction scan2d [OPTIONS]

  2D distance scan with harmonic restraints.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure file (.pdb, .cif, .mmcif,
                                  .xyz, _trj.xyz, ...).  [required]
  -s, --scan-lists TEXT           Scan targets: inline Python literal or a
                                  YAML/JSON spec file path. scan2d expects
                                  EXACTLY 2 quadruples (i, j, low, high) — one
                                  per scanned bond axis — e.g.
                                  '[(12,45,1.30,3.10),(10,55,1.20,3.20)]'. Atom
                                  indices may also be strings like 'CE SAM 216';
                                  use positional
                                  CHAIN:RESNAME:RESSEQ[ICODE]:ATOM when chain
                                  qualification is needed. Step count per axis
                                  is set via --max-step-size, NOT inside the
                                  tuple (scan2d does not accept a 5th element;
                                  the scan command's 3-tuple (i, j, target) form
                                  is rejected here too).  [required]
  -q, --charge INTEGER            Total charge. Required for non-.gjf inputs
                                  unless --ligand-charge is provided (PDB/mmCIF
                                  inputs or XYZ/GJF with --ref-pdb).
  --workers INTEGER               MLIP predictor workers; >1 spawns a parallel
                                  predictor. NOTE: with UMA, workers>1 plus an
                                  explicit Analytical Hessian request is an
                                  error; use workers=1 or FiniteDifference.
                                  [default: 1]
  --workers-per-node INTEGER      Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB/mmCIF input or --ref-
                                  pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1).  [default: (1)]
  --one-based / --zero-based      Interpret (i,j) indices in --scan-lists as
                                  1-based or 0-based.  [default: one-based]
  --max-step-size FLOAT           Maximum step size per scanned distance [Å].
                                  [default: 0.2]
  --bias-k FLOAT                  Harmonic well strength k [eV/Å^2]. YAML bias.k
                                  applies when this option is omitted; explicit
                                  CLI wins.  [default: (300.0)]
  --relax-max-cycles INTEGER RANGE
                                  Maximum optimizer cycles per grid relaxation.
                                  An explicitly provided value overrides YAML
                                  opt.max_cycles.  [default: (100000); x>=1]
  --opt-mode [grad|hess]          Relaxation mode: grad (=LBFGS) or hess (=RFO).
                                  [default: grad]
  --freeze-links / --no-freeze-links
                                  Freeze parent atoms of cap hydrogens
                                  (PDB/mmCIF input or XYZ/GJF with --ref-pdb).
                                  [default: freeze-links]
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., '1,3,5').
  --dump / --no-dump              Write inner scan trajectories per d1-step as
                                  TRJ under result_scan2d/grid/.  [default: no-
                                  dump]
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB/CIF/GJF
                                  companions based on the input
                                  topology/template.  [default: convert-files]
  --ref-pdb FILE                  Reference PDB/mmCIF topology to use when the
                                  input is XYZ/GJF (keeps XYZ coordinates).
  -o, --out-dir TEXT              Base output directory.  [default:
                                  ./result_scan2d/]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset (gau_loose|gau|gau_tight|ga
                                  u_vtight|baker|never).  [default: (baker)]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --preopt / --no-preopt          Pre-optimize the initial structure without
                                  bias before the scan.  [default: no-preopt]
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Experimental, computationally expensive xTB
                                  solvent delta correction. Examples: water,
                                  methanol, acetonitrile, dmso, thf, toluene.
                                  'none' disables it.  [default: none]
  --solvent-model [alpb|cpcmx]    xTB solvent model.  [default: alpb]
  --baseline [min|first]          Reference for relative energy (kcal/mol):
                                  'min' or 'first' (i=0,j=0).  [default: min]
  --zmin FLOAT                    Lower bound of color scale for plots
                                  (kcal/mol).  [default: (the surface minimum)]
  --zmax FLOAT                    Upper bound of color scale for plots
                                  (kcal/mol).  [default: (the surface maximum)]
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
                                  (cart|redund|dlc|tric).  [default: (cart)]
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
