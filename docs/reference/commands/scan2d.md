# `pdb2reaction scan2d`

```text
Usage: pdb2reaction scan2d [OPTIONS]

  2D internal-coordinate scan with harmonic restraints.

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
  -s, --scan-lists TEXT           Two scan ranges as an inline literal or
                                  YAML/JSON file: distance (i,j,low,high), angle
                                  (i,j,k,low,high), or dihedral
                                  (i,j,k,l,low,high). Atom indices may also be
                                  strings like 'CE SAM 216'; use positional
                                  CHAIN:RESNAME:RESSEQ[ICODE]:ATOM when chain
                                  qualification is needed. Distances use Å;
                                  angles and dihedrals use degrees.  [required]
  --solvent-xtb-cmd TEXT          Command for the xTB solvent correction,
                                  including optional xTB arguments. If xTB SCC
                                  convergence is poor, increasing --etemp may
                                  help (for example: 'xtb --etemp 1000').
  -q, --charge INTEGER            Total charge. Required for non-.gjf inputs
                                  unless --ligand-charge is provided (PDB/mmCIF
                                  inputs or XYZ/GJF with --ref-pdb).
  --uma-workers, --workers INTEGER
                                  MLIP predictor workers; >1 spawns a parallel
                                  predictor. NOTE: with UMA, workers>1 plus an
                                  explicit Analytical Hessian request is an
                                  error; use workers=1 or FiniteDifference.
                                  [default: 1]
  --uma-workers-per-node, --workers-per-node INTEGER
                                  Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB/mmCIF input or --ref-
                                  pdb).
  -m, --multiplicity INTEGER RANGE
                                  Spin multiplicity (2S+1).  [default: (1);
                                  x>=1]
  --one-based / --zero-based      Interpret atom indices in --scan-lists as
                                  1-based or 0-based.  [default: one-based]
  --max-step-size FLOAT           Maximum scanned distance change per step [Å].
                                  [default: 0.2]
  --max-angle-step-size FLOAT RANGE
                                  Maximum scanned angle change per step
                                  [degree].  [default: 5.0; x>0.0]
  --max-dihedral-step-size FLOAT RANGE
                                  Maximum scanned dihedral change per step
                                  [degree].  [default: 10.0; x>0.0]
  --restraint-k, --bias-k FLOAT   Harmonic well strength k [eV/Å^2 for
                                  distances; eV/rad^2 for angles]. YAML bias.k
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
  --baseline [min|first]          Reference for relative energy (kcal/mol):
                                  'min' or 'first' (i=0,j=0).  [default: min]
  --zmin FLOAT                    Lower bound of color scale for plots
                                  (kcal/mol).  [default: (the surface minimum)]
  --zmax FLOAT                    Upper bound of color scale for plots
                                  (kcal/mol).  [default: (the surface maximum)]
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
  --calc-factory, --calc-file-func-name TEXT
                                  Name of the callable in --calc-file that
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
  --func-basis TEXT               DFT method as FUNCTIONAL/BASIS; HF/BASIS is
                                  also accepted.  [default: (wb97m-v/def2-svp)]
  --dft-engine, --engine [gpu|cpu]
                                  PySCF execution engine used by --backend dft.
                                  [default: (gpu)]
  --save-scf-checkpoint / --no-save-scf-checkpoint
                                  Persist a structure-bound PySCF checkpoint
                                  (default: disabled).  [default: (disabled)]
  --scf-checkpoint FILE           Load/save the optional structure-bound PySCF
                                  checkpoint at PATH.
  --scf-stepwise-grid / --no-scf-stepwise-grid
                                  Converge the first SCF on a coarse grid, then
                                  on the final grid (later SCFs reuse the
                                  previous density as usual).  [default:
                                  (disabled)]
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
