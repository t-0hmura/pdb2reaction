# `pdb2reaction freq`

```text
Usage: pdb2reaction freq [OPTIONS]

  Vibrational frequency analysis and mode writer (+ default thermochemistry
  summary).

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+detailed step logging
                                  and deliverable paths; 3=everything (full
                                  config blocks, per-file paths, DEBUG logging).
                                  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure file (.pdb, .cif, .mmcif,
                                  .xyz, .gjf, _trj.xyz, ...).  [required]
  --workers INTEGER               MLIP predictor workers; >1 spawns a parallel
                                  predictor. NOTE: with UMA, workers>1 plus an
                                  explicit Analytical Hessian request is an
                                  error; use workers=1 or FiniteDifference.
                                  [default: 1]
  --workers-per-node INTEGER      Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  --freeze-links / --no-freeze-links
                                  Freeze parent atoms of cap hydrogens
                                  (PDB/mmCIF input or XYZ/GJF with --ref-pdb).
                                  [default: freeze-links]
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., '1,3,5').
  --tr-projection [constrained|legacy-active]
                                  Rigid-mode treatment for PHVA. 'constrained'
                                  removes only full-system rigid motions
                                  compatible with frozen anchors (default);
                                  'legacy-active' is deprecated comparison-only
                                  behavior and must not be used for pass/HOSP
                                  transition-state certification.
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB companions
                                  when a PDB template is available.  [default:
                                  convert-files]
  --ref-pdb FILE                  Reference PDB/mmCIF topology to use when the
                                  input is XYZ/GJF (keeps XYZ coordinates).
  --max-write INTEGER             How many modes to export (after sorting per
                                  --sort).  [default: 10]
  --amplitude-ang FLOAT           Animation amplitude (Å) used for both _trj.xyz
                                  and .pdb.  [default: 0.8]
  --n-frames INTEGER              Number of frames per mode animation.
                                  [default: 20]
  --sort [value|abs]              Sort modes by 'value' (cm^-1) or by absolute
                                  value.  [default: value]
  -o, --out-dir TEXT              Output directory.  [default: ./result_freq/]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --temperature FLOAT             Temperature (K) for thermochemistry summary.
                                  [default: 298.15]
  --pressure FLOAT                Pressure (atm) for thermochemistry summary.
                                  [default: 1.0]
  --symmetry-number INTEGER RANGE
                                  External rotational symmetry number used in
                                  the thermochemistry partition function.
                                  [default: 1; x>=1]
  --dump / --no-dump              When True, write 'thermoanalysis.yaml' under
                                  out-dir.  [default: no-dump]
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running frequency analysis.  [default:
                                  no-dry-run]
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
  --hessian-calc-mode [finitedifference|analytical]
                                  How the ML backend computes Hessian. Defaults
                                  to 'FiniteDifference' (can also be set via
                                  YAML).
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
                                  omitted (requires PDB/mmCIF input or --ref-
                                  pdb).
  -m, --multiplicity INTEGER RANGE
                                  Spin multiplicity (2S+1).  [x>=1]
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
                                  Calculator instance). CLI overrides config
                                  YAML; otherwise defaults to get_calculator.
  --deterministic / --no-deterministic
                                  Request strict same-stack PyTorch determinism
                                  (deterministic algorithms + index_reduce_
                                  shim). Slower; raises for detected unsupported
                                  ops; custom calculators are outside its scope.
  --allow-charge-mult-mismatch    Skip the cluster charge/multiplicity electron-
                                  parity check (logs that it was skipped). For
                                  an intentional open-shell or modified-residue
                                  cluster.
  --coord-type [cart|redund|dlc|tric]
                                  Optimization coordinate system
                                  (cart|redund|dlc|tric). cart is the reliable
                                  default used for the published results; dlc
                                  speeds up torsion-rich optimizations.
  -h, --help                      Show this message and exit.
```
