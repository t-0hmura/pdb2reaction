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
  --uma-workers, --workers INTEGER
                                  MLIP predictor workers; >1 spawns a parallel
                                  predictor. NOTE: with UMA, workers>1 plus an
                                  explicit Analytical Hessian request is an
                                  error; use workers=1 or FiniteDifference.
                                  [default: 1]
  --uma-workers-per-node, --workers-per-node INTEGER
                                  Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  --freeze-links / --no-freeze-links
                                  Freeze parent atoms of cap hydrogens
                                  (PDB/mmCIF input or XYZ/GJF with --ref-pdb).
                                  [default: freeze-links]
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., '1,3,5').
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB companions
                                  when a PDB template is available.  [default:
                                  convert-files]
  --ref-pdb FILE                  Reference PDB/mmCIF topology to use when the
                                  input is XYZ/GJF (keeps XYZ coordinates).
  --max-write INTEGER             How many modes to export (after sorting per
                                  --sort).  [default: 10]
  --amplitude-ang FLOAT           Mode-trajectory amplitude (Å) used for both
                                  _trj.xyz and .pdb.  [default: 0.8]
  --n-frames INTEGER              Number of frames per mode trajectory.
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
  --dump / --no-dump              Write 'thermoanalysis.yaml' under out-dir.
                                  [default: no-dump]
  --show-config / --no-show-config
                                  Print the loaded YAML file and its top-level
                                  keys, then continue.  [default: no-show-
                                  config]
  --dry-run / --no-dry-run        Validate options and inputs without running
                                  frequency analysis.  [default: no-dry-run]
  --read-hess FILE                Use the Hessian in this .npy file (e.g. from
                                  freq or tsopt --dump-hess) instead of
                                  computing it: the Cartesian Hessian of the
                                  input geometry in Hartree/bohr^2, for all
                                  atoms or only the movable ones.
  --dump-hess FILE                Save the Hessian as a NumPy .npy array
                                  (Cartesian, Hartree/bohr^2; movable atoms only
                                  when atoms are frozen) for '--read-hess' in
                                  freq, tsopt, or irc, or for other programs.
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
  --hessian-calc-mode [finitedifference|analytical]
                                  How the ML backend computes the Hessian (can
                                  also be set via YAML).  [default:
                                  (FiniteDifference)]
  --hess-device [auto|cuda|cpu]   Device for post-evaluation Hessian placement
                                  and diagonalization (auto/cuda/cpu). Use 'cpu'
                                  to move the evaluated Hessian off GPU before
                                  diagonalization. The calculator itself still
                                  runs on its own device.  [default: auto]
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
  -q, --charge INTEGER            Total charge. Required for non-.gjf inputs
                                  unless --ligand-charge is provided (.gjf
                                  templates inherit the charge automatically).
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB/mmCIF input or --ref-
                                  pdb).
  -m, --multiplicity INTEGER RANGE
                                  Spin multiplicity (2S+1).  [default: (1);
                                  x>=1]
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
