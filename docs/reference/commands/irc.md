# `pdb2reaction irc`

```text
Usage: pdb2reaction irc [OPTIONS]

  Run an IRC calculation with EulerPC. Only the documented CLI options are
  accepted; all other settings come from YAML.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure file (.pdb, .cif, .mmcif,
                                  .xyz, .gjf, _trj.xyz, etc.).  [required]
  --workers INTEGER               MLIP predictor workers; >1 spawns a parallel
                                  predictor. NOTE: with UMA, workers>1 plus an
                                  explicit Analytical Hessian request is an
                                  error; use workers=1 or FiniteDifference.
                                  [default: 1]
  --workers-per-node INTEGER      Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  --max-cycles INTEGER            Maximum number of IRC steps; an explicit value
                                  overrides YAML irc.max_cycles. Defaults to 125
                                  when not provided.  [default: (125)]
  --step-size FLOAT               Step length in Bohr (unweighted Cartesian
                                  coordinates); an explicit value overrides YAML
                                  irc.step_length. Default: 0.10 Bohr.
                                  [default: (0.10)]
  --never-stop / --no-never-stop  Ignore RMS-gradient, hard-gradient, energy-
                                  rise, and energy-change stops and trace until
                                  --max-cycles. Numerical/integration failures
                                  and external interruption still stop the run.
                                  An explicit toggle overrides YAML
                                  irc.never_stop; default off.  [default: (no-
                                  never-stop)]
  --root INTEGER                  Imaginary mode index used for the initial
                                  displacement; an explicit value overrides YAML
                                  irc.root. Defaults to 0.  [default: (0)]
  --forward / --no-forward        Run the forward IRC; an explicit toggle
                                  overrides YAML irc.forward. Defaults to True.
                                  [default: (forward)]
  --backward / --no-backward      Run the backward IRC; an explicit toggle
                                  overrides YAML irc.backward. Defaults to True.
                                  [default: (backward)]
  --freeze-links / --no-freeze-links
                                  Freeze parent atoms of cap hydrogens
                                  (PDB/mmCIF input or XYZ/GJF with --ref-pdb).
                                  [default: freeze-links]
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., '1,3,5').
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB/CIF
                                  companions based on the input format.
                                  [default: convert-files]
  --ref-pdb FILE                  Reference PDB/mmCIF topology to use when the
                                  input is XYZ/GJF (keeps XYZ coordinates).
  -o, --out-dir TEXT              Output directory.  [default: ./result_irc/]
  --hessian-calc-mode [finitedifference|analytical]
                                  How the ML backend builds the Hessian
                                  (Analytical or FiniteDifference); an explicit
                                  value overrides YAML calc.hessian_calc_mode.
                                  Defaults to 'FiniteDifference'.  [default:
                                  (FiniteDifference)]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running IRC.  [default: no-dry-run]
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
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
                                  / off:small for mace). Default: the backend's
                                  built-in model.  [default: (the selected
                                  backend's own model)]
  --calc-file FILE                Python file exposing get_calculator(...) -> an
                                  ASE Calculator, used as the energy/gradient
                                  backend (overrides --backend). Couples GFN-xTB
                                  / DFTB+ / any ASE engine. See --calc-file-
                                  func-name.
  --calc-file-func-name, --calc-factory TEXT
                                  Name of the callable in --calc-file that
                                  returns an ASE Calculator (or a module-level
                                  Calculator instance). CLI overrides config
                                  YAML; otherwise defaults to get_calculator.
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
  --irc-pos-def / --no-irc-pos-def
                                  Require pos-def Hessian at IRC convergence
                                  (blocks shoulder false-convergence).
                                  [default: (no-irc-pos-def (rms-only
                                  criterion))]
  -h, --help                      Show this message and exit.
```
