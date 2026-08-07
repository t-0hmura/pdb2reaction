# `pdb2reaction sp`

```text
Usage: pdb2reaction sp [OPTIONS]

  Compute a single-point MLIP energy + forces (and optionally Hessian).

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+detailed step logging
                                  and deliverable paths; 3=everything (full
                                  config blocks, per-file paths, DEBUG logging).
                                  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure (PDB / mmCIF / XYZ / GJF).
                                  [required]
  --workers INTEGER               MLIP predictor workers; >1 spawns a parallel
                                  predictor. NOTE: with UMA, workers>1 plus an
                                  explicit Analytical Hessian request is an
                                  error; use workers=1 or FiniteDifference.
                                  [default: 1]
  --workers-per-node INTEGER      Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  -o, --out-dir TEXT              Output directory.  [default: ./result_sp/]
  --hess / --no-hess              Also compute a Hessian and save it to
                                  hessian.npy (active block when YAML
                                  geom.freeze_atoms is non-empty).  [default:
                                  no-hess]
  --hessian-calc-mode [analytical|finitedifference]
                                  Hessian backend when --hess is set.
                                  FiniteDifference is the default for every
                                  backend; UMA, ORB, MACE, and AIMNet2 also
                                  support an explicit Analytical request.
                                  [default: (FiniteDifference)]
  --convert-files / --no-convert-files
                                  Compatibility toggle; sp writes no structure
                                  companions.  [default: convert-files]
  --config FILE                   YAML config file with sections (calc:, geom:,
                                  …).
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the plan without
                                  running the calculation.  [default: no-dry-
                                  run]
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
                                  omitted (requires PDB/mmCIF input).
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
  -h, --help                      Show this message and exit.
```
