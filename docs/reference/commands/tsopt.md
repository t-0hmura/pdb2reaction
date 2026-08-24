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
  -i, --input FILE                Single-geometry input (.pdb, .cif, .mmcif,
                                  .xyz, or .gjf). Extract a trajectory frame to
                                  .xyz before use.  [required]
  --ref-mode FILE                 Advanced/internal Cartesian reference
                                  direction(s) for Hessian TS root selection and
                                  overlap tracking. Accepts .npz path-mode
                                  caches, .npy arrays, or whitespace text
                                  containing one 3N vector or a 2-D candidate
                                  table. This guides mode identity; it does not
                                  replace the Hessian and is not supported by
                                  Dimer. The all workflow supplies it from the
                                  MEP; standalone tsopt users normally leave it
                                  unset.
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
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB/CIF/GJF
                                  companions based on the input format.
                                  [default: convert-files]
  --ref-pdb FILE                  Reference PDB/mmCIF topology to use when the
                                  input is XYZ/GJF (keeps XYZ coordinates).
  --max-cycles INTEGER RANGE      Maximum number of optimization cycles.
                                  [default: (100000); x>=1]
  --flatten / --no-flatten        Enable the extra-imaginary-mode flattening
                                  loop (grad: dimer loop, hess: post-RS-P-RFO).
                                  [default: no-flatten]
  --opt-mode [grad|hess|dimer|rsirfo|trim|rsprfo]
                                  TS optimizer: 'grad'/'dimer' → Hessian Guided
                                  Dimer; 'hess'/'rsprfo' → RS-P-RFO (Banerjee,
                                  default); 'trim' → TRIM (Helgaker); 'rsirfo' →
                                  RS-I-RFO.  [default: hess]
  --dump / --no-dump              Write the per-cycle trajectory
                                  ('optimization_trj.xyz' for RS-P-RFO/RS-I-
                                  RFO/TRIM, 'optimization_all_trj.xyz' for
                                  Dimer).  [default: no-dump]
  -o, --out-dir TEXT              Output directory.  [default: ./result_tsopt/]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for the active optimizer (g
                                  au_loose|gau|gau_tight|gau_vtight|baker|never)
                                  .  [default: (baker)]
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
                                  Choose MLIP Hessian evaluation mode. YAML
                                  supplies the value when this option is
                                  omitted; explicit CLI wins.  [default:
                                  (FiniteDifference)]
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Experimental, computationally expensive xTB
                                  solvent delta correction. Examples: water,
                                  methanol, acetonitrile, dmso, thf, toluene.
                                  'none' disables it.  [default: none]
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
  --coord-type [cart|redund|dlc|tric]
                                  Optimization coordinate system
                                  (cart|redund|dlc|tric).  [default: (cart)]
  --print-every INTEGER RANGE     Print optimizer status every N cycles.
                                  [default: (100); x>=1]
  --stop-plateau / --no-stop-plateau
                                  Stop when the energy stops changing while the
                                  convergence criteria are still unmet, and
                                  report the run as stalled. It never signals
                                  convergence; --max-cycles remains the real
                                  bound.  [default: no-stop-plateau]
  --stop-plateau-thresh FLOAT     Energy range (hartree) below which --stop-
                                  plateau treats the window as flat.  [default:
                                  (1e-4)]
  --stop-plateau-window INTEGER   Number of consecutive cycles --stop-plateau
                                  inspects.  [default: (50)]
  -h, --help                      Show this message and exit.
```
