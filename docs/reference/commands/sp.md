# `pdb2reaction sp`

```text
Usage: pdb2reaction sp [OPTIONS]

  Compute single-point energy and forces (and optionally the Hessian).

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+detailed step logging
                                  and output paths; 3=everything (full config
                                  blocks, per-file paths, DEBUG logging).
                                  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure (PDB / mmCIF / XYZ / GJF).
                                  [required]
  --uma-workers, --workers INTEGER
                                  MLIP predictor workers; >1 spawns a parallel
                                  predictor (needs fairchem-core[ray]). With
                                  UMA, workers>1 plus an explicit Analytical
                                  Hessian request is an error; use workers=1 or
                                  FiniteDifference.  [default: 1]
  --uma-workers-per-node, --workers-per-node INTEGER
                                  Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  -o, --out-dir TEXT              Output directory.  [default: ./result_sp/]
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., '1,3,5').
  --hess / --no-hess              Also compute a Hessian and save it to
                                  hessian.npy (movable atoms only when atoms are
                                  frozen).  [default: no-hess]
  --hessian-calc-mode [analytical|finitedifference]
                                  Hessian calculation method when --hess is set.
                                  FiniteDifference is the default for every
                                  backend; UMA, ORB, MACE, AIMNet2, and DFT also
                                  support Analytical.  [default:
                                  (FiniteDifference)]
  --config FILE                   YAML config file with sections (calc:, geom:,
                                  …).
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the plan without
                                  running the calculation.  [default: no-dry-
                                  run]
  --out-json / --no-out-json      Write result.json to the output directory.
                                  [default: no-out-json]
  -b, --backend [uma|orb|mace|aimnet2|dft]
                                  Energy/force calculator backend.  [default:
                                  uma]
  --solvent TEXT                  Optional xTB or AnySolv solvent correction for
                                  MLIP backends; dft uses native PySCF PCM/SMD.
                                  'none' disables it. AnySolv MoE is recommended
                                  for accuracy-sensitive reaction calculations;
                                  see docs for model access and setup.
                                  [default: none]
  --solvent-model [alpb|cpcmx|anysolv|pcm|smd]
                                  Solvent model: ALPB/CPCMx or AnySolv for MLIP
                                  backends; PCM/SMD for dft. For AnySolv, MoE is
                                  recommended; select its weights with
                                  calc.anysolv_checkpoint (see docs for model
                                  access).  [default: alpb]
  --solvent-xtb-cmd TEXT          Command for the xTB solvent correction,
                                  including optional xTB arguments. If xTB SCC
                                  convergence is poor, increasing --etemp may
                                  help (for example: 'xtb --etemp 1000').
  -q, --charge INTEGER            Total charge. Required unless --ligand-charge,
                                  YAML calc.charge, or a .gjf input gives it.
  -l, --ligand-charge TEXT        Total ligand charge or per-resname mapping
                                  (e.g., GPP:-3,SAM:1) used to derive charge
                                  when -q is omitted (requires PDB/mmCIF input).
  -m, --multiplicity INTEGER RANGE
                                  Spin multiplicity (2S+1).  [default: (1);
                                  x>=1]
  --precision [tf32|fp32|fp64]    MLIP backend precision: tf32, fp32 (IEEE) or
                                  fp64. tf32 uses TF32 matmul for speed (uma,
                                  orb only); avoid it for Hessians. aimnet2:
                                  fp32 only.  [default: (per backend: uma fp32;
                                  orb, mace fp64)]
  --compile / --no-compile        torch.compile the UMA or ORB model: slower
                                  first call, slightly faster afterwards; worth
                                  it only for long runs.  [default: (off)]
  --backend-model TEXT            Model variant for the selected --backend (e.g.
                                  uma-s-1p2p1 / uma-s-1p2p1-ft-d4rev / uma-m-1p1
                                  for uma, orbmol-v2 / orb_v3_conservative_omol
                                  for orb, MACE-OMOL-0 / MACE-POLAR-1-M for
                                  mace, aimnet2 / aimnet2-nse for aimnet2).
                                  [default: (the selected backend's own model)]
  --calc-file FILE                Python file exposing get_calculator(...) -> an
                                  ASE Calculator, used as the energy/gradient
                                  backend (overrides --backend). Works with GFN-
                                  xTB, DFTB+ or any ASE engine. See --calc-file-
                                  func-name.
  --calc-factory, --calc-file-func-name TEXT
                                  Name of the callable in --calc-file that
                                  returns an ASE Calculator (or a module-level
                                  Calculator instance).  [default:
                                  (get_calculator)]
  --deterministic / --no-deterministic
                                  Make repeated runs on the same software
                                  versions and hardware give identical results
                                  (PyTorch deterministic algorithms). Slower;
                                  raises an error on unsupported ops; custom
                                  calculators are outside its scope.  [default:
                                  no-deterministic]
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
                                  Save a structure-bound PySCF checkpoint.
                                  [default: (disabled)]
  --scf-checkpoint FILE           Load/save the optional structure-bound PySCF
                                  checkpoint.
  --scf-stepwise-grid / --no-scf-stepwise-grid
                                  Converge the first SCF on a coarse grid, then
                                  on the final grid (later SCFs reuse the
                                  previous density as usual).  [default:
                                  (enabled)]
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
