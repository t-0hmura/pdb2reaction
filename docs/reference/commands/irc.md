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
  -i, --input FILE                Input structure file (.pdb, .xyz, _trj.xyz,
                                  etc.).  [required]
  --workers INTEGER               MLIP predictor workers; >1 spawns a parallel
                                  predictor. NOTE: analytical Hessian raises a
                                  RuntimeError when workers>1; pass --hessian-
                                  calc-mode FiniteDifference explicitly.
                                  [default: 1]
  --workers-per-node INTEGER      Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  --max-cycles INTEGER            Maximum number of IRC steps; used unless YAML
                                  sets irc.max_cycles. Defaults to 125 when not
                                  provided.
  --step-size FLOAT               Step length in Bohr (unweighted Cartesian
                                  coordinates); used unless YAML sets
                                  irc.step_length. Default: 0.10 Bohr.
  --root INTEGER                  Imaginary mode index used for the initial
                                  displacement; used unless YAML sets irc.root.
                                  Defaults to 0.
  --forward / --no-forward        Run the forward IRC; used unless YAML sets
                                  irc.forward. Defaults to True.
  --backward / --no-backward      Run the backward IRC; used unless YAML sets
                                  irc.backward. Defaults to True.
  --freeze-links / --no-freeze-links
                                  Freeze parent atoms of link hydrogens (PDB
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
  -o, --out-dir TEXT              Output directory.  [default: ./result_irc/]
  --hessian-calc-mode [finitedifference|analytical]
                                  How the ML backend builds the Hessian
                                  (Analytical or FiniteDifference); used unless
                                  YAML sets calc.hessian_calc_mode. Defaults to
                                  'FiniteDifference'.
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
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1).
  --precision [fp32|fp64]         MLIP backend precision: fp32 (default) or
                                  fp64. Routed to backend-specific kwargs (UMA
                                  precision / ORB precision / MACE
                                  default_dtype). aimnet2: fp32 no-op; fp64
                                  rejected.
  --deterministic / --no-deterministic
                                  Strict bit-reproducible GPU runs
                                  (deterministic algorithms + index_reduce_
                                  shim). Slower; raises if unsupported. Default
                                  off.
  --allow-charge-mult-mismatch    Skip the cluster charge/multiplicity electron-
                                  parity check (logs that it was skipped). For
                                  an intentional open-shell or modified-residue
                                  cluster.
  --irc-pos-def / --no-irc-pos-def
                                  Require pos-def Hessian at IRC convergence
                                  (blocks shoulder false-convergence).
  -h, --help                      Show this message and exit.
```
