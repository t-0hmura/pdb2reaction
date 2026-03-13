# `pdb2reaction opt`

```text

Usage: pdb2reaction opt [OPTIONS]

  Single-structure geometry optimization using LBFGS or RFO.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure file (.pdb, .xyz, _trj.xyz,
                                  ...).  [required]
  -q, --charge INTEGER            Total charge. Required for non-.gjf inputs
                                  unless --ligand-charge is provided (PDB inputs
                                  or XYZ/GJF with --ref-pdb).
  --workers INTEGER               MLIP predictor workers; >1 spawns a parallel
                                  predictor (disables analytic Hessian).
                                  [default: 1]
  --workers-per-node INTEGER      Workers per node when using a parallel UMA
                                  predictor (workers>1).  [default: 1]
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1) for the ML region.
  --dist-freeze TEXT              Python-like list(s) of (i,j,target_A) to
                                  restrain distances (target optional).
  --one-based / --zero-based      Interpret --dist-freeze indices as 1-based
                                  (default) or 0-based.  [default: one-based]
  --bias-k FLOAT                  Harmonic restraint strength k [eV/Å^2] for
                                  --dist-freeze.  [default: 10.0]
  --freeze-links / --no-freeze-links
                                  Freeze parent atoms of link hydrogens (PDB
                                  only).  [default: freeze-links]
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB/GJF
                                  companions based on the input format.
                                  [default: convert-files]
  --ref-pdb FILE                  Reference PDB topology to use when the input
                                  is XYZ/GJF (keeps XYZ coordinates).
  --max-cycles INTEGER            Maximum number of optimization cycles.
                                  [default: 10000]
  --opt-mode [grad|hess|lbfgs|rfo]
                                  Optimization mode: grad (lbfgs) or hess (rfo).
                                  Aliases lbfgs/rfo are accepted.  [default:
                                  grad]
  --flatten / --no-flatten        Enable/disable imaginary-mode flatten loop
                                  after optimization.  [default: no-flatten]
  --dump / --no-dump              Write optimization trajectory to
                                  'optimization_trj.xyz'.  [default: no-dump]
  -o, --out-dir TEXT              Output directory.  [default: ./result_opt/]
  --thresh TEXT                   Convergence preset (gau_loose|gau|gau_tight|ga
                                  u_vtight|baker|never). Defaults to 'gau' when
                                  not provided.
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running optimization.  [default: no-
                                  dry-run]
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  --solvent-model [alpb|cpcmx]    xTB solvent model.  [default: alpb]
  -h, --help                      Show this message and exit.
```
