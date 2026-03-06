# `pdb2reaction tsopt`

```text

Usage: pdb2reaction tsopt [OPTIONS]

  Transition state optimization (Dimer or RS-I-RFO).

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure file (.pdb, .xyz, _trj.xyz,
                                  ...).  [required]
  -q, --charge INTEGER            Total charge. Required for non-.gjf inputs
                                  unless --ligand-charge is provided (PDB inputs
                                  or XYZ/GJF with --ref-pdb).
  --workers INTEGER               UMA predictor workers; >1 spawns a parallel
                                  predictor (disables analytic Hessian).
                                  [default: 1]
  --workers-per-node INTEGER      Workers per node when using a parallel UMA
                                  predictor (workers>1).  [default: 1]
  --ligand-charge TEXT            Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1) for the ML region.
  --freeze-links / --no-freeze-links
                                  Freeze parent atoms of link hydrogens (PDB
                                  only).  [default: freeze-links]
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB/GJF
                                  companions based on the input format.
                                  [default: convert-files]
  --ref-pdb FILE                  Reference PDB topology to use when the input
                                  is XYZ/GJF (keeps XYZ coordinates).
  --max-cycles INTEGER            Max cycles / steps cap.  [default: 10000]
  --flatten / --no-flatten        Enable the extra-imaginary-mode flattening
                                  loop (grad: dimer loop, hess: post-RSIRFO).
                                  [default: no-flatten]
  --opt-mode [grad|hess|dimer|rsirfo]
                                  grad (dimer) or hess (rsirfo). Aliases
                                  dimer/rsirfo are accepted.  [default: hess]
  --dump / --no-dump              Write optimization trajectory to the output
                                  directory.  [default: no-dump]
  --out-dir TEXT                  Output directory.  [default: ./result_tsopt/]
  --thresh TEXT                   Convergence preset for the active optimizer (g
                                  au_loose|gau|gau_tight|gau_vtight|baker|never)
                                  . Defaults to 'baker' when not provided.
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running TS optimization.  [default:
                                  no-dry-run]
  --hessian-calc-mode [finitedifference|analytical]
                                  Choose UMA Hessian evaluation mode (used
                                  unless YAML sets calc.hessian_calc_mode).
                                  Defaults to 'FiniteDifference'.
  --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  --solvent-model [alpb|cpcmx]    xTB solvent model.  [default: alpb]
  -h, --help                      Show this message and exit.
```
