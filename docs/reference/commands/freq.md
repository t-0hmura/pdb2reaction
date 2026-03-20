# `pdb2reaction freq`

```text

Usage: pdb2reaction freq [OPTIONS]

  Vibrational frequency analysis and mode writer (+ default thermochemistry
  summary).

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
  --workers-per-node INTEGER      Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1).
  --freeze-links / --no-freeze-links
                                  Freeze parent atoms of link hydrogens (PDB
                                  only).  [default: freeze-links]
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB companions
                                  when a PDB template is available.  [default:
                                  convert-files]
  --ref-pdb FILE                  Reference PDB topology to use when the input
                                  is XYZ/GJF (keeps XYZ coordinates).
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
  --dump / --no-dump              When True, write 'thermoanalysis.yaml' under
                                  out-dir.  [default: no-dump]
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running frequency analysis.  [default:
                                  no-dry-run]
  --hessian-calc-mode [finitedifference|analytical]
                                  How the ML backend computes Hessian. Defaults
                                  to 'FiniteDifference' (can also be set via
                                  YAML).
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  --solvent-model [alpb|cpcmx]    xTB solvent model.  [default: alpb]
  -h, --help                      Show this message and exit.
```
