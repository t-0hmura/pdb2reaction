# pdb2reaction tsopt

```
pdb2reaction ver. 0.2.1.dev116+gff86b942a

Usage: cli tsopt [OPTIONS]

  Transition state optimization (Dimer or RS-I-RFO).

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure file (.pdb, .xyz, _trj.xyz,
                                  ...).  [required]
  -q, --charge INTEGER            Total charge. Required for non-.gjf inputs
                                  unless --ligand-charge is provided (PDB inputs
                                  or XYZ/GJF with --ref-pdb).
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1).
  --max-cycles INTEGER            Maximum number of optimization cycles.
                                  [default: 10000]
  --opt-mode [grad|hess|dimer|rsirfo]
                                  TS optimizer: 'grad'/'dimer' → Hessian Guided
                                  Dimer; 'hess'/'rsirfo' → RS-I-RFO.  [default:
                                  hess]
  -o, --out-dir TEXT              Output directory.  [default: ./result_tsopt/]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  -h, --help                      Show this message and exit.
```
