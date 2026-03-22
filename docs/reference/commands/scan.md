# pdb2reaction scan

```
pdb2reaction ver. 0.2.1.dev116+gff86b942a

Usage: cli scan [OPTIONS]

  Bond-length driven scan with staged harmonic restraints and relaxation.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure file (.pdb, .xyz, _trj.xyz,
                                  ...).  [required]
  -s, --scan-lists TEXT           Scan targets: inline Python literal (e.g.
                                  '[(1,5,1.4)]') or a YAML/JSON spec file path.
                                  Multiple inline literals define sequential
                                  stages.
  -q, --charge INTEGER            Total charge. Required for non-.gjf inputs
                                  unless --ligand-charge is provided (PDB inputs
                                  or XYZ/GJF with --ref-pdb).
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1).
  -o, --out-dir TEXT              Base output directory.  [default:
                                  ./result_scan/]
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  -h, --help                      Show this message and exit.
```
