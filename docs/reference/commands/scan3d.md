# pdb2reaction scan3d

```
pdb2reaction ver. 0.2.1.dev116+gff86b942a

Usage: cli scan3d [OPTIONS]

  3D distance scan with harmonic restraints.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure file (.pdb, .xyz, _trj.xyz,
                                  ...). Required unless --csv is provided.
  -s, --scan-lists TEXT           Scan targets: inline Python literal or a
                                  YAML/JSON spec file path.
  --csv FILE                      If provided, skip the 3D scan and read a
                                  precomputed surface.csv from this path. When
                                  used, -i/--input and --scan-lists are
                                  optional.
  -q, --charge INTEGER            Total charge. Required for non-.gjf inputs
                                  unless --ligand-charge is provided (PDB inputs
                                  or XYZ/GJF with --ref-pdb).
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1).
  -o, --out-dir TEXT              Base output directory.  [default:
                                  ./result_scan3d/]
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  -h, --help                      Show this message and exit.
```
