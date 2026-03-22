# pdb2reaction irc

```
pdb2reaction ver. 0.2.1.dev116+gff86b942a

Usage: cli irc [OPTIONS]

  Run an IRC calculation with EulerPC. Only the documented CLI options are
  accepted; all other settings come from YAML.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure file (.pdb, .xyz, _trj.xyz,
                                  etc.).  [required]
  -q, --charge INTEGER            Total charge. Required for non-.gjf inputs
                                  unless --ligand-charge is provided (PDB inputs
                                  or XYZ/GJF with --ref-pdb).
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1).
  --max-cycles INTEGER            Maximum number of IRC steps; used unless YAML
                                  sets irc.max_cycles. Defaults to 125 when not
                                  provided.
  --step-size FLOAT               Step length in Bohr (unweighted Cartesian
                                  coordinates); used unless YAML sets
                                  irc.step_length. Default: 0.10 Bohr.
  --forward / --no-forward        Run the forward IRC; used unless YAML sets
                                  irc.forward. Defaults to True.
  --backward / --no-backward      Run the backward IRC; used unless YAML sets
                                  irc.backward. Defaults to True.
  -o, --out-dir TEXT              Output directory.  [default: ./result_irc/]
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  -h, --help                      Show this message and exit.
```
