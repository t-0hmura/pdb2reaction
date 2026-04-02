# pdb2reaction path-opt

```
pdb2reaction ver. 0.2.1.dev116+gff86b942a

Usage: cli path-opt [OPTIONS]

  MEP optimization via GSM or DMF.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE...             Two endpoint structures (reactant and
                                  product); accepts .pdb or .xyz.  [required]
  --mep-mode [gsm|dmf]            MEP optimizer: Growing String Method (gsm) or
                                  Direct Max Flux (dmf).  [default: gsm]
  -q, --charge INTEGER            Total charge. Required unless a .gjf template
                                  provides charge metadata or --ligand-charge is
                                  supplied for PDB inputs.
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1).
  --max-nodes INTEGER             Maximum number of internal nodes (string has
                                  up to max_nodes+2 images including endpoints).
                                  [default: 20]
  -o, --out-dir TEXT              Output directory.  [default:
                                  ./result_path_opt/]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  -h, --help                      Show this message and exit.
```
