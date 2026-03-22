# pdb2reaction path-search

```
pdb2reaction ver. 0.2.1.dev116+gff86b942a

Usage: cli path-search [OPTIONS]

  Multistep MEP search via recursive GSM/DMF segmentation.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Two or more structures in reaction order.
                                  Repeat -i/--input for each path.  [required]
  --mep-mode [gsm|dmf]            MEP optimizer: Growing String Method (gsm) or
                                  Direct Max Flux (dmf).  [default: gsm]
  --refine-mode [peak|minima]     Refinement seed selection around the highest-
                                  energy image: 'peak' uses HEI±1, 'minima' uses
                                  the nearest local minima in each direction.
                                  Defaults to peak for gsm and minima for dmf
                                  when omitted.
  -q, --charge INTEGER            Total charge. Required for non-.gjf inputs
                                  unless --ligand-charge derives it from PDB
                                  inputs.
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (PDB inputs only).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1; defaults from a .gjf
                                  template when available, otherwise 1).
  --max-nodes INTEGER             Number of internal nodes (string has
                                  max_nodes+2 images including endpoints). Used
                                  for *segment* GSM unless overridden by YAML
                                  search.max_nodes_segment.  [default: 20]
  -o, --out-dir TEXT              Output directory.  [default:
                                  ./result_path_search/]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  -h, --help                      Show this message and exit.
```
