# `pdb2reaction path-search`

```text

Usage: pdb2reaction path-search [OPTIONS]

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
  --workers INTEGER               UMA predictor workers; >1 spawns a parallel
                                  predictor (disables analytic Hessian).
                                  [default: 1]
  --workers-per-node INTEGER      Workers per node when using a parallel UMA
                                  predictor (workers>1).  [default: 1]
  --ligand-charge TEXT            Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (PDB inputs only).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1) for the ML region
                                  (defaults from a .gjf template when available,
                                  otherwise 1).
  --freeze-links / --no-freeze-links
                                  Freeze parent atoms of link hydrogens (PDB
                                  only).  [default: freeze-links]
  --max-nodes INTEGER             Number of internal nodes (string has
                                  max_nodes+2 images including endpoints). Used
                                  for *segment* GSM unless overridden by YAML
                                  search.max_nodes_segment.  [default: 20]
  --max-cycles INTEGER            Maximum GSM optimization cycles.  [default:
                                  300]
  --climb / --no-climb            Enable climbing image for standard GSM
                                  segments (bridge segments always disable
                                  climbing).  [default: climb]
  --opt-mode [grad|hess]          Single-structure optimizer: grad (=LBFGS) or
                                  hess (=RFO).  [default: grad]
  --dump / --no-dump              Write GSM/single-optimization trajectories
                                  during the run.  [default: no-dump]
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB/GJF
                                  companions based on the input format.
                                  [default: convert-files]
  --out-dir TEXT                  Output directory.  [default:
                                  ./result_path_search/]
  --thresh TEXT                   Convergence preset for single-structure
                                  optimizations only (gau_loose|gau|gau_tight|ga
                                  u_vtight|baker|never). Defaults to 'gau' when
                                  not provided.
  --thresh-stopt TEXT             Convergence preset for the string optimizer
                                  (stopt) (gau_loose|gau|gau_tight|gau_vtight|ba
                                  ker|never). Defaults to 'gau_loose' when not
                                  provided.
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running path search.  [default: no-
                                  dry-run]
  --preopt / --no-preopt          If False, skip initial single-structure
                                  optimizations of inputs.  [default: preopt]
  --align / --no-align            After preoptimization, align all inputs to the
                                  *first* input and match freeze_atoms using the
                                  align_freeze_atoms API. When --align is True
                                  and --ref-full-pdb is provided, the first
                                  reference PDB will be used for all pairs in
                                  the final merge.  [default: align]
  --ref-full-pdb FILE             Full-size template PDBs in the same reaction
                                  order as --input. With --align True, only the
                                  *first* provided reference PDB is used for all
                                  pairs in the final merge (you may pass just
                                  one).
  --ref-pdb FILE                  Pocket reference PDBs used only for the final
                                  full-system merge. Useful when --input uses
                                  XYZ/GJF intermediates but PDB snapshots exist
                                  for merging. Must match the number and order
                                  of --input.
  --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  --solvent-model [alpb|cpcmx]    xTB solvent model.  [default: alpb]
  -h, --help                      Show this message and exit.
```
