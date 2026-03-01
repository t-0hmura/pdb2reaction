# `pdb2reaction path-opt`

```text
pdb2reaction ver. 0.2.1.dev39+gd7f50241d

Usage: pdb2reaction path-opt [OPTIONS]

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
  --max-nodes INTEGER             Number of internal nodes (string has
                                  max_nodes+2 images including endpoints).
                                  [default: 10]
  --max-cycles INTEGER            Maximum optimization cycles.  [default: 300]
  --climb / --no-climb            Search for a transition state (climbing image)
                                  after path growth.  [default: climb]
  --opt-mode [grad|hess]          Single-structure optimizer for endpoint
                                  preoptimization: grad (=LBFGS) or hess (=RFO).
                                  [default: grad]
  --dump / --no-dump              Write optimizer trajectory and restarts during
                                  the run.  [default: no-dump]
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB/GJF
                                  companions based on the input format.
                                  [default: convert-files]
  --ref-pdb FILE                  Reference PDB topology to use when the input
                                  is XYZ/GJF (keeps XYZ coordinates).
  --out-dir TEXT                  Output directory.  [default:
                                  ./result_path_opt/]
  --thresh TEXT                   Convergence preset for endpoint
                                  preoptimization only (gau_loose|gau|gau_tight|
                                  gau_vtight|baker|never). Defaults to 'gau'
                                  when not provided.
  --thresh-stopt TEXT             Convergence preset for the string optimizer
                                  (stopt) (gau_loose|gau|gau_tight|gau_vtight|ba
                                  ker|never).  [default: gau]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running path optimization.  [default:
                                  no-dry-run]
  --preopt / --no-preopt          If True, preoptimize each endpoint via the
                                  selected single-structure optimizer
                                  (LBFGS/RFO) before alignment and GSM.
                                  [default: no-preopt]
  --preopt-max-cycles INTEGER     Maximum cycles for endpoint preoptimization
                                  (applies to the chosen optimizer; only used
                                  when --preopt True).  [default: 10000]
  --fix-ends / --no-fix-ends      Fix structures of input endpoints during GSM.
                                  [default: no-fix-ends]
  -h, --help                      Show this message and exit.
```
