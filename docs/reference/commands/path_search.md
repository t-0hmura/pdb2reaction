# `pdb2reaction path-search`

```text
Usage: pdb2reaction path-search [OPTIONS]

  Multistep MEP search via recursive GSM/DMF segmentation.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Two or more PDB, mmCIF, XYZ, or GJF structures
                                  in reaction order. Repeat -i/--input for each
                                  path.  [required]
  --mep-mode [gsm|dmf]            MEP optimizer: Growing String Method (gsm) or
                                  Direct Max Flux (dmf).  [default: gsm]
  --dmf-backend [cpu|gpu]         DMF compute backend (--mep-mode dmf only): gpu
                                  (dmf.torch / CUDA) or cpu (dmf / NumPy). On a
                                  GPU out-of-memory error, retry with cpu.
                                  [default: gpu]
  --refine-mode [peak|minima]     Refinement seed selection around the highest-
                                  energy image: 'peak' uses HEI±1, 'minima' uses
                                  the nearest local minima in each direction.
                                  Defaults to peak for gsm and minima for dmf
                                  when omitted.  [default: (peak for gsm, minima
                                  for dmf)]
  -q, --charge INTEGER            Total charge. Required for non-.gjf inputs
                                  unless --ligand-charge derives it from
                                  PDB/mmCIF inputs.
  --workers INTEGER               MLIP predictor workers; >1 spawns a parallel
                                  predictor. NOTE: with UMA, workers>1 plus an
                                  explicit Analytical Hessian request is an
                                  error; use workers=1 or FiniteDifference.
                                  [default: 1]
  --workers-per-node INTEGER      Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (PDB/mmCIF inputs only).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1; inherits from a .gjf
                                  template when available).  [default: (1 (or
                                  the .gjf template value))]
  --freeze-links / --no-freeze-links
                                  Freeze parent atoms of cap hydrogens
                                  (PDB/mmCIF input only).  [default: freeze-
                                  links]
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., '1,3,5').
  --max-nodes INTEGER             Number of movable internal images per GSM/DMF
                                  segment; the complete segment has max_nodes+2
                                  images including endpoints. When not given,
                                  YAML search.max_nodes_segment applies.
                                  [default: 20]
  --gsm-param [equi|energy]       GSM node parameterization after string growth.
                                  The energy scheme concentrates nodes in high-
                                  energy regions and may be tried when an
                                  equidistant path skips the reaction-coordinate
                                  region near the HEI.  [default: (equi)]
  --max-cycles-gsm INTEGER RANGE  Maximum GSM string-optimizer cycles for the
                                  MEP stage.  [default: (300); x>=1]
  --max-cycles-dmf INTEGER RANGE  Maximum IPOPT iterations for the DMF MEP
                                  stage. This is a solver iteration count, not a
                                  string-optimizer cycle count.  [default:
                                  (300); x>=1]
  --climb / --no-climb            Enable climbing image for standard GSM
                                  segments (bridge segments always disable
                                  climbing).  [default: climb]
  --opt-mode [grad|hess]          Single-structure optimizer: grad (=LBFGS) or
                                  hess (=RFO).  [default: grad]
  --dump / --no-dump              Write GSM/single-optimization trajectories
                                  during the run.  [default: no-dump]
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB/CIF/GJF
                                  companions based on the input format.
                                  [default: convert-files]
  -o, --out-dir TEXT              Output directory.  [default:
                                  ./result_path_search/]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for single-structure
                                  optimizations only (gau_loose|gau|gau_tight|ga
                                  u_vtight|baker|never).  [default: (gau)]
  --thresh-gsm [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for the GSM string
                                  optimizer (gau_loose|gau|gau_tight|gau_vtight|
                                  baker|never).  [default: (gau_loose)]
  --thresh-dmf TEXT               IPOPT dual-infeasibility tolerance for the DMF
                                  path optimizer: tight (0.04) | middle (0.10) |
                                  loose (0.20) or a positive float. This is not
                                  a Gaussian preset.  [default: (tight)]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running path search.  [default: no-
                                  dry-run]
  --preopt / --no-preopt          If True, run initial single-structure
                                  optimizations of inputs.  [default: preopt]
  --align / --no-align            After preoptimization, align all inputs to the
                                  *first* input and match freeze_atoms using the
                                  align_freeze_atoms API.  [default: align]
  --ref-full-pdb FILE             Static full-size PDB/mmCIF template for
                                  coordinate merging used for inspection. Use
                                  the template corresponding to the first -i
                                  input.
  --write-ref-merge / --no-write-ref-merge
                                  Write mep_w_ref/hei_w_ref coordinate
                                  composites for inspection. Requires --align
                                  and --ref-full-pdb.  [default: no-write-ref-
                                  merge]
  --ref-pdb FILE                  Pocket reference PDB/mmCIF files used only for
                                  the final full-system merge. Useful when
                                  --input uses XYZ/GJF intermediates but PDB
                                  snapshots exist for merging. Must match the
                                  number and order of --input.
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Experimental, computationally expensive xTB
                                  solvent delta correction. Examples: water,
                                  methanol, acetonitrile, dmso, thf, toluene.
                                  'none' disables it.  [default: none]
  --solvent-model [alpb|cpcmx]    xTB solvent model.  [default: alpb]
  --coord-type [cart|dlc]         Optimization coordinate system (cart|dlc).
                                  [default: (cart)]
  --precision [fp32|fp64]         MLIP backend precision: fp32 or fp64. Unset
                                  defaults per backend (uma: fp32; orb, mace:
                                  fp64). Routed to backend-specific kwargs (UMA
                                  precision / ORB precision / MACE
                                  default_dtype). aimnet2: fp32 no-op; fp64
                                  rejected.  [default: (per backend: uma fp32;
                                  orb, mace fp64)]
  --backend-model TEXT            Model variant for the selected --backend (e.g.
                                  uma-s-1p2 / uma-m-1p1 for uma,
                                  orb_v3_conservative_omol for orb, MACE-OMOL-0
                                  / off:small for mace).  [default: (the
                                  selected backend's own model)]
  --calc-file FILE                Python file exposing get_calculator(...) -> an
                                  ASE Calculator, used as the energy/gradient
                                  backend (overrides --backend). Couples GFN-xTB
                                  / DFTB+ / any ASE engine. See --calc-file-
                                  func-name.
  --calc-file-func-name TEXT      Name of the callable in --calc-file that
                                  returns an ASE Calculator (or a module-level
                                  Calculator instance). CLI overrides config
                                  YAML; otherwise defaults to get_calculator.
                                  [default: (get_calculator)]
  --deterministic / --no-deterministic
                                  Request strict same-stack PyTorch determinism
                                  (deterministic algorithms + index_reduce_
                                  shim). Slower; raises for detected unsupported
                                  ops; custom calculators are outside its scope.
                                  [default: no-deterministic]
  --allow-charge-mult-mismatch    Skip the cluster charge/multiplicity electron-
                                  parity check (logs that it was skipped). Open-
                                  shell clusters need a matching multiplicity
                                  instead; use this only for an intentionally
                                  nonstandard electron count.
  -h, --help                      Show this message and exit.
```
