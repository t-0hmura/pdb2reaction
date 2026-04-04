# `pdb2reaction scan`

```text

Usage: pdb2reaction scan [OPTIONS]

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
  --workers INTEGER               MLIP predictor workers; >1 spawns a parallel
                                  predictor (Hessian computation not supported
                                  with workers>1).  [default: 1]
  --workers-per-node INTEGER      Workers per node when using a parallel MLIP
                                  predictor (workers>1).  [default: 1]
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1).
  --one-based / --zero-based      Interpret (i,j) indices in --scan-lists as
                                  1-based (default) or 0-based.  [default: one-
                                  based]
  --max-step-size FLOAT           Maximum change in any scanned bond length per
                                  step [Å].  [default: 0.2]
  --bias-k FLOAT                  Harmonic well strength k [eV/Å^2].  [default:
                                  300.0]
  --relax-max-cycles INTEGER      Maximum optimizer cycles per grid relaxation.
                                  When explicitly provided, used unless YAML
                                  sets opt.max_cycles.  [default: 10000]
  --opt-mode [grad|hess]          Relaxation mode: grad (=LBFGS) or hess (=RFO).
                                  [default: grad]
  --freeze-links / --no-freeze-links
                                  Freeze parent atoms of link hydrogens (PDB
                                  input or XYZ/GJF with --ref-pdb).  [default:
                                  freeze-links]
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., '1,3,5').
  --dump / --no-dump              Write per-step optimizer trajectory files.
                                  scan_trj.xyz and scan.pdb are always written
                                  to out-dir regardless of this flag.  [default:
                                  no-dump]
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB/GJF
                                  companions based on the input format.
                                  [default: convert-files]
  --ref-pdb FILE                  Reference PDB topology to use when the input
                                  is XYZ/GJF (keeps XYZ coordinates).
  -o, --out-dir TEXT              Base output directory.  [default:
                                  ./result_scan/]
  --thresh TEXT                   Convergence preset (gau_loose|gau|gau_tight|ga
                                  u_vtight|baker|never).
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --preopt / --no-preopt          Pre-optimize the initial structure without
                                  bias before the scan.  [default: no-preopt]
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --solvent TEXT                  Implicit solvent name for xTB correction (e.g.
                                  'water'). 'none' to disable.  [default: none]
  --solvent-model [alpb|cpcmx]    xTB solvent model.  [default: alpb]
  --print-parsed / --no-print-parsed
                                  Print parsed scan targets after resolving
                                  --scan-lists.  [default: no-print-parsed]
  --endopt / --no-endopt          After each stage, run an additional unbiased
                                  optimization of the stage result.  [default:
                                  no-endopt]
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
  -h, --help                      Show this message and exit.
```
