# `pdb2reaction dft`

```text
Usage: pdb2reaction dft [OPTIONS]

  Single-point DFT using GPU4PySCF (CPU PySCF backend).

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
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
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1; inherits from .gjf
                                  when available; otherwise defaults to 1).
  --convert-files / --no-convert-files
                                  Accepted for interface consistency; dft does
                                  not emit PDB/GJF outputs.  [default: convert-
                                  files]
  --ref-pdb FILE                  Reference PDB topology to use when the input
                                  is XYZ/GJF (keeps XYZ coordinates).
  --func-basis TEXT               Exchange–correlation functional and basis set
                                  as 'FUNC/BASIS' (e.g., 'wb97m-v/6-31g**',
                                  'wb97m-v/def2-tzvpd').  [default:
                                  wb97m-v/def2-tzvpd]
  --max-cycle INTEGER             Maximum SCF iterations.  [default: 100]
  --conv-tol FLOAT                SCF convergence tolerance (Eh).  [default:
                                  1e-09]
  --grid-level INTEGER            Numerical integration grid level (PySCF
                                  grids.level).  [default: 3]
  -o, --out-dir TEXT              Output directory.  [default: ./result_dft/]
  --engine [gpu|cpu]              SCF backend: gpu (GPU4PySCF, raises error if
                                  unavailable) or cpu (PySCF).  [default: gpu]
  --lowmem / --no-lowmem          Use gpu4pyscf rks_lowmem.RKS for closed-shell
                                  GPU runs (memory-efficient direct JK; mutually
                                  exclusive with density fitting). Open-shell,
                                  CPU, or pre-rks_lowmem GPU4PySCF installs
                                  auto-fall back to standard RKS/UKS with
                                  density fitting.  [default: lowmem]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running DFT.  [default: no-dry-run]
  -h, --help                      Show this message and exit.
```
