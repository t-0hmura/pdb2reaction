# pdb2reaction dft

```
pdb2reaction ver. 0.2.1.dev116+gff86b942a

Usage: cli dft [OPTIONS]

  Single-point DFT using GPU4PySCF (CPU PySCF backend).

Options:
  --help-advanced             Show all options (including advanced settings) and
                              exit.
  -i, --input FILE            Input structure file (.pdb, .xyz, _trj.xyz, etc.;
                              loaded via pysisyphus.helpers.geom_loader).
                              [required]
  -q, --charge INTEGER        Total charge. Required for non-.gjf inputs unless
                              --ligand-charge is provided (PDB inputs or XYZ/GJF
                              with --ref-pdb).
  -l, --ligand-charge TEXT    Total charge or per-resname mapping (e.g.,
                              GPP:-3,SAM:1) used to derive charge when -q is
                              omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER  Spin multiplicity (2S+1; inherits from .gjf when
                              available; otherwise defaults to 1).
  --func-basis TEXT           Exchange–correlation functional and basis set as
                              'FUNC/BASIS' (e.g., 'wb97m-v/6-31g**',
                              'wb97m-v/def2-tzvpd').  [default:
                              wb97m-v/def2-tzvpd]
  -o, --out-dir TEXT          Output directory.  [default: ./result_dft/]
  --engine [gpu|cpu|auto]     Preferred SCF backend: GPU (strict), CPU, or auto
                              (try GPU then CPU if unavailable).  [default: gpu]
  --config FILE               Base YAML configuration file applied before
                              explicit CLI options.
  -h, --help                  Show this message and exit.
```
