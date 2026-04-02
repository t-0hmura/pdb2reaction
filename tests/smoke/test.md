## Test Data and Smoke Tests (Developer-Oriented)

The repository includes `tests/smoke/` with sample inputs and a convenience
script for quick smoke checks:

```
tests/smoke/
├── run.sh          # Batch of CLI smoke tests (redirects output to test*.out)
├── *.pdb/*.xyz/*.gjf
├── h2.gjf          # Minimal H2 molecule for DFT test
└── input.yaml      # Example YAML overrides
```

### Running the smoke tests

From the repository root:

```bash
bash tests/smoke/run.sh
```

Or from inside the `tests/smoke/` directory:

```bash
cd tests/smoke
bash run.sh
```

### Speed optimizations

All iterative tests use `--thresh gau_loose` and reduced `--max-cycles` (3-5)
to keep runtime short. Scan tests use `--preopt False --endopt False`, and
path commands use `--preopt False --climb False --max-nodes 5`.

### Test coverage (35 tests)

| # | Category | Command | Notes |
|---|----------|---------|-------|
| 1 | Subcommand | `opt` (grad/lbfgs) | `--dump True` for trj2fig |
| 2 | Subcommand | `opt` (hess/rfo) | |
| 3 | Subcommand | `tsopt` (hess) | |
| 4 | Subcommand | `freq` | |
| 5 | Subcommand | `irc` | max-cycles 3 |
| 6 | Subcommand | `dft` (hf/sto-3g) | `--engine cpu` |
| 7 | Subcommand | `scan` (1D) | preopt/endopt False |
| 8 | Subcommand | `scan2d` | extract model first |
| 9 | Subcommand | `scan3d` | |
| 10 | Subcommand | `path-opt` (gsm) | max-nodes 5 |
| 11 | Subcommand | `path-search` | max-nodes 5 |
| 12 | Input format | `all` (pdb+pdb) | |
| 13 | Input format | `all` (xyz+xyz, `--ref-pdb`) | XYZ with PDB topology |
| 14 | Input format | `all` (gjf+gjf) | |
| 15 | Input format | `all` (scan-lists, pdb) | |
| 16 | Input format | `all` (scan-lists, xyz) | |
| 17 | Complex | `opt` (complex, ligand-charge) | |
| 18 | Complex | `all` (complex, extract+scan) | |
| 19 | Complex | `all` (complex, multi-input) | |
| 20 | TSOPT-only | `all` (ts.pdb, `--tsopt --thermo`) | |
| 21 | TSOPT-only | `all` (ts.pdb, `--opt-mode-post hess`) | |
| 22 | MEP mode | `all` (pdb+pdb, `--mep-mode dmf`) | |
| 23 | Complex TS | `tsopt` (complex, hess) | |
| 24 | Complex TS | `tsopt` (complex, grad) | |
| 25-31 | Dry-run | opt/tsopt/freq/scan/dft/path-search/irc | `--dry-run` |
| 32 | Utility | `extract` | |
| 33 | Utility | `add-elem-info` | |
| 34 | Utility | `trj2fig` | conditional on test1 output |
| 35 | Utility | `energy-diagram` | |

Notes:
- Outputs (`test1/`, `test2/`, ...) and log files (`test1.out`, ...) are created in
  `tests/smoke/` (the script changes into its own directory).
- Some tests use GPU UMA models; make sure CUDA, PyTorch, and UMA are available
  and that you have authenticated to Hugging Face (see Getting Started).
- The DFT test uses `--engine cpu` to avoid GPU4PySCF requirements.
