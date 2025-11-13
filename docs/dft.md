# `dft` subcommand

## Purpose
Runs single-point DFT calculations using GPU4PySCF when available (falling back to CPU PySCF), reporting total energies and atomic charges.

## Usage
```bash
pdb2reaction dft -i INPUT -q CHARGE -s SPIN
                 [--func-basis FUNC/BASIS]
                 [--max-cycle N] [--conv-tol Eh] [--grid-level L]
                 [--out-dir DIR] [--args-yaml FILE]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge supplied to PySCF. | Required |
| `-s, --spin INT` | Spin multiplicity (2S+1). Converted to `2S` for PySCF. | Required |
| `--func-basis TEXT` | Functional and basis in `FUNC/BASIS` form. | `wb97m-v/6-31g**` |
| `--max-cycle INT` | Maximum SCF iterations (`dft.max_cycle`). | `100` |
| `--conv-tol FLOAT` | SCF convergence tolerance in Hartree (`dft.conv_tol`). | `1e-9` |
| `--grid-level INT` | PySCF numerical integration grid level (`dft.grid_level`). | `3` |
| `--out-dir TEXT` | Output directory. | `./result_dft/` |
| `--args-yaml FILE` | YAML overrides (see below). | _None_ |

### Section `dft`
- `conv_tol` (`1e-9`): SCF convergence threshold (Hartree).
- `max_cycle` (`100`): Maximum SCF iterations.
- `grid_level` (`3`): PySCF `grids.level`.
- `verbose` (`4`): PySCF verbosity (0–9).
- `out_dir` (`"./result_dft/"`): Output directory.

_Functional/basis selection and molecular charge/spin must be supplied on the CLI._

## Outputs
- `<out-dir>/input_geometry.xyz`: Geometry snapshot passed to PySCF.
- `<out-dir>/result.yaml`: Total energy (Hartree and kcal·mol⁻¹), convergence status, engine metadata, SCF wall time, and Mulliken/Löwdin/IAO charges.
- Console pretty block summarising charge, multiplicity, functional, basis, convergence settings, and output directory.

## Notes
- GPU4PySCF is used when available; otherwise a CPU SCF object is built. Nonlocal VV10 is enabled automatically for `-v` functionals.
- The YAML file must contain a mapping root; non-mapping roots raise an error via `load_yaml_dict`.
- Exit codes: `0` (converged), `3` (not converged), `2` (PySCF import failure), `1` (other errors), `130` (interrupt).

## YAML configuration (`--args-yaml`)
Accepts a mapping with top-level key `dft`. CLI overrides YAML values.

```yaml
dft:
  conv_tol: 1.0e-09
  max_cycle: 100
  grid_level: 3
  verbose: 4
  out_dir: ./result_dft/
```