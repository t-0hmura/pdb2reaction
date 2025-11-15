# `dft` subcommand

## Purpose
Runs single-point DFT calculations using GPU4PySCF when available (falling back to CPU PySCF), reporting total energies and atomic charges.

## Usage
```bash
pdb2reaction dft -i INPUT -q CHARGE [-s SPIN] \
                 [--func-basis "FUNC/BASIS"] \
                 [--max-cycle N] [--conv-tol Eh] [--grid-level L] \
                 [--out-dir DIR] [--args-yaml FILE]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge supplied to PySCF. | Required |
| `-s, --spin INT` | Spin multiplicity (2S+1). Converted to `2S` for PySCF. | `1` |
| `--func-basis TEXT` | Functional and basis in `FUNC/BASIS` form (quotes recommended when using `*`). | `wb97m-v/6-31g**` |
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

_Functional/basis selection and molecular charge must be supplied on the CLI. Spin defaults to `1` (singlet) but should be set explicitly for other states._

## Outputs
- `<out-dir>/input_geometry.xyz`: Geometry snapshot passed to PySCF (identical coordinates to the input file).
- `<out-dir>/result.yaml`:
  - `energy` block with Hartree/kcal·mol⁻¹ values, convergence flag, wall time, and engine metadata (`gpu4pyscf` vs `pyscf(cpu)`, `used_gpu`).
  - `charges`: Mulliken, meta-Löwdin, and IAO atomic charges (IAO may be `null` if unavailable).
  - `spin_densities`: Mulliken, meta-Löwdin, and IAO atomic spin densities (restricted cases report zero/`null` as appropriate).
- Console pretty block summarising charge, multiplicity, spin (2S), functional, basis, convergence knobs, and resolved output directory.

## Notes
- GPU4PySCF is used when available; otherwise a CPU SCF object is built. Nonlocal VV10 is enabled automatically for functionals ending with `-v` or containing `vv10`.
- The YAML file must contain a mapping root with top-level key `dft`; non-mapping roots raise an error via `load_yaml_dict`.
- `-q/--charge` is required. `-s/--spin` defaults to `1` (RKS) but should be specified explicitly for non-singlet states to ensure the correct (RKS vs UKS) solver.
- Exit codes: `0` (converged), `3` (not converged), `2` (PySCF import failure), `1` (other errors), `130` (interrupt).
- IAO spin/charge analysis may fail for challenging systems; in that case the IAO column values in `result.yaml` are `null` and a warning is printed.

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