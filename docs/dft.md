# `dft` subcommand

## Overview
Run single-point DFT calculations with a GPU-first policy (GPU4PySCF when available, CPU PySCF otherwise). In addition to total energies the command reports Mulliken, meta-Löwdin, and IAO atomic charges/spin densities so users can reuse the results downstream without reprocessing the PySCF objects.

## Usage
```bash
pdb2reaction dft -i INPUT.{pdb|xyz|gjf|...} -q CHARGE [-m 2S+1] \
                 --func-basis "FUNC/BASIS" \
                 [--max-cycle N] [--conv-tol Eh] [--grid-level L] \
                 [--out-dir DIR] [--engine gpu|cpu|auto] [--convert-files/--no-convert-files] \
                 [--args-yaml FILE]
```

### Examples
```bash
# Default GPU-first policy with explicit functional/basis
pdb2reaction dft -i input.pdb -q 0 -m 1 --func-basis "wb97m-v/6-31g**"

# Tighter controls, larger basis, CPU-only backend
pdb2reaction dft -i input.pdb -q 0 -m 2 --func-basis "wb97m-v/def2-tzvpd" \
                --max-cycle 150 --conv-tol 1e-9 --engine cpu
```

## Workflow
1. **Input handling** – Any file loadable by `geom_loader` (.pdb/.xyz/.trj/…) is accepted. Coordinates are re-exported as `input_geometry.xyz`, and when the source was a Gaussian template a matching `.gjf` snapshot is written for reference.
2. **Configuration merge** – Defaults → YAML (`dft` block) → CLI. Charge/multiplicity inherit `.gjf` metadata when present; otherwise `-q/--charge` is required and multiplicity defaults to `1`. Always set them explicitly so the correct RKS/UKS solver is selected.
3. **SCF build** – `--func-basis` is parsed into XC functional and orbital basis. Density fitting is enabled automatically and a practical JKFIT auxiliary basis is guessed from the orbital basis (def2, cc-pVXZ, Pople families). `--engine` controls GPU/CPU preference (`gpu` tries GPU4PySCF before falling back; `cpu` forces CPU; `auto` tries GPU then CPU). Nonlocal VV10 is activated for functionals whose names end in `-v` or contain `vv10`.
4. **Population analysis & outputs** – After convergence (or failure) the command writes `result.yaml` summarising energy (Hartree/kcal·mol⁻¹), convergence metadata, timing, backend info, and per-atom Mulliken/meta-Löwdin/IAO charges and spin densities (UKS only for spins). Any failed analysis column is set to `null` with a warning.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge supplied to PySCF (`calc.charge`). Required unless the input is a `.gjf` template that already stores charge. | Required when not in template |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). Converted to `2S` for PySCF. | `.gjf` template value or `1` |
| `--func-basis TEXT` | Functional/basis pair in `FUNC/BASIS` form (quote strings with `*`). | `wb97m-v/def2-tzvpd` |
| `--max-cycle INT` | Maximum SCF iterations (`dft.max_cycle`). | `100` |
| `--conv-tol FLOAT` | SCF convergence tolerance in Hartree (`dft.conv_tol`). | `1e-9` |
| `--grid-level INT` | PySCF numerical integration grid level (`dft.grid_level`). | `3` |
| `--out-dir TEXT` | Output directory (`dft.out_dir`). | `./result_dft/` |
| `--engine [gpu|cpu|auto]` | Backend policy: GPU4PySCF first, CPU only, or auto. | `gpu` |
| `--convert-files/--no-convert-files` | Toggle XYZ → PDB/GJF companions for PDB/Gaussian inputs. | `--convert-files` |
| `--args-yaml FILE` | YAML overrides (see below). | _None_ |

## Outputs
- `<out-dir>/input_geometry.xyz`: Geometry snapshot passed to PySCF (identical coordinates to the input file). `.gjf` is emitted when the input had Gaussian metadata and conversion is enabled.
- `<out-dir>/result.yaml`:
  - `energy` block with Hartree/kcal·mol⁻¹ values, convergence flag, wall time, and engine metadata (`gpu4pyscf` vs `pyscf(cpu)`, `used_gpu`).
  - `charges`: Mulliken, meta-Löwdin, and IAO atomic charges (IAO may be `null` if unavailable).
  - `spin_densities`: Mulliken, meta-Löwdin, and IAO atomic spin densities (restricted cases report zero/`null` as appropriate).
- Console pretty block summarising charge, multiplicity, spin (2S), functional, basis, convergence knobs, and resolved output directory.

## Notes
- GPU4PySCF is used whenever available; CPU PySCF is built otherwise (unless `--engine cpu` forces CPU). `--engine auto` mirrors the CLI docstring behaviour (GPU attempt, then CPU). **Meta-GGA functionals (e.g., `wb97m-v`, `scan`, `m06-2x`) are forced onto the CPU backend to avoid GPU4PySCF symbol errors such as `work_mgga` missing.**
- Density fitting is always attempted; auxiliary bases are auto-selected for def2/cc-pVXZ/Pople families and silently skipped when unsupported.
- The YAML file must contain a mapping root with top-level key `dft`; non-mapping roots raise an error via `load_yaml_dict`.
- Exit codes: `0` (converged), `3` (not converged), `2` (PySCF import failure), `1` (other errors), `130` (interrupt).
- IAO spin/charge analysis may fail for challenging systems; corresponding columns in `result.yaml` become `null` and a warning is printed.

## YAML configuration (`--args-yaml`)
Accepts a mapping with top-level key `dft`. CLI overrides YAML values.

`dft` keys (defaults in parentheses):
- `conv_tol` (`1e-9`): SCF convergence threshold (Hartree).
- `max_cycle` (`100`): Maximum SCF iterations.
- `grid_level` (`3`): PySCF `grids.level`.
- `verbose` (`4`): PySCF verbosity (0–9).
- `out_dir` (`"./result_dft/"`): Output directory root.

_Functional/basis selection defaults to `wb97m-v/def2-tzvpd` but can be overridden on the CLI. Charge/spin inherit `.gjf` template metadata when present; otherwise `-q/--charge` is required and spin defaults to `1`. Set them explicitly for non-default states._

```yaml
dft:
  conv_tol: 1.0e-09
  max_cycle: 100
  grid_level: 3
  verbose: 4
  out_dir: ./result_dft/
```
