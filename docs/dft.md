# `dft` subcommand

## Overview
Run single-point DFT calculations with a GPU (GPU4PySCF when available, CPU PySCF otherwise). In addition to total energies, the command reports Mulliken, meta-Löwdin, and IAO atomic charges/spin densities.

## Usage
```bash
pdb2reaction dft -i INPUT.{pdb|xyz|gjf|...} -q CHARGE [-m MULTIPLICITY] \
                 --func-basis 'FUNC/BASIS' \
                 [--max-cycle N] [--conv-tol Eh] [--grid-level L] \
                 [--out-dir DIR] [--engine gpu|cpu|auto] [--convert-files {True|False}] \
                 [--args-yaml FILE]
```

### Examples
```bash
# Default GPU-first policy with explicit functional/basis
pdb2reaction dft -i input.pdb -q 0 -m 1 --func-basis 'wb97m-v/6-31g**'

# Tighter controls, larger basis, CPU-only backend
pdb2reaction dft -i input.pdb -q 1 -m 2 --func-basis 'wb97m-v/def2-tzvpd' --max-cycle 150 --conv-tol 1e-9 --engine cpu
```

## Workflow
1. **Input handling** – Any file loadable by `geom_loader` (.pdb/.xyz/.trj/…) is accepted. Coordinates are re-exported as `input_geometry.xyz`.
2. **Configuration merge** – Defaults → CLI → YAML (`dft` block). YAML overrides take precedence over CLI flags. Charge/multiplicity inherit `.gjf` metadata when present. If `-q` is omitted but `--ligand-charge` is provided, the structure is treated as an enzyme–substrate complex and `extract.py`’s charge summary derives the total charge; explicit `-q` still overrides. Otherwise charge defaults to `0` and multiplicity to `1`.
3. **SCF build** – `--func-basis` is parsed into functional and basis. Density fitting is enabled automatically with PySCF defaults. `--engine` controls GPU/CPU preference (`gpu` tries GPU4PySCF before falling back; `cpu` forces CPU; `auto` tries GPU then CPU). Nonlocal corrections (e.g., VV10) are not configured explicitly beyond the backend defaults.
4. **Population analysis & outputs** – After convergence (or failure) the command writes `result.yaml` summarizing energy (Hartree/kcal·mol⁻¹), convergence metadata, timing, backend info, and per-atom Mulliken/meta-Löwdin/IAO charges and spin densities (UKS only for spins). Any failed analysis column is set to `null` with a warning.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge supplied to PySCF (`calc.charge`). Required unless the input is a `.gjf` template that already stores charge. Overrides `--ligand-charge` when both are set. | Required when not in template |
| `--ligand-charge TEXT` | Total charge or per-resname mapping used when `-q` is omitted. Triggers extract-style charge derivation on the full complex. | `None` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). Converted to `2S` for PySCF. | `.gjf` template value or `1` |
| `--func-basis TEXT` | Functional/basis pair in `FUNC/BASIS` form (quote strings with `*`). | `wb97m-v/def2-tzvpd` |
| `--max-cycle INT` | Maximum SCF iterations (`dft.max_cycle`). | `100` |
| `--conv-tol FLOAT` | SCF convergence tolerance in Hartree (`dft.conv_tol`). | `1e-9` |
| `--grid-level INT` | PySCF numerical integration grid level (`dft.grid_level`). | `3` |
| `--out-dir TEXT` | Output directory (`dft.out_dir`). | `./result_dft/` |
| `--engine [gpu\|cpu\|auto]` | Backend policy: GPU4PySCF first, CPU only, or auto. | `gpu` |
| `--convert-files {True|False}` | Toggle XYZ → PDB/GJF companions for PDB/Gaussian inputs. | `True` |
| `--args-yaml FILE` | YAML overrides (see below). | _None_ |

## Outputs
```
out_dir/ (default: ./result_dft/)
├─ input_geometry.xyz   # Geometry snapshot sent to PySCF
└─ result.yaml          # Energy/charge/spin summaries with convergence/engine metadata
```
- `result.yaml` expands to:
  - `energy`: Hartree/kcal·mol⁻¹ values, convergence flag, wall time, engine metadata
    (`gpu4pyscf` vs `pyscf(cpu)`, `used_gpu`).
  - `charges`: Mulliken, meta-Löwdin, and IAO atomic charges (`null` when a method fails).
  - `spin_densities`: Mulliken, meta-Löwdin, and IAO spin densities (UKS-only for spins).
- It also summarizes charge, multiplicity, spin (2S), functional, basis,
  convergence knobs, and resolved output directory.

## Notes
- GPU4PySCF is used whenever available; CPU PySCF is built otherwise (unless `--engine cpu` forces CPU). `--engine auto` mirrors the GPU-first fallback logic, automatically retrying on the CPU backend when GPU import/runtime errors occur.
- If **Blackwell architecture** GPUs are detected, a warning is emitted because current GPU4PySCF may be unsupported.
- Compiled GPU4PySCF wheels may not support Blackwell-architecture GPUs, and non-x86 systems require compiling from source; we recommend using the CPU backend or building GPU4PySCF yourself in these situations. (see https://github.com/pyscf/gpu4pyscf)
- Density fitting is always attempted with PySCF defaults (no auxiliary basis guessing is implemented).
- The YAML input file must contain a mapping root with top-level key `dft`; non-mapping roots raise an error via `load_yaml_dict`.
- IAO spin/charge analysis may fail for challenging systems; corresponding columns in `result.yaml` become `null` and a warning is printed.

## YAML configuration (`--args-yaml`)
Accepts a mapping with top-level key `dft`. YAML values override CLI values.

`dft` keys (defaults in parentheses):
- `func` (`"wb97m-v"`): Exchange–correlation functional.
- `basis` (`"def2-tzvpd"`): Basis set name.
- `func_basis` (_None_): Optional combined `FUNC/BASIS` string that overrides `func`/`basis` when provided.
- `conv_tol` (`1e-9`): SCF convergence threshold (Hartree).
- `max_cycle` (`100`): Maximum SCF iterations.
- `grid_level` (`3`): PySCF `grids.level`.
- `verbose` (`0`): PySCF verbosity (0–9). The CLI constructs the configuration with this quiet default unless overridden.
- `out_dir` (`"./result_dft/"`): Output directory root.

_Functional/basis selection defaults to `wb97m-v/def2-tzvpd` but can be overridden on the CLI. Charge/spin inherit `.gjf` template metadata when present. If `-q` is omitted but `--ligand-charge` is provided, the input is treated as an enzyme–substrate complex and `extract.py`’s charge summary computes the total charge; explicit `-q` still overrides. Otherwise charge defaults to `0` and spin to `1`. Set them explicitly for non-default states._

```yaml
dft:
  func: wb97m-v         # exchange–correlation functional
  basis: def2-tzvpd     # basis set name (alternatively use func_basis: "FUNC/BASIS")
  conv_tol: 1.0e-09     # SCF convergence tolerance (Hartree)
  max_cycle: 100        # maximum SCF iterations
  grid_level: 3         # PySCF grid level
  verbose: 0            # PySCF verbosity (0-9)
  out_dir: ./result_dft/  # output directory root
```
