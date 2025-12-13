# `dft` subcommand

## Overview
Run single-point DFT calculations with a GPU (GPU4PySCF when available, CPU PySCF otherwise). In addition to total energies the command reports Mulliken, meta-Löwdin, and IAO atomic charges/spin densities.

## Usage
```bash
pdb2reaction dft -i INPUT.{pdb|xyz|gjf|...} -q CHARGE [-m MULTIPLICITY] \
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
1. **Input handling** – Any file loadable by `geom_loader` (.pdb/.xyz/.trj/…) is accepted. Coordinates are re-exported as `input_geometry.xyz`.
2. **Configuration merge** – Defaults → CLI → YAML (`dft` block). YAML values take final precedence and can override `charge`, `multiplicity`, `func_basis`, and `engine` in addition to SCF knobs. Charge/multiplicity inherit `.gjf` metadata when present; otherwise `-q/--charge` is required and multiplicity defaults to `1`.
3. **SCF build** – `--func-basis` is parsed into functional and basis. Density fitting is enabled automatically with PySCF defaults. `--engine` controls GPU/CPU preference (`gpu` tries GPU4PySCF before falling back; `cpu` forces CPU; `auto` tries GPU then CPU). Nonlocal VV10 corrections are not configured explicitly.
4. **Population analysis & outputs** – After convergence (or failure) the command writes `result.yaml` summarising energy (Hartree/kcal·mol⁻¹), convergence metadata, backend info, and per-atom Mulliken/meta-Löwdin/IAO charges and spin densities (UKS only for spins). Any failed analysis column is set to `null` with a warning.

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
```
out_dir/ (default: ./result_dft/)
├─ input_geometry.xyz   # Geometry snapshot sent to PySCF
└─ result.yaml          # Energy/charge/spin summaries with convergence/engine metadata
```
- `result.yaml` expands to:
  - `energy`: Hartree/kcal·mol⁻¹ values, convergence flag, and engine metadata
    (`gpu4pyscf` vs `pyscf(cpu)`, `used_gpu`).
  - `charges`: Mulliken, meta-Löwdin, and IAO atomic charges (`null` when a method fails).
  - `spin_densities`: Mulliken, meta-Löwdin, and IAO spin densities (UKS-only for spins).
- Console pretty block summarising charge, multiplicity, spin (2S), functional, basis,
  convergence knobs, and resolved output directory.

## Notes
- GPU4PySCF is used whenever available; CPU PySCF is built otherwise (unless `--engine cpu` forces CPU). `--engine auto` mirrors the GPU-first fallback logic, automatically retrying on the CPU backend when GPU import/runtime errors occur. **Blackwell architecture** GPUs are detected and forced to CPU with a warning to avoid unsupported GPU4PySCF configurations.
- GPU4PySCF is required to by compiled from source when you do not use **x86** archtecture. (See <https://github.com/pyscf/gpu4pyscf>)
- Density fitting is always attempted with PySCF defaults (no auxiliary basis guessing is implemented).
- YAML overrides are read from a mapping root; `dft` entries override CLI/defaults (including `charge`, `multiplicity`, `func_basis`, and `engine`) and an optional `geom` block is passed to `geom_loader`.
- Exit codes: `0` (converged), `3` (not converged), `2` (PySCF import failure), `1` (other errors), `130` (interrupt).
- IAO spin/charge analysis may fail for challenging systems; corresponding columns in `result.yaml` become `null` and a warning is printed.

## YAML configuration (`--args-yaml`)
Accepts a mapping with optional top-level key `dft` (and `geom`). YAML values override CLI values.

`dft` keys (defaults in parentheses):
- `charge` (`null` → template or `-q` required): Total charge.
- `multiplicity` (`1`): Spin multiplicity (2S+1).
- `func_basis` (`"wb97m-v/def2-tzvpd"`): Functional/basis pair in `FUNC/BASIS` form.
- `engine` (`"gpu"`): Backend policy; accepts `gpu`, `cpu`, or `auto`.
- `conv_tol` (`1e-9`): SCF convergence threshold (Hartree).
- `max_cycle` (`100`): Maximum SCF iterations.
- `grid_level` (`3`): PySCF `grids.level`.
- `verbose` (`0`): PySCF verbosity (0–9). The CLI constructs the configuration with this quiet default unless overridden.
- `out_dir` (`"./result_dft/"`): Output directory root.

Charge/spin inherit `.gjf` template metadata when present; otherwise `-q/--charge` is required and spin defaults to `1`. YAML values override both CLI inputs and template metadata.

```yaml
dft:
  charge: 0             # total charge (overrides CLI/template)
  multiplicity: 1       # spin multiplicity (2S+1)
  func_basis: wb97m-v/def2-tzvpd  # functional/basis
  engine: gpu           # backend policy: gpu, cpu, or auto
  conv_tol: 1.0e-09     # SCF convergence tolerance (Hartree)
  max_cycle: 100        # maximum SCF iterations
  grid_level: 3         # PySCF grid level
  verbose: 0            # PySCF verbosity (0-9)
  out_dir: ./result_dft/  # output directory root
```
