# `dft`

## Overview

> **Summary:** Runs single-point DFT with GPU4PySCF or CPU PySCF. The default functional/basis is ωB97M-V/def2-TZVPD. Results include energies and population analysis (Mulliken, meta-Löwdin, IAO charges).

`pdb2reaction dft` runs single-point DFT calculations using PySCF (CPU) or GPU4PySCF (GPU). The default functional/basis is ωB97M-V/def2-TZVPD. Results include energy and population analysis (Mulliken, meta-Löwdin, IAO charges).

The backend is controlled by `--engine`:
- `gpu` (default): Uses GPU4PySCF. **Raises an error if GPU is unavailable.**
- `cpu`: Forces CPU PySCF.
- `auto` (recommended for portability): Attempts GPU4PySCF first, falls back to CPU PySCF if GPU is unavailable.

In addition to total energies, the command reports Mulliken, meta-Löwdin, and IAO atomic charges and spin densities.

## Minimal example

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 --engine auto --out-dir ./result_dft
```

## Output checklist

- `result_dft/summary.md`
- `result_dft/key_input_geometry.xyz`, `result_dft/key_result.yaml` (symlink/copy shortcuts)
- `result_dft/input_geometry.xyz`
- `result_dft/result.yaml`
- Engine metadata (`gpu4pyscf` / `pyscf(cpu)`) in `result.yaml`

## Common examples

1. Run with a larger basis and tighter SCF settings.

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 \
  --func-basis 'wb97m-v/def2-tzvpd' --conv-tol 1e-10 --max-cycle 200 \
  --engine auto --out-dir ./result_dft_tight
```

2. Force CPU backend for portability.

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 --engine cpu --out-dir ./result_dft_cpu
```

3. Derive total charge from ligand mapping when `-q` is omitted.

```bash
pdb2reaction dft -i input.pdb --ligand-charge 'LIG:0' -m 1 \
  --engine auto --out-dir ./result_dft_ligand
```

## Usage
```bash
pdb2reaction dft -i INPUT.{pdb|xyz|gjf|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m MULTIPLICITY] \
                 [--func-basis 'FUNC/BASIS'] \
                 [--max-cycle N] [--conv-tol Eh] [--grid-level L] \
                 [--out-dir DIR] [--engine gpu|cpu|auto] [--convert-files/--no-convert-files] \
                 [--ref-pdb FILE] [--config FILE] [--override-yaml FILE] \
                 [--show-config] [--dry-run] [--args-yaml FILE]
```

### Examples
```bash
# Default GPU-first policy with explicit functional/basis
pdb2reaction dft -i input.pdb -q 0 -m 1 --func-basis 'wb97m-v/6-31g**'

# Tighter controls, larger basis, CPU-only backend
pdb2reaction dft -i input.pdb -q 1 -m 2 --func-basis 'wb97m-v/def2-tzvpd' --max-cycle 150 --conv-tol 1e-9 --engine cpu
```

## Workflow
1. **Input handling** – Any file loadable by `geom_loader` (.pdb/.xyz/.trj/…) is accepted. Coordinates are re-exported as `input_geometry.xyz`. For XYZ/GJF inputs, `--ref-pdb` supplies a reference PDB topology for atom-count validation and (if you also use `--ligand-charge`) charge derivation; the DFT stage itself does **not** emit PDB/GJF outputs.
2. **Configuration merge** – Defaults → `--config` → explicit CLI options → `--override-yaml` (`dft` block). `--args-yaml` is retained as a legacy alias of `--override-yaml`. Charge/multiplicity inherit `.gjf` metadata when present. If `-q` is omitted but `--ligand-charge` is provided, the structure is treated as an enzyme–substrate complex and `extract.py`’s charge summary derives the total charge; explicit `-q` still overrides. For non-`.gjf` inputs, omitting `-q` without `--ligand-charge` aborts; multiplicity defaults to `1` when omitted.
3. **SCF build** – `--func-basis` is parsed into functional and basis. Density fitting is enabled automatically with PySCF defaults. `--engine` controls GPU/CPU preference (`gpu` requires GPU4PySCF; `cpu` forces CPU; `auto` tries GPU then CPU). Nonlocal corrections (e.g., VV10) are not configured explicitly beyond the backend defaults.
4. **Population analysis & outputs** – After convergence (or failure) the command writes `result.yaml` summarizing energy (Hartree/kcal·mol⁻¹), convergence metadata, timing, backend info, and per-atom Mulliken/meta-Löwdin/IAO charges and spin densities (UKS only for spins). Any failed analysis column is set to `null` with a warning.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge supplied to PySCF (`calc.charge`). Required unless a `.gjf` template or `--ligand-charge` (PDB inputs or XYZ/GJF with `--ref-pdb`) supplies it. Overrides `--ligand-charge` when both are set. | Required unless template/derivation applies |
| `--ligand-charge TEXT` | Total charge or per-resname mapping used when `-q` is omitted. Triggers extract-style charge derivation on the full complex (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). Converted to `2S` for PySCF. | `.gjf` template value or `1` |
| `--func-basis TEXT` | Functional/basis pair in `FUNC/BASIS` form (quote strings with `*`). | `wb97m-v/def2-tzvpd` |
| `--max-cycle INT` | Maximum SCF iterations (`dft.max_cycle`). | `100` |
| `--conv-tol FLOAT` | SCF convergence tolerance in Hartree (`dft.conv_tol`). | `1e-9` |
| `--grid-level INT` | PySCF numerical integration grid level (`dft.grid_level`). | `3` |
| `--out-dir TEXT` | Output directory (`dft.out_dir`). | `./result_dft/` |
| `--engine [gpu\|cpu\|auto]` | Backend policy: GPU4PySCF first, CPU only, or auto. | `gpu` |
| `--convert-files/--no-convert-files` | Accepted for interface consistency; no PDB/GJF outputs are produced by `dft`. | `True` |
| `--ref-pdb FILE` | Reference PDB topology to validate atom counts and enable ligand-charge derivation for XYZ/GJF inputs (no output conversion). | _None_ |
| `--config FILE` | Base YAML configuration file applied before explicit CLI options. | _None_ |
| `--override-yaml FILE` | Final YAML override file (highest priority YAML layer). | _None_ |
| `--args-yaml FILE` | _Legacy alias_ of `--override-yaml` (deprecated). | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration and continue execution. | `False` |
| `--dry-run/--no-dry-run` | Validate options and print execution plan without running DFT. | `False` |

## Outputs
```
out_dir/ (default: ./result_dft/)
├─ input_geometry.xyz   # Geometry snapshot sent to PySCF
├─ result.yaml          # Energy/charge/spin summaries with convergence/engine metadata
├─ summary.md           # Quick guide to key outputs
├─ key_input_geometry.xyz  # Shortcut to geometry snapshot (symlink/copy)
└─ key_result.yaml      # Shortcut to result summary (symlink/copy)
```
- `result.yaml` expands to:
  - `energy`: Hartree/kcal·mol⁻¹ values, convergence flag, wall time, engine metadata
    (`gpu4pyscf` vs `pyscf(cpu)`, `used_gpu`).
  - `charges`: Mulliken, meta-Löwdin, and IAO atomic charges (`null` when a method fails).
  - `spin_densities`: Mulliken, meta-Löwdin, and IAO spin densities (UKS-only for spins).
- It also summarizes charge, multiplicity, spin (2S), functional, basis,
  convergence knobs, and resolved output directory.

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- `--engine gpu` (default) requires GPU4PySCF and **raises an error** if a GPU is unavailable. Use `--engine auto` for automatic fallback to CPU PySCF when GPU resources are not detected, or `--engine cpu` to force CPU-only execution.
- If **Blackwell architecture** GPUs are detected, a warning is emitted because current GPU4PySCF may be unsupported.
- Compiled GPU4PySCF wheels may not support Blackwell-architecture GPUs, and non-x86 systems require compiling from source; we recommend using the CPU backend or building GPU4PySCF yourself in these situations. (see https://github.com/pyscf/gpu4pyscf)
- Density fitting is always attempted with PySCF defaults (no auxiliary basis guessing is implemented).
- The YAML input file must have a mapping root; the `dft` section is optional. Non-mapping roots raise an error via `load_yaml_dict`.
- IAO spin/charge analysis may fail for challenging systems; corresponding columns in `result.yaml` become `null` and a warning is printed.

## YAML configuration (`--config` / `--override-yaml` / `--args-yaml`)
Accepts a mapping root; the `dft` section (and optional `geom`) is applied when present. Merge order is:
- defaults
- `--config`
- explicit CLI options
- `--override-yaml` (or `--args-yaml`, legacy alias)

`dft` keys (defaults in parentheses):
- `func` (`"wb97m-v"`): Exchange–correlation functional.
- `basis` (`"def2-tzvpd"`): Basis set name.
- `func_basis` (_None_): Optional combined `FUNC/BASIS` string that overrides `func`/`basis` when provided.
- `conv_tol` (`1e-9`): SCF convergence threshold (Hartree).
- `max_cycle` (`100`): Maximum SCF iterations.
- `grid_level` (`3`): PySCF `grids.level`.
- `verbose` (`0`): PySCF verbosity (0–9). The CLI constructs the configuration with this quiet default unless overridden.
- `out_dir` (`"./result_dft/"`): Output directory root.

_Functional/basis selection defaults to `wb97m-v/def2-tzvpd` but can be overridden on the CLI. Charge/spin inherit `.gjf` template metadata when present. If `-q` is omitted but `--ligand-charge` is provided, the input is treated as an enzyme–substrate complex and `extract.py`’s charge summary computes the total charge; explicit `-q` still overrides. For non-`.gjf` inputs, omitting `-q` without `--ligand-charge` aborts; multiplicity defaults to `1` when omitted. Set them explicitly for non-default states._

```yaml
geom:
  coord_type: cart       # optional geom_loader settings
dft:
  func: wb97m-v         # exchange–correlation functional
  basis: def2-tzvpd     # basis set name (alternatively use func_basis: "FUNC/BASIS")
  conv_tol: 1.0e-09     # SCF convergence tolerance (Hartree)
  max_cycle: 100        # maximum SCF iterations
  grid_level: 3         # PySCF grid level
  verbose: 0            # PySCF verbosity (0-9)
  out_dir: ./result_dft/  # output directory root
```

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing

- [freq](freq.md) — UMA-based vibrational analysis (often precedes DFT refinement)
- [all](all.md) — End-to-end workflow with `--dft`
- [YAML Reference](yaml_reference.md) — Full `dft` configuration options
- [Glossary](glossary.md) — Definitions of DFT, SP (Single Point)
