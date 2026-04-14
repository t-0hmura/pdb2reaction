# `dft`

## Overview

> **Summary:** Runs single-point DFT with GPU4PySCF or CPU PySCF. The default functional/basis is ωB97M-V/def2-tzvpd. Results include energy and population analysis (Mulliken, meta-Löwdin, IAO charges).

`pdb2reaction dft` runs single-point DFT calculations using PySCF (CPU) or GPU4PySCF (GPU). The default functional/basis is ωB97M-V/def2-tzvpd. Results include energy and population analysis (Mulliken, meta-Löwdin, IAO charges).

The backend is controlled by `--engine`:
- `gpu` (default): Uses GPU4PySCF. **Raises an error if GPU is unavailable.** Best for production runs on GPU-equipped nodes where you want to guarantee GPU acceleration.
- `cpu`: Forces CPU PySCF. Use when no GPU is available or when you need deterministic CPU-only execution (e.g., portability or debugging).
> **Note:** The default basis `def2-tzvpd` is a triple-zeta diffuse-augmented set and is computationally expensive for large systems. Consider a smaller basis (e.g., `6-31g**` or `def2-svp`) for exploratory calculations.

> **Prerequisites:** DFT dependencies (PySCF, GPU4PySCF) are **not** included in the default install. Install them with `pip install "pdb2reaction[dft]"`.

> **System size limit:** DFT single-point calculations are practical only for systems up to **~500 atoms**. Larger systems require excessive compute time and memory. For such systems, HPC clusters with high-end GPUs (e.g. A100, H200) are typically required. For enzyme systems, extract a small active site model (binding pocket) before running DFT.

In addition to total energies, the command reports Mulliken, meta-Löwdin, and IAO atomic charges and spin densities.

## Minimal example

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 --engine gpu --out-dir ./result_dft
```

## Output checklist

- `result_dft/input_geometry.xyz`
- `result_dft/result.yaml`
- Engine metadata (`gpu4pyscf` / `pyscf(cpu)`) in `result.yaml`

## Common examples

1. Run with a larger basis and tighter SCF settings.

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 \
 --func-basis 'wb97m-v/def2-tzvpd' --conv-tol 1e-10 --max-cycle 200 \
 --engine gpu --out-dir ./result_dft_tight
```

> **Caveat:** The tight `def2-tzvpd` setting above is appropriate only for **small** active-site models (≲150 atoms) on a GPU with sufficient VRAM. For larger systems on 16–24 GB GPUs this combination will OOM (see the "Out-of-memory with def2-TZVPD" note below); switch to `--func-basis 'wb97m-v/def2-svp'` or use an external DFT program (ORCA, Gaussian) for production work on full systems.

2. Force CPU backend for portability.

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 --engine cpu --out-dir ./result_dft_cpu
```

3. Derive total charge from ligand mapping when `-q` is omitted.

```bash
pdb2reaction dft -i input.pdb -l 'LIG:0' -m 1 \
 --engine gpu --out-dir ./result_dft_ligand
```

## Usage
```bash
pdb2reaction dft -i INPUT.{pdb|xyz|gjf|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULTIPLICITY] \
 [--func-basis 'FUNC/BASIS'] \
 [--max-cycle N] [--conv-tol Eh] [--grid-level L] \
 [--out-dir DIR] [--engine gpu|cpu] [--convert-files/--no-convert-files] \
 [--ref-pdb FILE] [--config FILE] [--show-config] [--dry-run]
```

### Examples
```bash
# Default GPU-first policy with explicit functional/basis
pdb2reaction dft -i input.pdb -q 0 -m 1 --func-basis 'wb97m-v/6-31g**'

# Tighter controls, larger basis, CPU-only backend
pdb2reaction dft -i input.pdb -q 1 -m 2 \
 --func-basis 'wb97m-v/def2-tzvpd' --max-cycle 150 \
 --conv-tol 1e-9 --engine cpu
```

## Workflow
1. **Input handling** – Any file loadable by `geom_loader` (.pdb/.xyz/_trj.xyz/…) is accepted. Coordinates are re-exported as `input_geometry.xyz`. For XYZ/GJF inputs, `--ref-pdb` supplies a reference PDB topology for atom-count validation and (if you also use `--ligand-charge`) charge derivation; the DFT stage itself does **not** emit PDB/GJF outputs.
2. **SCF build** – `--func-basis` is parsed into functional and basis. Density fitting is enabled automatically with PySCF defaults. `--engine` controls GPU/CPU preference (`gpu` requires GPU4PySCF and raises an error if unavailable; `cpu` forces CPU). Nonlocal corrections (e.g., VV10) are not configured explicitly beyond the backend defaults.
3. **Population analysis & outputs** – After convergence (or failure) the command writes `result.yaml` summarizing the energy (in hartree and kcal/mol), convergence metadata, timing, backend info, and per-atom Mulliken/meta-Löwdin/IAO charges and spin densities (UKS only for spins). Any failed analysis column is set to `null` with a warning.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge supplied to PySCF (`calc.charge`). Required unless a `.gjf` template or `--ligand-charge` (PDB inputs or XYZ/GJF with `--ref-pdb`) supplies it. Overrides `--ligand-charge` when both are set. | Required unless template/derivation applies |
| `-l, --ligand-charge TEXT` | Per-residue charge mapping (e.g., `GPP:-3,SAM:1`). Automatically derives the total system charge from PDB residue charges — no manual counting needed. Used when `-q` is omitted (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). Converted to `2S` for PySCF. | `.gjf` template value or `1` |
| `--func-basis TEXT` | Functional/basis pair in `FUNC/BASIS` form (quote strings with `*`). | `wb97m-v/def2-tzvpd` |
| `--max-cycle INT` | Maximum SCF iterations (`dft.max_cycle`). | `100` |
| `--conv-tol FLOAT` | SCF convergence tolerance in hartree (`dft.conv_tol`). | `1e-9` |
| `--grid-level INT` | PySCF numerical integration grid level (`dft.grid_level`). | `3` |
| `-o, --out-dir TEXT` | Output directory (`dft.out_dir`). | `./result_dft/` |
| `--engine [gpu\|cpu]` | SCF backend: gpu (GPU4PySCF) or cpu (PySCF). On the standalone `dft` subcommand this option is named `--engine`; the same knob is named `--dft-engine` when forwarded through `pdb2reaction all`. | `gpu` |
| `--convert-files/--no-convert-files` | **No-op on `dft`.** Accepted purely for interface consistency with the other subcommands; `dft` never produces PDB or GJF outputs (only `input_geometry.xyz` + `result.yaml`). The flag's value is ignored. | `True` |
| `--ref-pdb FILE` | Reference PDB topology to validate atom counts and enable ligand-charge derivation for XYZ/GJF inputs (no output conversion). | _None_ |
| `--config FILE` | Base YAML configuration file applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration and continue execution. | `False` |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` to `out_dir`. See [JSON Output Schema](json-output.md) for the schema. | `False` |
| `--dry-run/--no-dry-run` | Validate options and print execution plan without running DFT. | `False` |

## Outputs
```
out_dir/ (default:./result_dft/)
├─ input_geometry.xyz # Geometry snapshot sent to PySCF
├─ result.yaml # Energy/charge/spin summaries with convergence/engine metadata
```
- `result.yaml` expands to:
 - `energy`: energy in hartree and kcal/mol, convergence flag, wall time, engine metadata
  (`gpu4pyscf` vs `pyscf(cpu)`, `used_gpu`).
 - `charges`: Mulliken, meta-Löwdin, and IAO atomic charges (`null` when a method fails).
 - `spin_densities`: Mulliken, meta-Löwdin, and IAO spin densities (UKS-only for spins).
- It also summarizes charge, multiplicity, spin (2S), functional, basis,
  convergence knobs, and resolved output directory.

## Exit codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Unexpected error |
| 2 | Import failure (PySCF or GPU4PySCF not installed) |
| 3 | SCF not converged |
| 130 | Keyboard interrupt |

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes-common-errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- `--engine gpu` (default) requires GPU4PySCF and **raises an error** if GPU is unavailable. Use `--engine cpu` to force CPU-only execution.
- **Blackwell-architecture GPUs** (RTX 50xx): GPU4PySCF may fail with out-of-memory errors even for small systems (~100 atoms). Use `--engine cpu` or an external DFT program (ORCA, Gaussian) for production calculations on these GPUs. Note that `--engine cpu` is only practical for small active-site models (≲150 atoms) and small basis sets (e.g. `def2-svp`); larger systems on CPU become prohibitively slow, so an external DFT program is the recommended path for full systems.
- **Out-of-memory with def2-TZVPD**: The default basis set `def2-tzvpd` is large and may cause OOM for systems with >150 atoms on 16–24 GB GPUs. Use `--func-basis 'wb97m-v/def2-svp'` as a practical alternative; barrier height errors between def2-SVP and def2-TZVPD are typically 1–3 kcal/mol.
- Compiled GPU4PySCF wheels may not support non-x86 systems; build from source in that case (see https://github.com/pyscf/gpu4pyscf).
- Density fitting is always attempted with PySCF defaults (no auxiliary basis guessing is implemented).
- The YAML input file must have a mapping root; the `dft` section is optional. Non-mapping roots raise an error via `load_yaml_dict`.
- IAO spin/charge analysis may fail for challenging systems; corresponding columns in `result.yaml` become `null` and a warning is printed.

Accepts a mapping root; the `dft` section (and optional `geom`) is applied when present. Merge order is:
- defaults
- `--config`
- explicit CLI options

`dft` keys (defaults in parentheses):
- `func` (`"wb97m-v"`): Exchange-correlation functional.
- `basis` (`"def2-tzvpd"`): Basis set name.
- `func_basis` (_None_): Optional combined `FUNC/BASIS` string that overrides `func`/`basis` when provided.
- `conv_tol` (`1e-9`): SCF convergence threshold (hartree).
- `max_cycle` (`100`): Maximum SCF iterations.
- `grid_level` (`3`): PySCF `grids.level`.
- `verbose` (`0`): PySCF verbosity (0–9). The CLI constructs the configuration with this quiet default unless overridden.
- `out_dir` (`"./result_dft/"`): Output directory root.

_Functional/basis selection defaults to `wb97m-v/def2-tzvpd` but can be overridden on the CLI. Charge/spin inherit `.gjf` template metadata when present. If `-q` is omitted but `--ligand-charge` is provided, the input is treated as an enzyme–substrate complex and `extract.py`’s charge summary computes the total charge; explicit `-q` still overrides. For non-`.gjf` inputs, omitting `-q` without `--ligand-charge` aborts; multiplicity defaults to `1` when omitted. Set them explicitly for non-default states._

```yaml
geom:
 coord_type: cart # optional geom_loader settings
dft:
 func: wb97m-v # exchange–correlation functional
 basis: def2-tzvpd # basis set name (alternatively use func_basis: "FUNC/BASIS")
 conv_tol: 1.0e-09 # SCF convergence tolerance (hartree)
 max_cycle: 100 # maximum SCF iterations
 grid_level: 3 # PySCF grid level
 verbose: 0 # PySCF verbosity (0-9)
 out_dir: ./result_dft/ # output directory root
```

---

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing

- [freq](freq.md) — MLIP-based vibrational analysis (often precedes DFT refinement)
- [all](all.md) — End-to-end workflow with `--dft`
- [YAML Reference](yaml-reference.md) — Full `dft` configuration options
- [Glossary](glossary.md) — Definitions of DFT, SP (Single Point)
