# `dft`

Runs single-point DFT with GPU4PySCF or CPU PySCF, reporting energy and population analysis (Mulliken, meta-Löwdin, IAO charges). The default functional/basis is ωB97M-V/def2-tzvpd. Use it for single-point DFT energy (and population analysis) on a small active-site model — typically to refine MLIP-optimized R/TS/P structures — selecting the backend with `--engine` (default `gpu`); use `cpu` when no GPU is available, or for portable/debug runs.

> See {ref}`engine-vs-dft-engine` for the `--engine` (standalone `dft`) vs `--dft-engine` (forwarded through `pdb2reaction all`) naming convention.

> **Prerequisites:** DFT dependencies (PySCF, GPU4PySCF) are **not** included in the default install. Install them with `pip install "pdb2reaction[dft]"`.

## Examples

Command form:

```bash
pdb2reaction dft -i INPUT.{pdb|xyz|gjf|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULTIPLICITY] \
 [--func-basis 'FUNC/BASIS'] \
 [--max-cycle N] [--conv-tol Eh] [--grid-level L] \
 [--out-dir DIR] [--engine gpu|cpu] [--convert-files/--no-convert-files] \
 [--ref-pdb FILE] [--config FILE] [--show-config] [--dry-run]
```

Basic GPU single point.

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 --engine gpu --out-dir ./result_dft
```

Run with a larger basis and tighter SCF settings.

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 \
 --func-basis 'wb97m-v/def2-tzvpd' --conv-tol 1e-10 --max-cycle 200 \
 --engine gpu --out-dir ./result_dft_tight
```

> **Caveat:** The tight `def2-tzvpd` setting above suits only small active-site models (≲150 atoms on a GPU with sufficient VRAM). See [Notes](#notes) for size / basis / backend thresholds.

Force CPU backend for portability.

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 --engine cpu --out-dir ./result_dft_cpu
```

Derive total charge from ligand mapping when `-q` is omitted.

```bash
pdb2reaction dft -i input.pdb -l 'LIG:0' -m 1 \
 --engine gpu --out-dir ./result_dft_ligand
```

When `-q` is omitted but `--ligand-charge/-l` is provided, the input is treated as an enzyme–substrate complex and `extract.py`’s charge summary computes the total charge; an explicit `-q` still overrides. For non-`.gjf` inputs, omitting `-q` without `--ligand-charge/-l` aborts.

## Workflow

1. **Input handling** – Any file loadable by `geom_loader` (.pdb/.xyz/_trj.xyz/…) is accepted. Coordinates are re-exported as `input_geometry.xyz`. For XYZ/GJF inputs, `--ref-pdb` supplies a reference PDB topology for atom-count validation and (if you also use `--ligand-charge/-l`) charge derivation; the DFT stage itself does **not** emit PDB/GJF outputs.
2. **SCF build** – `--func-basis` is parsed into functional and basis. `--engine` controls GPU/CPU preference (`gpu` requires GPU4PySCF and raises an error if unavailable; `cpu` forces CPU). On the closed-shell GPU path with `--lowmem` (default), the SCF object is `gpu4pyscf.dft.rks_lowmem.RKS`, which uses a memory-efficient direct-JK pipeline (no density fitting); on the open-shell GPU, CPU, or `--no-lowmem` paths, density fitting is enabled automatically with PySCF defaults. Nonlocal corrections (e.g., VV10) are not configured explicitly beyond the backend defaults.
3. **Population analysis & outputs** – After convergence (or failure) the command writes `result.yaml` summarizing the energy (in hartree and kcal/mol), convergence metadata, backend info, and per-atom Mulliken/meta-Löwdin/IAO charges and spin densities (UKS only for spins). Any failed analysis column is set to `null` with a warning.

## Outputs

```text
out_dir/ (default:./result_dft/)
├─ input_geometry.xyz # Geometry snapshot sent to PySCF
├─ result.yaml # Energy/charge/spin summaries with convergence/engine metadata
```

- `result.yaml` expands to:
 - `energy`: energy in hartree and kcal/mol, convergence flag, engine metadata
  (`engine`: `gpu4pyscf(rks_lowmem)` / `gpu4pyscf` / `pyscf(cpu)`; `used_gpu`; `used_lowmem`).
 - `charges`: Mulliken, meta-Löwdin, and IAO atomic charges (`null` when a method fails).
 - `spin_densities`: Mulliken, meta-Löwdin, and IAO spin densities (UKS-only for spins).
- It also summarizes charge, multiplicity, spin (2S), functional, basis,
  convergence knobs, and resolved output directory.

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge supplied to PySCF (`calc.charge`). Required unless a `.gjf` template or `--ligand-charge/-l` (PDB inputs or XYZ/GJF with `--ref-pdb`) supplies it. Overrides `--ligand-charge/-l` when both are set. | Required unless template/derivation applies |
| `-l, --ligand-charge TEXT` | Either a scalar integer (e.g., `-1`) for the total ligand charge, or a per-residue mapping (e.g., `GPP:-3,SAM:1`) that derives the total from PDB residue charges. Used when `-q` is omitted (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). Converted to `2S` for PySCF. | `.gjf` template value or `1` |
| `--func-basis TEXT` | Functional/basis pair in `FUNC/BASIS` form (quote strings with `*`). | `wb97m-v/def2-tzvpd` |
| `--max-cycle INT` | Maximum SCF iterations (`dft.max_cycle`). | `100` |
| `--conv-tol FLOAT` | SCF convergence tolerance in hartree (`dft.conv_tol`). | `1e-9` |
| `--grid-level INT` | PySCF numerical integration grid level (`dft.grid_level`). | `3` |
| `-o, --out-dir TEXT` | Output directory (`dft.out_dir`). | `./result_dft/` |
| `--engine [gpu\|cpu]` | SCF backend: gpu (GPU4PySCF) or cpu (PySCF). See {ref}`engine-vs-dft-engine` for the `--engine` vs `--dft-engine` naming convention. | `gpu` |
| `--lowmem/--no-lowmem` | Use `gpu4pyscf.dft.rks_lowmem.RKS` for closed-shell GPU runs (skips density fitting in favor of memory-efficient direct JK). Open-shell, CPU, or pre-`rks_lowmem` GPU4PySCF installs fall back to standard RKS/UKS automatically. | `True` |
| `--convert-files/--no-convert-files` | **No-op on `dft`.** Accepted purely for interface consistency with the other subcommands; `dft` never produces PDB or GJF outputs (only `input_geometry.xyz` + `result.yaml`). The flag's value is ignored. | `True` |
| `--ref-pdb FILE` | Reference PDB topology to validate atom counts and enable ligand-charge derivation for XYZ/GJF inputs (no output conversion). | _None_ |
| `--config FILE` | Base YAML configuration file applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration and continue execution. | `False` |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` to `out_dir`. See [JSON Output Schema](json-output.md) for the schema. | `False` |
| `--dry-run/--no-dry-run` | Validate options and print execution plan without running DFT. | `False` |

## YAML configuration

Accepts a mapping root; the `dft` section (and optional `geom`) is applied when present. Merge order is:

- defaults
- `--config`
- explicit CLI options

```yaml
geom:
 coord_type: cart # optional geom_loader settings
dft:
 func: wb97m-v # exchange–correlation functional
 basis: def2-tzvpd # basis set name (alternatively use func_basis: "FUNC/BASIS")
 conv_tol: 1.0e-09 # SCF convergence tolerance (hartree)
 max_cycle: 100 # maximum SCF iterations
 grid_level: 3 # PySCF grid level
 verbose: 0 # PySCF verbosity (0-9); CLI -v 2/3 raises runtime PySCF verbosity to >=4
 out_dir: ./result_dft/ # output directory root
```

Full schema (every key and default): [YAML Reference](yaml-reference.md).

## Exit codes

See {ref}`exit-codes` in CLI Conventions.

## Notes

- **System size / basis cost:** `def2-tzvpd` is expensive; on 16-24 GB GPUs it OOMs above ~150 atoms, and DFT single points are practical only up to ~300 atoms. Use `--func-basis 'wb97m-v/def2-svp'` (1-3 kcal/mol barrier-height error) or an external program (ORCA, Gaussian) for larger systems, and extract a small active-site model first.
- **Blackwell GPUs (RTX 50xx):** GPU4PySCF may OOM even at ~100 atoms; use `--engine cpu` or an external DFT program.
- **CPU backend:** `--engine cpu` is practical only for small models (<=150 atoms) and small basis sets; larger systems are prohibitively slow.
- **HPC scratch:** PySCF / GPU4PySCF write to `$PYSCF_TMPDIR` (then `$TMPDIR`, `/tmp`); on nodes with a small or tmpfs `/tmp`, set `PYSCF_TMPDIR` to the job filesystem (e.g. `export PYSCF_TMPDIR="$PBS_O_WORKDIR"`) before launching.
- Compiled GPU4PySCF wheels may not support non-x86 systems; build from source in that case (see https://github.com/pyscf/gpu4pyscf).
- No auxiliary basis guessing is implemented; density-fitting behavior is described under Workflow (SCF build) and the `--lowmem` CLI option.
- The YAML input file must have a mapping root; the `dft` section is optional. Non-mapping roots raise an error via `load_yaml_dict`.
- IAO spin/charge analysis may fail for challenging systems; corresponding columns in `result.yaml` become `null` and a warning is printed.

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed fixes for common failure modes
- [freq](freq.md) — MLIP-based vibrational analysis (often precedes DFT refinement)
- [all](all.md) — End-to-end workflow with `--dft`
- [YAML Reference](yaml-reference.md) — Full `dft` configuration options
- [Glossary](glossary.md) — Definitions of DFT, SP (Single Point)
