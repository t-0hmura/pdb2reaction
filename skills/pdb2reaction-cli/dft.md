# `pdb2reaction dft`

## Purpose

Single-point DFT energy on an arbitrary geometry, via PySCF (CPU) or
GPU4PySCF (CUDA, x86_64). Use as a post-MLIP refinement on R / TS / P
geometries from `irc` / `tsopt`, or as a standalone DFT driver on any
input.

## Synopsis

```bash
pdb2reaction dft -i geom.{pdb,cif,mmcif,xyz,gjf} \
    [-q 0 -m 1] [-l 'RES:Q,...'] \
    [--func-basis 'wb97m-v/def2-tzvpd'] \
    [--engine gpu|cpu] \
    [-o ./result_dft/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | `.pdb` / `.cif` / `.mmcif` / `.xyz` / `.gjf` |
| `-q` / `-l` / `-m` | — | — | Charge / spin (required for `.xyz` without `--ref-pdb`) |
| `--ref-pdb` | path | none | Reference PDB/mmCIF so `-l` works on `.xyz` input |
| `--func-basis` | str | `wb97m-v/def2-tzvpd` | `'FUNC/BASIS'` |
| `--engine` | str | `gpu` | `gpu` (GPU4PySCF) or `cpu` (PySCF) |
| `--lowmem / --no-lowmem` | toggle | `--lowmem` | Closed-shell GPU uses `rks_lowmem.RKS` (no DF); open-shell / CPU / older gpu4pyscf automatically fall back to RKS/UKS+DF |
| `--config` | path | none | YAML config file |
| `-o, --out-dir` | path | `./result_dft/` | Output directory |
| `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

## Examples

### Default DFT//MLIP on a TS

```bash
pdb2reaction dft -i seg_01/ts.pdb \
    -l 'SAM:1,GPP:-3' \
    --func-basis 'wb97m-v/def2-tzvpd' \
    --engine gpu --out-json
```

`--out-json` enables the `result.json` example below; omit it for `result.yaml` only.

### CPU PySCF (aarch64 / no GPU)

```bash
pdb2reaction dft -i ts.xyz -q 0 -m 1 \
    --func-basis 'wb97m-v/def2-svp' \
    --engine cpu \
    -o result_dft_cpu
```

## Output

| Path | When | Content |
|---|---|---|
| `<out_dir>/result.yaml` | successful DFT calculation | energy + per-atom Mulliken / Loewdin / IAO charges & spin densities |
| `<out_dir>/result.json` | `--out-json` | machine-readable result |
| `<out_dir>/input_geometry.xyz` | input preparation succeeds | geometry snapshot sent to PySCF |

`result.json` keys (written only when `--out-json` is passed):

```python
import json
d = json.load(open("result_dft/result.json"))
print(d["energy_hartree"])
print(d["xc_functional"], d["basis_set"])  # e.g. "wb97m-v", "def2-tzvpd"
print(d["engine"])           # "gpu4pyscf(rks_lowmem)", "gpu4pyscf", or "pyscf(cpu)"
print(d["used_gpu"], d["used_lowmem"])  # bool, bool (lowmem False on open-shell / CPU / --no-lowmem)
print(d["converged"])        # True / False (exit code 3 if False)
```

`result.yaml` carries the run info: grid_level, Mulliken / Loewdin /
IAO charges, spin densities.

## Common errors

| Symptom | Fix |
|---|---|
| `OSError: libcusolver.so.11 not found` | Capture `pip check` and compare the clean-environment library-loading test in [`env-cuda.md`](../pdb2reaction-install-backends/env-cuda.md); do not guess a library path. |
| `cupy ... invalid device ordinal` | Keep the scheduler-provided `CUDA_VISIBLE_DEVICES` and select a valid local ordinal (usually device 0 in a one-GPU allocation). Do not unset the scheduler's isolation variable. |
| `RuntimeError: CUDA out of memory` | First try the same method with `--engine cpu` or a larger-memory GPU. Lowering the grid or basis changes the scientific method, so do it only as an explicit new calculation and label it accordingly. |
| aarch64 `--engine gpu` raises `ClickException` ("GPU backend failed...") | PyPI wheel is x86_64-only; re-submit with `--engine cpu` or build `gpu4pyscf` from source (https://github.com/pyscf/gpu4pyscf) |

## Caveats

- `pdb2reaction dft` runs only **single points**, not optimization.
  The built-in `tsopt` / `opt` backend choices are MLIPs
  (`-b uma|orb|mace|aimnet2`), but those optimizers can also use an
  arbitrary ASE Calculator through `--calc-file`. There is no built-in
  PySCF/GPU4PySCF geometry optimizer; DFT-level refinement therefore needs
  either a suitable ASE calculator adapter or a separate QM program.
- `--func-basis` follows PySCF naming. Test a basis name directly, for example
  `python -c "from pyscf import gto; print(len(gto.basis.load('def2-tzvpd', 'C')))"`.
- The standalone `dft` subcommand does not accept `--solvent` /
  `--solvent-model` / `-b/--backend`. xTB-ALPB is an MLIP-stage correction.
  A DFT single point on a solvent-corrected-MLIP geometry still reports the
  configured PySCF energy only; it does not add an implicit-solvent or xTB
  correction to the DFT result.

## See also

- [`pdb2reaction-install-backends/dft.md`](../pdb2reaction-install-backends/dft.md) — install + aarch64 handling.
- `tsopt.md`, `irc.md` — produce the geometries used for DFT single points.
- `pdb2reaction-workflows-output/SKILL.md` — DFT//MLIP recipe.
- Defaults: `import pdb2reaction.core.defaults as d; print(d.GEOM_KW_DEFAULT)`
