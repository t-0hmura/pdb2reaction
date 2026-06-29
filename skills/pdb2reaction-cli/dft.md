# `pdb2reaction dft`

## Purpose

Single-point DFT energy on an arbitrary geometry, via PySCF (CPU) or
GPU4PySCF (CUDA, x86_64). Use as a post-MLIP refinement on R / TS / P
geometries from `irc` / `tsopt`, or as a standalone DFT driver on any
input.

## Synopsis

```bash
pdb2reaction dft -i geom.{pdb,xyz,gjf} \
    [-q 0 -m 1] [-l 'RES:Q,...'] \
    [--func-basis 'wb97m-v/def2-tzvpd'] \
    [--engine gpu|cpu] \
    [-o ./result_dft/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | `.pdb` / `.xyz` / `.gjf` |
| `-q` / `-l` / `-m` | — | — | Charge / spin (required for `.xyz` without `--ref-pdb`) |
| `--ref-pdb` | path | none | Reference PDB so `-l` works on `.xyz` input |
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
| `<out_dir>/result.yaml` | always | energy + per-atom Mulliken / Loewdin / IAO charges & spin densities |
| `<out_dir>/result.json` | `--out-json` | machine-readable result |
| `<out_dir>/input_geometry.xyz` | always | geometry snapshot sent to PySCF |

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
| `OSError: libcusolver.so.11 not found` | [`pdb2reaction-install-backends/env-cuda.md`](../pdb2reaction-install-backends/env-cuda.md) (LD_LIBRARY_PATH order) |
| `cupy ... invalid device ordinal` | `unset CUDA_VISIBLE_DEVICES` |
| `RuntimeError: CUDA out of memory` | Lower `grid_level`, switch to `def2-svp`, or `--engine cpu` |
| aarch64 `--engine gpu` raises `ClickException` ("GPU backend failed...") | PyPI wheel is x86_64-only; re-submit with `--engine cpu` or build `gpu4pyscf` from source (https://github.com/pyscf/gpu4pyscf) |

## Caveats

- `pdb2reaction dft` runs only **single points**, not optimization.
  `tsopt` / `opt` accept only MLIP backends (`-b uma|orb|mace|aimnet2`),
  so DFT-level geometry refinement requires a separate QM code (e.g.
  Gaussian, ORCA, PySCF).
- `--func-basis` follows PySCF naming; cross-check with
  `python -c "from pyscf import gto; print(gto.basis._BASIS_DEFAULT)"`.
- The standalone `dft` subcommand does not accept `--solvent` / `--solvent-model` / `-b/--backend`. xTB-ALPB solvent corrections are MLIP-stage flags (`scan` / `scan2d` / `scan3d`, `path-search`, `path-opt`, `tsopt`, `freq`, `irc`, `opt`, `all`); to combine with DFT, run them at the MLIP stage and then the `dft` single point on the MLIP-optimized geometry.

## See also

- [`pdb2reaction-install-backends/dft.md`](../pdb2reaction-install-backends/dft.md) — install + aarch64 handling.
- `tsopt.md`, `irc.md` — produce the geometries you DFT-refine.
- `pdb2reaction-workflows-output/SKILL.md` — DFT//MLIP recipe.
- Defaults: `import pdb2reaction.core.defaults as d; print(d.GEOM_KW_DEFAULT)`
