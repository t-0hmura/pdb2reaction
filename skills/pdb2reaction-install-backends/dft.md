# DFT backend — PySCF + GPU4PySCF (dft.md)

`pdb2reaction dft` is a single-point DFT driver that re-evaluates
stationary-point energies (R / TS / IM / P) at a higher level of theory
than the MLIP used for the geometry. It runs through PySCF on CPU or
GPU4PySCF on CUDA-enabled x86_64.

The DFT backend is **optional** — it ships in the `[dft]` extras and is
not pulled by the default install.

## Install

```bash
pip install 'pdb2reaction[dft]'
```

This pulls (canonical pin in `pyproject.toml`):

| Package | Purpose | Platform |
|---|---|---|
| `pyscf>=2.13.0` | Reference SCF / DFT engine | All |
| `gpu4pyscf-cuda12x>=1.7.0` | CUDA acceleration of PySCF | **x86_64 only** |
| `cupy-cuda12x>=13.0,!=13.4.0` | Tensor backend for GPU4PySCF | x86_64 only |
| `basis-set-exchange>=0.11` | Programmatic basis-set lookup | All |

`cutensor-cu12` (cuTENSOR / RIJCOSX kernels, x86_64 only) is a separate
optional installation; neither the `[dft]` extra nor GPU4PySCF metadata pulls
it in automatically. Install it only when the selected GPU4PySCF path requires it.

On `aarch64` (`uname -m`), the `gpu4pyscf-cuda12x` PyPI wheel is
x86_64-only. The extras install will succeed for `pyscf` and
`basis-set-exchange` but skip GPU4PySCF, leaving you on CPU PySCF.
Build `gpu4pyscf` from source
(https://github.com/pyscf/gpu4pyscf) to use `--engine gpu` on aarch64.

Verify:

```bash
python -c "import pyscf; print('pyscf       :', pyscf.__version__)"
python -c "import gpu4pyscf; print('gpu4pyscf   :', gpu4pyscf.__version__)"   # wheel: x86_64; source build possible
python -c "import cupy; print('cupy        :', cupy.__version__)"
```

## CPU vs GPU choice

| `--engine` | When to pick | Approximate cost |
|---|---|---|
| `gpu` (default) | x86_64 + a working GPU4PySCF/CuPy stack. **Raises `ClickException` if GPU unavailable** — does **not** auto-fallback to CPU | Benchmark a representative structure; functional, basis, grid, and hardware dominate |
| `cpu` | aarch64, no compatible GPU wheel, or an explicit CPU calculation | Usually slower for large hybrid-DFT jobs; do not assume a fixed factor |

## CLI usage

```bash
pdb2reaction dft -i ts.pdb \
    -l 'SAM:1,GPP:-3' \
    --func-basis 'wb97m-v/def2-svp' \
    --engine gpu                  # default; use 'cpu' to force PySCF CPU
```

Common flag set:

| Flag | Purpose | Default |
|---|---|---|
| `-i, --input` | `.pdb`, `.cif` / `.mmcif`, `.xyz`, or `.gjf` input | required |
| `-q, --charge` / `-l, --ligand-charge` | Total charge or per-residue mapping | required for `.xyz` without `--ref-pdb` |
| `-m, --multiplicity` | Spin multiplicity (2S+1) | 1 |
| `--func-basis` | `'FUNC/BASIS'` like `'wb97m-v/def2-tzvpd'` | `wb97m-v/def2-tzvpd` |
| `--engine` | `gpu` / `cpu` | `gpu` |
| `-o, --out-dir` | Output directory | `./result_dft/` |

Inspect the live default kwargs:

```bash
python -c "from pdb2reaction.workflows.dft import DFT_KW; print(DFT_KW)"
```

## Failure diagnosis

| Symptom | Likely cause | Fix |
|---|---|---|
| `OSError: libcusolver.so.11 not found` | Missing/mixed CUDA wheel dependency or environment library collision | Run the clean-environment diagnostics in `env-cuda.md`; avoid a guessed hard-coded library path |
| `cupy.cuda.runtime.CUDARuntimeError: invalid device ordinal` | Requested local device index is outside the scheduler-visible set | Keep the scheduler's `CUDA_VISIBLE_DEVICES`; select a valid **local** ordinal (usually 0 in a one-GPU job) |
| `RuntimeError: CUDA out of memory` mid-SCF | The selected method/system exceeds available VRAM | Try the same method with `--engine cpu` or a larger-memory GPU. A smaller basis/grid is a different scientific method and must be labeled/revalidated. |
| `gpu4pyscf` imports but SCF stalls/fails near startup | The message alone does not identify cuTENSOR or another dependency | Capture the full traceback/log, run `pip check`, and compare with the installed GPU4PySCF version's official requirements before installing extra CUDA libraries. |
| aarch64: `--engine gpu` requested but no `gpu4pyscf` | x86_64-only wheel | Re-run with `--engine cpu` or build from source (see Install). |

## Resource planning

Atom count alone does not determine memory or wall time. Basis size, element
types, functional, integration grid, density fitting, and GPU model all matter.
Run one representative single point with scheduler memory/VRAM monitoring
before submitting a batch. This command does not expose a general multi-GPU
SCF mode, so do not plan capacity by summing VRAM across requested GPUs.

## See also

- `env-cuda.md` — `LD_LIBRARY_PATH` and torch CUDA pairing.
- [`pdb2reaction-cli/dft.md`](../pdb2reaction-cli/dft.md) — full subcommand flag reference.
- `pdb2reaction-workflows-output/SKILL.md` — DFT//MLIP single-point
  workflow (run `pdb2reaction dft` after `pdb2reaction all`).
