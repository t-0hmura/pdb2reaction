# Orb backend (orb.md)

Orb-v3 (Orbital Materials) is a fast MLIP useful for **screening** large
candidate sets where you can tolerate slightly lower TS-region accuracy
than UMA / MACE.

## Install

```bash
pip install 'pdb2reaction[orb]'         # pulls orb-models
```

> **torch_scatter has no PyPI binary wheel** (only an sdist), so `[orb]` source-builds it and fails
> under PEP517 isolation (`No module named 'torch'`). Install from PyG's prebuilt-wheel index matching
> your torch+CUDA tag (single command, no compiler needed):
> ```bash
> pip install 'pdb2reaction[orb]' -f https://data.pyg.org/whl/torch-2.8.0+cu129.html
> ```
> Fallback (CUDA toolchain present): `pip install torch` → `pip install torch_scatter --no-build-isolation` → `[orb]`.

Or, if `pdb2reaction` is already installed:

```bash
pip install orb-models
```

Confirm:

```bash
python -c "import orb_models; print('orb_models:', orb_models.__version__)"
python -c "from pdb2reaction.backends import create_calculator; create_calculator(backend='orb', charge=0, spin=1)"
```

Orb model weights are downloaded on first use; no separate auth required.

## CLI usage

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -b orb
```

Default model: `orb_v3_conservative_omol`. The `_conservative_`
checkpoint uses energy-conservative forces (forces = ∇E); the `_omol`
suffix indicates training on the OMol25 dataset (organics +
biomolecules + transition-metal complexes).

Inspect the default kwarg dict:

```bash
python -c "import pdb2reaction.core.defaults as d; print(d.ORB_BACKEND_DEFAULTS)"
```

## Backend-specific flags

Orb accepts (canonical list in
`backends/__init__.py:_BACKEND_ACCEPTED_KEYS['orb']`):

| Key | Purpose |
|---|---|
| `charge`, `spin` | Total charge and spin multiplicity |
| `device` | `'cuda'`, `'cpu'`, `'auto'` |
| `model` | Override the default Orb checkpoint |
| `precision` | `'float64'` is the default: an unset `--precision` resolves per backend (`ORB_BACKEND_DEFAULTS["precision"]`, `_BACKEND_DEFAULT_PRECISION["orb"] = "fp64"`). `'float32-high'` selects ORB's TF32 matmul mode — fast, but its force noise inflates finite-difference Hessians into spurious imaginary modes, so keep it for screening only. pdb2reaction normalizes `'float32'`/`'fp32'` → `'float32-high'` and `'fp64'` → `'float64'` before the loader (raw orb_models rejects `'float32'` and demotes to a slow path). |
| `compile_model` | `True` to torch-compile (faster after first call, slower start) |
| `freeze_atoms`, `hessian_calc_mode`, `return_partial_hessian`, `hessian_double` | Same as UMA |

## Strengths and weaknesses

| Strength | Weakness |
|---|---|
| Faster per call than UMA on large systems | TS curvature less accurate than UMA / MACE; many TS searches end up with multiple imaginary modes |
| Trained on broad organic chemistry | Not the right tool for fine-grained `wB97M-V` benchmarking |
| Easy install, no auth gate | Smaller user community than UMA |

In practice: use Orb-v3 to **filter** down a list of candidates, then
re-run survivors with UMA or MACE for the final TS / IRC.

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `RuntimeError: ... mat1 and mat2 shapes ... ` during Hessian | Appears once the run has been downgraded to `float32-high` (TF32) — e.g. `--precision fp32`, or a `--config` YAML carrying `calc.precision: fp32`. Keep the `'float64'` default for every Hessian / freq step. |
| `compile_model=True` adds torch-compile overhead on first call | Subsequent calls are faster. Disable for short jobs. |
| TS converges with > 1 imaginary mode | Common with Orb on aromatic or metalloenzyme systems. Re-run that step with UMA/MACE. |

## See also

- `env-cuda.md` — torch / CUDA prerequisites.
- `uma.md` — recommended primary backend for production runs.
- `mace.md` — alternative high-accuracy backend (separate env).
- [`pdb2reaction-cli/tsopt.md`](../pdb2reaction-cli/tsopt.md) — diagnosing TS convergence problems.