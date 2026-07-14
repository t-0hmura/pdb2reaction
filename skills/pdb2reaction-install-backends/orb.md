# Orb backend (orb.md)

Orb-v3 (Orbital Materials) is available as an alternative MLIP integration.
pdb2reaction uses fp64 by default because its own ORB benchmark found noisy
finite-difference Hessians in the explicit reduced-fp32 mode. Runtime and
domain suitability must be measured for the actual system, and TS results
need the same frequency/IRC validation as every backend.

## Install

```bash
pip install 'pdb2reaction[orb]'         # pulls orb-models
```

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
checkpoint uses energy-derived conservative forces; the `_omol`
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
| `compile_model` | Request torch compilation; startup/runtime benefit is version, model, and workload dependent, so benchmark it |
| `freeze_atoms`, `hessian_calc_mode`, `return_partial_hessian`, `hessian_double` | Same as UMA |

## Integration constraints

| Property | Interpretation |
|---|---|
| No gated model login in the standard install path | Weight availability can still depend on the installed orb-models release/network cache |
| Explicit fp32 maps to `float32-high` (reduced CUDA matmul, commonly TF32) | Use for screening only when its numerical effect is acceptable; pdb2reaction defaults ORB to fp64 |
| Conservative OMol checkpoint selected by default | This identifies model construction/training family, not guaranteed accuracy for a particular reaction |

For final results, keep fp64 and require exactly one independently
recomputed imaginary mode plus the expected IRC connectivity. Cross-check
with another backend when the chemistry or mode assignment remains ambiguous.

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| Extra low-magnitude imaginary modes after a finite-difference Hessian | Check the effective precision and keep the backend-default `'float64'`; explicit `--precision fp32` selects `float32-high`/TF32 and can amplify force noise. Recompute independently before classifying the stationary point. |
| `compile_model=True` changes startup/runtime behavior | Benchmark compiled and uncompiled modes with the installed torch/orb versions; disable it if compilation fails or does not amortize. |
| TS result has more than one imaginary mode | Treat it as not converged: inspect modes, retry from a better MEP seed and/or `--flatten`, then cross-check with UMA/MACE if ambiguity remains. |

## See also

- `env-cuda.md` — torch / CUDA prerequisites.
- `uma.md` — default backend integration.
- `mace.md` — alternative OMol backend integration (separate env).
- [`pdb2reaction-cli/tsopt.md`](../pdb2reaction-cli/tsopt.md) — diagnosing TS convergence problems.
