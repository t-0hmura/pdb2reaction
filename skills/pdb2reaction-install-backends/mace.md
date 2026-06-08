# MACE backend (mace.md)

MACE backend uses `MACE-OMOL-0` (OMol25; 83 elements). Requires a
dedicated env (e3nn pin conflict with UMA).

## Critical: separate environment required

`mace-torch` pins a **different `e3nn` version** than `fairchem-core`
(UMA). The two cannot coexist. Plan: keep UMA in your default env and
put MACE in a separate env, e.g. `<your_mace_env>`.

```bash
conda create -n <your_mace_env> python=3.11 -y
conda activate <your_mace_env>

# torch matching your CUDA driver (see env-cuda.md)
pip install torch --index-url https://download.pytorch.org/whl/<cu_index>

# Install pdb2reaction (pulls in fairchem-core / UMA), then swap UMA
# out for MACE. mace-torch will reinstall e3nn at the version it
# requires.
pip install pdb2reaction
pip uninstall -y fairchem-core
pip install mace-torch
```

If you accidentally keep both UMA and MACE in one env, you'll see
errors like:

```
ImportError: e3nn 0.5.x requires ... but mace-torch installed e3nn 0.4.x
```

The fix is to remove the env and recreate it (`conda env remove -n <env>`), then re-run the install steps above.

## Confirm install

```bash
python -c "import mace; print('mace:', mace.__version__)"
python -c "from pdb2reaction.backends import create_calculator; create_calculator(backend='mace', charge=0, spin=1)"
```

## CLI usage

```bash
conda activate <your_mace_env>
pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -b mace
```

Default model: `MACE-OMOL-0`. Inspect:

```bash
python -c "import pdb2reaction.core.defaults as d; print(d.MACE_BACKEND_DEFAULTS)"
```

## Backend-specific flags

MACE accepts (from `backends/__init__.py:_BACKEND_ACCEPTED_KEYS['mace']`):

| Key | Purpose |
|---|---|
| `charge`, `spin` | Total charge and spin multiplicity |
| `device` | `'cuda'`, `'cpu'`, `'auto'` |
| `model` | Override the default MACE checkpoint |
| `default_dtype` | `'float64'` (default; from `MACE_BACKEND_DEFAULTS`) or `'float32'` for faster, lower-precision inference |
| `freeze_atoms`, `hessian_calc_mode`, `return_partial_hessian`, `hessian_double`, `out_hess_torch`, `print_timing` | Standard cross-backend |

## Strengths and weaknesses

| Strength | Weakness |
|---|---|
| Mature, widely benchmarked | Separate env required (e3nn pin conflict with UMA) |
| — | No multi-GPU sharding API |
| — | First call slow (model load + torch graph build) |

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `e3nn` import error | UMA + MACE in the same env. Use a fresh env. |
| `RuntimeError: Expected all tensors to be on the same device` | Mixed `cpu`/`cuda` tensors after a `.to()` round-trip. Restart Python and ensure `device='cuda'` consistently. |
| Slow Hessian on the default `default_dtype='float64'` | Expected: float64 Hessian is materially slower than float32. Switch to `default_dtype='float32'` if precision allows. |

## See also

- `env-cuda.md` — torch + CUDA prereq.
- `core.md` — `pdb2reaction` install (do this **inside** `<your_mace_env>`).
- `uma.md` — primary backend; keep it in a separate env.
- [`pdb2reaction-cli/tsopt.md`](../pdb2reaction-cli/tsopt.md) — TS solver choice (RS-I-RFO default; Dimer alternative) interacts with backend.