# MACE backend (mace.md)

The MACE backend defaults to `MACE-OMOL-0`. It requires a dedicated
environment because its current `e3nn` requirement conflicts with UMA's
`fairchem-core` stack.

## Critical: separate environment required

`mace-torch` pins a **different `e3nn` version** than `fairchem-core`
(UMA). The two cannot coexist. Plan: keep UMA in your default env and
put MACE in a separate env, e.g. `<your_mace_env>`.

```bash
conda create -n <your_mace_env> python=3.11 -y
conda activate <your_mace_env>

# torch matching your CUDA driver (see env-cuda.md)
pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/<cu_index>

# Install pdb2reaction (pulls in fairchem-core / UMA), then swap UMA
# out for MACE. mace-torch will reinstall e3nn at the version it
# requires.
pip install pdb2reaction
pip uninstall -y fairchem-core
pip install mace-torch
```

This workaround deliberately leaves the installed `pdb2reaction` distribution's
declared `fairchem-core` requirement unsatisfied, so `pip check` reports it as
missing even though the MACE backend can run. Do not treat that report as proof
that fairchem should be reinstalled in this env. A later
`pip install --upgrade pdb2reaction` may pull fairchem back in; repeat the
uninstall/install swap and the smoke check afterward. This packaging limitation
is why the environment must remain MACE-only.

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
| `default_dtype` | `'float64'` (default; from `MACE_BACKEND_DEFAULTS`) or `'float32'`; benchmark the performance/numerical tradeoff on the target calculation |
| `freeze_atoms`, `hessian_calc_mode`, `return_partial_hessian`, `hessian_double`, `out_hess_torch`, `print_timing` | Standard cross-backend |

## Strengths and weaknesses

| Strength | Weakness |
|---|---|
| OMol model is available through the ASE-compatible MACE runtime | Separate env required (e3nn pin conflict with UMA) |
| — | No multi-GPU sharding API |
| — | First-call time may include download, model loading, and runtime initialization |

Model availability is not evidence that a reaction lies in its reliable domain.
Validate stationary points with an independent frequency calculation, the
intended imaginary displacement, and IRC connectivity.

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `e3nn` import error | UMA + MACE in the same env. Use a fresh env. |
| `RuntimeError: Expected all tensors to be on the same device` | Capture the full traceback and verify that the requested device and every custom tensor agree; the message alone does not identify a pdb2reaction cause. Reproduce in a fresh process before changing the environment. |
| Hessian runtime is high with the default `default_dtype='float64'` | Benchmark supported precisions on the target system and retain independent frequency/IRC validation. |

## See also

- `env-cuda.md` — torch + CUDA prereq.
- `core.md` — `pdb2reaction` install (do this **inside** `<your_mace_env>`).
- `uma.md` — alternate backend; keep it in a separate env.
- [`pdb2reaction-cli/tsopt.md`](../pdb2reaction-cli/tsopt.md) — TS solver choice (RS-P-RFO default; Dimer alternative) interacts with backend.
