# AIMNet2 backend (aimnet2.md)

AIMNet2 (IsayevLab) is a lightweight MLIP for small organic molecules
with narrow element coverage (see list below).

## Element coverage

AIMNet2 supports H, B, C, N, O, F, Si, P, S, Cl, As, Se, Br, I.
**No metals** — neither alkaline-earth (Mg, Ca) nor first-row
transition metals (Mn, Fe, Co, Ni, Cu, Zn) are covered.
Metalloenzymes need UMA or MACE.

## Install

```bash
pip install 'pdb2reaction[aimnet]'         # pulls aimnet>=0.2.0
```

Or, if `pdb2reaction` is already installed:

```bash
pip install 'aimnet>=0.2.0'
```

Confirm:

```bash
python -c "import aimnet; print('aimnet:', aimnet.__version__)"
python -c "from pdb2reaction.backends import create_calculator; create_calculator(backend='aimnet2', charge=0, spin=1)"
```

## CLI usage

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c 'GLU' -l 'GLU:-1' \
    --tsopt --thermo \
    -b aimnet2
```

Default model: `aimnet2`. Inspect:

```bash
python -c "import pdb2reaction.core.defaults as d; print(d.AIMNET2_BACKEND_DEFAULTS)"
```

## Backend-specific flags

AIMNet2 accepts (from `backends/__init__.py:_BACKEND_ACCEPTED_KEYS['aimnet2']`):

| Key | Purpose |
|---|---|
| `charge`, `spin` | Total charge and spin multiplicity |
| `device` | `'cuda'`, `'cpu'`, `'auto'` |
| `model` | Override the default checkpoint |
| `freeze_atoms`, `hessian_calc_mode`, `return_partial_hessian`, `hessian_double`, `out_hess_torch`, `print_timing` | Standard cross-backend |

## When to use AIMNet2

| Use it when | Don't use it when |
|---|---|
| Substrate is small organics (≤ 100 atoms, no metals) | Active site contains Zn, Mg, Mn, Fe, etc. |
| You want a cheap CPU baseline run | High accuracy on TS curvature is important |
| Pre-screening many transition states | Reproducing published `wB97M-V` benchmark numbers |

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `KeyError` on element during atom-type lookup | Unsupported element; switch to UMA or MACE. |
| `RuntimeError: charge mismatch` | AIMNet2 charge is per-atom-network output; supply `-q TOTAL` matching the cluster. |
| Convergence failures on radicals | AIMNet2 was trained only on closed-shell systems; use `-m 1` only. |

## See also

- `env-cuda.md` — torch / CUDA prereq.
- `core.md` — `pdb2reaction` install.
- `uma.md` — recommended for anything containing metals.
- [`pdb2reaction-structure-io/charge-multiplicity.md`](../pdb2reaction-structure-io/charge-multiplicity.md) — figuring out
  `-q` and `-m` for an unfamiliar substrate.