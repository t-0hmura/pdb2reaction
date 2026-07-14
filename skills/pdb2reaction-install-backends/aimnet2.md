# AIMNet2 backend (aimnet2.md)

AIMNet2 (IsayevLab) is an MLIP family for organic and elemental-organic
systems. Element and electronic-state coverage depends on the selected model.

## Element coverage

The default `aimnet2` model supports H, B, C, N, O, F, Si, P, S, Cl,
As, Se, Br, and I. It does not support Mg or first-row transition metals.
The separately selected `aimnet2-pd` model adds Pd (with a different training
domain); that does not imply coverage of other metals. Check the installed
`aimnet` model table before choosing a model for a new element set.

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
    -c 'LIG' -l 'LIG:-1' \
    --tsopt --thermo \
    -b aimnet2
```

Default model: `aimnet2`. Inspect:

```bash
python -c "import pdb2reaction.core.defaults as d; print(d.AIMNET2_BACKEND_DEFAULTS)"
```

Useful upstream model names in `aimnet>=0.2.0` include `aimnet2-2025`
(general organic model), `aimnet2-nse` (open-shell/radical chemistry),
`aimnet2-pd` (Pd domain), and `aimnet2-rxn` (H/C/N/O reactive chemistry;
closed-shell, net-neutral systems only, for relative energies within the same
composition and model family).
Pass one with `--backend-model`; model availability is controlled by the
installed `aimnet` version and is not inferred automatically from multiplicity.

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
| Elements and electronic state are covered by the selected model | Active site contains an unsupported element such as Zn, Mg, Mn, or Fe |
| You need a CPU-capable baseline and the selected model covers the system | You need to assume TS curvature quality without an independent frequency/IRC check |
| Pre-screening candidates within the model domain | Reproducing a different electronic-structure reference without validation |

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `KeyError` on element during atom-type lookup | The selected model does not support that element; choose a model/backend whose published element domain includes the whole cluster. |
| Charge/multiplicity result looks inconsistent | Supply the **total cluster** `-q` and multiplicity `-m`; these are model inputs, not per-atom charges. |
| Radical/open-shell result is unreliable | Select an open-shell-capable model such as `aimnet2-nse`, pass the actual multiplicity, and validate independently. The default model is not silently replaced. |
| `aimnet2-rxn` is considered for a charged/open-shell cluster | Do not use it outside its closed-shell, net-neutral H/C/N/O domain; choose a compatible model/backend and validate independently. |

## See also

- `env-cuda.md` — torch / CUDA prereq.
- `core.md` — `pdb2reaction` install.
- `uma.md` / `mace.md` — alternative integrations; verify the chosen model's element/state domain rather than selecting only by backend name.
- [`pdb2reaction-structure-io/charge-multiplicity.md`](../pdb2reaction-structure-io/charge-multiplicity.md) — figuring out
  `-q` and `-m` for an unfamiliar substrate.
