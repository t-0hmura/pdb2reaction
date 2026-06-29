# MLIP Backends

pdb2reaction dispatches every workflow stage (`opt`, `scan`, `tsopt`, `freq`,
`irc`, `path-search`,...) through a single `MLIPCalculator` adapter. This page
documents how to select a backend, the per-backend kwargs, and how to add a new
backend.

## Backend dispatcher pattern

```python
from pdb2reaction.backends import create_calculator, create_ase_calculator

calc = create_calculator(
 backend="uma", # one of: "uma", "orb", "mace", "aimnet2", "auto"
 charge=0, spin=1,
 device="cuda", workers=1,
 model="uma-s-1p1",
)
# calc is a pysisyphus-compatible MLIPCalculator; pass it to any stage runner.

# ASE-based stages (e.g. DMF path optimization) use the ASE factory:
ase_calc = create_ase_calculator(backend="uma", model="uma-s-1p1", device="cuda")
```

`create_calculator(...)` filters `**kwargs` against each backend's
`_BACKEND_ACCEPTED_KEYS` set, so workflow code can pass a superset and unknown
keys are silently dropped. The same pattern applies to `create_ase_calculator()`.

`backend="auto"` resolves UMA → Orb → MACE → AIMNet2 in that order and picks the
first one whose import succeeds.

## File map

| file | role |
|------|------|
| `pdb2reaction/backends/__init__.py` | `BACKEND_REGISTRY` dict + `create_calculator()` / `create_ase_calculator()` factories + `resolve_backend('auto')` UMA-first fallback |
| `pdb2reaction/backends/base.py` | `MLIPCalculator(pysisyphus.calculators.Calculator)` ABC + FD-Hessian assembly + unit conversion + `BackendError(RuntimeError)` |
| `pdb2reaction/backends/uma.py` | UMA (Meta FAIR fairchem-core) — autograd Hessian path |
| `pdb2reaction/backends/orb.py` | Orb (Orbital Materials) — precision / compile_model |
| `pdb2reaction/backends/mace.py` | MACE — default_dtype |
| `pdb2reaction/backends/aimnet2.py` | AIMNet2 — charge-aware (excluded from 5-backend benchmark for the p2r paper) |
| `pdb2reaction/backends/solvent.py` | xTB ALPB wrap — wraps any base calculator with a solvent ΔE correction |
| `pdb2reaction/backends/xtb_alpb_correction.py` | Standalone xTB ALPB correction (distinct API from the wrapper above; do not confuse) |

## Per-backend characteristics

| backend | install | model identifier | precision option |
|---------|---------|------------------|------------------|
| `uma` | `pip install fairchem-core` + HF auth | `uma-s-1p1` / `uma-s-1p2` | `precision="fp32" \| "fp64"` |
| `orb` | `pip install orb-models` | `orb_v3_conservative_omol` | `precision="float32-high" \|...` |
| `mace` | dedicated env: `pip uninstall -y fairchem-core && pip install mace-torch` (`e3nn` pin conflicts with UMA) | `MACE-OMOL-0` | `default_dtype="float64"` |
| `aimnet2` | `pip install aimnet` | `aimnet2` | n/a |

### UMA fp64

Switching OMol-trained UMA from default fp32 to fp64 can have non-trivial impact
on TSopt and Hessian results. Enable via:

```bash
pdb2reaction tsopt -i ts.pdb -q 0 --precision fp64 ...
pdb2reaction freq -i opt.pdb -q 0 --precision fp64 ...
pdb2reaction irc -i ts.pdb -q 0 --precision fp64 ...
```

`--precision` is backend-agnostic: it routes to UMA `precision` /
ORB `precision` / MACE `default_dtype` automatically; AIMNet2 treats fp32 as a no-op and rejects fp64 (its model inputs are cast to float32 upstream).

`--backend-model NAME` overrides the model variant for the selected `--backend`
(e.g. `--backend uma --backend-model uma-s-1p2`); unset keeps the backend's
built-in default model.

Or via YAML config:

```yaml
calc:
 precision: fp64
```

Requires `fairchem-core ≥ 2.0` for the `InferenceSettings` API.

## Add-a-backend recipe (5 steps)

To add a new backend `XYZModel` exposed as `--backend xyz`:

1. **Create `pdb2reaction/backends/xyz.py`**: implement
 `XYZCalculator(MLIPCalculator)` (pysisyphus path) and `XYZASECalculator(...)`
 (ASE path). Both must accept the common kwargs `charge / spin / device /
 freeze_atoms / hessian_calc_mode / return_partial_hessian / hessian_double /
 print_timing / model` and any backend-specific kwargs (`precision`,
 `default_dtype`, etc.).
2. **Inherit `MLIPCalculator`** (`backends/base.py`) and implement the
 subclass hook `_compute_energy_forces_ev(elem, coord_ang) -> (energy_eV,
 forces_eV_Ang)`. Optionally override `_compute_analytical_hessian_ev(elem,
 coord_ang) -> hessian_eV_Ang2` if the backend exposes an analytical Hessian;
 otherwise the base class supplies the FD-Hessian assembly + unit conversion
 (eV/Å → Hartree/Bohr) for free.
3. **Register in `BACKEND_REGISTRY`** (`backends/__init__.py`): add
 `"xyz": {"module": "pdb2reaction.backends.xyz", "pysis_cls": "XYZCalculator",
 "ase_cls": "XYZASECalculator"}` next to existing entries.
4. **Declare accepted kwargs**: add the set to `_BACKEND_ACCEPTED_KEYS["xyz"]`
 and `_ASE_ACCEPTED_KEYS["xyz"]`. Add `"xyz"` to the `resolve_backend('auto')`
 fallback order if it should participate.
5. **Document + smoke**: add an entry to this page's file map / per-backend
 table, document model identifiers + install command, and add an `xyz` line
 to `tests/smoke/run.sh` so the new backend is exercised end-to-end.

## See Also

- [Architecture](architecture.md) — 6-layer directory map + dependency direction.
- [CONTRIBUTING](https://github.com/t-0hmura/pdb2reaction/blob/main/CONTRIBUTING.md) — Recipe 3.2 "Add an MLIP backend" with full gate cycle references.
