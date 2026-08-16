# MLIP Backends

pdb2reaction uses `MLIPCalculator` for pysisyphus-based geometry and path
stages, and a separate ASE calculator factory for ASE-based stages such as DMF.
The DFT stage uses PySCF/GPU4PySCF directly. This page documents how to select a
backend, the per-backend kwargs, and how to add a new backend.

## Backend dispatcher pattern

```python
from pdb2reaction.backends import create_calculator, create_ase_calculator

calc = create_calculator(
 backend="uma", # one of: "uma", "orb", "mace", "aimnet2", "auto"
 charge=0, spin=1,
 device="cuda", workers=1,
 model="uma-s-1p2",
)
# calc is a pysisyphus-compatible MLIPCalculator.

# ASE-based stages (e.g. DMF path optimization) use the ASE factory:
ase_calc = create_ase_calculator(backend="uma", model="uma-s-1p2", device="cuda")
```

`create_calculator(...)` filters `**kwargs` against each backend's
`_BACKEND_ACCEPTED_KEYS` set. Unknown keys are dropped with a warning.
`create_ase_calculator()` uses its own `_ASE_ACCEPTED_KEYS` set and silently
filters unknown keys.

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

## Per-backend characteristics

| backend | install | model identifier | precision option |
|---------|---------|------------------|------------------|
| `uma` | `pip install fairchem-core` + HF auth | `uma-s-1p2` / `uma-s-1p1` | `precision="fp32" \| "fp64"` |
| `orb` | `pip install orb-models` | `orb_v3_conservative_omol` | `precision="float32-high" \| "float32-highest" \| "float64"` (`fp32` / `float32` are normalized aliases) |
| `mace` | dedicated conda env: install pdb2reaction, then `pip uninstall -y fairchem-core && pip install 'mace-torch>=0.3.8'` (`mace-torch` pins `e3nn==0.4.4`, `fairchem-core` requires `e3nn>=0.5`) | `MACE-OMOL-0` | `default_dtype="float64"` |
| `aimnet2` | `pip install aimnet` | `aimnet2` | n/a |

### Precision

`--precision` is backend-agnostic: it routes to UMA `precision` /
ORB `precision` / MACE `default_dtype` automatically; AIMNet2 treats fp32 as a no-op and rejects fp64 (its model inputs are cast to float32 upstream).

When `--precision` is not given, each backend takes its own default:

| backend | default | why |
|---------|---------|-----|
| `uma` | fp32 | The upstream fairchem baseline. |
| `orb` | fp64 | Use `--precision fp32` to select ORB's reduced `float32-high` mode explicitly. |
| `mace` | fp64 | MACE ships `default_dtype="float64"` upstream. |
| `aimnet2` | fp32 | No precision knob. |

Switching OMol-trained UMA from its default fp32 to fp64 can have non-trivial
impact on TSopt and Hessian results. Enable via:

```bash
pdb2reaction tsopt -i ts.pdb -q 0 --precision fp64 ...
pdb2reaction freq -i opt.pdb -q 0 --precision fp64 ...
pdb2reaction irc -i ts.pdb -q 0 --precision fp64 ...
```

Conversely `--precision fp32` explicitly lowers ORB / MACE precision for throughput.
This is not recommended outside screening. If you do use it, expect noisier Hessians
and check the imaginary-mode count.

`--backend-model NAME` overrides the model variant for the selected `--backend`
(e.g. `--backend uma --backend-model uma-s-1p2`); unset keeps the backend's
built-in default model.

Or via YAML config:

```yaml
calc:
 precision: fp64
```

Requires `fairchem-core ≥ 2.0` for the `InferenceSettings` API.

## Custom backend — bring your own ASE Calculator (`--calc-file`)

Beyond the built-in MLIP backends, any [ASE](https://wiki.fysik.dtu.dk/ase/)
Calculator can be supplied at run time with `--calc-file`, without modifying
pdb2reaction. This couples the pipeline to GFN-xTB (via `tblite` / `xtb-python`),
DFTB+, ORCA, Psi4, or any ASE-compatible engine — the boundary is the standard
ASE Calculator interface (energy in eV, forces in eV/Å).

Write a Python file exposing a `get_calculator` factory that returns an ASE
Calculator:

```python
# my_calc.py  (minimal illustrative example)
from ase.calculators.emt import EMT

def get_calculator(charge=0, spin=1, device="auto", **kwargs):
    return EMT()
```

Swap `EMT()` for the engine you want — e.g. `tblite.ase.TBLite(...)` for
GFN-xTB, the DFTB+ ASE calculator, or `ase.calculators.orca.ORCA(...)`. Then
pass the file to a stage or to `all` (it selects the `custom` backend,
overriding `--backend`):

```bash
pdb2reaction sp     -i model.xyz --calc-file my_calc.py -q 0 -m 1
pdb2reaction opt    -i model.xyz --calc-file my_calc.py -q 0 -m 1
pdb2reaction tsopt  -i ts.xyz    --calc-file my_calc.py -q 0 -m 1
pdb2reaction freq   -i ts.xyz    --calc-file my_calc.py -q 0 -m 1
pdb2reaction all    -i R.pdb P.pdb -c 'LIG' --calc-file my_calc.py -q 0 -m 1
```

Notes:

- The factory receives `charge`, `spin` (multiplicity; also offered as `mult` /
  `multiplicity`), and `device` when its signature accepts them, or
  unconditionally if it declares `**kwargs`, so engines that need the total
  charge (e.g. xTB) can be configured. Use a different factory name with
  `--calc-file-func-name NAME`; a module-level Calculator instance is also accepted.
- Hessians use the finite-difference path inherited from `MLIPCalculator`, so
  `freq` and `tsopt --opt-mode hess` work with any engine. Frozen atoms
  (`--freeze-links` / `--freeze-atoms`) are honored as usual.
- Available on `all` and the standalone subcommands (`sp`, `opt`, `tsopt`,
  `freq`, `irc`, `scan` / `scan2d` / `scan3d`, `path-opt`, `path-search`).
  `all` forwards the same factory to its calculator-backed child stages. For
  a permanent, installable backend with its own `--backend` name, see below.

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
 and `_ASE_ACCEPTED_KEYS["xyz"]`. If it should participate in
 `resolve_backend('auto')`, also register its import probe in
 `_BACKEND_AVAILABILITY_MODULES` and add `"xyz"` to the fallback tuple.
5. **Document + smoke**: add an entry to this page's file map / per-backend
 table, document model identifiers + install command, and add an `xyz` line
 to `tests/smoke/run.sh` so the new backend is exercised end-to-end.

## See Also

- [Architecture](architecture.md) — 6-layer directory map + dependency direction.
- [CONTRIBUTING](https://github.com/t-0hmura/pdb2reaction/blob/main/CONTRIBUTING.md) — Recipe 3.2 "Add an MLIP backend" with full gate cycle references.
