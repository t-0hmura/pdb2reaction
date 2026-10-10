# MLIP Backends

This page explains how to choose a calculation backend and lists, for each
backend, the install command, model names, precision, reproducibility settings, and Hessian
evaluation mode. The default backend is **UMA** (Meta's Universal Models for
Atoms). `-b/--backend` also selects **ORB**, **MACE**, and **AIMNet2**. All four
are machine-learning interatomic potentials (MLIPs).

## Per-backend characteristics

Select a backend with `-b/--backend` on any calculation command:

```bash
# UMA (default)
pdb2reaction opt -i input.pdb -q 0

# ORB
pdb2reaction opt -i input.pdb -q 0 -b orb

# MACE
pdb2reaction opt -i input.pdb -q 0 -b mace

# AIMNet2
pdb2reaction opt -i input.pdb -q 0 -b aimnet2
```

| backend | install | model identifier | `--precision` | analytical Hessian | multiple workers |
|---------|---------|------------------|------------------|--------------------|------------------|
| `uma` | included (`fairchem-core` ≥ 2.22 is a core dependency) + [Hugging Face login](installation.md#required) | `uma-s-1p2` (default) / `uma-m-1p1` | `fp32` / `fp64` | yes (autograd) | yes |
| `orb` | `pip install "pdb2reaction[orb]"` | `orb_v3_conservative_omol` (conservative models only) | `fp32` / `fp64` | yes (autograd) | no |
| `mace` | dedicated conda env: install pdb2reaction, then `pip uninstall -y fairchem-core && pip install 'mace-torch>=0.3.8'` (UMA does not run in this env) | `MACE-OMOL-0` | `fp32` / `fp64` | yes | no |
| `aimnet2` | `pip install "pdb2reaction[aimnet]"` | `aimnet2` | `fp32` only | yes | no |

`--backend-model NAME` overrides the model variant for the selected `--backend`
(e.g. `--backend uma --backend-model uma-m-1p1`).

The run prints the backend and model it loads as
`[backend] Preparing MLIP model (<backend> / <model>)...`. With the default
models the parentheses read `UMA / UMA-S-1.2 (OMol)`,
`ORB / ORB-v3-conservative-OMol`, `MACE / MACE-OMOL-0`, or
`AIMNet2 / aimnet2`.

### Precision

`--precision fp32|fp64` sets the floating-point precision of MLIP inference for
every backend: the value is passed to UMA `precision`, ORB `precision`, and MACE
`default_dtype`. AIMNet2 has no precision setting, so `fp32` changes nothing.

When `--precision` is not given, each backend takes its own default:

| backend | default | why |
|---------|---------|-----|
| `uma` | fp32 | The upstream fairchem baseline. |
| `orb` | fp64 | Use `--precision fp32` to select ORB's reduced `float32-high` mode explicitly. |
| `mace` | fp64 | MACE ships `default_dtype="float64"` upstream. |
| `aimnet2` | fp32 | No precision knob. |

Which value to choose depends on the purpose:

| Purpose | Recommended | Why |
| --- | --- | --- |
| Routine run | Leave unset | Keeps the defaults above: UMA/AIMNet2 fp32, ORB/MACE fp64. |
| Speed screening | `--precision fp32` only when needed | This lowers ORB/MACE precision (see [Notes](#notes)). |
| Final TS/Hessian | Leave unset; with UMA, compare `--precision fp64` when n_imag ≥ 2 ([tsopt](tsopt.md#wrong-imaginary-mode-count-after-optimization)) | Whatever the precision, check n_imag from the final Hessian of `tsopt` and confirm with IRC and the endpoint optimizations that the TS connects the intended R and P. |

Enable fp64 with:

```bash
pdb2reaction tsopt -i ts.pdb -q 0 --precision fp64 ...
pdb2reaction freq -i opt.pdb -q 0 --precision fp64 ...
pdb2reaction irc -i ts.pdb -q 0 --precision fp64 ...
```

Or via YAML config:

```yaml
calc:
 precision: fp64
```

## Determinism and reproducibility

`--deterministic` makes repeated runs with the same input give the same result
on the same software and GPU. Without it, two GPU runs with identical inputs can
differ in the last digits.

`--deterministic` turns on PyTorch's
deterministic algorithms (`torch.use_deterministic_algorithms`) and replaces one
PyTorch operation that has no deterministic GPU version. It does not control DFT
calculations, a custom ASE calculator, or third-party GPU code outside PyTorch.

```bash
pdb2reaction opt -i input.pdb -q 0 --deterministic
pdb2reaction all -i r.pdb p.pdb -q -1 --tsopt --deterministic
```

- It applies to the whole process: set on `all`, it covers every MLIP stage that
  `all` runs, so you do not pass it per stage.
- It is slower: the deterministic GPU operations have lower throughput. Use it
  only when you need repeated runs to match.
- It stops with an error when PyTorch has no deterministic version of an
  operation in the run, instead of silently giving non-reproducible output.

| Backend | `--deterministic` |
|---|---|
| `uma` | Supported |
| `orb` / `mace` | PyTorch's deterministic mode is turned on; check that two runs match for the installed backend version |
| `aimnet2` | **Not supported**: the run stops with an error (see [Notes](#notes)) |
| `custom` | Up to the user-supplied ASE calculator; the flag cannot guarantee it |

## Workers and Hessian mode

`--uma-workers N` (default 1) runs N parallel UMA predictors, and
`--uma-workers-per-node` (default 1) sets how many of them run on each node.
Both flags exist on
`opt`, `tsopt`, `freq`, `irc`, `sp`, `all`, `path-opt`, `path-search`, `scan`,
`scan2d`, and `scan3d`. ORB, MACE, and AIMNet2 ignore them with a warning. Whether more workers
shorten a run is covered in [HPC example › Walltime budgeting](hpc-example.md#walltime-budgeting).

(hessian-evaluation)=
### Hessian evaluation mode

`--hessian-calc-mode` chooses how the Hessian is computed. `FiniteDifference` (default) takes central differences of the forces;
`Analytical` uses second-order autograd on the selected device. UMA, ORB, MACE,
AIMNet2, and DFT can compute analytical Hessians, while a custom calculator
supports only `FiniteDifference`.

## xTB solvent correction

For MLIP backends, `--solvent NAME` adds the xTB solvation energy,
`E_xTB(solvent) - E_xTB(vacuum)`, and the matching force and Hessian differences
to the MLIP surface. It is intended mainly for small-molecule solution
calculations. `--solvent-model` selects `alpb` (default) or `cpcmx` for MLIP
backends, and `--solvent-xtb-cmd` passes the
xTB command with extra arguments, for example `'xtb --etemp 1000'` when the xTB
self-consistent charge (SCC) iterations do not converge.

Each correction runs xTB twice (in solvent and in vacuum), so about 200–300
atoms is a practical upper size rather than a hard limit; benchmark the actual
system and hardware. The option models a named bulk solvent, so on an enzyme
cluster use it only when that environment is justified, and keep the reacting
species, charge, multiplicity, backend and model, and energy references the same
when you compare a solution barrier with a cluster barrier.

## DFT backend

`sp`, `opt`, `tsopt`, `irc`, `freq`, `scan`, `scan2d`, `scan3d`, `path-opt`,
`path-search`, and `all` accept `-b dft --func-basis FUNCTIONAL/BASIS
--dft-engine gpu|cpu` (defaults `wb97m-v/def2-svp` and `gpu`), which computes
every energy and force with PySCF/GPU4PySCF. The separate `pdb2reaction dft`
command gives single points with population analysis. Low-memory mode, CPU
threads and host RAM, and SCF checkpoints are described in
[Refine an MLIP TS with DFT](dft-backend.md).

(backends-custom-calculator)=
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
pass the file to a stage or to `all`; it selects the `custom` backend and
overrides `--backend`:

```bash
pdb2reaction sp     -i model.xyz --calc-file my_calc.py -q 0 -m 1
pdb2reaction opt    -i model.xyz --calc-file my_calc.py -q 0 -m 1
pdb2reaction tsopt  -i ts.xyz    --calc-file my_calc.py -q 0 -m 1
pdb2reaction freq   -i ts.xyz    --calc-file my_calc.py -q 0 -m 1
pdb2reaction all    -i R.pdb P.pdb -c 'A:LIG' --calc-file my_calc.py -q 0 -m 1
```

- The factory receives `charge`, `spin` (also as `mult` / `multiplicity`),
  and `device` when its signature accepts them, or
  unconditionally if it declares `**kwargs`, so engines that need the total
  charge (e.g. xTB) can be configured. Use a different factory name with
  `--calc-file-func-name NAME`; a Calculator instance assigned to that name is
  also accepted.
- Hessians are computed by finite differences of the forces, so `freq` and
  `tsopt --opt-mode hess` work with any engine. `--freeze-links` /
  `--freeze-atoms` are honored as usual.
- Available on `all`, `sp`, `opt`, `tsopt`, `freq`, `irc`, `scan` / `scan2d` /
  `scan3d`, `path-opt`, and `path-search`. `all` passes the same factory to every stage that uses a calculator. For a
  permanent, installable backend with its own `--backend` name, see
  [For developers](#for-developers).

## Python API

### Quick start

```python
import numpy as np
from pdb2reaction.backends.uma import UMACalculator

# Example: a neutral singlet diatomic on GPU when available
calc = UMACalculator(charge=0, spin=1, model="uma-s-1p2", device="auto")

# UMACalculator expects coordinates in Bohr (shape: [n_atoms, 3])
coords_bohr = np.array([
 [0.0, 0.0, 0.0],
 [2.2, 0.0, 0.0], # ~1.16 Å
])

symbols = ["C", "O"]

# NOTE: These methods return dicts; extract values with the appropriate key
energy_h = calc.get_energy(symbols, coords_bohr)["energy"] # float (hartree)
forces_h_bohr = calc.get_forces(symbols, coords_bohr)["forces"] # ndarray (hartree/bohr)
hessian_h_bohr2 = calc.get_hessian(symbols, coords_bohr)["hessian"] # ndarray (hartree/bohr²)
```

- Coordinates are supplied in **bohr**; the wrapper converts to Angstrom for UMA
  and converts energies/derivatives back to hartree / hartree bohr⁻¹ /
  hartree bohr⁻².
- `device="auto"` selects CUDA when it is available, otherwise CPU.
- Attach the calculator to a geometry object of `pysisyphus` (the bundled
  optimization library) or call it directly as above.

### Calculator factory

The `backends` module provides a factory for creating MLIP calculators
programmatically:

```python
from pdb2reaction.backends import create_calculator, create_ase_calculator
```

| Function | Description |
|----------|-------------|
| `create_calculator(backend="uma", **kwargs)` | Create a pysisyphus-compatible MLIP calculator. Accepted kwargs are backend-specific; keys the selected backend does not accept are dropped with a warning, except UMA keys (such as `task_name` or `max_neigh`), which other backends drop without one. Direct Python `freeze_atoms` indices are 0-based (CLI/YAML convert their 1-based values). |
| `create_ase_calculator(backend="uma", **kwargs)` | Create an ASE-compatible calculator. Accepted kwargs are backend-specific, and unsupported keys are dropped without a warning. The UMA/ORB/MACE calculators read charge and spin for each frame from `atoms.info`; AIMNet2 takes `charge`/`spin` as constructor arguments instead. |

```python
from pdb2reaction.backends import create_calculator

# UMA calculator with analytical Hessians
calc = create_calculator(
    backend="uma",
    charge=0,
    spin=1,
    device="auto",
    hessian_calc_mode="Analytical",
)

```

The returned calculator implements the pysisyphus calculator interface:
`get_energy`, `get_forces`, and `get_hessian` accept
`(atoms: List[str], coords: np.ndarray)` where coords are in **Bohr** and return
dicts with `"energy"` (hartree), `"forces"` (hartree/bohr), and `"hessian"`
(hartree/bohr²). Frozen atoms get zero forces, and the Hessian is either the
block of the movable atoms (`return_partial_hessian=True`) or the full matrix
with the frozen rows and columns set to zero.

(configuration-reference)=
### Configuration reference

Common calculator keywords, which are also the keys of the YAML `calc` section.
All `calc` keys, including the xTB
solvent keys, are listed in [YAML Reference › calc](yaml-reference.md#calc).

| Option | Description | Default |
| --- | --- | --- |
| `backend` | MLIP backend engine. | `"uma"` |
| `charge` | Total system charge; used only when written in YAML, and `-q`/`-l` override it. | none |
| `spin` | Spin multiplicity (2S+1). | `1` |
| `model` | Model of the selected backend (`--backend-model`). Left at the UMA default, ORB / MACE / AIMNet2 use `orb_v3_conservative_omol` / `MACE-OMOL-0` / `aimnet2`. | `"uma-s-1p2"` |
| `precision` | MLIP numerical precision (`"fp32"` or `"fp64"`); `"auto"` gives fp32 for UMA/AIMNet2 and fp64 for ORB/MACE. | `"auto"` |
| `task_name` | Task tag recorded in UMA batches. | `"omol"` |
| `device` | "cuda", "cpu", or automatic selection. | `"auto"` |
| `workers` / `workers_per_node` | Parallel UMA predictors (UMA backend; ORB / MACE / AIMNet2 ignore them with a warning). | `1` / `1` |
| `max_neigh`, `radius`, `r_edges` | Optional overrides for UMA neighborhood construction. | `None`, `None`, `False` |
| `freeze_atoms` | List of 0-based atom indices for the direct Python API (CLI/YAML are 1-based). | _None_ |
| `hessian_calc_mode` | "Analytical" or "FiniteDifference" for Hessian evaluation. | `"FiniteDifference"` |
| `return_partial_hessian` | Return the Hessian of the movable atoms only instead of the full matrix. | `True` |
| `hessian_double` | Assemble and return the Hessian in float64 precision. | `True` |
| `out_hess_torch` | Return Hessians as `torch.Tensor` objects. | `True` |
| `print_timing` | Print Hessian computation timing breakdown. | `True` |
| `print_vram` | Print CUDA VRAM usage during Hessian evaluation (UMA backend only). | `True` |

## For developers

### Backend dispatcher pattern

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

pysisyphus-based geometry and path stages use `create_calculator(...)`;
ASE-based stages such as Direct Max Flux (DMF) use `create_ase_calculator(...)`.
`backend="auto"` tries UMA, ORB, MACE, and AIMNet2 in that order and uses the
first one that imports; YAML `calc.backend: auto` does the same, while `-b`
does not take `auto`.

### File map

| file | role |
|------|------|
| `pdb2reaction/backends/__init__.py` | `BACKEND_REGISTRY` dict + `create_calculator()` / `create_ase_calculator()` factories + `resolve_backend('auto')` UMA-first fallback |
| `pdb2reaction/backends/base.py` | `MLIPCalculator(pysisyphus.calculators.Calculator)` base class: frozen-atom handling, finite-difference Hessian assembly, unit conversion, and the error type for backend failures |
| `pdb2reaction/backends/uma.py` | UMA (Meta FAIR fairchem-core) — autograd Hessian path |
| `pdb2reaction/backends/orb.py` | Orb (Orbital Materials) — precision / compile_model |
| `pdb2reaction/backends/mace.py` | MACE — default_dtype |
| `pdb2reaction/backends/aimnet2.py` | AIMNet2 — charge-aware |
| `pdb2reaction/backends/pyscf_dft.py` | PySCF/GPU4PySCF DFT/HF calculator that reuses the SCF state between steps, caches results for identical coordinates, and writes optional SCF checkpoints |

To add a built-in backend with its own `--backend` name, follow recipe 3.2
"Add an MLIP backend" in
[CONTRIBUTING](https://github.com/t-0hmura/pdb2reaction/blob/main/CONTRIBUTING.md).

## Notes

- `--precision fp32` on ORB or MACE is for screening only: expect noisier
  finite-difference Hessians and check n_imag before you use the result.
- AIMNet2 supports neither `--precision fp64` nor `--deterministic`; both stop
  the run with an error. AIMNet2 casts its model inputs to float32, and it
  computes forces with its own CUDA code outside PyTorch's deterministic mode,
  so its forces are not bit-reproducible (the energy is). When you need
  repeatable runs, use UMA, ORB, or MACE with `--deterministic` and run twice in
  the same environment to compare.

(workers-analytical-error)=
- With UMA, `--uma-workers` above 1 cannot be combined with
  `--hessian-calc-mode Analytical`: the run stops with an error because the
  parallel predictor has no autograd model. Use `--uma-workers 1` for an
  analytical Hessian, or `FiniteDifference` with several workers.
- Model precision and Hessian precision are separate settings. Energies and forces are always returned in float64, and the Hessian is assembled in float64 by default; `calc.hessian_double: false` returns it in the model's native dtype (typically float32). With `--precision fp64`, the Hessian is always float64: a `hessian_double: false` in the config is overridden with a warning.
- For CI jobs or the Python API `create_calculator`, the environment variable `PDB2REACTION_STRICT_DETERMINISTIC=1` turns on the same mode as `--deterministic`.

## See Also

- [Architecture](architecture.md) — directory map and dependency direction.
- [HPC example](hpc-example.md) — PBS + Open MPI + Ray template for scaling `workers` / `workers_per_node` across nodes.
- [Refine an MLIP TS with DFT](dft-backend.md) — DFT settings, memory, and checkpoints.
- [Troubleshooting](troubleshooting.md) — detailed troubleshooting guide.
- [opt](opt.md), [path-opt](path-opt.md), [all](all.md) — commands that run on the selected backend.
