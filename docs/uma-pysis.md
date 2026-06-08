# MLIP Calculator

## Overview
`pdb2reaction` supports multiple machine-learning interatomic potentials (MLIPs) as calculator backends for pysisyphus. The default backend is **UMA** (Meta's Universal Models for Atoms), but **ORB**, **MACE**, and **AIMNet2** are also available. Each backend returns energies, forces, and Hessian matrices in hartree-based atomic units while handling GPU/CPU dispatch and bohr↔Å conversion internally. For cluster-sized systems (hundreds of atoms), the 3N × 3N Hessian tensors processed by augmented-Hessian eigensolves (RFO / RS-I-RFO) and IRC propagation would otherwise dominate wall-clock through repeated host–device synchronization, so the Hessian is kept on-device through pysisyphus (energies and forces flow as scalars / numpy arrays at negligible cost). The calculator is used throughout `pdb2reaction` for optimization, path searches, thermochemistry, and trajectory post-processing.

`pdb2reaction` ships a GPU-accelerated `pysisyphus` fork: the RFO / RS-I-RFO single-structure optimizers, the EulerPC IRC integrator, and the vibrational-mode (Hessian) diagonalization run on CUDA when `device="cuda"` (or `"auto"` on a GPU host). On CPU hosts the same routines fall back transparently to the upstream NumPy/SciPy paths.

## Quick start
```python
import numpy as np
from pdb2reaction.backends.uma import UMACalculator

# Example: a neutral singlet diatomic on GPU when available
calc = UMACalculator(charge=0, spin=1, model="uma-s-1p1", device="auto")

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

- Coordinates are supplied in **bohr**; the wrapper converts to Angstrom for UMA and converts energies/derivatives back to hartree / hartree bohr⁻¹ / hartree bohr⁻².
- Attach the calculator to a `pysisyphus` geometry object or call it directly as above.

## Python API: Calculator Factory

The `backends` module provides a factory for creating MLIP calculators programmatically:

```python
from pdb2reaction.backends import create_calculator, create_ase_calculator
```

| Function | Description |
|----------|-------------|
| `create_calculator(backend="uma", **kwargs)` | Create a pysisyphus-compatible MLIP calculator. Accepts `charge`, `spin`, `model`, `device`, `solvent`, `solvent_model`, `hessian_calc_mode`, `freeze_atoms`, and other backend-specific kwargs. Unknown keys are silently filtered per-backend. |
| `create_ase_calculator(backend="uma", **kwargs)` | Create an ASE-compatible MLIP calculator (used for DMF workflows and ASE-based tools). Same kwargs as `create_calculator`. |

### Example

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

# ORB calculator with implicit solvent
calc_orb = create_calculator(
    backend="orb",
    charge=-1,
    spin=1,
    solvent="water",
    solvent_model="alpb",
)
```

The returned calculator implements the pysisyphus calculator interface (`get_energy`, `get_forces`, `get_hessian`). All methods accept `(atoms: List[str], coords: np.ndarray)` where coords are in **Bohr** and return dicts with `"energy"` (hartree), `"forces"` (hartree/bohr), and `"hessian"` (hartree/bohr²).

## Backend selection

Select a backend with `-b/--backend` on any command, or set `calc.backend` in YAML:

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

| Backend | Install | Analytical Hessian | Multi-worker | Notes |
|---------|---------|-------------------|-------------|-------|
| **UMA** | included | Yes (autograd) | Yes | Full feature set via `fairchem-core` |
| **ORB** | `pip install "pdb2reaction[orb]"` | Yes (autograd) | No | orb-models (conservative models only) |
| **MACE** | `pip uninstall -y fairchem-core && pip install mace-torch` | Yes (`calc.get_hessian`) | No | mace-torch >= 0.3.8 |
| **AIMNet2** | `pip install "pdb2reaction[aimnet]"` | Yes (native) | No | aimnet |

### Implicit solvent correction

All backends support xTB-based implicit solvent corrections via `--solvent`:

```bash
pdb2reaction opt -i input.pdb -q 0 --solvent water
pdb2reaction opt -i input.pdb -q 0 -b orb --solvent water --solvent-model cpcmx
```

The correction uses a delta approach: ΔE = E_xTB(solvent) - E_xTB(vacuum), added to the MLIP energy/forces/Hessian. Requires `xtb` to be installed and accessible on `PATH`.

## Key features

- **MLIP backends** – the default UMA backend loads pretrained UMA checkpoints via `fairchem-core`'s `pretrained_mlip` helpers and forwards charge/spin metadata in the AtomicData batch. Alternative backends (ORB, MACE, AIMNet2) are available via `-b/--backend`.
- **Device handling** – `device="auto"` selects CUDA when available, otherwise CPU. Graph construction happens on the chosen device; when `workers>1`, the parallel predictor manages device transfers.
- **Freeze atoms** – provide 1-based indices via `freeze_atoms`; frozen atoms receive zeroed forces. The Hessian matrix either drops frozen degrees of freedom (`return_partial_hessian=True`) or zeroes the corresponding rows and columns in the full matrix.
- **Precision control** – energies and forces are always returned as float64. Set `hessian_double=False` to obtain the Hessian matrix in the model's native dtype (typically float32).
- **Multi-worker inference** – `workers>1` spawns `fairchem-core`'s `ParallelMLIPPredictUnit` with `workers_per_node` workers per node, useful for batch throughput.

(workers-fd-downgrade)=
### `workers > 1` silent FD downgrade (UMA backend)

```{warning}
When the UMA backend is used with `workers > 1`, analytical Hessians are silently downgraded to finite differences even if `hessian_calc_mode="Analytical"` is explicitly set (an internal flag, not user-configurable). **No log marker is emitted when this downgrade occurs.** The only way to diagnose it is to notice that the Hessian evaluation time grows to the FD-equivalent cost relative to an analytical-Hessian reference run of the same atom count. Use `workers = 1` if you need analytical Hessians. This downgrade applies to every subcommand that exposes `--workers` / `--workers-per-node` (`opt`, `tsopt`, `freq`, `irc`, `all`, `path-opt`, `path-search`, and the scan family). On non-UMA backends (ORB, MACE, AIMNet2) `workers` / `workers_per_node` are silently filtered out per `_BACKEND_ACCEPTED_KEYS`, so this rule does not apply.
```

(hessian-evaluation)=
### Hessian evaluation mode

`hessian_calc_mode="Analytical"` uses second-order autograd on the selected device; `"FiniteDifference"` (default) computes central differences of forces. Analytical mode is silently downgraded to finite differences when multiple inference workers are requested (see the note above); pick `workers = 1` if you need analytical Hessians, or stay on FiniteDifference.

## HPC example: PBS + Open MPI + Ray

For multi-node deployment with `workers` / `workers_per_node` scaled across a Ray cluster under PBS + Open MPI, see the dedicated [HPC example page](hpc-example.md). The template there is copy-paste-ready after adjusting module names, conda path, and PBS resource requests to your site.

(configuration-reference)=
## Configuration reference
Common constructor keywords (defaults shown in the rightmost column):

| Option | Description | Default |
| --- | --- | --- |
| `backend` | MLIP backend engine. | `"uma"` |
| `charge` | Total system charge. | `0` |
| `spin` | Spin multiplicity (2S+1). | `1` |
| `model` | UMA pretrained model name (`uma-s-1p1`, `uma-m-1p1`). | `"uma-s-1p1"` |
| `task_name` | Task tag recorded in UMA batches. | `"omol"` |
| `device` | "cuda", "cpu", or automatic selection. | `"auto"` |
| `workers` / `workers_per_node` | Parallel UMA predictors (UMA backend; ignored by ORB / MACE / AIMNet2); requesting `workers > 1` silently downgrades analytical Hessians to finite differences (no warning). FD is the default `hessian_calc_mode` anyway, so this usually matters only when `Analytical` was explicitly requested. | `1` / `1` |
| `max_neigh`, `radius`, `r_edges` | Optional overrides for UMA neighborhood construction. | `None`, `None`, `False` |
| `freeze_atoms` | List of 1-based atom indices to freeze. | _None_ |
| `hessian_calc_mode` | "Analytical" or "FiniteDifference" for Hessian evaluation. | `"FiniteDifference"` |
| `return_partial_hessian` | Return only the active-DOF Hessian block instead of the full matrix. | `True` |
| `hessian_double` | Assemble and return the Hessian in float64 precision. | `True` |
| `out_hess_torch` | Return Hessians as `torch.Tensor` objects. | `True` |
| `print_timing` | Print Hessian computation timing breakdown. | `True` |
| `print_vram` | Print CUDA VRAM usage during Hessian evaluation (UMA backend only). | `True` |
| `solvent` | Implicit solvent name (e.g. `"water"`) or `"none"`. | `"none"` |
| `solvent_model` | xTB solvent model: `"alpb"` or `"cpcmx"`. | `"alpb"` |
| `xtb_cmd` | Command used to invoke xTB for the solvent correction. | `"xtb"` |
| `xtb_acc` | xTB accuracy setting passed to the solvent-correction run. | `0.2` |


## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [opt](opt.md) -- Single-structure optimization using an MLIP backend
- [path-opt](path-opt.md) -- MEP optimization with MLIP backend
- [all](all.md) -- End-to-end workflow using MLIP across stages
