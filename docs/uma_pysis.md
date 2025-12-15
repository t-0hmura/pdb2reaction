# `uma_pysis` calculator

## Overview
`uma_pysis` exposes Meta's UMA machine-learning interatomic potentials to PySisyphus as an ASE-compatible calculator. It returns energies, forces, and Hessians (via analytical autograd or finite differences) in Hartree units while handling device placement, graph construction, and unit conversions internally. The calculator is used throughout `pdb2reaction` for optimization, path searches, thermochemistry, and trajectory post-processing.

## Quick start
```python
from pdb2reaction.uma_pysis import uma_pysis

# Build a calculator for a neutral singlet system on GPU when available
calc = uma_pysis(charge=0, spin=1, model="uma-s-1p1", device="auto")

energy_only = calc.get_energy(["C", "O"], coords_bohr)
forces = calc.get_forces(["C", "O"], coords_bohr)
hessian = calc.get_hessian(["C", "O"], coords_bohr)
```
- Coordinates are supplied in **Bohr**; the wrapper converts to Å for UMA and converts energies/derivatives back to Hartree / Hartree·Bohr⁻¹ / Hartree·Bohr⁻².
- Attach the calculator to a `pysisyphus` geometry object or call it directly as above.

## Key features
- **UMA backend** – loads pretrained UMA checkpoints via FAIR-Chem's `pretrained_mlip` helpers and forwards charge/spin metadata in the AtomicData batch.
- **Device handling** – `device="auto"` selects CUDA when available, otherwise CPU. Graph construction happens on the chosen device; when `workers>1`, the parallel predictor manages device transfers.
- **Hessian modes** – `hessian_calc_mode="Analytical"` uses second-order autograd on the selected device; `"FiniteDifference"` (default) computes central differences of forces. Analytical mode is automatically disabled when multiple inference workers are requested.
- **Freeze atoms** – provide 0-based indices via `freeze_atoms`; frozen atoms receive zeroed forces. Hessians either drop frozen degrees of freedom (`return_partial_hessian=True`) or zero corresponding columns in the full matrix.
- **Precision control** – energies and forces are always returned as float64. Set `hessian_double=False` to obtain the Hessian in the model's native dtype (typically float32).
- **Multi-worker inference** – `workers>1` spawns FAIR-Chem's `ParallelMLIPPredictUnit` with `workers_per_nodes` workers per node, useful for batch throughput. Analytical Hessians are skipped in this mode.

## Configuration reference
Common constructor keywords (defaults shown in parentheses):

| Option | Description |
| --- | --- |
| `charge` (0) | Total system charge. |
| `spin` (1) | Spin multiplicity (2S+1). |
| `model` ("uma-s-1p1") | UMA pretrained model name. |
| `task_name` ("omol") | Task tag recorded in UMA batches. |
| `device` ("auto") | "cuda", "cpu", or automatic selection. |
| `workers` (1) / `workers_per_nodes` (1) | Parallel UMA predictors; when `workers>1`, analytical Hessians are disabled. |
| `max_neigh`, `radius`, `r_edges` (all `None`/`False`) | Optional overrides for UMA neighborhood construction. |
| `freeze_atoms` (`None`) | List of 0-based atom indices to freeze. |
| `hessian_calc_mode` ("FiniteDifference") | "Analytical" or "FiniteDifference" for Hessian evaluation. |
| `return_partial_hessian` (`False`) | Return only the active-DOF Hessian block instead of the full matrix. |
| `hessian_double` (`True`) | Assemble and return the Hessian in float64 precision. |
| `out_hess_torch` (`True`) | Return Hessians as `torch.Tensor` objects. |

## CLI and YAML usage
`uma_pysis` is registered as a PySisyphus calculator entry point. With a YAML input file matching these keywords you can run:

```bash
uma_pysis input.yaml
```

Within `pdb2reaction` commands (e.g., `all`, `opt`, `path-opt`), calculator settings can be supplied via `--args-yaml` under `calc.kwargs` to reuse the same UMA configuration across stages.
