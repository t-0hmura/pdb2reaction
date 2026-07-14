# `pysisyphus/` (bundled fork)

> **This is a repo-internal fork of [pysisyphus](https://github.com/eljost/pysisyphus), NOT the upstream PyPI package. Do not `pip install pysisyphus` alongside this package — it will silently overwrite the bundled copy and break long-running jobs.**

The fork is shipped inside this repository; treat it as part of `pdb2reaction` rather than a swappable upstream.

## Why a fork?

The upstream `pysisyphus` does not natively handle:

- **MLIP backends with autograd Hessians** evaluated on a GPU, with the optimizer iterating on CPU coordinates
- **VRAM-aware stage handoff** — explicit `del` between IRC / tsopt / freq stages to free CUDA memory before the next stage loads its model
- **Initial-displacement memory hygiene** in IRC for large systems (16 GB+ Hessians)
- **bofill_update advanced-indexing scatter** when only a subset of internal coordinates is updated (GPU OOM on the naive path)
- **Atomic trial rollback and first-order-saddle validation** for MLIP optimizations, including coordinate-bound exact PHVA checks and path-mode tracking
- **Frozen-boundary PHVA** that removes only full-system rigid motions compatible with fixed anchors

The bundled fork keeps these divergences explicit and limits changes outside the files listed below, so future upstream improvements remain reviewable.

## Divergent files (do NOT replace with upstream)

| file | divergence | rule |
|------|------------|------|
| `pysisyphus/Geometry.py`, `pysisyphus/tr_projection.py` | ordered full/active PHVA metadata, non-mutating normal modes, and constrained rigid-null projection | frozen-boundary vibrational physics |
| `pysisyphus/irc/IRC.py` | initial-displacement memory hygiene for large-system IRC; constrained rigid-null treatment; opt-in `require_pos_def_hessian` PSD convergence guard | chemistry-rule adjacent, freq-stage VRAM invariant |
| `pysisyphus/optimizers/Optimizer.py`, `LBFGS.py`, `RFOptimizer.py` | atomic coordinate/history rollback and uphill-trial rejection for minimizers | optimizer state integrity |
| `pysisyphus/optimizers/HessianOptimizer.py` | opt-in `trust_band` rho-band trust update, `hessian_update_window` multistep TS-BFGS, `weighted_trust` per-coord-type L_inf trust; coordinate/history rollback for rejected minimizer trials; `get_xp`-based torch/numpy dispatch where shared API allows | optimizer step control / trust radius |
| `pysisyphus/optimizers/hessian_updates.py` | `bofill_update` advanced-indexing scatter (CHEMISTRY-RULE:7); `multistep_ts_bfgs_update` helper; re-exports `_outer / _dot` from `_array` shim | scatter on subset of internals |
| `pysisyphus/optimizers/restrict_step.py` | `per_coord_type_weights` + `weighted_max_internal_step` helpers | weighted L∞ trust check |
| `pysisyphus/optimizers/gdiis.py` | `get_xp`-based torch/numpy dispatch (xp.linalg.norm / xp.sum) | torch/numpy backend share |
| `pysisyphus/tsoptimizers/{TSHessianOptimizer,RSIRFOptimizer,RSPRFOptimizer,TRIM}.py` | mode-loss rollback, exact PHVA saddle-order validation, path-mode identity, and bounded recovery for MLIP-driven TS search | tsopt convergence |
| `pysisyphus/_array.py` | torch/numpy backend dispatch shim (`get_xp`, `_outer`, `_dot`, `_eigh`, `as_numpy`, `to_xp`) | used by `hessian_updates.py` + `HessianOptimizer.py` + `gdiis.py` |


## Scope vs upstream

This fork ships only the modules that pdb2reaction needs:
`Calculator` / `Dimer` (calculators), `cos.GrowingString`, `irc.EulerPC`,
`optimizers.{LBFGS, RFOptimizer, StringOptimizer, HessianOptimizer, gdiis,
restrict_step, hessian_updates,...}`, `tsoptimizers.{RSIRFOptimizer,
RSPRFOptimizer, TRIM, TSHessianOptimizer}`, `intcoords/`, `io/`, `cos/`,
`helpers`, `constants`, `elem_data`, `TablePrinter`, `xyzloader`,
`exceptions`, and the `_array` shim.

The upstream `pysisyphus run` CLI driver (`run.py`, `trj.py`, `drivers/`,
`dynamics/`, `color.py`), the wavefunction-analysis subtree
(`wavefunction/`), and the 30+ QM calculator backends (Gaussian / ORCA /
MOPAC / Psi4 / Turbomole /...) are not part of this fork. Any external
caller that imported `from pysisyphus import run` to register a calculator
into `CALC_DICT` must use pdb2reaction's own subcommand layer
(`pdb2reaction/cli/app.py`).

## release scope

During routine polish, **only annotation edits are allowed** on this directory:

- docstring additions / improvements
- type hints
- section banners (`# ===... ===`)
- per-file module docstring

**Forbidden** during polish:

- any change to numerical behaviour, control flow, or function signatures of the divergent files listed above
- any new external dependency
- any rename of public symbols (would break `pdb2reaction/tsopt.py`, `irc.py`, `path_opt.py` callers)

Logic edits require a demonstrated defect or approved numerical feature, a
focused regression test, and the relevant HEAVY/GPU benchmark before merge.
The v0.4.12 optimizer-safeguard changes follow that policy; their regression
and GPU benchmark evidence is recorded in the release validation and changelog.

## Upstream compatibility

If you `pip install pysisyphus` into the same environment as `pdb2reaction`, Python's import machinery may resolve to the bundled copy or the upstream copy depending on `sys.path` order. The flat-top placement at `<repo-top>/pysisyphus/` is required for the bundled copy to take precedence in the editable-install path; do not move it under `pdb2reaction/`.

## See also

- [`../docs/architecture.md`](../docs/architecture.md) §5.3, §6 — repo-internal fork policy
- [`../CONTRIBUTING.md`](../CONTRIBUTING.md) §4.3 — validation policy for logic edits in bundled forks
- `THIRD_PARTY_NOTICES.txt` — pysisyphus upstream attribution
