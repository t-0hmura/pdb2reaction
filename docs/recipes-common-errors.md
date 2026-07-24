# Common Error Recipes

Use this page when you know the symptom but do not know which subcommand page to open first.
For full details, keep [Troubleshooting](troubleshooting.md) open in parallel.

## Quick routing

Each row deep-links into the relevant reference section (mostly [Troubleshooting](troubleshooting.md)) via the section title in the **Then read** column.

| Symptom | Start here | Then read (reference section) |
| --- | --- | --- |
| **Input & extraction** | | |
| Missing element columns / extraction aborts | `add-elem-info` on the original PDB | {ref}`Input / extraction problems <input-extraction-problems>` |
| `[multi] Atom count mismatch` / `[multi] Atom order mismatch` | Regenerate all PDBs with the same prep tool + settings; never reorder atoms once the initial atom order is fixed | {ref}`Input / extraction problems <input-extraction-problems>` |
| **Charge & spin** | | |
| `-q/--charge is required` errors | Set `-q/--charge` or `--ligand-charge/-l` explicitly | {ref}`Charge / spin problems <charge-spin-problems>` |
| Energies/states look wrong after a run | Re-check charge/multiplicity policy in CLI conventions | {ref}`Charge / spin problems <charge-spin-problems>` |
| **Calculation & convergence** | | |
| UMA raises `BackendError` for `--workers > 1` plus `--hessian-calc-mode Analytical` | Use `--workers 1` for an analytical Hessian, or select FiniteDifference | {ref}`workers > 1 Hessian restriction <workers-analytical-error>` |
| CUDA out-of-memory at runtime (`torch.cuda.OutOfMemoryError`) | Re-extract with smaller `--radius` (extract / all only), switch to `--opt-mode grad`, keep default FD Hessian, or move to a larger GPU | {ref}`Calculation / convergence problems <calculation-convergence-problems>` |
| TS converged but extra small imaginary modes remain | Add `--flatten` (available on `tsopt`, `opt`, and `pdb2reaction all`) | {ref}`Calculation / convergence problems <calculation-convergence-problems>` |
| TSOPT does not converge | For the gradient-based optimizers (Limited-memory BFGS (L-BFGS) / Dimer): reduce `max_step`. For the Hessian-based optimizers (Rational Function Optimization (RFO) / Restricted-Step Partitioned RFO (RS-P-RFO, `tsopt` default) / Restricted-Step Image RFO (RS-I-RFO)): reduce `trust_radius`/`trust_min`/`trust_max`. Increase cycles, validate TS quality first | {ref}`Calculation / convergence problems <calculation-convergence-problems>` |
| IRC does not terminate | Reduce `--step-size`, increase `--max-cycles`, confirm a single imaginary mode | {ref}`Calculation / convergence problems <calculation-convergence-problems>` |
| Opt/TSOPT stops with `stalled` after an energy plateau | Treat it as non-converged; inspect the final geometry and force/step criteria, then retry with an appropriate threshold or optimizer setting | {ref}`Calculation / convergence problems <calculation-convergence-problems>` |
| Minimum-energy-path (MEP) search (Growing String Method (GSM) / Direct Max Flux (DMF)) fails | Increase `--max-nodes` above default 20, enable `--preopt` (default `True` for `all`/`path-search`/`path-opt`; `False` for `scan*`), try the alternative `--mep-mode` | {ref}`Calculation / convergence problems <calculation-convergence-problems>` |
| **Installation & environment** | | |
| DMF mode import error (`cyipopt`), or `No module named 'dmf'` | `conda install -c conda-forge cyipopt` (`pydmf` ships with `pdb2reaction`) | {ref}`Installation / environment problems <installation-environment-problems>` |
| UMA model 401/403 / gated-repo error (`huggingface_hub.errors.GatedRepoError`) | Run `hf auth login` and accept the UMA model license | {ref}`Installation / environment problems <installation-environment-problems>` |
| `e3nn` / `fairchem-core` import conflict (MACE in UMA env) | Use a dedicated conda env for MACE (`mace-torch` pins `e3nn==0.4.4`, incompatible with `fairchem-core`'s `e3nn>=0.5`): `pip uninstall -y fairchem-core && pip install mace-torch` | {ref}`Installation / environment problems <installation-environment-problems>` |
| `ImportError: orb-models is required` (or similar for AIMNet2 / MACE) | Install the backend extras: `pip install "pdb2reaction[orb]"`; MACE goes in a separate env | {ref}`Installation / environment problems <installation-environment-problems>` |
| CUDA/GPU runtime mismatch | Verify `torch.cuda.is_available()` and CUDA build pairing | {ref}`Installation / environment problems <installation-environment-problems>` |
| Plot export failures | Run `plotly_get_chrome -y` to install headless Chrome | {ref}`Installation / environment problems <installation-environment-problems>` |

## Recipe 1: Extraction fails before MEP starts

Signal:
- Errors mention missing element symbols, atom-count mismatch, or empty active site models (binding pockets).

First checks:
- Confirm all inputs are prepared by the same workflow and atom ordering is consistent.
- Ensure element columns are present before running `extract` or `all`.

Typical fix path:
- Repair elements with `pdb2reaction add-elem-info -i input.pdb -o input_fixed.pdb` — rerun extraction — confirm active site model size (`--radius`) and residue inclusion (`--selected-resn`). See {ref}`selected-resn-takes-ids` in CLI Conventions for accepted residue names, IDs, and chain-qualified selector forms.

## Recipe 2: Charge/spin validation fails

Signal:
- Errors mention unresolved charge, especially for non-`.gjf` inputs.

First checks:
- Ensure net charge and multiplicity are physically correct for the target state.
- If using residue maps, validate each residue key in `--ligand-charge/-l`.
- Verify the resolution rules in [CLI Conventions](cli-conventions.md) when results look physically inconsistent.

Typical fix path:
- Prefer explicit `-q/--charge` or `--ligand-charge/-l` and `-m` for critical runs, then retry scan/path/tsopt.

## Recipe 3: Environment blockers

Signal:  
- DMF import failures, CUDA mismatch, or unavailable plotting backends.

First checks:  
- Confirm optional dependencies are installed in the active env.
- Validate GPU visibility and PyTorch CUDA compatibility.

Typical fix path:
- Repair environment first, verify with `pdb2reaction --version` and `python -c "import torch; print(torch.cuda.is_available())"`, then rerun with `--dry-run` before full execution.

## Recipe 4: Convergence and post-processing failures

Signal:
- TSOPT stalls, IRC branches look unstable, or MEP refinement stops unexpectedly.

First checks:
- Confirm TS candidate quality: exactly one imaginary frequency, and the corresponding imaginary mode shows displacement along the reaction coordinate. The detection cutoff is `hessian_dimer.neg_freq_thresh_cm` (default 5 cm⁻¹).
- For a wrong imaginary-mode count (a spurious second small mode, or no dominant reaction mode), raise precision with `--precision fp64` and/or switch to `--coord-type dlc`, then add `--flatten` for residual small modes — these levers are complementary. See [`tsopt` → Wrong imaginary-mode count after optimization](tsopt.md#wrong-imaginary-mode-count-after-optimization).
- Tune step sizes / trust radii (YAML knobs `max_step`, `trust_radius`/`trust_min`/`trust_max`) and optimizer mode / flattening (CLI flags `--opt-mode`, `--flatten`); these are complementary. For YAML section layout see [YAML Reference](yaml-reference.md); for the canonical fix path see {ref}`Calculation / convergence problems <calculation-convergence-problems>`.
- If an energy plateau stops the run as `stalled`, keep it classified as non-converged and follow {ref}`Calculation / convergence problems <calculation-convergence-problems>` before retrying.

Typical fix path:
- Run a smaller diagnostic case, tune thresholds/step sizes, then scale back up.

(recipe-cpcmx-build)=
## Recipe 5: Building xTB with CPCM-X support

`conda install -c conda-forge xtb` ships with the Analytical Linearized Poisson-Boltzmann (ALPB) solvent model (default `--solvent-model alpb`) but does **not** include CPCM-X (the Conductor-like Polarizable Continuum Model, extended).

To use `--solvent-model cpcmx` you have to build xTB from source with the CPCM-X feature enabled.

Build (requires GCC ≥ 8 per upstream xtb; if CMake cannot autodetect BLAS/LAPACK, pass `-DBLAS_LIBRARIES=...` and `-DLAPACK_LIBRARIES=...`):

```bash
git clone --depth 1 https://github.com/grimme-lab/xtb.git
cd xtb
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DWITH_CPCMX=ON
# Two-step build avoids a parallel race in tblite:
make -C build tblite-lib -j8
make -C build xtb-exe -j8
```

Set `CPXHOME` to `<xtb-checkout>/build/_deps/cpcmx-src/` at runtime (the directory must contain `cpcmx.toml` and the `DB/` parameter folder), and either prepend `<xtb-checkout>/build` to `$PATH` or set `calc.xtb_cmd` (YAML) / `--xtb-cmd` (where exposed) to the custom binary.
