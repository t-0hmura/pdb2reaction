# Troubleshooting

Find your symptom in the quick table, then read the fix in the section it points to.

(troubleshooting-quick-table)=
## Quick routing

| Symptom | Start here | Then read |
| --- | --- | --- |
| **Input & extraction** | | |
| Blank element columns stop `extract` (`Element symbols are missing in '...'`); `all` fills them itself and stops when some atoms cannot be assigned | Run `add-elem-info` on the original PDB | {ref}`Input / extraction <input-extraction-problems>` |
| `[multi] Atom count mismatch` / `[multi] Atom order mismatch` | Regenerate all PDBs with the same preparation tool and settings; never reorder atoms once the atom order is fixed | {ref}`Input / extraction <input-extraction-problems>` |
| A utility command (`add-elem-info`, `fix-altloc`, `bond-summary`, `energy-diagram`, `trj2fig`) stops with an error or prints a warning | Find the message on that command's page: errors are in its Notes, and `add-elem-info`'s `[WARN] Could not confidently assign` is in its Output files | [add-elem-info](add-elem-info.md#output-files), [fix-altloc](fix-altloc.md#notes), [bond-summary](bond-summary.md#notes), [energy-diagram](energy-diagram.md#notes), [trj2fig](trj2fig.md#notes) |
| `bond-summary` prints `ERROR: Atom types and ordering must be identical.` for a pair | Compare structures that have the same atoms in the same order | [bond-summary](bond-summary.md#notes), {ref}`Input / extraction <input-extraction-problems>` |
| `fix-altloc` stops with `Output exists: <path> (use --overwrite to overwrite)` | Add `--overwrite`, or choose another output path with `-o` | [fix-altloc](fix-altloc.md#main-options) |
| **Charge & spin** | | |
| `-q/--charge is required` / `Total charge could not be resolved` | Set `-q/--charge` or `-l/--ligand-charge` explicitly | {ref}`Charge / spin <charge-spin-problems>` |
| `Cluster electron count inconsistent` (in `all --dry-run`: `--dry-run parity check failed`) | The charge and multiplicity do not fit the electron count. Fix `-q`/`-l`, or set `-m` (for example `-m 2` for an odd electron count) | {ref}`Charge / spin <charge-spin-problems>` |
| Energies or states look wrong after a run | Re-check the charge and multiplicity you passed | {ref}`Charge / spin <charge-spin-problems>` |
| **Calculation & convergence** | | |
| `all` ends with a `Scientific status:` other than `success` and prints `RESULT WARNING:` lines | Each line names the segment and the stage that failed (TS optimization, MEP, IRC, or endpoint optimization); find that stage in the rows below. `summary.json` lists the same reasons in `scientific_status_reasons` (for example `all:segment_1:tsopt:ts_optimization_not_converged`) | [JSON Output](json-output.md#execution-and-requested-stage-completion) |
| UMA raises `BackendError` for `--uma-workers` above 1 with `--hessian-calc-mode Analytical` | Use `--uma-workers 1` for an analytical Hessian, or select `FiniteDifference` | {ref}`Performance <troubleshooting-performance>`, {ref}`Workers and analytical Hessians <workers-analytical-error>` |
| CUDA out of memory (`torch.cuda.OutOfMemoryError`) | Keep the default `FiniteDifference` Hessian, lower `--max-nodes` or use a smaller MLIP model (`--backend-model`), or move to a larger GPU; re-extract with a smaller `--radius` only as the last step | {ref}`GPU memory <troubleshooting-gpu-memory>` |
| TS optimization converged, but n_imag is not 1 | n_imag ≥ 2 (`TS imaginary-mode validation found n_imag=…`): add `--flatten` (`tsopt`, `opt`, and `all`). n_imag = 0 (`[tsopt] No imaginary mode detected.`): start from another candidate, or use `--refine-path` in `all` | {ref}`TS optimization <troubleshooting-ts>`, {ref}`When the TS search fails <ts-search-fails>` |
| TS optimization does not converge (`TS optimization did not converge`) | Check the TS candidate first, then switch the optimizer (`tsopt --opt-mode` / `all --opt-mode-post`), then reduce the step size or trust radius in YAML | {ref}`TS optimization <troubleshooting-ts>`, {ref}`When the TS search fails <ts-search-fails>` |
| IRC does not terminate | Standalone `irc`: reduce `--step-size`, raise `--max-cycles`. In `all`: `--irc-step-size` / `--irc-max-cycles`. Check the endpoints first | {ref}`IRC <troubleshooting-irc>` |
| Opt / TS optimization stops as `stalled` after an energy plateau (`TS optimization status is stalled`) | Treat it as not converged; inspect the final geometry and the force / step criteria, then retry with another threshold or optimizer setting | {ref}`max_cycles and plateau stops <troubleshooting-max-cycles>` |
| Minimum energy path (MEP) search (GSM / DMF) fails (`MEP optimization did not converge`) | Raise `--max-nodes` above the default 20, keep `--preopt` (on by default in `all`, `path-search`, and `path-opt`; off in `scan`, `scan2d`, and `scan3d`), try the other `--mep-mode` | {ref}`MEP search <troubleshooting-mep>` |
| `freq` stops with an error | Leave at least one atom movable | {ref}`freq errors <troubleshooting-freq>` |
| DFT SCF does not converge | Try `--no-dft-low-memory` (density fitting) or a level shift in YAML: `mf: {level_shift: 0.2}` under `dft.pyscf` (`dft` and `--dft`) or `calc.dft.pyscf` (`-b dft`) | [dft notes](dft.md#notes), [DFT backend notes](dft-backend.md#notes) |
| DFT runs out of GPU memory | Use a smaller basis, trim the model, or move to a larger GPU; if `--dft` runs out, run `pdb2reaction dft` separately | [DFT backend notes](dft-backend.md#notes) |
| **Installation & environment** | | |
| DMF import error (`cyipopt`), or `No module named 'dmf'` | `conda install -c conda-forge cyipopt` (`pydmf` is installed with `pdb2reaction`); if the import still fails with the default GPU backend, `pip install 'pydmf[torch]'` | {ref}`Installation / environment <installation-environment-problems>` |
| UMA model 401 / 403 or gated-repo error (`huggingface_hub.errors.GatedRepoError`) | Run `hf auth login` and accept the UMA model license | {ref}`Installation / environment <installation-environment-problems>` |
| `e3nn` / `fairchem-core` import conflict (MACE in the UMA env) | Use a dedicated environment for MACE | {ref}`Installation / environment <installation-environment-problems>` |
| `ORB backend requires orb-models and torch` (or the same for AIMNet2 / MACE) | Install the backend extra: `pip install "pdb2reaction[orb]"`; MACE goes in a separate env | {ref}`Installation / environment <installation-environment-problems>` |
| CUDA / GPU runtime mismatch | Check the GPU, the PyTorch build, and the driver together | {ref}`Installation / environment <installation-environment-problems>` |
| Plot export fails | Run `plotly_get_chrome -y` to install headless Chrome | {ref}`Installation / environment <installation-environment-problems>` |

## Preflight checklist

Before a long run, check that:

- A Hugging Face login is set up on this machine (needed for the default UMA model).
- Your input PDB/mmCIF structures contain hydrogens and element symbols.
- When you pass several PDBs, they share the same atoms in the same order.

---

(input-extraction-problems)=
## Input / extraction

### `Element symbols are missing in '...'`

- **Symptom**: `extract` stops with `Element symbols are missing in '...'. For PDB input, run pdb2reaction add-elem-info -i ... before extract`. `all` fills blank element columns itself before extraction, and stops with the same message when some atoms cannot be assigned.
- **Cause**: many PDBs leave the element column (columns 77–78) blank, and `extract` needs the elements to type the atoms. mmCIF input must provide `_atom_site.type_symbol`.
- **Fix**: fill the column with `add-elem-info` and rerun with the new file. For each atom listed under `[WARN] Could not confidently assign N atoms`, write its element symbol by hand, right-justified in columns 77–78.

  ```bash
  pdb2reaction add-elem-info -i input.pdb -o input_with_elem.pdb
  ```

### `[multi] Atom count mismatch` / `[multi] Atom order mismatch`

- **Symptom**: a run with several inputs stops with `[multi] Atom count mismatch between input #1 and input #2: ...` or `[multi] Atom order mismatch between input #1 and input #2.`
- **Cause**: the structures were prepared with different tools or settings, or the atom order changed after re-protonation or re-parametrization.
- **Fix**: regenerate **all** structures with the same protonation tool and settings. For MD snapshots, take every frame from the same trajectory and topology. Never reorder PDB atoms after the topology is built.

### The active-site model is empty or misses catalytic residues

- **Symptom**: the extracted model is smaller than expected, or catalytic residues are missing.
- **Cause**: the radius (`-r/--radius`, default 2.6 Å) is too small for this site, or `--exclude-backbone` removed too much.
- **Fix**: raise `--radius` (for example 2.6 → 3.5 Å), or add the residue. `--selected-resn 'A:TYR:44'` adds at least its side chain, and adding it to `-c` keeps it whole. See [Make the model larger](model-setup.md#make-the-model-larger) and [Residue selectors](cli-conventions.md#residue-selectors). If you passed `--exclude-backbone` and it trims too much, pass `--no-exclude-backbone`.

### Energies or barriers shift with model size

- **Symptom**: energies or barriers look unreasonable, or change a lot when the model grows.
- **Cause**: the extracted model is too small.
- **Fix**: enlarge the radius and check how the result depends on the model size and boundary.

  ```bash
  pdb2reaction extract -i complex.pdb -c 'SUB' -o model.pdb -r 4.0
  ```

### A modified residue is not truncated

- **Symptom**: an unlisted modified amino acid keeps its full backbone and gets no cap hydrogens.
- **Cause**: backbone truncation and cap-hydrogen placement need an entry in the amino-acid catalog. SEP, TPO, and MLY are already built in.
- **Fix**: register only the unlisted residue and give its nominal charge, for example `--modified-residue "XAA:0"`. A bare name assigns charge 0, so do not use it to re-register a charged built-in residue. If the backbone topology is unusual, [build the active-site model by hand](model-setup.md) and pass it directly to the downstream commands.

---

(charge-spin-problems)=
## Charge / spin

First check that the total charge and the multiplicity are correct for the target state, and that each residue key in `-l/--ligand-charge` exists in the structure. For important runs, give `-q/--charge` or `-l/--ligand-charge` and `-m` explicitly; the rules are in {ref}`Charge specification <charge-specification>`.

### `-q/--charge is required` / `Total charge could not be resolved`

- **Symptom**: a run on non-`.gjf` input stops with `-q/--charge is required unless the input is a .gjf template with charge metadata.` or, in `all`, `[all] Total charge could not be resolved.`
- **Cause**: without `-q/--charge`, the workflow resolves the charge from `-l/--ligand-charge` with PDB/mmCIF input (or XYZ/GJF with `--ref-pdb`), from YAML `calc.charge`, or from a `.gjf` template. None of them applied.
- **Fix**: give the charge and multiplicity, or a per-residue charge map with extraction.

  ```bash
  pdb2reaction path-search -i R.pdb P.pdb -q 0 -m 1
  pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
  ```

---

(installation-environment-problems)=
## Installation / environment

First confirm that the optional packages are installed in the active environment and that PyTorch sees the GPU. After a repair, check with `pdb2reaction --version` and `python -c "import torch; print(torch.cuda.is_available())"`, then run the command once with `--dry-run` to check the options and input before the full run.

| Symptom | Cause | Fix |
| --- | --- | --- |
| UMA download fails (`huggingface_hub.errors.GatedRepoError`, `401`, `403`) | No Hugging Face login, or the UMA model license is not accepted | Run `hf auth login` once per environment and machine, and accept the UMA model license on its Hugging Face page. On HPC, make sure compute nodes can write to the Hugging Face cache directory |
| `ORB backend requires orb-models and torch`, `AIMNet2 backend requires torch and aimnet`, or `Could not import mace.calculators because mace-torch is not installed` | The backend package is not installed in this environment | ORB: `pip install "pdb2reaction[orb]"`. AIMNet2: `pip install "pdb2reaction[aimnet]"`. MACE: a separate environment (next row) |
| `e3nn` / `fairchem-core` import conflict | MACE was installed into the UMA environment: `mace-torch` pins `e3nn==0.4.4`, while `fairchem-core` needs `e3nn>=0.5` | Use a dedicated conda environment for MACE: `pip uninstall -y fairchem-core && pip install 'mace-torch>=0.3.8'` |
| ORB import still fails after installing the extra | A package in the environment conflicts with `orb-models` | Run `python -m pip check`. Install PyG / `torch_scatter` only when the error names it; current `orb-models` does not require it |
| `torch.cuda.is_available()` returns `False`, or a CUDA runtime error at import | The PyTorch build does not match the GPU or driver of the node | Check the assigned GPU, the installed wheel, and the driver with `nvidia-smi`, `python -m torch.utils.collect_env`, and `python -m pip check`; the boolean alone does not tell the cause. The `CUDA Version` that `nvidia-smi` shows is the newest CUDA the driver supports; install a wheel for that CUDA or an older one (`cu126`, `cu130`, `cu132`). A local CUDA toolkit is not needed |
| `--mep-mode dmf` fails with `DMF mode (--mep-mode dmf) requires ase, cyipopt, and pydmf>=1.2`, or `No module named 'dmf'` | `cyipopt` is missing, or `pydmf` lacks its GPU part | Run `conda install -c conda-forge cyipopt`, preferably before installing `pdb2reaction`. `pydmf` is installed with `pdb2reaction`; if the import still fails with the default GPU DMF backend, run `pip install 'pydmf[torch]'`, as the error message suggests |
| Plot export fails (Plotly / Chrome) | Headless Chrome is missing | Run `plotly_get_chrome -y` once |

### DMF is unusually slow inside IPOPT

If IPOPT/MUMPS uses multithreaded BLIS, nested threads can cause long waits.
Set `BLIS_NUM_THREADS=1` before starting Python or the CLI, for example in
the job script. Leave the other thread settings, such as `OMP_NUM_THREADS`, unchanged.
Manual `BLIS_JC_NT`, `BLIS_PC_NT`, `BLIS_IC_NT`, `BLIS_JR_NT`, or
`BLIS_IR_NT` settings override this limit; remove them from that job's
configuration. Restart an existing notebook kernel after changing the settings.
See [BLIS thread controls](https://github.com/flame/blis/blob/2.0/docs/Multithreading.md).

---

(calculation-convergence-problems)=
## Calculation / convergence

First check the TS candidate: a successful TS optimization gives one imaginary mode along the reaction coordinate. A mode counts as imaginary when ν < −5.00 cm⁻¹; YAML `freq.zero_cutoff_cm` changes the cutoff. When the TS optimization does not converge, read {ref}`TS optimization <troubleshooting-ts>`; when it converges but n_imag is not 1, read {ref}`When the TS search fails <ts-search-fails>`.

(troubleshooting-max-cycles)=
### Optimizer reaches `max_cycles` with `max(force)` slightly above threshold

- **Symptom**: the optimizer runs to `max_cycles`, and the final summary shows `max(force)` or `rms(force)` just above the selected threshold while the energy has stopped changing.
- **Cause**: MLIP force noise and flatness can keep a force threshold out of reach; the level depends on the backend, model, precision, system, and hardware.
- **Fix**: inspect the final geometry and the forces, then rerun:
  - Loosen the convergence preset with `--thresh gau_loose`. The defaults are `gau` for `opt` and [`baker`](tsopt.md#how-it-works) for `tsopt`.
  - To end such a run early instead of at `--max-cycles`, add `--stop-plateau`. The run then stops as `stalled` (**not converged**) once the energy range over the last `--stop-plateau-window` steps (default 50) falls below `--stop-plateau-thresh` (default `1×10⁻⁴ au`).

(troubleshooting-ts)=
### TS optimization does not converge / multiple imaginary modes remain

- **Symptom**: the TS optimization runs many cycles without converging, or n_imag is 2 or more after the optimization.
- **What is reported**: a TS optimization (RS-P-RFO, RS-I-RFO, TRIM, or Dimer) that reaches max cycles without converging does not compute the Hessian, so no n_imag is reported. A run stopped on an energy plateau always computes the Hessian and reports n_imag.
- **Fix when the optimization does not converge**: try the following in order.
  1. Switch the optimizer between RS-P-RFO (the default) and the Dimer method: `tsopt --opt-mode hess` / `dimer`, or `all --opt-mode-post hess` / `grad` (Dimer).
  2. Reduce the step size in YAML: `rsirfo.trust_radius` / `trust_min` / `trust_max` for RS-P-RFO, RS-I-RFO, and TRIM, or `hessian_dimer.lbfgs.max_step` for Dimer; see [YAML Reference](yaml-reference.md#ts-optimization-sections).
  3. Start from another candidate, such as a better HEI (highest-energy image) of the path.
- **Fix when n_imag ≥ 2 remains**: re-optimize with `--flatten`, or tighten the convergence preset from the default `baker` to `gau_tight` or `gau_vtight` (`tsopt --thresh` or `all --thresh-post`). {ref}`When the TS search fails <ts-search-fails>` lists other moves: `--refine-path`, splitting the reaction into stages, a new starting structure.

(troubleshooting-irc)=
### IRC does not terminate properly

An IRC that stops before it converges is still usable when the endpoint optimizations reach the intended R and P, so check the optimized endpoints first.

- **Symptom**: the IRC stops before it reaches a clear minimum, or the energy oscillates and the gradient norm stays large.
- **Cause**: the step is too large for this surface, the cycle limit is too low, or the starting structure has more than one imaginary mode.
- **Fix**:
  - Standalone `irc`: `--step-size 0.05` (default 0.10 bohr) and, if needed, `--max-cycles 200` (default 125).
  - `all`: `--irc-step-size 0.05` and, if needed, `--irc-max-cycles 200`.
  - Confirm that the starting structure has n_imag = 1.
  - To ignore the physical stop criteria and trace to the cycle limit, use `irc --never-stop` or `all --irc-never-stop`, then inspect the trajectory and endpoints.

(troubleshooting-mep)=
### MEP search (GSM / DMF) fails or misses bonds

- **Symptom**: the MEP search ends without a usable path, or misses an expected bond change.
- **Fix**:
  - Raise `--max-nodes` (default 20) to 30 or 40 for complex reactions.
  - Keep endpoint preoptimization on (the default); remove `--no-preopt` if you passed it.
  - Try the other method: `--mep-mode dmf` ↔ `gsm`.
  - Tune bond-change detection with YAML `bond.bond_factor` and `bond.delta_fraction`.

(troubleshooting-freq)=
### `freq` stops with an error

- **Every atom is frozen**: with nothing movable there is no vibration to analyze, and `freq` stops with an error. Check `--freeze-atoms`, YAML `geom.freeze_atoms`, and the cap-hydrogen parents frozen by `--freeze-links`; see {ref}`Freeze atoms and restrain distances <freeze-atoms-and-restraints>`.
- **`--uma-workers` above 1 with `--hessian-calc-mode Analytical`**: see {ref}`Performance <troubleshooting-performance>`.

---

(troubleshooting-performance)=
## Performance / stability tips

- **`workers > 1`** — may improve UMA throughput, depending on the hardware and workload, but the parallel predictor has no analytical Hessian. An explicit `Analytical` request raises `BackendError` (a `RuntimeError` subclass) with `Analytical Hessian cannot be combined with UMA workers>1`; use `--uma-workers 1` for an analytical Hessian, or select `FiniteDifference`. To spread workers over several nodes, use the job script in the [HPC example](hpc-example.md).
- **Large systems** — make a chemically justified smaller active-site model and run radius/boundary sensitivity checks; multi-GPU support is backend- and workflow-specific, so do not assume that increasing GPU count reduces memory per worker.
- **DFT scratch on HPC** — if PySCF/GPU4PySCF uses temporary disk for the chosen calculation, point `PYSCF_TMPDIR` at a filesystem with verified capacity and performance. Do not assume that node-local `/tmp`, `$PBS_O_WORKDIR`, or another shared path is suitable at every site.

(troubleshooting-gpu-memory)=
## GPU memory (VRAM) requirements

VRAM depends on the MLIP model, precision, Hessian mode, and frozen atoms, not on the atom count alone, so measure a representative pilot with the same setup. On `torch.cuda.OutOfMemoryError`, try the following in order:

1. Keep or switch to the default `--hessian-calc-mode FiniteDifference`; an analytical Hessian normally has a larger memory peak.
2. Lower `--max-nodes`, or use a smaller MLIP model (`--backend-model`). In `opt` and `scan`, keep {ref}`--opt-mode grad <opt-mode-semantics>` (L-BFGS, no Hessian) instead of `hess`.
3. Move to a GPU with more memory.
4. Reduce the cluster model, only after checking that the required residues remain and the boundary is well placed.

## How to report an issue

Include the exact command, `summary.log` (or console output), the smallest reproducing inputs, and your env (OS / Python / CUDA / PyTorch).

## See also

- [Tips for studying reaction mechanisms](mechanism-tips.md) — what to try when the TS search fails
- [Installation](installation.md) — environment setup and optional backends
- [MLIP Backends](backends.md) — choosing a backend
- [Building the cluster model](model-setup.md) — check, trim, or enlarge the active-site model
