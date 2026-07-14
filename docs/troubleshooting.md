# Troubleshooting

For a symptom-first index see [Common Error Recipes](recipes-common-errors.md).

## Preflight checklist

Before a long run, verify:

- A Hugging Face token is set up on this machine (required for the default UMA backend to download model weights).
- Your input PDB/mmCIF structures contain hydrogens and element symbols.
- When you pass multiple PDBs, they share the same atoms in the same order.

---

(input-extraction-problems)=
## Input / extraction

| Symptom | Cause | Fix |
|---|---|---|
| `Element symbols are missing in '...'. Please run pdb2reaction add-elem-info` | Many PDBs leave the element column (cols 77–78) blank; `extract` needs them for atom typing. | `pdb2reaction add-elem-info -i input.pdb -o input_with_elem.pdb`, then re-run. |
| `[multi] Atom count mismatch` / `[multi] Atom order mismatch` | Inputs prepared by different tools / settings, or atom order changed after re-protonation / re-parametrization. | Regenerate **all** structures with the same protonation tool + settings. For MD ensembles, extract frames from the same trajectory + topology. Never reorder PDB atoms after building topology. |
| Active-site model empty / catalytic residues missing | Default radius too small. | Increase `--radius` (e.g. 2.6 → 3.5 Å); force-include with `--selected-resn 'A:123,B:456'` (see {ref}`selected-resn-takes-ids`); if `--exclude-backbone` over-trims, pass `--no-exclude-backbone`. |
| Unreliable energies / barriers shifting with model size | Extracted model too small. | Increase `-r` (e.g. `pdb2reaction extract -i complex.pdb -c 'SUB' -o model.pdb -r 4.0`). |
| Non-standard residues not truncated (SEP / TPO / MLY / D-amino acids) | Backbone truncation + cap-H placement only apply to known three-letter codes. | `--modified-residue "SEP,TPO,MLY"`. If insufficient (unusual backbone topology), build the active-site model manually and pass it (`-i model.pdb`) directly to downstream subcommands. |

---

(charge-spin-problems)=
## Charge / spin

Most stages need a net charge when the input is not `.gjf`. If you omit `-q/--charge`, the workflow tries `--ligand-charge/-l` (PDB only) or a `.gjf` template; if neither resolves, it errors out.

```bash
pdb2reaction path-search -i R.pdb P.pdb -q 0 -m 1
pdb2reaction         -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'    # extraction route
```

---

(installation-environment-problems)=
## Installation / environment

| Symptom | Fix |
|---|---|
| UMA download fails / HF auth missing (`huggingface_hub.errors.GatedRepoError`, `401`, `403`) | `hf auth login` once per env / machine; accept the UMA model license on the HF page. On HPC, ensure HF cache dir is writable from compute nodes. |
| `ImportError: orb-models is required` (or similar for AIMNet2 / MACE) | For ORB: `pip install "pdb2reaction[orb]"`. For AIMNet2: `pip install "pdb2reaction[aimnet]"`. MACE installs into a separate env. |
| ORB backend import fails | Install the current extra with `pip install "pdb2reaction[orb]"`, then run `python -m pip check`. Do not add PyG/`torch_scatter` instructions unless the actual error names that independently installed package; current `orb-models` does not require it. |
| `torch.cuda.is_available()` returns `False` | Verify the assigned GPU, installed wheel, driver and environment with `nvidia-smi`, `python -m torch.utils.collect_env`, and `python -m pip check`; do not infer the cause from the boolean alone. |
| DMF fails (`--mep-mode dmf`: `cyipopt` missing or `No module named 'dmf'`) | `conda install -c conda-forge cyipopt` before installing `pdb2reaction`. `pydmf` ships as a dep; if missing, `pip install --force-reinstall pdb2reaction`. |
| Plot export fails (Plotly / Chrome) | `plotly_get_chrome -y`. |

---

(calculation-convergence-problems)=
## Calculation / convergence

### Optimizer reaches `max_cycles` with `max(force)` slightly above threshold

Measured MLIP force noise/flatness can prevent a selected force threshold from
being reached; the level depends on backend, model, precision, system, and
hardware stack. The **energy-plateau fallback** handles this case:
`opt.energy_plateau: true` declares minimum optimization convergence when the
energy range over the last `opt.energy_plateau_window` (default 50) steps falls
below `opt.energy_plateau_thresh` (default `1×10⁻⁴ au`). Hessian-family TS
optimizers additionally require exact first-order-saddle validation.

To override, do one of:

- Loosen the force threshold (`--thresh gau` default / `--thresh gau_loose`).
- Tune `opt.energy_plateau_thresh` / `opt.energy_plateau_window` in YAML.
- Disable the fallback with `opt.energy_plateau: false` in YAML.

The fallback is **skipped for chain-of-states optimizers** (optimizers that move a whole chain of path images together), namely the `path-opt` and `path-search` Growing String Method (GSM) / Direct Max Flux (DMF) stages.

### TS optimization does not converge / multiple imaginary modes remain

Try the following, in order:

1. Switch the optimizer mode: `--opt-mode grad` (Dimer Method) ↔ `--opt-mode hess` (Restricted-Step Partitioned-RFO, RS-P-RFO).
2. Add `--flatten` (available on standalone `tsopt` / `opt` / `pdb2reaction all`).
3. If the HEI came from a coarse MEP, retry `all` with `--refine-path`. This may split a poor path into unnecessary segments and multiply cost, so it is off by default and should follow inspection of the coarse MEP.
4. Raise the cycle limit: `--max-cycles 20000` (standalone `tsopt`) or `--tsopt-max-cycles 20000` (`all`).
5. Tighten the force threshold: `--thresh baker` / `gau_tight`.
6. Reduce step sizes / trust radii via YAML: `lbfgs.max_step`, `hessian_dimer.lbfgs.max_step`, `rfo.trust_radius` / `trust_min` / `trust_max`, the `rsirfo` block — see [YAML Reference](yaml-reference.md).

### IRC does not terminate properly

Reduce `--step-size 0.05` (default 0.10 bohr), especially when a branch stops after only a few frames; raise `--max-cycles 200`; confirm the TS candidate has exactly one imaginary frequency before IRC (detection cutoff `hessian_dimer.neg_freq_thresh_cm`, default 5 cm⁻¹). If a small uphill/flat shoulder is intentional, add opt-in `--never-stop`; it ignores only energy-based stops, so inspect the resulting trajectory and endpoints.

### MEP search (GSM / DMF) fails or misses bonds

The minimum energy path (MEP) search can stall or skip an expected bond change. Try the following:

- Raise `--max-nodes 30` / `40` for complex reactions.
- Add `--preopt`.
- Try the alternate method: `--mep-mode dmf` ↔ `gsm`.
- Tune `bond.bond_factor` / `bond.delta_fraction` in YAML.

---

## Performance / stability tips

- **OOM** — reduce the active-site model only after checking that required residues remain, lower `--max-nodes`, or use `--opt-mode grad` when a full-Hessian optimizer is not required.
- **Analytical Hessian** — keep the portable `FiniteDifference` default until the selected backend/model and representative atom count have been piloted. Analytical autograd normally has a larger memory peak, but no universal atom-count or VRAM cutoff applies.
- **`workers > 1`** — may improve UMA throughput, depending on the hardware and workload, but the parallel predictor has no analytical Hessian. An explicit `Analytical` request raises `BackendError` (a `RuntimeError` subclass); use `--workers 1` for an analytical Hessian, or select `FiniteDifference`.
- **Large systems** — make a chemically justified smaller active-site model and run radius/boundary sensitivity checks; multi-GPU support is backend- and workflow-specific, so do not assume that increasing GPU count reduces memory per worker.
- **DFT scratch on HPC** — if PySCF/GPU4PySCF uses temporary disk for the chosen calculation, point `PYSCF_TMPDIR` at a filesystem with verified capacity and performance. Do not assume that node-local `/tmp`, `$PBS_O_WORKDIR`, or another shared path is suitable at every site.

## Choosing a backend

Do not infer speed, memory use, or chemical accuracy from the backend name alone.
Record the backend, model, precision, package versions, GPU, atom count, and
workflow options, then pilot the actual calculation. UMA is the default and is
the only built-in backend with multi-worker inference. ORB defaults to fp64
because explicit fp32 uses reduced `float32-high`/TF32 matmul and can make
finite-difference Hessians noisy. MACE currently needs a dedicated environment
because its supported package stack conflicts with `fairchem-core`. AIMNet2
charge/spin coverage depends on the selected model. For every backend, accept a
TS only after optimizer convergence, one significant intended imaginary mode,
and correct IRC endpoints; use cross-backend or DFT checks when the scientific
claim requires them.

## GPU memory (VRAM) requirements

VRAM does not scale from atom count alone: model architecture, neighbor count,
precision, Hessian mode, frozen degrees of freedom, and software versions all
matter. Measure a representative pilot with the same setup and leave headroom
for transient peaks. On `torch.cuda.OutOfMemoryError`, keep or switch to
`--hessian-calc-mode FiniteDifference`, use a smaller supported model, or reduce
the cluster only after checking its chemistry and boundary placement.

## How to report an issue

Include the exact command, `summary.log` (or console output), the smallest reproducing inputs, and your env (OS / Python / CUDA / PyTorch).
