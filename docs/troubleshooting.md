# Troubleshooting

For a symptom-first index see [Common Error Recipes](recipes-common-errors.md).

## Preflight checklist

Before a long run, verify:

- A Hugging Face token is set up on this machine (required for the default UMA backend to download model weights).
- Your input PDB(s) contain hydrogens and element symbols.
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
| `[orb]` install fails building **torch_scatter** (`No module named 'torch'`) | torch_scatter ships no PyPI binary wheel (only an sdist) → source-build fails under PEP517 build isolation. Install from PyG's prebuilt-wheel index matching your torch+CUDA tag: `pip install "pdb2reaction[orb]" -f https://data.pyg.org/whl/torch-2.8.0+cu129.html`. Fallback (CUDA toolchain present): `pip install torch_scatter --no-build-isolation`. |
| `torch.cuda.is_available()` returns `False` | Install PyTorch matching your cluster CUDA runtime; verify `nvidia-smi` + `python -c "import torch; print(torch.version.cuda, torch.cuda.is_available())"`. |
| DMF fails (`--mep-mode dmf`: `cyipopt` missing or `No module named 'dmf'`) | `conda install -c conda-forge cyipopt` before installing `pdb2reaction`. `pydmf` ships as a dep; if missing, `pip install --force-reinstall pdb2reaction`. |
| Plot export fails (Plotly / Chrome) | `plotly_get_chrome -y`. |

---

(calculation-convergence-problems)=
## Calculation / convergence

### Optimizer reaches `max_cycles` with `max(force)` slightly above threshold

MLIP gradients carry a stochastic noise floor (~10⁻⁴ Ha/Bohr) that can exceed the `baker` force-convergence criterion (the `baker` preset's max-force threshold, 3×10⁻⁴ au) even after the geometry has effectively converged. The **energy-plateau fallback** handles this automatically: `opt.energy_plateau: true` declares convergence when the energy range over the last `opt.energy_plateau_window` (default 50) steps falls below `opt.energy_plateau_thresh` (default `1×10⁻⁴ au ≈ 0.06 kcal/mol`).

To override, do one of:

- Loosen the force threshold (`--thresh gau` default / `--thresh gau_loose`).
- Tune `opt.energy_plateau_thresh` / `opt.energy_plateau_window` in YAML.
- Disable the fallback with `opt.energy_plateau: false` in YAML.

The fallback is **skipped for chain-of-states optimizers** (optimizers that move a whole chain of path images together), namely the `path-opt` and `path-search` Growing String Method (GSM) / Direct Max Flux (DMF) stages.

### TS optimization does not converge / multiple imaginary modes remain

Try the following, in order:

1. Switch the optimizer mode: `--opt-mode grad` (Dimer Method) ↔ `--opt-mode hess` (Restricted-Step Partitioned-RFO, RS-P-RFO).
2. Add `--flatten` (available on standalone `tsopt` / `opt` / `pdb2reaction all`).
3. Raise the cycle limit: `--max-cycles 20000` (standalone `tsopt`) or `--tsopt-max-cycles 20000` (`all`).
4. Tighten the force threshold: `--thresh baker` / `gau_tight`.
5. Reduce step sizes / trust radii via YAML: `lbfgs.max_step`, `hessian_dimer.lbfgs.max_step`, `rfo.trust_radius` / `trust_min` / `trust_max`, the `rsirfo` block — see [YAML Reference](yaml-reference.md).

### IRC does not terminate properly

Reduce `--step-size 0.05` (default 0.10 bohr); raise `--max-cycles 200`; confirm the TS candidate has exactly one imaginary frequency before IRC (detection cutoff `hessian_dimer.neg_freq_thresh_cm`, default 5 cm⁻¹).

### MEP search (GSM / DMF) fails or misses bonds

The minimum energy path (MEP) search can stall or skip an expected bond change. Try the following:

- Raise `--max-nodes 30` / `40` for complex reactions.
- Add `--preopt`.
- Try the alternate method: `--mep-mode dmf` ↔ `gsm`.
- Tune `bond.bond_factor` / `bond.delta_fraction` in YAML.

---

## Performance / stability tips

- **OOM** — shrink active-site model (`--radius`), lower `--max-nodes`, use lighter `--opt-mode grad`.
- **Analytical Hessian** — keep the default `FiniteDifference`; only set `--hessian-calc-mode Analytical` if you have 16 GB+ VRAM, and note that on 16 GB it is feasible only up to ~200 atoms (see the VRAM table below).
- **`workers > 1`** — improves UMA throughput on HPC, but is incompatible with the analytical Hessian (it raises a `RuntimeError`); pass `--hessian-calc-mode FiniteDifference` explicitly, or run with `--workers 1` for Hessian-based modes.
- **Large systems (1000+ atoms)** — extract a smaller active-site model (`--radius 2.5`) or run multi-GPU.
- **DFT scratch on HPC (~hundreds of atoms)** — PySCF / GPU4PySCF spills integrals to `$PYSCF_TMPDIR` (or `$TMPDIR` / `/tmp` if unset). Node-local `/tmp` is often small / `tmpfs`-backed and can fill up mid-SCF. Set `export PYSCF_TMPDIR="$PBS_O_WORKDIR"` before launching `dft`.

## Choosing a backend

Order-of-magnitude per-step L-BFGS cost on small-to-medium cluster models, measured on a 16 GB consumer GPU (not a rigorous benchmark):

| Backend (`-b/--backend`, model id) | s/step | VRAM | Notes |
|---|---|---|---|
| `uma` (`uma-s-1p2`, default) | 0.03 | ~2 GB | Fast, good for exploration. |
| `uma` (`uma-m-1p1`) | 0.22 | ~8 GB | Medium model, higher VRAM. |
| `mace` (`MACE-OMOL-0`) | 0.37 | ~4 GB | Separate env (`e3nn` conflict with fairchem-core). |
| `orb` (`orb_v3_conservative_omol`) | 0.02 | ~2 GB | Fastest; see caveat. |

Start with the default UMA model for rapid screening, then cross-check key results with MACE or `uma-m-1p1`. For SAM-dependent S~N~2 / methyltransfer chemistries MACE and default UMA are complementary — try both when one produces a suspect TS. The Orb backend often identifies the correct reaction coordinate but tends to converge to TS geometries with extra small imaginary modes — good for initial *mechanism recovery*, but for quantitative kinetics or freq analysis re-score Orb geometries with UMA / MACE / DFT.

## GPU memory (VRAM) requirements

Approximate VRAM by system size:

| Atoms | L-BFGS opt | Hessian (analytical) | Hessian (finite diff) |
|---|---|---|---|
| 50 | ~2 GB | ~3 GB | ~2 GB |
| 100 | ~3 GB | ~6 GB | ~3 GB |
| 200 | ~4 GB | ~12 GB | ~4 GB |
| 500 | ~6 GB | OOM on 16 GB | ~6 GB |

On `torch.cuda.OutOfMemoryError`: switch to `--hessian-calc-mode FiniteDifference`, reduce `--radius`, or pick a smaller model (`calc.model: uma-s-1p2` instead of `uma-m-1p1` in YAML).

## How to report an issue

Include the exact command, `summary.log` (or console output), the smallest reproducing inputs, and your env (OS / Python / CUDA / PyTorch).
