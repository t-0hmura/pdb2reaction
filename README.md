# `pdb2reaction`: End-to-End Reaction-Path Elucidation from PDB Structures Using Machine-Learning Interatomic Potentials

## Overview

<img src="https://raw.githubusercontent.com/t-0hmura/pdb2reaction/main/docs/overview.png" alt="pdb2reaction workflow overview" width="90%">

`pdb2reaction` is a Python CLI for elucidating **enzymatic reaction pathways** from **PDB or mmCIF structures** using machine-learning interatomic potentials (MLIPs). Given (i) two or more reaction-ordered structures, (ii) one structure with `--scan-lists`, or (iii) one TS candidate with `--tsopt`, it can run an **MEP search** and optionally chain **TS optimization → IRC → thermochemical correction → DFT single-point**. Active-site extraction is performed only when `-c/--center` is supplied; otherwise the PDB/mmCIF/XYZ/GJF model is used as-is. Each stage is also exposed as an [individual subcommand](#cli-subcommands).

Test a reaction mechanism in a single command:

```bash
# Multi-PDB mode (R + P endpoints → MEP, with TS optimization + thermo)
pdb2reaction all -i R.pdb P.pdb -c 'LIG' -l 'LIG:-1' --tsopt --thermo
```

Inputs are not limited to full enzyme PDBs: mmCIF is accepted directly, including multi-character chains and large residue IDs. You can also pass a small molecule as `.xyz` / `.gjf`, or a cluster model you built yourself as PDB/mmCIF, and omit `--center/-c` to skip extraction.

> **Prerequisites:** PDB/mmCIF inputs must already contain hydrogens; reaction-ordered structures must share the same atom identities and order (only coordinates differ). Small-molecule `.xyz` / `.gjf` inputs work when `--center/-c` and `--ligand-charge/-l` are omitted.

## Related tools

| Tool | Use case |
|---|---|
| [**mlmm-toolkit**](https://github.com/t-0hmura/mlmm_toolkit) | **ML/MM ONIOM** with the full protein environment; automates MM parameterization and ML-region assignment from a single PDB. |
| [**uma_pysis**](https://github.com/t-0hmura/uma_pysis) | Lightweight **YAML-driven UMA–pysisyphus interface** for quick/exploratory reaction-mechanism studies (GS / TS / IRC / ΔG). |

> `pdb2reaction` bundles a GPU-optimized pysisyphus fork that is **not** compatible with upstream pysisyphus — do not install it into an environment that already has upstream pysisyphus.

## Documentation

- [Getting Started](docs/getting-started.md) · [mmCIF and large structures](docs/cif.md) · [Installation](docs/installation.md) · [Examples](examples/) · [Troubleshooting](docs/troubleshooting.md)
- [YAML Reference](docs/yaml-reference.md) · [JSON Output Schema](docs/json-output.md)
- Full site: <https://t-0hmura.github.io/pdb2reaction/>

## Colab GUI

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/t-0hmura/pdb2reaction/blob/main/examples/pdb2reaction_colab.ipynb)

[`examples/pdb2reaction_colab.ipynb`](examples/pdb2reaction_colab.ipynb) provides
an ordered PDB/mmCIF/XYZ/GJF input queue and a compact
**Input → Viewer → Options (optional) → Results** workspace. The Viewer places
the workflow and relevant selectors beside a 4:3 molecular view; an exact pick
shows the residue as thick element-colored sticks and keeps the atom distinct
with a dark halo. Key and advanced options stay collapsed until needed, each
control has click-to-open help, and inapplicable CLI options remain hidden.
Validation covers the exact editable command and its current input files.
Results are scoped to the current run and preview supported structures,
diagrams, tables, documents, and trajectories. Setup installs the pinned PyPI
release and matching-tag examples into the user's GPU runtime; MACE, UMA, and
ORB are installed one backend per runtime.

## System requirements

| Component | Requirement |
|---|---|
| OS / Python | Linux recommended. Python >= 3.11. |
| GPU / CUDA / VRAM | CUDA-capable NVIDIA GPU recommended for production, with a compatible driver and official PyTorch 2.8 CUDA wheel (`cu126`, `cu128`, or `cu129`). Required VRAM is backend/model/system/workflow dependent; pilot the real calculation. |
| RAM / Disk | Size for the selected environment, model cache, structures, trajectories, and DFT scratch; no atom-count-only minimum is reliable. |

CPU-only execution works but is usually much slower; benchmark the selected backend/model. Full requirement and tuning details: [docs/installation.md](docs/installation.md).

## Installation

```bash
# 1. CUDA-enabled PyTorch (choose the official 2.8 wheel for your driver/GPU)
pip install 'torch==2.8.0' --index-url https://download.pytorch.org/whl/cu126

# 2. pdb2reaction (editable from a local clone, or `pip install pdb2reaction`)
pip install -e .

# 3. Authenticate Hugging Face once (only required for the default UMA backend)
#    Accept the FAIR Chemistry License v1 at https://huggingface.co/facebook/UMA, then:
hf auth login                               # interactive
# OR: export HF_TOKEN=hf_xxx && hf auth login --token "$HF_TOKEN" --add-to-git-credential   # CI / HPC
```

**Optional extras** (install only what you need):

| Extra | Adds |
|---|---|
| `[orb]` / `[aimnet]` | Orb / AIMNet2 MLIP backend (`-b orb` / `-b aimnet2`) — *not* HF-gated |
| `[dft]` | PySCF + GPU4PySCF single-point DFT (`--dft` / `pdb2reaction dft`) |
| `[mcp]` | Model Context Protocol server for agent clients |

The MACE backend (`-b mace`) is **not** a pip extra: `mace-torch` pins `e3nn==0.4.4`, which conflicts with `fairchem-core`'s `e3nn>=0.5` (UMA), so it needs a dedicated environment — `pip uninstall -y fairchem-core && pip install mace-torch` (see [docs/installation.md](docs/installation.md)).

CUDA module loads, alternative-backend recipes, DMF/`cyipopt` setup, Plotly Chromium, and HPC job-script templates: [docs/installation.md](docs/installation.md) and [docs/hpc-example.md](docs/hpc-example.md).

## Quick Examples

Examples use GPP C6-methyltransferase BezA ([Tsutsumi et al., *Angew. Chem. Int. Ed.* 2022, 61, e202111217](https://doi.org/10.1002/anie.202111217)) — runnable MEP and scan commands are in [`examples/run.sh`](examples/run.sh).

```bash
# Multi-structure MEP (R + P → MEP, with TS + thermochemistry)
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --out-dir result_mep

# Scan mode (single structure → staged bond scan → MEP)
pdb2reaction -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","GPP 321 C7",1.60)]' --tsopt --thermo --out-dir result_scan

# TS-only validation (single TS candidate → tsopt → IRC → freq)
pdb2reaction -i TS_candidate.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo --out-dir result_tsonly
```

`pdb2reaction` can also be used to investigate reaction mechanisms of **small molecules** and **user-defined cluster models**.

```bash
# Small molecule (gas-phase): .xyz / .gjf input — omit -c, set charge with -q
pdb2reaction -i reactant.xyz product.xyz -q 0 --tsopt --thermo --out-dir result_small

# Your own cluster model (already-trimmed PDB): omit -c to use it as-is
pdb2reaction -i cluster_R.pdb cluster_P.pdb -q 0 --tsopt --thermo --out-dir result_cluster
```

For hand-built clusters, standardize backbone ends at Cα (`CA`), place other
boundaries on aliphatic C–C single bonds whenever possible, avoid cutting
peptide/polar/conjugated/metal bonds, and use the identical atom order and cap
topology for every state. See [the cluster-boundary checklist](docs/extract.md#building-or-auditing-a-cluster-model-manually).

Per-stage walkthrough (`extract` → `opt` → `path-opt` → `tsopt` → `freq` → `irc` → `dft`): [docs/getting-started.md](docs/getting-started.md) and [docs/quickstart-all.md](docs/quickstart-all.md).

## Output

A non-dry `all` run writes the deliverables reached by its enabled stages to
`--out-dir` (default `./result_all/`):

- `segments/seg_NN/{reactant,ts,product}.*` — the canonical R / TS / P structures to cite
- `mep_trj.xyz` (plus `mep.pdb` when topology is available and `mep.cif` for bridged inputs) — the merged reaction path in MEP/scan-list modes
- `energy_diagram_MEP.png` — MEP diagram when MEP construction and static-image export succeed
- `summary.log` (human-readable) / `summary.json` (machine-readable)

Pipeline scratch lives under `_work/` (safe to delete). Full layout and filename conventions: [docs/output-layout.md](docs/output-layout.md).

## CLI Subcommands

| Subcommand | Role | Doc |
|---|---|---|
| `all` (default) | End-to-end: extract → MEP → TS → IRC → freq → DFT | [all](docs/all.md) |
| `extract` | Build active-site cluster model | [extract](docs/extract.md) |
| `fix-altloc` | Resolve PDB alternate conformations | [fix-altloc](docs/fix-altloc.md) |
| `add-elem-info` | Repair PDB element columns (77–78) | [add-elem-info](docs/add-elem-info.md) |
| `opt` | Geometry optimization (L-BFGS / RFO) | [opt](docs/opt.md) |
| `tsopt` | TS optimization (Dimer / RS-P-RFO) | [tsopt](docs/tsopt.md) |
| `path-opt` | MEP via GSM or DMF | [path-opt](docs/path-opt.md) |
| `path-search` | Recursive MEP search with refinement | [path-search](docs/path-search.md) |
| `scan` / `scan2d` / `scan3d` | 1D / 2D / 3D bond-distance scans | [scan](docs/scan.md) · [scan2d](docs/scan2d.md) · [scan3d](docs/scan3d.md) |
| `freq` | Vibrational analysis + thermochemistry | [freq](docs/freq.md) |
| `irc` | IRC (EulerPC) | [irc](docs/irc.md) |
| `dft` | Single-point DFT (GPU4PySCF / PySCF) | [dft](docs/dft.md) |
| `sp` | Single-point MLIP energy / forces / Hessian | [sp](docs/sp.md) |
| `bond-summary` | Compare structures, report bond changes | [bond-summary](docs/bond-summary.md) |
| `trj2fig` / `energy-diagram` | Energy plot / R→TS→P diagram | [trj2fig](docs/trj2fig.md) · [energy-diagram](docs/energy-diagram.md) |

## Getting Help

```bash
pdb2reaction --help                       # top-level
pdb2reaction <subcmd> --help              # core options
pdb2reaction <subcmd> --help-advanced     # full option set
```

Issues: <https://github.com/t-0hmura/pdb2reaction/issues>.

## Citation

```bibtex
@misc{ohmura2026pdb2reaction,
  author = {Ohmura, Takuto and Sato, Hajime and Terada, Tohru},
  title  = {pdb2reaction: End-to-End Reaction-Path Elucidation from PDB Structures Using Machine-Learning Interatomic Potentials},
  year   = {2026}, doi = {10.26434/chemrxiv.15003538/v1}, note = {ChemRxiv preprint}
}
```

## Agent Skills

Agent Skills for Claude Code / Codex / Cursor etc. in [`skills/`](skills/) — copy into your project's skill location (e.g. `.claude/skills/`) to let an agent drive `pdb2reaction` workflows and subcommands.

## Known limitations

- **MACE + UMA cannot coexist** (`e3nn` version conflict). Use separate conda envs.
- **DFT single-point cost** depends strongly on basis, functional, grid, elements, and hardware; pilot one representative structure before batching.
- **Every backend's TS** requires an independent frequency calculation and IRC connectivity check. ORB `fp32`/TF32 is unsuitable for trusted finite-difference Hessians; the pdb2reaction ORB default is fp64.
- **CPU-only execution** is supported but usually much slower than GPU.

## Contributing

Issues and pull requests are welcome — see [CONTRIBUTING.md](CONTRIBUTING.md).

## License

GNU General Public License v3 (GPL-3.0).
