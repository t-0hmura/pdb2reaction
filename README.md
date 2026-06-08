# `pdb2reaction`: End-to-End Reaction-Path Elucidation from PDB Structures Using Machine-Learning Interatomic Potentials

## Overview

<img src="./docs/overview.png" alt="pdb2reaction workflow overview" width="90%">

`pdb2reaction` is a Python CLI for elucidating **enzymatic reaction pathways** from **PDB structures** using machine-learning interatomic potentials (MLIPs). Given (i) two or more PDB files (R → ... → P), (ii) one PDB with `--scan-lists`, or (iii) one TS candidate with `--tsopt`, it extracts an **active-site cluster model**, runs an **MEP search**, and optionally chains **TS optimization → IRC → frequencies → DFT single-point**. Each stage is also exposed as an [individual subcommand](#cli-subcommands).

An initial reaction path is one command:

```bash
# Multi-PDB mode (R + P endpoints → MEP, with TS optimisation + thermo)
pdb2reaction all -i R.pdb P.pdb -q 0 --tsopt --thermo
```

> **Prerequisites:** input PDBs must already contain hydrogens; multiple PDBs must share the same atoms in the same order (only coordinates differ). Small-molecule `.xyz` / `.gjf` inputs work when `-c` / `--ligand-charge` are omitted.

## Related tools

| Tool | Use case |
|---|---|
| **`pdb2reaction`** (this repo) | Pure-MLIP **cluster-model** reaction paths from PDB / XYZ / GJF — no MM force field required. |
| [**mlmm-toolkit**](https://github.com/t-0hmura/mlmm_toolkit) | **ML/MM ONIOM** with the full protein environment; automates MM parameterisation and ML-region assignment from a single PDB. |
| [**uma_pysis**](https://github.com/t-0hmura/uma_pysis) | YAML-input reaction-mechanism analysis for **small molecules**. |

`pdb2reaction` and `mlmm-toolkit` bundle the same GPU-optimised pysisyphus fork; it is **not** compatible with upstream pysisyphus — do not install them side by side.

## Documentation

- [Getting Started](docs/getting-started.md) · [Installation](docs/installation.md) · [Examples](examples/) · [Troubleshooting](docs/troubleshooting.md)
- [YAML Reference](docs/yaml-reference.md) · [JSON Output Schema](docs/json-output.md)
- Full site: <https://t-0hmura.github.io/pdb2reaction/>

## System requirements

| Component | Requirement |
|---|---|
| OS / Python | Linux x86_64 (validated); macOS / WSL 2 for CPU-only smoke tests. Python >= 3.11 (3.12 tested). |
| GPU / CUDA / VRAM | NVIDIA GPU, CUDA >= 12.6 (12.9 recommended; required for RTX 50-series, matched to the PyTorch wheel). 8 GB VRAM minimum, 16 GB recommended (24 GB for analytical Hessian on 500+-atom regions). |
| RAM / Disk | 32 GB RAM minimum (60 GB recommended); 20 GB free disk for the conda env, UMA cache, and artefacts. |

CPU-only execution works but is 10–100× slower; not recommended for full TS / IRC / Hessian workflows. Full requirement and tuning details: [docs/installation.md](docs/installation.md).

## Installation

```bash
# 1. CUDA-enabled PyTorch (match your CUDA runtime)
pip install torch --index-url https://download.pytorch.org/whl/cu129

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

Examples use GPP C6-methyltransferase BezA ([Tsutsumi et al., *Angew. Chem. Int. Ed.* 2022, 61, e202111217](https://doi.org/10.1002/anie.202111217)) — full commands in [`examples/run.sh`](examples/).

```bash
# Multi-structure MEP (R + P → MEP, with TS + thermochemistry)
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --out-dir result_mep

# Scan mode (single structure → staged bond scan → MEP)
pdb2reaction -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","GPP 321 C7",1.60)]' --tsopt --thermo --out-dir result_scan

# TS-only validation (single TS candidate → tsopt → IRC → freq)
pdb2reaction -i TS_candidate.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt
```

Per-stage walkthrough (`extract` → `opt` → `path-search` → `tsopt` → `freq` → `irc` → `dft`): [docs/getting-started.md](docs/getting-started.md) and [docs/quickstart-all.md](docs/quickstart-all.md).

## Output

A run writes its deliverables to `--out-dir` (default `./result_all/`):

- `segments/seg_NN/{reactant,ts,product}.*` — the canonical R / TS / P structures to cite
- `mep.pdb` / `mep_trj.xyz` — the merged reaction path
- `energy_diagram_MEP.png` — barrier diagram across all segments
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
| `tsopt` | TS optimization (Dimer / RS-I-RFO) | [tsopt](docs/tsopt.md) |
| `path-opt` | MEP via GSM or DMF | [path-opt](docs/path-opt.md) |
| `path-search` | Recursive MEP search with refinement | [path-search](docs/path-search.md) |
| `scan` / `scan2d` / `scan3d` | 1D / 2D / 3D bond-distance scans | [scan](docs/scan.md) · [scan2d](docs/scan2d.md) · [scan3d](docs/scan3d.md) |
| `freq` | Vibrational analysis + thermochemistry | [freq](docs/freq.md) |
| `irc` | IRC (EulerPC) | [irc](docs/irc.md) |
| `dft` | Single-point DFT (GPU4PySCF / PySCF) | [dft](docs/dft.md) |
| `sp` | Single-point MLIP energy / forces / Hessian | [sp](docs/sp.md) |
| `bond-summary` | Compare structures, report bond changes | [bond-summary](docs/bond-summary.md) |
| `trj2fig` / `energy-diagram` | Energy plot / R→TS→P diagram | [trj2fig](docs/trj2fig.md) · [energy-diagram](docs/energy-diagram.md) |

> `tsopt`, `freq`, `irc`: set `--hessian-calc-mode Analytical` when VRAM allows.

## Getting Help

```bash
pdb2reaction --help                       # top-level
pdb2reaction <subcmd> --help              # core options
pdb2reaction <subcmd> --help-advanced     # full option set
p2r --help                                # short alias
python -m pdb2reaction --help             # equivalent module invocation
```

Issues: <https://github.com/t-0hmura/pdb2reaction/issues>.

## Citation

```bibtex
@misc{ohmura2026pdb2reaction,
  author = {Ohmura, Takuto and Sato, Hajime and Terada, Tohru},
  title  = {pdb2reaction: End-to-End Reaction-Path Elucidation from PDB Structures Using Machine-Learning Interatomic Potentials},
  year   = {2026}, doi = {10.26434/chemrxiv.15003538}, note = {ChemRxiv preprint}
}
```

## Agent Skills

Agent instructions for Claude Code / Codex / Cursor live in [`skills/`](skills/) — copy into your project's skill location (e.g. `.claude/skills/`) to let an agent drive `extract` / `path-search` / `tsopt` / `irc` / `dft` end-to-end.

## Known limitations

- **MACE + UMA cannot coexist** (`e3nn` version conflict). Use separate conda envs.
- **DFT single-point** is practical up to ~300 atoms; larger systems need fragmentation.
- **ORB backend** sometimes converges TS with extra soft imaginary modes — mechanism recovery is fine, but for clean single-saddle spectra prefer UMA / MACE or re-score with DFT.
- **CPU-only execution** is 10–100× slower than GPU.

## Contributing

Issues and pull requests are welcome — see [CONTRIBUTING.md](CONTRIBUTING.md).

## License

GNU General Public License v3 (GPL-3.0).
