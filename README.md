# `pdb2reaction`: End-to-End Reaction-Path Elucidation from PDB Structures Using Machine-Learning Interatomic Potentials

[![PyPI](https://img.shields.io/pypi/v/pdb2reaction.svg)](https://pypi.org/project/pdb2reaction/) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/t-0hmura/pdb2reaction/blob/main/examples/pdb2reaction_colab.ipynb)

`pdb2reaction` is a Python CLI that uses machine-learning interatomic potentials (MLIPs) to search for candidate **enzymatic reaction pathways** from **PDB or mmCIF structures**.

## Overview

<img src="https://raw.githubusercontent.com/t-0hmura/pdb2reaction/main/docs/overview.png" alt="pdb2reaction workflow overview" width="90%">

Test a reaction mechanism in a single command:

```bash
# Multi-PDB mode (R + P endpoints → MEP, with TS optimization + thermo)
pdb2reaction all -i R.pdb P.pdb -c 'LIG' -l 'LIG:-1' --tsopt --thermo
```

The run writes the R / TS / P structures, energy diagrams, and `summary.log` / `summary.json` to `result_all/`.

### What it is for

* **Trial and error on reaction mechanisms**: screen large systems for which DFT alone would take too long
* **Starting structures for quantum chemistry**: build cluster models of the reactant (R), transition state (TS), and product (P)
* **High-throughput calculations over many systems**: explore reaction pathways systematically across substrate variants and enzyme mutants

Once MLIP finds a plausible path, pdb2reaction can take that TS straight into DFT TS optimization: the TS optimization → IRC → endpoint optimization → frequency workflow runs as GPU-accelerated DFT with GPU4PySCF. See [DFT backend](docs/dft-backend.md).

## Installation

Requirements: Linux, Python 3.11 or later (3.12 recommended; ORB requires 3.11 or 3.12), and an NVIDIA GPU with the official PyTorch 2.13 CUDA wheel that matches your driver/GPU (`cu130` recommended). Details: [docs/installation.md](docs/installation.md).

```bash
# 1. CUDA-enabled PyTorch (choose the official 2.13 wheel for your driver/GPU)
pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130

# 2. Install pdb2reaction
pip install pdb2reaction

# 3. Authenticate Hugging Face once (only required for the default UMA backend)
#    Accept the FAIR Chemistry License v1 at https://huggingface.co/facebook/UMA, then:
hf auth login                               # interactive
# OR, for non-interactive CI/HPC jobs: export HF_TOKEN=hf_xxx
```

**Optional extras** (install only what you need):

| Extra | Adds |
|---|---|
| `[orb]` / `[aimnet]` | Orb / AIMNet2 MLIP backend (`-b orb` / `-b aimnet2`) — *not* HF-gated |
| `[dft]` / `[dft-cuda12]` | DFT calculator and standalone command with native CUDA 13 / CUDA 12 GPU4PySCF |
| `[mcp]` | Model Context Protocol server for agent clients |

The MACE backend (`-b mace`) does not install into the same environment as UMA; create a dedicated environment as described in [docs/installation.md](docs/installation.md).

CUDA module loads, alternative-backend recipes, DMF/`cyipopt` setup, Plotly Chromium, and HPC job-script templates: [docs/installation.md](docs/installation.md) and [docs/hpc-example.md](docs/hpc-example.md).

## Quick Examples

> **Before you start:**
>
> - PDB/mmCIF inputs must already contain hydrogens.
> - Reaction-ordered structures must share the same atom identities and order (only coordinates differ).
> - mmCIF files (including multi-character chains and large residue IDs) and large structures: see [mmCIF and large structures](docs/cli-conventions.md#mmcif-and-large-structures).
> - Small-molecule `.xyz` / `.gjf` inputs work when `--center/-c` and `--ligand-charge/-l` are omitted.

Examples use GPP C6-methyltransferase BezA ([Tsutsumi et al., *Angew. Chem. Int. Ed.* 2022, 61, e202111217](https://doi.org/10.1002/anie.202111217)). Run the commands below from the repository root (`git clone https://github.com/t-0hmura/pdb2reaction && cd pdb2reaction`); the complete MEP and scan examples are in [`examples/run.sh`](examples/run.sh).

```bash
# Multi-structure MEP (R + P → MEP, with TS + thermochemistry)
pdb2reaction -i examples/1.R.pdb examples/3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --out-dir result_mep

# Scan mode (single structure → staged bond scan → MEP)
pdb2reaction -i examples/1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","GPP 321 C7",1.60)]' --tsopt --thermo --out-dir result_scan

# TS-only validation (your own TS candidate → tsopt → IRC → freq)
pdb2reaction -i TS_candidate.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo --out-dir result_tsonly
```

`pdb2reaction` can also be used to investigate reaction mechanisms of **small molecules** and **user-defined cluster models**.

```bash
# Small molecule (gas-phase): .xyz / .gjf input — omit -c, set charge with -q
pdb2reaction -i reactant.xyz product.xyz -q 0 --tsopt --thermo --out-dir result_small

# Your own cluster model (already-trimmed PDB): omit -c to use it as-is
pdb2reaction -i cluster_R.pdb cluster_P.pdb -q 0 --tsopt --thermo --out-dir result_cluster
```

For a hand-built cluster, check its boundaries with [the cluster-boundary checklist](docs/model-setup.md#building-or-auditing-a-cluster-model-manually).

Each stage (`extract` → `opt` → `path-opt` → `tsopt` → `irc` → `freq` → `dft`) also runs as its own subcommand; see [CLI Subcommands](#cli-subcommands) for the per-stage pages.

## Output

A non-dry `all` run writes the deliverables reached by its enabled stages to
`--out-dir` (default `./result_all/`):

- `segments/seg_NN/{reactant,ts,product}.*` — the canonical R / TS / P structures to cite
- `mep_trj.xyz` (plus `mep_trj.pdb` when topology is available and `mep_trj.cif` for mmCIF input or very large PDB) — the merged reaction path in MEP/scan-list modes
- `energy_diagram_MEP.png` — MEP energy diagram
- `summary.log` / `summary.json`

Pipeline scratch lives under `_work/`; keep it if you may redo the post-processing with `--resume-segment`. Full layout and filename conventions: [docs/output-layout.md](docs/output-layout.md).

## Colab GUI workspace

**An interactive GUI workspace is available in Google Colab.** It brings ordered structure input, Mol* visualization and atom picking, controls generated from the live CLI, execution, and linked MEP/IRC/result inspection into one notebook. Choose a GPU runtime and [open the Colab GUI workspace](https://colab.research.google.com/github/t-0hmura/pdb2reaction/blob/main/examples/pdb2reaction_colab.ipynb).

<img src="https://raw.githubusercontent.com/t-0hmura/pdb2reaction/main/docs/colab_workspace.png" alt="pdb2reaction Colab GUI workspace showing Mol* structure setup and active-site selection controls" width="90%">

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
| `sp` | Single-point calculator energy / forces / Hessian (MLIP or `-b dft`) | [sp](docs/sp.md) |
| `bond-summary` | Compare structures, report bond changes | [bond-summary](docs/bond-summary.md) |
| `trj2fig` / `energy-diagram` | Energy plot / R→TS→P diagram | [trj2fig](docs/trj2fig.md) · [energy-diagram](docs/energy-diagram.md) |

## Documentation

- [Getting Started](docs/getting-started.md) · [Installation](docs/installation.md) · [Quickstart: all](docs/quickstart-all.md) · [Building the cluster model](docs/model-setup.md) · [DFT backend](docs/dft-backend.md) · [Troubleshooting](docs/troubleshooting.md)
- Full site: <https://t-0hmura.github.io/pdb2reaction/>

## Agent Skills

`skills/` holds Agent Skills that let an AI coding agent run `pdb2reaction` workflows and subcommands. Copy the skill folders into `.claude/skills/` or `~/.claude/skills/` for Claude Code, or into `.agents/skills/` or `~/.agents/skills/` for Codex. The list and an example copy command are in [`skills/README.md`](skills/README.md).

## Getting Help

```bash
pdb2reaction --help                       # top-level
pdb2reaction <subcmd> --help              # core options
pdb2reaction <subcmd> --help-advanced     # full option set
```

Issues: <https://github.com/t-0hmura/pdb2reaction/issues>.

## Related tools

| Tool | Use case |
|---|---|
| [**mlmm-toolkit**](https://github.com/t-0hmura/mlmm_toolkit) | **ML/MM ONIOM** with the full protein environment; automates MM parameterization and ML-region assignment from a single PDB. |
| [**uma_pysis**](https://github.com/t-0hmura/uma_pysis) | Lightweight **YAML-driven UMA–pysisyphus interface** for quick/exploratory reaction-mechanism studies (GS / TS / IRC / ΔG). |

## Known limitations

- **MACE + UMA cannot coexist** (`e3nn` version conflict). Use separate conda envs.
- **DFT single-point cost** depends strongly on basis, functional, grid, elements, and hardware; pilot one representative structure before batching.
- **Check every TS.** When tsopt converges or stops on an energy plateau, it computes the Hessian of the final geometry and reports n_imag; a successful TS optimization shows one imaginary mode along the reaction. IRC and endpoint optimization then confirm that it connects R and P (reaching the right R and P counts even if IRC stops early). ORB runs in fp64 by default; validate frequencies and IRC yourself if you switch to fp32/TF32.
- **CPU-only execution** is supported but usually much slower than GPU.
- `pdb2reaction` bundles a GPU-optimized pysisyphus fork that is **not** compatible with upstream pysisyphus — do not install it into an environment that already has upstream pysisyphus.

## Citation

```bibtex
@misc{ohmura2026pdb2reaction,
  author = {Ohmura, Takuto and Sato, Hajime and Terada, Tohru},
  title  = {pdb2reaction: End-to-End Reaction-Path Elucidation from PDB Structures Using Machine-Learning Interatomic Potentials},
  year   = {2026}, doi = {10.26434/chemrxiv.15003538/v1}, note = {ChemRxiv preprint}
}
```

## Contributing

Issues and pull requests are welcome — see [CONTRIBUTING.md](CONTRIBUTING.md).

## License

GNU General Public License v3 (GPL-3.0).
