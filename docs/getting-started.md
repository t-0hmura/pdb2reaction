# Getting Started

## Overview

`pdb2reaction` is a Python CLI toolkit for **modeling enzymatic reaction pathways from PDB structures** using machine-learning interatomic potentials (MLIPs).

In many workflows, a **single command** like the one below is enough to generate a useful initial reaction path:
```bash
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3'
```

---
You can also run **Minimum Energy Path (MEP) search → Transition State (TS) optimization → Intrinsic Reaction Coordinate (IRC) → thermochemistry → single-point DFT** in a single run by adding `--tsopt --thermo --dft`:
```bash
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo --dft
```
---

Given **(i) two or more full protein–ligand PDB files** (R → … → P), **or (ii) one PDB with `--scan-lists`**, **or (iii) one TS candidate with `--tsopt`**, `pdb2reaction` automatically:

- extracts an **active site model (binding pocket)** around user-defined substrates to build a **cluster model**,
- explores **minimum-energy paths (MEPs)** with path optimization methods such as the Growing String Method (GSM) and Direct Max Flux (DMF),
- _optionally_ optimizes **transition states**, runs **vibrational analysis**, **IRC calculations**, and **single-point DFT calculations**.

Calculations use machine-learning interatomic potentials (MLIPs). The default backend is Meta's **UMA**, but **ORB**, **MACE**, and **AIMNet2** are also supported via `-b/--backend`. Typical use cases include:

- **Trial-and-error exploration of reaction mechanisms** at a scale where DFT-level verification would be prohibitively slow
- **Generating initial geometries** (reactant/TS/product cluster models) for subsequent quantum-chemistry refinement
- **High-throughput screening** of reaction pathways across substrate variants or enzyme mutants

The CLI is designed to generate **multi-step enzymatic reaction mechanisms** with minimal manual intervention. The same workflow also works for small-molecule systems. When you skip active site model extraction (omit `--center/-c` and `--ligand-charge`), you can also use `.xyz` or `.gjf` inputs.

On **HPC clusters or multi-GPU workstations**, `pdb2reaction` can scale to large cluster models (and optionally **full protein–ligand complexes**) by parallelizing UMA inference across nodes. Set `workers` and `workers_per_node` to enable multi-worker inference; see [MLIP Calculator](uma-pysis.md) for configuration details. Alternative backends (ORB, MACE, AIMNet2) can be selected with `-b/--backend`.

```{important}
- Input PDB files must already contain **hydrogen atoms**.
- When you provide multiple PDBs, they must contain **the same atoms in the same order** (only coordinates may differ); otherwise an error is raised.
```

```{tip}
For symptom-first diagnosis, start with [Common Error Recipes](recipes-common-errors.md).
If you encounter an error during setup or runtime, refer to [Troubleshooting](troubleshooting.md).
```

### CLI conventions

| Convention | Example | Notes |
|------------|---------|-------|
| **Residue selectors** | `'SAM,GPP'` or `'A:123,B:456'` | Quote multi-value strings to prevent shell expansion |
| **Charge mapping** | `-l 'SAM:1,GPP:-3'` | Colon separates name and charge; comma separates entries |
| **Atom selectors** | `'TYR,285,CA'` or `'TYR 285 CA'` | Delimiters: space, comma, slash, backtick, backslash |

For full details, see [CLI Conventions](cli-conventions.md).

### Recommended tools for hydrogen addition

If your PDB lacks hydrogen atoms, use one of the following tools before running pdb2reaction:

| Tool | Example Command | Notes |
|------|-----------------|-------|
| **reduce** (Richardson Lab) | `reduce input.pdb > output.pdb` | Fast, widely used for crystallographic structures |
| **pdb2pqr** | `pdb2pqr --ff=AMBER input.pdb output.pqr` | Adds hydrogens and assigns partial charges |
| **Open Babel** | `obabel input.pdb -O output.pdb -h` | General-purpose cheminformatics toolkit |
| **PyMOL** | `h_add` (in PyMOL) | Molecular visualization tool with hydrogen addition |
| **tleap** (AmberTools) | `tleap -f leapin` | Amber force-field preparation tool |

To ensure identical atom ordering across multiple PDB inputs, apply the same hydrogen-addition tool with consistent settings to all structures.

---

## Quickstart routes (recommended)

For setup and dependency installation, see [Installation](installation.md).

- [Quickstart: run `pdb2reaction all`](quickstart-all.md)
- [Quickstart: run a single-structure staged scan with `pdb2reaction scan`](quickstart-scan.md)
- [Quickstart: TS optimization and validation with `pdb2reaction tsopt`](quickstart-tsopt-freq.md)

---

## Command line basics

The main entry point is the `pdb2reaction` command, installed via `pip`. A shorthand alias **`p2r`** is also available — all commands can be run with either name. Internally it uses the **Click** library, and the default subcommand is `all`.

That means:

```bash
pdb2reaction [OPTIONS]...
# is equivalent to
pdb2reaction all [OPTIONS]...
```

The [`all`](all.md) command runs the full pipeline—cluster extraction, MEP search, TS optimization, vibrational analysis, and optional DFT—in a single invocation.

All high-level workflows share two important options when you use cluster extraction:

- `-i/--input`: one or more **full structures** (reactant, intermediate(s), product).
- `-c/--center`: how to define the **substrate / extraction center** (e.g., residue names or residue IDs).

If you omit `--center/-c`, cluster extraction is skipped and the **full input structure** is used directly.

---

## Main workflow modes

### Multi-structure MEP workflow (reactant → product)

Use this when you already have several full PDB structures along a putative reaction coordinate (e.g., R → I1 → I2 → P).

**Minimal example**

```bash
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3'
```

**Richer example**

```bash
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --out-dir ./result_all --tsopt --thermo --dft
```

Behavior:

- takes two or more **full systems** in reaction order,
- extracts cluster models for each structure,
- performs a **recursive MEP search** via `path-search` by default (outputs under `path_search/`),
- optionally switches to a **single-pass** [`path-opt`](path-opt.md) run with `--no-refine-path`,
- when PDB templates are available, merges the cluster-model MEP back into the **full system**,
- optionally runs TS optimization, vibrational analysis, and single-point DFT calculations for each segment.

This is the recommended mode when you can generate reasonably spaced intermediates (e.g., from docking, MD, or manual modeling).

```{important}
`pdb2reaction` assumes that multiple input PDBs contain **exactly the same atoms in the same order** (only coordinates may differ). If any non-coordinate fields differ across inputs, an error is raised. Input PDB files must also contain **hydrogen atoms**.
```

---

### Single-structure + staged scan (feeds MEP refinement)

Use this when you only have **one PDB structure**, but you know which inter-atomic distances should change along the reaction.

Provide a single `-i` together with `--scan-lists`:

**Minimal example**

```bash
pdb2reaction -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 -s '[("CS1 SAM 320","GPP 321 C7",1.60)]' \
 '[("GPP 321 H11","GLU 186 OE2",0.90)]'
```

**Richer example**

```bash
pdb2reaction -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 -s '[("CS1 SAM 320","GPP 321 C7",1.60)]' \
 '[("GPP 321 H11","GLU 186 OE2",0.90)]' \
 --multiplicity 1 --out-dir ./result_scan --tsopt --thermo --dft
```

Key points:

- `--scan-lists` describes **staged distance scans** on the extracted cluster model.
- Each tuple `(i, j, target_Å)` is:
 - a PDB atom selector string like `'TYR,285,CA'` (**delimiters can be: space/comma/slash/backtick/backslash ` ` `,` `/` `` ` `` `\`**) **or** a 1-based atom index,
 - automatically remapped to the cluster-model indices.
- Supplying one `-s/--scan-lists` literal runs a single scan stage; multiple literals run sequential stages (e.g. `-s '[(…)]' '[(…)]'`).
- Each stage writes a `stage_XX/result.pdb`, which is treated as a candidate intermediate or product.
- The default `all` workflow refines the concatenated stages with recursive `path-search`.
- With `--no-refine-path`, it instead performs a single-pass `path-opt` chain and skips the recursive refiner (no merged `mep_w_ref*.pdb`).

This mode is useful for building reaction paths starting from a single structure.

---

### Single-structure TSOPT-only mode

Use this when you already have a **transition-state candidate** and only want to optimize it and proceed to IRC calculations.

Provide exactly one PDB and enable `--tsopt`:

**Minimal example**

```bash
pdb2reaction -i TS_CANDIDATE.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt
```

**Richer example**

```bash
pdb2reaction -i TS_CANDIDATE.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft --out-dir ./result_tsopt_only
```

Behavior:

- skips the MEP/path search entirely,
- performs **transition-state optimization** on the cluster model,
- runs an **IRC** in both directions and optimizes the endpoints to obtain R and P minima,
- can then run vibrational analysis ([`freq`](freq.md)) and single-point DFT (`dft`) on the R/TS/P structures,
- produces MLIP, Gibbs, and DFT//MLIP (DFT single-point energies at MLIP-optimized geometries) energy diagrams.

Outputs such as `energy_diagram_*_all.png` and `irc_plot_all.png` are mirrored under the top-level `--out-dir`.

```{important}
Single-input runs require **either** `--scan-lists` (staged scan → GSM) **or** `--tsopt` (TSOPT-only). Supplying only a single `-i` without one of these will not trigger a full workflow.
```

---

## Important CLI options and behaviors

Below are the most commonly used options across workflows.

| Option | Description |
|--------|-------------|
| `-i, --input PATH...` | Input structures. **≥ 2 PDBs** → MEP search; **1 PDB + `--scan-lists`** → staged scan → GSM; **1 PDB + `--tsopt`** → TSOPT-only mode. |
| `-c, --center TEXT` | Defines the substrate / extraction center. Supports residue names (`'SAM,GPP'`), residue IDs (`A:123,B:456`), or PDB paths. |
| `-l, --ligand-charge TEXT` | Charge info: mapping (`'SAM:1,GPP:-3'`) or single integer. |
| `-q, --charge INT` | Hard override of net system charge. |
| `-m, --multiplicity INT` | Spin multiplicity (e.g., `1` for singlet). |
| `-s, --scan-lists TEXT...` | Staged distance scans for single-input runs. |
| `-o, --out-dir PATH` | Top-level output directory. |
| `--tsopt/--no-tsopt` | Enable TS optimization and IRC. |
| `--thermo/--no-thermo` | Run vibrational analysis and thermochemistry. |
| `--dft/--no-dft` | Perform single-point DFT calculations. |
| `--refine-path/--no-refine-path` | Recursive MEP refinement (default: enabled) vs single-pass. |
| `--opt-mode grad\|hess` | Workflow-level preset in `all` (`grad` -> LBFGS/Dimer, `hess` -> RFO/RS-I-RFO; default `grad`). For direct commands, prefer `opt --opt-mode grad|hess` and `tsopt --opt-mode grad|hess`. |
| `--mep-mode gsm\|dmf` | MEP method: Growing String Method or Direct Max Flux. |
| `--hessian-calc-mode Analytical\|FiniteDifference` | Hessian evaluation method. For Hessian evaluation modes, see [MLIP Calculator](uma-pysis.md#hessian-evaluation). |

For a full matrix of options and YAML schemas, see [all](all.md) and [YAML Reference](yaml-reference.md).

---

## Run summaries

Every `pdb2reaction all` run writes:

- `summary.log` – formatted summary for quick inspection, and
- `summary.yaml` – YAML-format summary.

They typically contain:

- the exact CLI command invoked,
- global MEP statistics (e.g. maximum barrier, path length),
- per-segment barrier heights and key bond changes,
- energies from the MLIP backend, thermochemistry, and DFT post-processing (where enabled).

Each segment directory under `path_search/` also gets its own `summary.log` and `summary.yaml`, so you can inspect local refinements independently.

---

## CLI commands

Most users will primarily call `pdb2reaction all`. The CLI also exposes individual subcommands—each supports `-h/--help`.
`pdb2reaction all --help` shows core options and `pdb2reaction all --help-advanced` shows the complete list.
`scan`, `scan2d`, `scan3d`, and the calculation commands (`opt`, `path-opt`, `path-search`, `tsopt`, `freq`, [`irc`](irc.md), `dft`) now follow the same progressive-help pattern (`--help` core, `--help-advanced` full). `add-elem-info`, `trj2fig`, and `energy-diagram` also use the same pattern. `extract` and `fix-altloc` also support progressive help (`--help` core, `--help-advanced` full parser options).

| Subcommand | Role | Documentation |
|------------|------|---------------|
| `all` | End-to-end workflow | [all](all.md) |
| `extract` | Extract active site model (cluster model) | [extract](extract.md) |
| `opt` | Geometry optimization | [opt](opt.md) |
| `tsopt` | Transition state optimization | [tsopt](tsopt.md) |
| `path-opt` | MEP optimization (GSM/DMF) | [path-opt](path-opt.md) |
| `path-search` | Recursive MEP search | [path-search](path-search.md) |
| `scan` | 1D bond-length scan | [scan](scan.md) |
| `scan2d` | 2D distance scan | [scan2d](scan2d.md) |
| `scan3d` | 3D distance scan | [scan3d](scan3d.md) |
| `irc` | IRC calculation | [irc](irc.md) |
| `freq` | Vibrational analysis | [freq](freq.md) |
| `dft` | Single-point DFT | [dft](dft.md) |
| `trj2fig` | Plot energy profiles | [trj2fig](trj2fig.md) |
| `energy-diagram` | Draw state energy diagram from numeric values | [energy-diagram](energy-diagram.md) |
| `fix-altloc` | Resolve PDB alternate-location indicators | [fix-altloc](fix-altloc.md) |
| `add-elem-info` | Repair PDB element columns | [add-elem-info](add-elem-info.md) |

```{tip}
For Hessian evaluation modes, see [MLIP Calculator](uma-pysis.md#hessian-evaluation).
```

---

## Quick reference

**Common command patterns:**

```bash
# Basic MEP search (2+ structures)
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3'

# Full workflow with post-processing
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft

# Single structure with staged scan
pdb2reaction -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 -s '[("CS1 SAM 320","GPP 321 C7",1.60)]' '[("GPP 321 H11","GLU 186 OE2",0.90)]'

# TS-only optimization
pdb2reaction -i TS.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo
```

**Essential options:**

| Option | Purpose |
|--------|---------|
| `-i` | Input structure(s) |
| `-c` | Substrate definition for active site model extraction |
| `-l, --ligand-charge` | Substrate charges (e.g., `'SAM:1,GPP:-3'`) |
| `--tsopt` | Enable TS optimization + IRC |
| `--thermo` | Run vibrational analysis |
| `--dft` | Run single-point DFT |
| `-o, --out-dir` | Output directory |

---

## Getting help

For any subcommand:

```bash
pdb2reaction <subcommand> --help
pdb2reaction <subcommand> --help-advanced
pdb2reaction all --help-advanced
```

For `all`, `--help` is intentionally short. Use `--help-advanced` to see every option.
For detailed MLIP backend options, see [MLIP Calculator](uma-pysis.md).

If you encounter any issues, please open an Issue on the [GitHub repository](https://github.com/t-0hmura/pdb2reaction).
