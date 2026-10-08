# Getting Started

## Overview

<img src="./overview.png" alt="pdb2reaction workflow overview" width="90%">

`pdb2reaction` is a Python command-line toolkit that uses machine-learning interatomic potentials (MLIPs) to **search automatically for candidate enzyme reaction pathways, starting from PDB / mmCIF structures**.

The MLIPs are neural networks trained on DFT data. They approximate a DFT-level potential energy surface at a tiny fraction of the cost, which makes pathway searches fast.

In many cases, a **single command** like this one gives a first draft of the reaction pathway:

```bash
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3'
```

`Scientific status: success` near the end of the terminal output means that every requested stage converged.

---

Add `--tsopt --thermo --dft` and the same run continues automatically through **minimum energy path (MEP) search → transition state (TS) optimization → intrinsic reaction coordinate (IRC) → vibrational analysis and thermochemical correction → DFT single points**.

```bash
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo --dft
```

---

> **Examples:** the [`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples) directory holds the structures used above (`1.R.pdb`, `3.P.pdb`) and a set of workflow scripts (MEP search and scan pipelines) built around the GPP C6-methyltransferase BezA ([Tsutsumi et al., *Angew. Chem. Int. Ed.* 2022, 61, e202111217](https://doi.org/10.1002/anie.202111217)). After [installation](installation.md), get it with `git clone https://github.com/t-0hmura/pdb2reaction && cd pdb2reaction/examples` and run the commands above there.

### What it is for

* **Trial and error on reaction mechanisms** at a scale that DFT or other quantum chemistry takes too long to check
* **Starting structures** for quantum chemistry (cluster models of the reactant, TS, and product)
* **Many reaction-path calculations** across substrate variants and enzyme mutants

### What it automates

Provide one of three inputs: (1) several PDB structures in reaction order (R → … → P), (2) one structure plus a distance scan, or (3) one structure plus TS optimization. `pdb2reaction` then handles the following automatically.

1. **Cluster model construction**: cuts out the active site (binding pocket) around the specified substrates
2. **Minimum energy path (MEP) search**: searches the pathway with the Growing String Method (GSM) or Direct Max Flux (DMF)
3. **High-accuracy checks**: TS optimization, IRC, vibrational analysis, and DFT single points

Once MLIP has found a reasonable pathway, `pdb2reaction` can take its TS straight into a DFT TS optimization. It runs the TS optimization → IRC → endpoint optimization → frequency workflow with GPU-accelerated DFT through GPU4PySCF. See [Refine an MLIP TS with DFT](dft-backend.md) for details.

---

## Workflow and pipeline

### The pipeline

The `all` subcommand (the default) runs the whole workflow in one go, stage by stage in this order:

```text
Input structure(s) (PDB / mmCIF)
  │
  ▼
[extract] extraction: cut out the active-site model (only with -c)
  │
  ▼
[scan] scan: staged distance-restrained scan (only with -s)
  │
  ▼
[path-opt / path-search] path search: find the MEP (minimum energy path); skipped in TS-only mode
  │
  ▼
[tsopt] TS optimization: refine the transition state (only with --tsopt)
  │
  ▼
[irc] IRC: follow the intrinsic reaction coordinate and optimize its endpoints (only with --tsopt)
  │
  ▼
[freq] vibrational analysis: compute the thermochemical correction (only with --tsopt --thermo)
  │
  ▼
[dft] DFT single points: compute DFT energies (only with --tsopt --dft)
```

Each stage also runs on its own as a subcommand ([`extract`](extract.md), [`tsopt`](tsopt.md), [`irc`](irc.md), and so on; see the [subcommand list](index.md#subcommands)).

---

## Where to start

For environment setup, see the [Installation guide](installation.md).

* **Try it in a web browser**: the [Colab GUI notebook](https://colab.research.google.com/github/t-0hmura/pdb2reaction/blob/main/examples/pdb2reaction_colab.ipynb) (pick residues in a 3D preview)
* **Start from several PDB structures**: [Quickstart: `pdb2reaction all`](quickstart-all.md)
* **Explore from one PDB structure with a scan**: [Quickstart: `pdb2reaction all --scan-lists`](quickstart-scan.md)
* **Optimize and check a TS candidate**: [Quickstart: TS-only mode](quickstart-tsopt.md)

---

## How the command works

Installation provides `pdb2reaction` and the short alias `p2r`. Without a subcommand, `all` runs.

```bash
# These two do the same thing
pdb2reaction [OPTIONS]...
pdb2reaction all [OPTIONS]...
```

### Choosing an input mode

| Mode | Input | What happens |
| --- | --- | --- |
| **Multi-structure MEP search** | Two or more PDBs (`-i R.pdb P.pdb`) | Extracts cluster models from the structures (with `-c`) and searches the MEP |
| **Single structure + scan** | One PDB + `--scan-lists` (`-s`) | Changes the chosen bond distances step by step to build the pathway |
| **TS-only mode** | One PDB + `--tsopt` | Skips the MEP search and goes straight to optimizing the TS candidate and running IRC |

> **Note:** a single structure without `--scan-lists/-s` or `--tsopt` stops with an error.

---

## Main CLI options

| Option | Example | Description |
| --- | --- | --- |
| `-i, --input` | `1.R.pdb 3.P.pdb` | Input structure files (PDB / mmCIF); accepts several |
| `-c, --center` | `'SAM,GPP'` / `'A:SAM:123'` | Extraction center (substrate residue names, residue IDs, or a PDB file of the substrate with the same coordinates as the input). Without it, no extraction runs and the whole structure is used |
| `-l, --ligand-charge` | `'SAM:1,GPP:-3'` | Formal charge of each ligand, as a mapping (standard residues and ions are counted automatically) |
| `-q, --charge` | `-2` | Total charge of the whole extracted model (set it to override the automatic value) |
| `-m, --multiplicity` | `1` | Spin multiplicity (default `1`, a singlet) |
| `--tsopt` | (flag) | Turns on TS optimization and IRC |
| `--thermo` | (flag) | Runs vibrational analysis and thermochemical correction with the QRRHO (quasi-rigid-rotor harmonic oscillator) model (with `--tsopt`) |
| `--dft` | (flag) | Runs DFT single points on the resulting structures (with `--tsopt`) |
| `-b, --backend` | `uma` / `orb` / `mace` | MLIP backend to use (default `uma`) |

`--dft` needs the DFT extra from step 7 of {ref}`Step-by-step installation <step-by-step-installation>`.

For the syntax rules, see [Common options and selectors](cli-conventions.md); for every option, see the [`all` CLI reference](reference/commands/all.md).

---

## Before you run: the input structures

### 1. Add hydrogens (required)

Input structures must contain **every hydrogen atom**. When a structure lacks hydrogens (a crystal structure, for example), add them beforehand with a tool such as these:

| Recommended tool | Example command | Notes |
| --- | --- | --- |
| **reduce** (Richardson Lab) | `reduce input.pdb > output.pdb` | Fast; widely used to add hydrogens to crystal structures |
| **Open Babel** | `obabel input.pdb -O output.pdb -h` | General-purpose cheminformatics toolkit |
| **PyMOL** | `h_add` at the command line | Add hydrogens while checking the structure visually |
| **tleap** (AmberTools) | `tleap -f leapin` | Careful placement based on the Amber force field |

`all` fills blank element columns (77–78) by itself; before a standalone command such as `extract`, fill them with [`add-elem-info`](add-elem-info.md). When reading a PDB, every command keeps the {ref}`alternate location (altLoc) with the highest mean occupancy <mmcif-input>` of each residue; [`fix-altloc`](fix-altloc.md) saves that choice to a file.

### 2. Keep the same atom order (multiple structures)

When the input has several structures, such as a reactant (R) and a product (P), **every structure must list the same atoms in the same order** (only the coordinates differ). Run the hydrogen tool on every structure with the same settings. An atom that moves to another residue keeps its residue and atom name from R: in the bundled example, the hydrogen that GPP passes to Glu186 is still `H11` of `GPP 321` in `3.P.pdb`.

mmCIF (`.cif`, `.mmcif`) and large structures (≥10,000 residues or ≥99,999 atoms) work with the same commands as PDB. See {ref}`mmCIF input <mmcif-input>` for details.

---

## Output files

When the run finishes, the output directory (`-o`, default `./result_all/`) contains the following files. [Output Directory Layout](output-layout.md) lists the main files, and {ref}`JSON Output Reference <summary-json-path-search-all>` every key of `summary.json`.

| File / folder | Contents |
| --- | --- |
| `summary.log` | Text summary (directory layout and progress of each stage) |
| `summary.json` | Machine-readable results (barriers, energies of each state, bond changes) |
| `energy_diagram_*.png` | Energy profile plots (electronic energy / Gibbs-corrected) |
| `mep_trj.pdb` / `mep_trj.cif` | Animated trajectory of the minimum energy path (MEP) |
| `segments/seg_NN/` | Detailed results for each reaction segment (optimized R/TS/P structures, IRC trajectories, and more; with `--tsopt`) |

At the end of the terminal output, the `Scientific status:` line under `====== Pipeline summary ======` (`scientific_status` in `summary.json`) is `success` when every requested stage converged and, with `--tsopt`, the TS has n_imag = 1. Otherwise it is `partial` or `failed`, with the reasons in `scientific_status_reasons`. Check yourself that the imaginary mode moves the intended bonds and that the IRC endpoints are the intended R and P (`segments[].bond_changes`); each quickstart lists the files to open.

---

## AI agent skills

`pdb2reaction` ships instructions for AI agents (Claude Code, Codex, Cursor, and others) in the `skills/` directory.

They define the CLI subcommands, the input/output rules for PDB/mmCIF/XYZ/GJF, backend installation steps, and best practices for HPC parallel runs. Load `skills/` into an agent, and it can run and analyze calculations from plain-language instructions. For where to place the files and the full list of skills, see [`skills/README.md`](https://github.com/t-0hmura/pdb2reaction/blob/main/skills/README.md). To call the commands as tools from an MCP client, see [MCP server](mcp_server.md).

---

## Troubleshooting and support

If an error occurs during a run, see these pages:

* {ref}`Troubleshooting <troubleshooting-quick-table>`: fixes by error symptom, and solutions for installation and environment problems
* [MLIP Backends](backends.md): details on GPU memory and parallel runs; [HPC example](hpc-example.md) for a job script over several GPU nodes

To see every option of a command, use the help options:

```bash
pdb2reaction <subcommand> --help
pdb2reaction all --help-advanced
```

Report unresolved problems and bugs on [GitHub Issues](https://github.com/t-0hmura/pdb2reaction/issues).
