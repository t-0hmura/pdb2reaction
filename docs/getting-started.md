# Getting Started

## Overview

<img src="./overview.png" alt="pdb2reaction workflow overview" width="90%">

`pdb2reaction` is a Python CLI for **elucidating enzymatic reaction pathways from PDB structures** using machine-learning interatomic potentials (MLIPs). The default backend is Meta's UMA; `orb`, `mace`, and `aimnet2` are also supported via `-b/--backend`. Foundation-model MLIPs make cluster-model TS optimisation, IRC verification, and QRRHO thermochemistry tractable on a single GPU — lowering the DFT-bound cost barrier that previously throttled mechanistic screening.

A single command generates a useful initial reaction path:

```bash
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3'                       # MEP only
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo --dft   # full
```

Given (i) ≥ 2 PDBs (R → ... → P), (ii) one PDB with `--scan-lists/-s`, or (iii) one TS candidate with `--tsopt`, `pdb2reaction` extracts an **active-site cluster model**, runs an **MEP search** (GSM / DMF), and optionally chains TS optimisation, IRC, frequencies, and single-point DFT. The same workflow also works for small-molecule systems — omit `--center/-c` and `--ligand-charge/-l` to use `.xyz` / `.gjf` inputs.

Working examples (BezA C6-methyltransferase, both multi-structure MEP and scan modes): [`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples). For setup see [Installation](installation.md); for symptom-first diagnosis see [Common Error Recipes](recipes-common-errors.md) and [Troubleshooting](troubleshooting.md).

### Pipeline (the `all` subcommand)

```text
PDB(s) → [extract] → [scan] (optional, --scan-lists) → [path-search] (MEP) → [tsopt] → [irc] → [freq] → [dft] (optional)
```

Each stage is also a standalone subcommand; `all` orchestrates them and writes unified `summary.json` + `summary.log`.

### Key output files

| File | Description |
|---|---|
| `summary.json` | Machine-readable results (barriers, energies, bond changes, environment) |
| `summary.log` | Human-readable text summary with directory tree |
| `seg_XX/` | IRC-optimised R/TS/P structures per reaction step |
| `mep.pdb` | Merged MEP trajectory (PyMOL / VMD) |
| `energy_diagram_*.png` | Energy profile plots (electronic / Gibbs-corrected) |

```{important}
- Input PDBs must already contain **hydrogen atoms**.
- When you provide multiple PDBs, they must contain the same atoms in the same order (only coordinates may differ).
```

### CLI conventions

| Convention | Example | Notes |
|---|---|---|
| Residue selectors | `'SAM,GPP'` or `'A:123,B:456'` | Quote multi-value strings. |
| Charge mapping | `-l 'SAM:1,GPP:-3'` | Colon separates name and charge; comma separates entries. |
| Atom selectors | `'TYR,285,CA'` or `'TYR 285 CA'` | Delimiters: space / comma / slash / backtick / backslash. |

Full table: [CLI Conventions](cli-conventions.md).

### Hydrogen addition (if your PDB lacks H)

`reduce input.pdb > out.pdb` (fast, crystallographic structures) · `pdb2pqr --ff=AMBER input.pdb out.pqr` (also assigns partial charges) · `obabel input.pdb -O out.pdb -h` (general cheminformatics) · PyMOL `h_add` · AmberTools `tleap` (Amber force-field prep). Apply the same tool with consistent settings to **every** input to keep atom order matched across structures.

## Quickstart routes

- [Quickstart: `pdb2reaction all`](quickstart-all.md) — multi-structure MEP
- [Quickstart: single-structure staged scan](quickstart-scan.md) — one PDB + `--scan-lists`
- [Quickstart: TS-only mode](quickstart-tsopt-freq.md) — `pdb2reaction all --tsopt`

## Command line basics

The CLI entry point is `pdb2reaction` (alias `p2r`; both register from the same setuptools entry point). The default subcommand is `all`:

```bash
pdb2reaction [OPTIONS]...    # equivalent to:  pdb2reaction all [OPTIONS]...
```

Two key options on the workflows that use cluster extraction:

- `-i/--input` — one or more full structures (reactant, intermediate(s), product).
- `-c/--center` — substrate / extraction center (residue names, residue IDs, or PDB paths). Omit to skip extraction and feed the full input structure directly.

## Main workflow modes

| Mode | Trigger | Use when | Quickstart |
|---|---|---|---|
| Multi-structure MEP | `-i R.pdb [I1.pdb ...] P.pdb` | You have ≥ 2 endpoints / intermediates. | [quickstart-all](quickstart-all.md) |
| Staged scan | `-i ONE.pdb --scan-lists '[...]' [ '[...]' ...]` | You'd rather define the reaction coordinates than provide endpoints. | [quickstart-scan](quickstart-scan.md) |
| TS-only | `-i TS_CANDIDATE.pdb --tsopt` | You already have a TS guess. | [quickstart-tsopt-freq](quickstart-tsopt-freq.md) |

```{important}
Single-input runs require **either** `--scan-lists/-s` or `--tsopt` — a bare `-i ONE.pdb` will not trigger a full workflow.
```

## Common options

| Option | Description |
|---|---|
| `-i, --input PATH...` | Input structures. ≥ 2 PDBs → MEP; 1 PDB + `--scan-lists` → staged scan; 1 PDB + `--tsopt` → TS-only. |
| `-c, --center TEXT` | Substrate / extraction center (residue names, residue IDs, or PDB paths). |
| `-l, --ligand-charge TEXT` | Charge mapping (`'SAM:1,GPP:-3'`) or single integer. |
| `-q, --charge INT` / `-m, --multiplicity INT` | Net system charge / spin multiplicity. |
| `--tsopt` / `--thermo` / `--dft` | TS optimisation + IRC / vibrational analysis / single-point DFT. |
| `-b, --backend uma\|orb\|mace\|aimnet2` | MLIP backend (default `uma`). |

Full option matrix: [CLI Conventions](cli-conventions.md) and the generated CLI reference under [reference/commands/index](reference/commands/index.md). Backend cost / VRAM comparison: see [Troubleshooting › Choosing a backend](troubleshooting.md#choosing-a-backend).

## Run summaries

Every `pdb2reaction all` run writes `summary.log` (human) + `summary.json` (machine) with the CLI command, global MEP statistics, per-segment barriers / bond changes, and MLIP / thermo / DFT energies (when enabled). Each `path_search/seg_NN/` (or `path_opt/` with `--refine-path False`) carries its own summaries.

## HPC / multi-GPU

`pdb2reaction` parallelises UMA inference across nodes — set `workers` and `workers_per_node` to enable multi-worker mode. Job-script templates: [docs/hpc-example.md](hpc-example.md). Backend configuration: [MLIP Calculator](uma-pysis.md).

## Agent Skills

`pdb2reaction` ships agent-readable instructions under [`skills/`](https://github.com/t-0hmura/pdb2reaction/tree/main/skills) — copy into your project as `.claude/skills/` (or `~/.claude/skills/` for user-global) to let Claude Code / Cursor / OpenCode drive the CLI end-to-end.

## Getting help

```bash
pdb2reaction <subcommand> --help               # core options
pdb2reaction <subcommand> --help-advanced      # full option set
```

`--help-advanced` is available on the calculation / scan / extract / utility commands; for the per-command index see the [documentation home](index.md#subcommands).

Issues: <https://github.com/t-0hmura/pdb2reaction/issues>.
