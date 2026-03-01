# Concepts & Workflow

This page explains the key terms in pdb2reaction—pockets, templates, segments, and images—and how the `all` command ties together the subcommands.

---

## Workflow at a glance

Most workflows follow this flow:

```text
Full system(s) (PDB/XYZ/GJF)
 │
 ├─ (optional) pocket extraction [extract] ← requires PDB when you use --center/-c
 │ ↓
 │ Pocket/cluster model(s) (PDB)
 │ │
 │ ├─ (optional) staged scan [scan] ← single-structure workflows
 │ │ ↓
 │ │ Ordered intermediates
 │ │ ↓
 │ └─ MEP search [path-search] or [path-opt]
 │ ↓
 │ MEP trajectory (mep_trj.xyz) + energy diagrams
 │ ↓
 └─ (optional) TS optimization + IRC [tsopt] → [irc]
 └─ (optional) thermo [freq]
 └─ (optional) single-point DFT [dft]
```

Each stage is available as an individual subcommand. The `pdb2reaction all` command runs many stages end-to-end.

```{important}
Transition states: treat HEI / `tsopt` outputs as **TS candidates** until validated via `freq` (a single imaginary mode) and `irc` (endpoints reach intended minima).
```

---

## Key objects and terms

### Full system vs. pocket (cluster model)
- **Full system**: your original structure(s). In enzyme use-cases this is typically a protein–ligand complex.
- **Pocket / cluster model**: a truncated structure around the substrate(s) used to reduce system size for MEP/TS search.

Pocket extraction is controlled by:
- `-c/--center`: how to locate the substrate (residue IDs, residue names, or a substrate-only PDB).
- `-r/--radius`, `--radius-het2het`, `--include-H2O`, `--exclude-backbone`, `--add-linkH`, `--selected-resn`.

### Images and segments
- **Image**: a single geometry (one “node”) along a chain-of-states path.
- **Segment**: an MEP between two adjacent endpoints (e.g., R → I1, I1 → I2, …). A multi-structure run is decomposed into segments.

### Templates and file conversion (`--convert-files`)
`pdb2reaction` often writes a **trajectory** (e.g., `mep_trj.xyz`, `irc_trj.xyz`). When you supply topology-aware inputs (PDB templates or Gaussian inputs), it can optionally write companion files:
- `.pdb` companions when a PDB template exists
- `.gjf` companions when a Gaussian template exists

This behavior is controlled globally by `--convert-files/--no-convert-files` (default: `True`).

---

## Three common workflow modes

### 1) Multi-structure MEP search (R → … → P)
Use this when you already have **two or more** full structures along a reaction coordinate.

Typical command:

```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

### 2) Single-structure staged scan → MEP
Use this when you only have **one** structure, but you can define a scan that generates endpoints.

Typical command:

```bash
pdb2reaction -i holo.pdb -c '308,309' \
 --scan-lists '[("TYR,285,CA","MMT,309,C10",2.20)]'
```

### 3) TSOPT-only mode (pocket TS optimization)
Use this when you already have a TS candidate (or want a quick TS optimization on one structure).

Typical command:

```bash
pdb2reaction -i ts_guess.pdb -c 'SAM,GPP' --tsopt
```

---

## When to use `all` vs individual subcommands

### Prefer `pdb2reaction all` when…
- You want an **end-to-end** run (extract → MEP → TSOPT/IRC → freq/DFT).
- You are still exploring the workflow and want a single command to manage outputs.

### Prefer subcommands when…
- You want to **debug** a specific stage (e.g., only `extract`, only `path-search`).
- You want to mix-and-match a custom workflow (e.g., your own endpoint preparation).

---

## A few CLI conventions worth knowing

```{important}
- Boolean options accept both `--flag` / `--no-flag` and value style `--flag True/False` (`yes/no`, `1/0` are also accepted). Prefer toggle style.
- With multiple PDB inputs, all files should have the **same atoms in the same order** (only coordinates differ).
- For enzyme use-cases, you usually want hydrogens present in the input PDB.
```

---

## Next steps

### Getting started
- [Getting Started](getting_started.md) — quick start and workflow overview
- [Installation](installation.md) — setup and dependencies
- [Common Error Recipes](recipes_common_errors.md) — symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — common errors and fixes

### Core subcommands
| Subcommand | Purpose | Documentation |
|------------|---------|---------------|
| `all` | End-to-end workflow | [all.md](all.md) |
| `extract` | Pocket extraction | [extract.md](extract.md) |
| `path-search` | Recursive MEP search | [path_search.md](path_search.md) |
| `tsopt` | TS optimization | [tsopt.md](tsopt.md) |
| `freq` | Vibrational analysis | [freq.md](freq.md) |
| `dft` | Single-point DFT | [dft.md](dft.md) |

### Reference
- [YAML Reference](yaml_reference.md) — complete YAML configuration options
- [Glossary](glossary.md) — terminology reference
