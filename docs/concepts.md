# Concepts & Workflow

This page explains the key terms in pdb2reaction—pockets, templates, segments, and images—and how the `all` command ties together the subcommands.

---

## Workflow at a glance

Most workflows follow this flow:

```text
Full system(s) (PDB/XYZ/GJF)
 │
 ├─ (optional) pocket extraction [`extract`](extract.md) ← requires PDB when you use --center/-c
 │ ↓
 │ Pocket/cluster model(s) (PDB)
 │ │
 │ ├─ (optional) staged scan [`scan`](scan.md) ← single-structure workflows
 │ │ ↓
 │ │ Ordered intermediates
 │ │ ↓
 │ └─ MEP search [`path-search`](path_search.md) or [`path-opt`](path_opt.md)
 │ ↓
 │ MEP trajectory (mep_trj.xyz) + energy diagrams
 │ ↓
 └─ (optional) TS optimization + IRC [`tsopt`](tsopt.md) → [`irc`](irc.md)
 └─ (optional) thermo [`freq`](freq.md)
 └─ (optional) single-point DFT [`dft`](dft.md)
```

Each stage is available as an individual subcommand. The `pdb2reaction all` command runs many stages end-to-end.

```{important}
Transition states (first-order saddle points): treat the Highest-Energy Image (HEI) / `tsopt` outputs as **TS candidates** until validated via `irc` (endpoints reach intended minima). `tsopt` already performs a final imaginary-frequency check internally — look for exactly one imaginary frequency (|ν| ≥ 100 cm⁻¹) in its output. If multiple imaginary frequencies remain, consider applying `--flatten`.
```

### MLIP backends

pdb2reaction supports multiple machine-learning interatomic potential (MLIP) backends. Pass `-b/--backend` to any calculation command to choose one:

| Backend | Flag | Install | Notes |
|---------|------|---------|-------|
| **UMA** (default) | `-b uma` | included | Full feature set including analytical Hessians and multi-worker inference |
| **ORB** | `-b orb` | `pip install 'pdb2reaction[orb]'` | orb-models; FD Hessians only |
| **MACE** | `-b mace` | separate env (see README) | mace-torch; conflicts with fairchem-core |
| **AIMNet2** | `-b aimnet2` | `pip install 'pdb2reaction[aimnet2]'` | aimnet; FD Hessians only |

All backends share the same `--solvent` option for xTB-based implicit solvent corrections.

---

## Key objects and terms

### Full system vs. pocket (cluster model)
- **Full system**: your original structure(s). In enzyme use cases this is typically a protein–ligand complex.
- **Pocket**: the extraction region around the substrate(s), defined by `-c/--center` and `-r/--radius`.
- **Cluster model**: the computational subsystem cut from the pocket. Severed bonds are capped with link hydrogens, and the model is used for MEP/TS search.

Pocket extraction is controlled by:
- `-c/--center`: how to locate the substrate (residue IDs, residue names, or a substrate-only PDB).
- `-r/--radius`, `--radius-het2het`, `--include-H2O`, `--exclude-backbone`, `--add-linkH`, `--selected-resn`.

For charge and spin specification, see [CLI Conventions: Charge specification](cli_conventions.md#charge-specification).

### Images and segments
- **Image**: a single geometry (one “node”) along a chain-of-states path.
- **Segment**: an MEP between two adjacent endpoints (e.g., R → I1, I1 → I2, …). A multi-structure run is decomposed into segments.

### Templates and file conversion (`--convert-files`)
`pdb2reaction` often writes a **trajectory** (e.g., `mep_trj.xyz`, `irc_trj.xyz`). When you supply topology-aware inputs (PDB templates or Gaussian inputs), it can optionally write companion files:
- `.pdb` companions when a PDB template exists
- `.gjf` companions when a Gaussian template exists

This behavior is controlled globally by `--convert-files/--no-convert-files` (default: `True`).

---

### Link hydrogen and frozen atoms

When pdb2reaction extracts a pocket from a larger structure, severed bonds are capped with **link hydrogens**. By default (`--freeze-links`), the parent atoms of these link hydrogens are frozen during optimization and path searches to prevent unphysical rearrangement at the boundary.

- **Forces**: frozen atoms receive zeroed forces.
- **Hessian**: frozen degrees of freedom are either removed (`return_partial_hessian: true`) or zeroed in the full matrix.
- **Vibrational analysis**: when frozen atoms are present, `freq` automatically performs Partial Hessian Vibrational Analysis (PHVA), diagonalizing only the active block of the Hessian.

Frozen atoms can also be set manually via the `geom.freeze_atoms` YAML key (1-based indices). CLI-detected link atoms are merged with YAML-specified atoms.

## Three common workflow modes

### 1) Multi-structure MEP search (R → … → P)
Use this when you already have **two or more** full structures along a reaction coordinate.

Typical command:

```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
```

### 2) Single-structure staged scan → MEP
Use this when you only have **one** structure, but you can define a scan that generates endpoints.

Typical command:

```bash
pdb2reaction -i holo.pdb -c '308,309' -l 'MMT:-1' \
 --scan-lists '[("TYR,285,CA","SAM,309,C10",2.20)]'
```

### 3) TSOPT-only mode (pocket TS optimization)
Use this when you already have a TS candidate (or want a quick TS optimization on one structure).

Typical command:

```bash
pdb2reaction -i ts_guess.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt
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
- For enzyme use cases, you usually want hydrogens present in the input PDB.
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
