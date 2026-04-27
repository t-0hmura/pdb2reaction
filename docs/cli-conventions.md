# CLI Conventions

This page documents the conventions used across all `pdb2reaction` commands. Understanding these conventions helps you write correct commands and avoid common errors.

---

## Boolean Options

Boolean options are normalized at the root CLI.
Both notations are accepted:

```bash
# Recommended
--tsopt --thermo --no-dft

# Also accepted
--tsopt True --thermo yes --dft 0
```

For options that are defined only as `--flag`, the root CLI also accepts `--no-flag` and `--flag False` as compatibility aliases.
All subcommands (including `extract` and `fix-altloc`) use Click as their CLI backend.

Common boolean options:
- `--tsopt`, `--thermo`, `--dft` — enable post-processing stages
- `--freeze-links` — freeze link-hydrogen parents (default: `True`)
- `--dump` — write trajectory files
- `--preopt`, `--endopt` — pre/post optimization toggles
- `--climb` — enable climbing image in MEP search
- `--convert-files` — generate PDB/GJF companion files

---

## Progressive Help (`all`)

`pdb2reaction all` uses two help levels:

```bash
pdb2reaction all --help # core options only
pdb2reaction all --help-advanced # full option list
```

The following subcommands support progressive help (`--help` for core options, `--help-advanced` for the full list):

`all`, `scan`, `scan2d`, `scan3d`, `opt`, `path-opt`, `path-search`, `tsopt`, `freq`, `irc`, `dft`, `add-elem-info`, `trj2fig`, `energy-diagram`, `extract`, `fix-altloc`.

---

## Residue Selectors

Residue selectors identify which residues to use as substrates or extraction centers.

### By residue name
```bash
-c 'SAM,GPP' # Select all residues named SAM or GPP
-c 'LIG' # Select all residues named LIG
```

### By residue ID
```bash
-c '123,456' # Residues 123 and 456
-c 'A:123,B:456' # Chain A residue 123, Chain B residue 456
-c '123A' # Residue 123 with insertion code A
-c 'A:123A' # Chain A, residue 123, insertion code A
```

### By PDB file
```bash
-c substrate.pdb # Use coordinates from a separate PDB to locate substrates
```

```{note}
When selecting by residue name, if multiple residues share the same name, **all** matches are included and a warning is logged.
```

```{warning}
**`--selected-resn` accepts residue IDs, not names.** Despite its name, the `--selected-resn` flag (on `extract` and `all`) accepts **residue IDs** (colon-separated integers with optional chains/insertion codes, e.g. `A:123A`), **not** 3-letter residue names. Use `-c/--center 'GPP,SAM'` if you want residue-name-based substrate selection. See the [`extract` CLI options table](extract.md) for the canonical description.
```

---

(charge-specification)=
## Charge Specification

For PDB inputs, `--ligand-charge/-l` lets you specify charges only for non-standard residues (substrates, cofactors, metal ions). The total system charge is then **automatically derived** by summing standard amino-acid charges, ion charges, and your ligand charges — no need to manually count atoms across the entire complex. This is especially useful for large enzyme–substrate systems where the total charge is not obvious.

### Per-residue mapping (recommended)
```bash
-l 'SAM:1,GPP:-3' # SAM has charge +1, GPP has charge -3
-l 'LIG:-2' # LIG has charge -2
```

### Total charge override
```bash
-q 0 # Force total system charge to 0
-q -1 # Force total system charge to -1
```

### Charge resolution order
1. `-q/--charge` (explicit CLI override) — highest priority
2. Active site model (binding pocket) extraction (sums amino acids, ions, `--ligand-charge/-l`) — only when `-c/--center` is passed and extraction actually runs (e.g. `all`, `extract`)
3. `--ligand-charge/-l` as fallback (when extraction skipped)
4. `.gjf` template metadata
5. Default: none (unresolved charge aborts; provide `-q` or `.gjf` charge metadata, or use PDB `--ligand-charge/-l`)

```{note}
Step 2 (extraction-based charge derivation) only fires for commands like `all` that take `-c/--center`. For standalone subcommands such as `opt`/`tsopt`/`freq`, or when `-c` is omitted, extraction is skipped and resolution becomes `1 → 3 → 4 → 5`.
```

```{note}
`--ligand-charge/-l` derivation is only applied for PDB inputs (including XYZ/GJF inputs when `--ref-pdb` is supplied) and only when charge is **not yet resolved**. In that unresolved case, ligand-derived charge is attempted before `.gjf` metadata fallback.
```

```{tip}
Always provide `--ligand-charge/-l` for non-standard residues (substrates, cofactors, unusual ligands) to ensure correct charge propagation.
```

---

## Spin Multiplicity

```bash
-m 1 # Singlet (default)
-m 2 # Doublet
-m 3 # Triplet
```

```{note}
Use `-m/--multiplicity` consistently in `all` and other subcommands.
```

---

## Atom Selectors

Atom selectors identify specific atoms for scans and restraints. They can be:

### Integer index (1-based by default)
```bash
--scan-lists '[(1, 5, 2.0)]' # Atoms 1 and 5, target distance 2.0 Å
```

### PDB-style selector string
```bash
--scan-lists '[("TYR,285,CA", "SAM,309,C10", 2.20)]'
```

Selector fields can be separated by:
- Space: `'TYR 285 CA'`
- Comma: `'TYR,285,CA'`
- Slash: `'TYR/285/CA'`
- Backtick: `` 'TYR`285`CA' ``
- Backslash: `'TYR\285\CA'`

The three tokens (residue name, residue number, atom name) can appear in any order—the parser uses a fallback heuristic if the order is non-standard.

---

## Input File Requirements

### PDB files
- Must contain **hydrogen atoms** (use `reduce`, `pdb2pqr`, or Open Babel to add them)
- Must have **element symbols** in columns 77–78 (use `pdb2reaction add-elem-info` if missing)
- Multiple PDBs must have **identical atoms in the same order** (only coordinates may differ)

### XYZ and GJF files
- Can be used when active site model extraction is skipped (omit `-c/--center`)
- `.gjf` files can provide charge/spin defaults from embedded metadata

---

(exit-codes)=
## Exit codes

`pdb2reaction` subcommands follow a largely shared exit-code convention, but the codes a given subcommand can actually emit differ (see the "Typical emitter" column below). For the exhaustive list per subcommand, consult that subcommand's page.

| Code | Meaning | Typical emitter |
|------|---------|-----------------|
| `0` | Success | every subcommand |
| `1` | Unexpected error (any unhandled exception) | every subcommand |
| `2` | Zero step length (step norm below minimum) **or** missing dependency on import | `opt`, `tsopt`, `path-opt`; `dft` (PySCF/GPU4PySCF not installed) |
| `3` | Optimizer failure **or** SCF not converged | `opt`, `tsopt`, `path-opt`; `dft` |
| `4` | Trajectory write error | `path-opt` |
| `5` | HEI export error | `path-opt` |
| `130` | Keyboard interrupt (SIGINT) | every subcommand |

Subcommands that only use `0 / 1 / 130` (such as `irc` and `freq`) still follow the same codes — they simply do not currently raise the optimizer-specific errors.

---

(opt-mode-semantics)=
## `--opt-mode` (subcommand-dependent)

```{warning}
The same `--opt-mode` token selects **different optimizer algorithms** depending on the subcommand, and its default is **not uniform**. Always check the per-subcommand table before copying a recipe.
```

| Subcommand | `grad` alias selects | `hess` alias selects | Default |
|------------|---------------------|----------------------|---------|
| `opt` | L-BFGS (`lbfgs`) | RFO (`rfo`) | `grad` (L-BFGS) |
| `tsopt` | Dimer (`dimer`) | RS-I-RFO (`rsirfo`) | `hess` (RS-I-RFO) |
| `path-opt` (endpoint preopt) | L-BFGS | RFO | `grad` |
| `path-search` (endpoint preopt) | L-BFGS | RFO | `grad` |
| `scan` / `scan2d` / `scan3d` (endpoint preopt) | L-BFGS | RFO | `grad` |
| `all` (pre-opt stage, `--opt-mode`) | L-BFGS | RFO | `grad` |
| `all` (post-opt — TSOPT preset, `--opt-mode-post`) | Dimer (`dimer`) | RS-I-RFO (`rsirfo`) | `hess` |
| `all` (post-opt — post-IRC endpoint optimizer, `--opt-mode-post`) | L-BFGS | RFO | `hess` |

**Accepted aliases** are subcommand-specific:

- `opt` accepts `grad` / `lbfgs` and `hess` / `rfo`.
- `tsopt` accepts `grad` / `dimer` and `hess` / `rsirfo`.

As a result, `--opt-mode grad` on `tsopt` is a **Dimer** TS search, not an L-BFGS minimization. Use the explicit algorithm alias (`--opt-mode lbfgs`, `--opt-mode rsirfo`, etc.) if you want to be unambiguous across subcommands.

---

## CLI ↔ YAML name mismatches

Some CLI flags use slightly different names than their YAML counterparts, and a few are renamed when wrapped in `all`. The full mapping table lives in {ref}`YAML Reference: Common CLI-to-YAML mapping <common-cli-to-yaml-mapping>`; the two most frequently misremembered cases are:

(pressure-vs-pressure-atm)=
### `--pressure` (CLI) vs `pressure_atm` (YAML)

- **CLI flag:** `--pressure FLOAT` (on `freq`; in `all` it is exposed as `--freq-pressure`).
- **YAML key:** `thermo.pressure_atm` (explicit unit suffix).
- Both carry values in **atm** and are converted to Pa internally.

(engine-vs-dft-engine)=
### `--engine` (standalone `dft`) vs `--dft-engine` (in `all`)

- On the **standalone `dft`** subcommand the backend selector is named **`--engine`** (values: `gpu`, `cpu`).
- Inside **`pdb2reaction all`** the exact same option is renamed **`--dft-engine`** (prefix-disambiguated so it does not collide with other engine-like flags under the `all` wrapper).
- In YAML both resolve to the same `dft` section setting; see {ref}`YAML Reference: dft section <dft-section>`.

Equivalent commands:

```bash
# Standalone dft
pdb2reaction dft -i ts.xyz -q 0 --engine gpu

# Same thing inside the all wrapper
pdb2reaction all -i r.pdb p.pdb -c SAM --dft --dft-engine gpu
```

---

## YAML Configuration

Advanced settings can be passed via layered YAML inputs:

```bash
pdb2reaction -i r.pdb p.pdb -q -1 --config my_settings.yaml --out-dir result/
```

See [YAML Reference](yaml-reference.md) for all available options.

(configuration-precedence)=
### Configuration precedence

Settings are resolved in the following order (later sources override earlier ones):

```
built-in defaults  <  --config (YAML)  <  CLI options
```

- **Built-in defaults** — hard-coded values for every parameter (see `pdb2reaction/defaults.py`).
- **`--config`** — a YAML file that overrides defaults. Useful for site-wide or project-wide settings.
- **CLI options** — explicit flags on the command line (e.g., `--backend orb`). Only *explicitly supplied* values override YAML; options left at their CLI default do not mask YAML values.

This precedence applies uniformly to `all`, `opt`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, and `dft`. See also {ref}`YAML Reference: Configuration precedence <yaml-configuration-precedence>`.

---

## Output Directory

Use `-o/--out-dir` to specify where results are saved:

```bash
-o ./my_results/ # Custom output directory
```

Default output directories:
- `all`: `./result_all/`
- `extract`: current directory or specified `-o`
- `opt`: `./result_opt/`
- `tsopt`: `./result_tsopt/`
- `path-opt`: `./result_path_opt/`
- `path-search`: `./result_path_search/`
- `scan`: `./result_scan/`
- `scan2d`: `./result_scan2d/`
- `scan3d`: `./result_scan3d/`
- `freq`: `./result_freq/`
- `irc`: `./result_irc/`
- `dft`: `./result_dft/`

---

## See Also

- [Installation](installation.md) — Setup and dependency installation
- [Getting Started](getting-started.md) — First run and workflow overview
- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Common errors and fixes
- [YAML Reference](yaml-reference.md) — Complete configuration options
