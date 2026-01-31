# CLI Conventions

This page documents the conventions used across all `pdb2reaction` commands. Understanding these conventions helps you write correct commands and avoid common errors.

---

## Boolean Options

All boolean CLI options require explicit `True` or `False` (case-sensitive):

```bash
# Correct
--tsopt True --thermo True --dft False

# Wrong (will not work)
--tsopt true    # lowercase
--tsopt 1       # numeric
--tsopt         # flag-style (no value)
```

Common boolean options:
- `--tsopt`, `--thermo`, `--dft` — enable post-processing stages
- `--freeze-links` — freeze link-hydrogen parents (default: `True`)
- `--dump` — write trajectory files
- `--preopt`, `--endopt` — pre/post optimization toggles
- `--climb` — enable climbing image in MEP search
- `--convert-files` — generate PDB/GJF companion files

---

## Residue Selectors

Residue selectors identify which residues to use as substrates or extraction centers.

### By residue name
```bash
-c 'SAM,GPP'          # Select all residues named SAM or GPP
-c 'LIG'              # Select all residues named LIG
```

### By residue ID
```bash
-c '123,456'          # Residues 123 and 456
-c 'A:123,B:456'      # Chain A residue 123, Chain B residue 456
-c '123A'             # Residue 123 with insertion code A
-c 'A:123A'           # Chain A, residue 123, insertion code A
```

### By PDB file
```bash
-c substrate.pdb      # Use coordinates from a separate PDB to locate substrates
```

```{note}
When selecting by residue name, if multiple residues share the same name, **all** matches are included and a warning is logged.
```

---

## Charge Specification

### Per-residue mapping (recommended)
```bash
--ligand-charge 'SAM:1,GPP:-3'    # SAM has charge +1, GPP has charge -3
--ligand-charge 'LIG:-2'          # LIG has charge -2
```

### Total charge override
```bash
-q 0                              # Force total system charge to 0
-q -1                             # Force total system charge to -1
```

### Charge resolution order
1. `-q/--charge` (explicit CLI override) — highest priority
2. Pocket extraction (sums amino acids, ions, `--ligand-charge`)
3. `--ligand-charge` as fallback (when extraction skipped)
4. `.gjf` template metadata
5. Default: 0

```{tip}
Always provide `--ligand-charge` for non-standard residues (substrates, cofactors, unusual ligands) to ensure correct charge propagation.
```

---

## Spin Multiplicity

```bash
-m 1      # Singlet (default)
-m 2      # Doublet
-m 3      # Triplet
```

```{note}
In the `all` command, use `-m/--mult`. In other subcommands, use `-m/--multiplicity`.
```

---

## Atom Selectors

Atom selectors identify specific atoms for scans and restraints. They can be:

### Integer index (1-based by default)
```bash
--scan-lists '[(1, 5, 2.0)]'      # Atoms 1 and 5, target distance 2.0 Å
```

### PDB-style selector string
```bash
--scan-lists '[("TYR,285,CA", "MMT,309,C10", 2.20)]'
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
- Must have **element symbols** in columns 77-78 (use `pdb2reaction add-elem-info` if missing)
- Multiple PDBs must have **identical atoms in the same order** (only coordinates may differ)

### XYZ and GJF files
- Can be used when pocket extraction is skipped (omit `-c/--center`)
- `.gjf` files can provide charge/spin defaults from embedded metadata

---

## YAML Configuration

Advanced settings can be passed via `--args-yaml`:

```bash
pdb2reaction all -i R.pdb P.pdb -c 'LIG' --args-yaml config.yaml
```

YAML values take **highest precedence**:
```
defaults → CLI options → YAML (wins)
```

See [YAML Reference](yaml-reference.md) for all available options.

---

## Output Directory

Use `--out-dir` to specify where results are saved:

```bash
--out-dir ./my_results/    # Custom output directory
```

Default output directories:
- `all`: `./result_all/`
- `extract`: current directory or specified `-o`
- `opt`: `./result_opt/`
- `tsopt`: `./result_tsopt/`
- `path-opt`: `./result_path_opt/`
- `path-search`: `./result_path_search/`
- `scan`: `./result_scan/`
- `freq`: `./result_freq/`
- `irc`: `./result_irc/`
- `dft`: `./result_dft/`

---

## See Also

- [Getting Started](getting-started.md) — Installation and first run
- [Troubleshooting](troubleshooting.md) — Common errors and fixes
- [YAML Reference](yaml-reference.md) — Complete configuration options
