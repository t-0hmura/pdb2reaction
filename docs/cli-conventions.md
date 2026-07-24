# CLI Conventions

Conventions shared across all `pdb2reaction` commands.

## Boolean options

Use paired flags in commands, documentation, and agent instructions:

| Form | Example |
|---|---|
| Positive flag | `--tsopt` |
| Negative flag | `--no-tsopt` |
| Enable | `--tsopt` |
| Disable | `--no-tsopt` |

```bash
--tsopt --thermo --no-dft
```

Older value-style invocations remain accepted for compatibility, but are not
the canonical syntax and must not be emitted by docs, skills, or generated
commands. The root CLI normalizes both forms and
`tests/test_bool_compat_cli.py` keeps the compatibility path covered.

Common toggles: `--tsopt` / `--thermo` / `--dft` (post-processing stages) · `--freeze-links` (freeze cap-H parents, default `True`) · `--dump` (write trajectory files) · `--preopt` / `--endopt` (pre/post optimization) · `--climb` (climbing-image MEP) · `--convert-files` (generate format-aware PDB / CIF / GJF companions).

### Contributing a new bool flag

When adding a boolean flag inside a subcommand, always route it through one of the `add_*_option()` factories in `pdb2reaction/cli/common_options.py` and register the long name in the matching `_COMMAND_BOOL_*_OPTIONS` table in `pdb2reaction/cli/app.py`. Avoid writing `@click.option("--foo/--no-foo", ...)` or `type=click.BOOL` directly in the subcommand body — that bypasses the registry and falls out of compatibility-test coverage.

## Progressive help

```bash
pdb2reaction <subcmd> --help               # core options
pdb2reaction <subcmd> --help-advanced      # full option set
```

Supported by `all`, `scan` / `scan2d` / `scan3d`, `opt`, `path-opt`, `path-search`, `tsopt`, `freq`, `irc`, `dft`, `sp`, `add-elem-info`, `trj2fig`, `energy-diagram`, `bond-summary`, `extract`, `fix-altloc`.

(verbosity-levels)=

## Verbosity levels

`-v/--verbose LEVEL` is an integer from 0 to 3 (**default 2**) that sets how much each command prints to the console. It is a per-command option, so write it with the subcommand, e.g. `pdb2reaction opt -v 1 ...`. The same four levels apply to every command; command pages describe only their command-specific payload (e.g. the `opt` cycle table or the `freq` thermochemistry summary).

| Level | What you see |
|---|---|
| `-v 0` | Silent. Confirm success from the exit code and the output artifacts. |
| `-v 1` | Milestones only: version, input summary, key settings, output location, dry-run / final status. No banner, `[command]`, `[mode]`, or config dump. |
| `-v 2` | Default. Adds the banner, `[command]`, `[mode]`, stage progress, the main optimizer cycle table, terminal status, the one-line Hessian summary, thermo / DFT summaries, and elapsed time. |
| `-v 3` | Debug: resolved config / dry-run plan, backend DEBUG, raw optimizer and internal-coordinate chatter, `[HessianTiming]`, and `[HessianVRAM]`. |

A semantic failure is a failure at any level: a `Traceback` that appears only at `-v 3` still means the run failed.

## Residue selectors

| Form | Example | Notes |
|---|---|---|
| By residue name | `-c 'SAM,GPP'` / `-c 'LIG'` | If multiple residues share a name, **all** matches are included (warning logged). |
| By residue ID | `-c '123,456'` / `-c 'A:123,B:456'` / `-c '123A'` / `-c 'A:123A'` | Optional chain prefix; trailing letter = insertion code. |
| By chain + name | `-c 'A:SAM'` / `-c 'A:SAM:123'` | First form selects all SAM in chain A; add resSeq to select one. |
| By structure file | `-c substrate.pdb` / `-c substrate.cif` | Use coordinates from a separate PDB/mmCIF to locate substrates. |

(selected-resn-takes-ids)=
### `--selected-resn` uses the same residue selectors

`--selected-resn` on `extract` and `all` accepts numeric IDs, residue names,
and chain-qualified names. For example, `A:123A` force-includes one insertion-
code residue, `A:SAM` includes every SAM in chain A, and `A:SAM:123` narrows
that selection to one residue. An unqualified name such as `TYR` includes all
matches and warns when more than one is present.

(charge-specification)=

## Charge specification

For PDB/mmCIF inputs, `--ligand-charge/-l` lets you specify charges only for non-standard residues (substrates, cofactors, metal ions). The total system charge is then **automatically derived** by summing standard amino-acid charges, ions, and your ligand charges.

```bash
-l 'SAM:1,GPP:-3'        # per-residue mapping (recommended; `=` separator also accepted)
-l 'LIG:-2'              # single mapping
-l -3                    # single integer = total ligand charge
-q 0                     # explicit total system charge override
```

**Resolution order** (highest priority first):

1. `-q/--charge` or `--ligand-charge/-l` supplied explicitly on the CLI.
2. `calc.charge` from `--config` (when neither CLI charge form is supplied).
3. Active-site model extraction summary (amino acids + ions + `--ligand-charge`; only when `-c/--center` is set and extraction runs — e.g. `all`, `extract`). In `all`, `-q` or `calc.charge` is checked as an assertion against this value.
4. `.gjf` template metadata.
5. Abort if unresolved.

For per-stage subcommands (`opt` / `tsopt` / `freq` …) or when `-c` is omitted, extraction is skipped. Ligand-derived or configured charge is resolved before `.gjf` metadata.

```{tip}
Always provide `--ligand-charge/-l` for non-standard residues (substrates, cofactors, unusual ligands) to ensure correct charge propagation.
```

## Spin multiplicity

```bash
-m 1    # singlet (default)
-m 2    # doublet
-m 3    # triplet
```

Use `-m/--multiplicity` consistently in `all` and per-stage subcommands.

## Atom selectors

```bash
--scan-lists '[(1, 5, 2.0)]'                                          # 1-based integer indices
--scan-lists '[("TYR,285,CA", "SAM,309,C10", 2.20)]'                  # PDB-style selector strings
--scan-lists '[("A:TYR:285:CA", "B:SAM:309:C10", 2.20)]'              # chain-qualified
```

Three-field selector delimiters are space · comma · slash · backtick ·
backslash, and their residue-name / residue-number / atom-name tokens may
appear in any order. To disambiguate repeated names or numbering, use the
positional four-field form `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`.

(scan-list-spec)=

### Scan-list spec

`--scan-lists/-s` (on `scan`, `scan2d`, `scan3d`, and `all`) accepts one or more inline Python literals. The standalone `scan` / `scan2d` / `scan3d` commands additionally accept a YAML / JSON spec file path; use a file for complex multi-stage runs, inline literals for short cases.

**YAML / JSON spec file** (root = mapping; key is `stages` for `scan`, `pairs` for `scan2d` / `scan3d`):

```yaml
one_based: true            # optional; defaults to CLI --one-based
stages:                    # scan
  - [[1, 5, 1.35]]
  - [[1, 5, 2.20], [2, 8, 1.80]]
```

```yaml
one_based: true
pairs:                     # scan2d (exactly 2 entries) / scan3d (exactly 3 entries)
  - [1, 5, 1.30, 3.10]
  - [2, 8, 1.20, 3.20]
```

Each `scan` stage is a list of `(i, j, target_Å)` triples; each `scan2d` / `scan3d` axis is `(i, j, low_Å, high_Å)`. Indices may be integers, three-field selectors, or positional `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` selectors.

**Inline literal**: wrap in **single quotes** so the shell does not interpret parens / spaces; use double-quoted PDB selectors inside.

```bash
-s '[(atom1, atom2, target_Å), ...]'             # scan: triples
-s '[(atom1, atom2, low_Å, high_Å), ...]'        # scan2d / scan3d: quadruples
-s '[("TYR,285,CA","SAM,309,C10",1.35)]'         # quoted selectors
-s "[(\"TYR,285,CA\",\"SAM,309,C10\",1.35)]"       # avoid: double-quoted outer literal requires escaping inner quotes
```

For `scan`, one literal = one **stage**; multiple stages → multiple literals after a single `--scan-lists` flag. For `scan2d` / `scan3d`, only one literal is accepted (no multi-stage support).

## Input file requirements

- **PDB** — must contain hydrogens (add via `reduce` / `pdb2pqr` / Open Babel) and element symbols in cols 77–78 (`pdb2reaction add-elem-info` if missing). Multiple PDBs must share identical atoms in the same order.
- **mmCIF** — use `.cif` / `.mmcif` for multi-character chains, large residue/atom identifiers, or 10,000+ residues. The common bridge calculates through an internal PDB and restores original IDs in CIF output. Reaction-ordered inputs still require identical atoms and order.
- **XYZ / GJF** — accepted when active-site extraction is skipped (omit `-c/--center`). `.gjf` files can provide charge / spin defaults from embedded metadata.

(exit-codes)=

## Exit codes

| Code | Meaning | Typical emitter |
|---|---|---|
| `0` | Success | every subcommand |
| `1` | Unexpected error (any unhandled exception) | every subcommand |
| `2` | Zero step length **or** missing import dependency | `opt`, `tsopt`; `dft` (PySCF / GPU4PySCF not installed) |
| `3` | Optimizer failure **or** SCF not converged | `opt`, `tsopt`, `path-opt`; `dft` |
| `4` | Trajectory write error | `path-opt` |
| `5` | HEI export error | `path-opt` |
| `130` | Keyboard interrupt (SIGINT) | every subcommand |

Subcommands that only use `0 / 1 / 130` (e.g. `irc`, `freq`) follow the same scheme; they simply don't currently raise the optimizer-specific errors.

(opt-mode-semantics)=

## `--opt-mode` (subcommand-dependent)

```{warning}
The same `--opt-mode` token selects **different algorithms** by subcommand, and defaults differ. Always check the table before copying a recipe.
```

| Subcommand | `grad` alias selects | `hess` alias selects | Default |
|---|---|---|---|
| `opt` | L-BFGS (`lbfgs`) | RFO (`rfo`) | `grad` (L-BFGS) |
| `tsopt` | Dimer (`dimer`) | RS-P-RFO (`rsprfo`) | `hess` (RS-P-RFO) |
| `path-opt` (endpoint preopt) | L-BFGS | RFO | `grad` |
| `path-search` (HEI±1 / kink-node single-structure) | L-BFGS | RFO | `grad` |
| `scan` / `scan2d` / `scan3d` (per-grid relaxation) | L-BFGS | RFO | `grad` |
| `all` (pre-opt, `--opt-mode`) | L-BFGS | RFO | `grad` |
| `all` (TSOPT preset, `--opt-mode-post`) | Dimer | RS-P-RFO | `hess` |
| `all` (post-IRC endpoint, `--opt-mode-post`) | L-BFGS | RFO | `hess` |

Algorithm aliases are accepted on `opt` (`lbfgs` / `rfo`) and `tsopt` (`dimer` / `rsirfo` / `trim` / `rsprfo`); all other subcommands accept only `grad` / `hess`. So `--opt-mode grad` on `tsopt` is a **Dimer** TS search, not L-BFGS minimization — use `--opt-mode dimer|rsirfo` on `tsopt` and `--opt-mode lbfgs|rfo` on `opt` to be unambiguous.

## CLI ↔ YAML name mismatches

A few CLI flags use slightly different names than their YAML counterparts, and a few are renamed when wrapped in `all`. Full mapping table: {ref}`YAML Reference › Common CLI-to-YAML mapping <common-cli-to-yaml-mapping>`. The two most-asked cases:

(pressure-vs-pressure-atm)=
- **`--pressure` (CLI) vs `pressure_atm` (YAML)** — on `freq` the flag is `--pressure FLOAT`; in `all` it is exposed as `--freq-pressure`. YAML key: `thermo.pressure_atm`. Both carry **atm** values (converted to Pa internally).

(engine-vs-dft-engine)=
### `--engine` (`dft`) vs `--dft-engine` (`all`)

- **`--engine` (`dft`) vs `--dft-engine` (`all`)** — standalone `dft` accepts `--engine gpu|cpu`; inside `pdb2reaction all` it is renamed `--dft-engine` (prefix-disambiguated). Both map to the same YAML `dft` section setting.

```bash
pdb2reaction dft -i ts.xyz -q 0 --engine gpu                                  # standalone
pdb2reaction all -i r.pdb p.pdb -c SAM --tsopt --dft --dft-engine gpu         # after TS optimization in `all`
```

## YAML configuration

```bash
pdb2reaction -i r.pdb p.pdb -q -1 --config my_settings.yaml --out-dir result/
```

(configuration-precedence)=

```
built-in defaults  <  --config (YAML)  <  CLI options
```

Built-in defaults are in `pdb2reaction/core/defaults.py`. Only *explicitly supplied* CLI values override YAML; options left at their CLI default do not mask YAML values. Applies uniformly to all calc subcommands. Full schema: [YAML Reference](yaml-reference.md).

- **Known default exception**: `flatten_max_iter` starts at 0 before YAML is
  applied. An omitted toggle therefore retains an explicit YAML value;
  `--flatten` enables the configured/built-in positive value and
  `--no-flatten` forces 0. See {ref}`flatten-precedence-caveat`.

## Output directory

`-o/--out-dir ./my_results/` overrides. Defaults: `all → ./result_all/`, per-stage subcommand → `./result_<subcmd>/` with hyphens in subcommand names mapped to underscores (e.g. `result_opt/`, `result_tsopt/`, `result_path_opt/`, `result_path_search/`). `extract` defaults to the current directory or the explicit `-o`.

## See Also

[Installation](installation.md) · [Getting Started](getting-started.md) · [Common Error Recipes](recipes-common-errors.md) · [Troubleshooting](troubleshooting.md) · [YAML Reference](yaml-reference.md).
