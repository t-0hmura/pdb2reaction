# Common options and selectors

This page collects the conventions shared by every `pdb2reaction` command: flags, residue and atom selectors, charge and multiplicity, exit codes, and configuration precedence.

## Boolean options

Turn a stage or behavior on or off with paired flags:

| Form | Example |
|---|---|
| Positive flag | `--tsopt` |
| Negative flag | `--no-tsopt` |

```bash
--tsopt --thermo --no-dft
```

Common toggles:

- `--tsopt` / `--thermo` / `--dft`: post-processing stages
- `--freeze-links`: freeze the parent atoms of the cap hydrogens (on by default)
- `--dump`: write trajectory files
- `--preopt` / `--endopt`: pre- / post-optimization
- `--climb`: climbing image in the minimum energy path search
- `--convert-files`: also write PDB / CIF / GJF copies of the outputs in the input format

## Progressive help

```bash
pdb2reaction <subcmd> --help               # core options
pdb2reaction <subcmd> --help-advanced      # full option set
```

(verbosity-levels)=

## Verbosity levels

`-v/--verbose LEVEL` is an integer from 0 to 3 (**default 2**) that sets how much each command prints to the console. It is a per-command option, so write it with the subcommand, e.g. `pdb2reaction opt -v 1 ...`. The same four levels apply to every command; command pages describe only what their own command adds.

| Level | What you see |
|---|---|
| `-v 0` | Silent. Confirm success from the exit code and the output files. |
| `-v 1` | Milestones only: version, input summary, key settings, output location, dry-run / final status. No banner, `[command]`, `[mode]`, or configuration printout. |
| `-v 2` | Default. Adds the banner, `[command]`, `[mode]`, stage progress, the main optimizer cycle table, terminal status, the one-line Hessian summary, thermo / DFT summaries, and elapsed time. |
| `-v 3` | Debug: the full configuration in effect, backend DEBUG, raw optimizer and internal-coordinate output, `[HessianTiming]`, and `[HessianVRAM]`. |

The level changes only what is printed, not the exit code; judge a run by its {ref}`exit code <exit-codes>`.

## Residue selectors

`-c/--center` (on `extract` and `all`) names the residues at the center of the model. The forms below run from the most specific to the broadest:

| Form | Example | What it selects |
|---|---|---|
| Chain + name + number (recommended) | `-c 'A:TYR:44'` / `-c 'A:TYR:44,A:SAM:123'` | Exactly one residue per entry. |
| Chain + name | `-c 'A:SAM'` | Every SAM in chain A; a warning is logged when more than one matches. |
| Chain + number | `-c 'A:123'` / `-c 'A:123,B:456'` / `-c 'A:123A'` | Residue 123 of chain A; a trailing letter is the insertion code. |
| Name only | `-c 'SAM,GPP'` / `-c 'LIG'` | Every residue with that name in any chain; a warning is logged when more than one matches. |
| Number only | `-c '123,456'` / `-c '123A'` | The residue with that number in every chain. |
| Structure file | `-c substrate.pdb` / `-c substrate.cif` | The residues whose coordinates match a separate PDB / mmCIF file. |

Long mmCIF chain IDs and residue numbers above 9999 use the same forms. Chain IDs are case-sensitive; residue names are not. A PDB with an empty chain column, such as the bundled example PDBs, takes only the name or number forms (`-c 'SAM,GPP,MG'`, `--selected-resn '44,63,186'`).

```bash
pdb2reaction extract -i complex.cif -c 'LONG_CHAIN:SAM' -o model.pdb        # every SAM in chain LONG_CHAIN
pdb2reaction extract -i complex.cif -c 'LONG_CHAIN:SAM:10001' -o model.pdb  # one SAM
pdb2reaction extract -i complex.cif -c 'LONG_CHAIN:10001' -o model.pdb      # chain + number
```

(selected-resn-takes-ids)=
### `--selected-resn` uses the same residue selectors

`--selected-resn` on `extract` and `all` force-includes residues in the model and accepts the same forms. For example, `A:TYR:44` includes one residue, `A:SAM` every SAM in chain A, and `A:123A` one insertion-code residue. A name without a chain, such as `TYR`, includes every match and warns when more than one is present.

(charge-specification)=

## Charge specification

For PDB/mmCIF inputs, `--ligand-charge/-l` lets you specify charges only for non-standard residues (substrates, cofactors, metal ions). The total system charge is then **automatically derived** by summing standard amino-acid charges, ions, and your ligand charges.

```bash
-l 'SAM:1,GPP:-3'        # per-residue mapping (recommended)
-l 'LIG:-2'              # single mapping
-l -3                    # single integer = total ligand charge
-q 0                     # explicit total system charge
```

**Resolution order** (highest priority first):

1. Explicit `-q/--charge`.
2. With `--ligand-charge/-l`: the total of the standard residues, ions, and your ligand charges in the PDB/mmCIF input (with `all -c`, in the extracted model).
3. `calc.charge` from `--config`.
4. Without `-l`: with `all -c`, the total of the standard residues and ions in the extracted model, other residues counted as 0; for a `.gjf` input, the charge in its template.
5. Otherwise, stop with an error.

The console prints the derived charge as `Total active site model charge`; [Check the model](model-setup.md#check-the-model) lists the lines to read after `extract`.

```{tip}
Always provide `--ligand-charge/-l` for non-standard residues to ensure correct charge propagation.
```

## Spin multiplicity

```bash
-m 1    # singlet (default)
-m 2    # doublet
-m 3    # triplet
```

Use `-m/--multiplicity` consistently in `all` and per-stage subcommands.

## Atom selectors

Atom selectors name single atoms in `--scan-lists` and in `--distance-restraint` of `opt`; `--freeze-atoms` takes only 1-based atom numbers (see [Freeze atoms and restrain distances](model-setup.md#freeze-atoms-and-restrain-distances)).

```bash
--scan-lists '[(1, 5, 2.0)]'                                          # 1-based integer indices
--scan-lists '[("SAM,320,CS1", "GPP,321,C7", 1.60)]'                  # residue name, number, atom name
--scan-lists '[("A:SAM:320:CS1", "A:GPP:321:C7", 1.60)]'              # with chain ID
```

A three-field selector gives the residue name, residue number, and atom name in any order, separated by spaces, commas, colons, slashes, backticks, or backslashes (`"SAM,320,CS1"`, `"SAM 320 CS1"`, and `"320,SAM,CS1"` select the same atom). Three fields never include a chain; to name the chain, use the four-field form `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` in this order, with any insertion code after the number (`A:SAM:12B:C1`).

(scan-list-spec)=

### Scan-list spec

On `scan` and `all`, `--scan-lists/-s` accepts one or more inline Python literals; `scan2d` and `scan3d` accept exactly one. The standalone `scan` / `scan2d` / `scan3d` commands additionally accept a YAML / JSON spec file path; use a file for complex multi-stage runs, inline literals for short cases.

**YAML / JSON spec file** (root = mapping; key is `stages` for `scan`, `pairs` for `scan2d` / `scan3d`):

```yaml
one_based: true            # optional; defaults to the command's --one-based/--zero-based (1-based)
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

Each `scan2d` / `scan3d` axis is `(i,j,low,high)`, `(i,j,k,low,high)`, or `(i,j,k,l,low,high)`. Indices may be integers, three-field selectors, or positional `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` selectors.

**Inline literal**: wrap in **single quotes** so the shell does not interpret parens / spaces; use double-quoted PDB selectors inside.

```bash
-s '[(atom1, atom2, target_Å), ...]'             # scan: triples
-s '[(atom1, atom2, low_Å, high_Å), ...]'        # distance range
-s '[(atom1, atom2, atom3, low_deg, high_deg)]'  # angle range
-s '[("SAM,320,CS1","GPP,321,C7",1.60)]'         # quoted selectors
-s "[(\"SAM,320,CS1\",\"GPP,321,C7\",1.60)]"       # avoid: double-quoted outer literal requires escaping inner quotes
```

For `scan`, one literal = one **stage**; multiple stages → multiple literals after a single `--scan-lists` flag. For `scan2d` / `scan3d`, only one literal is accepted.

| Command | Accepted scan specification |
| --- | --- |
| `scan` | Inline distance targets `(i,j,target)`, or ranges `(i,j,low,high)`, `(i,j,k,low,high)`, and `(i,j,k,l,low,high)` scanned in both directions from the input geometry ([Bidirectional scan](scan.md#bidirectional-scan-4-tuple)); YAML/JSON is also accepted |
| `all --scan-lists` | Inline targets only: distance `(i,j,target)`, angle `(i,j,k,deg)`, or dihedral `(i,j,k,l,deg)` (no ranges, no YAML/JSON) |
| `scan2d` | One literal/file containing exactly two distance, angle, or dihedral axes |
| `scan3d` | One literal/file containing exactly three distance, angle, or dihedral axes |

A four-element tuple is therefore a distance range in `scan` and an angle target in `all`.

## Input file requirements

- **PDB** — must contain hydrogens (add via `reduce` / `pdb2pqr` / Open Babel) and element symbols in cols 77–78 (`pdb2reaction add-elem-info` if missing). Multiple PDBs must share identical atoms in the same order.
- **mmCIF** — see {ref}`mmCIF and large structures <mmcif-input>` below.
- **XYZ / GJF** — accepted when active-site extraction is skipped (omit `-c/--center`). `.gjf` files can provide charge / spin defaults from embedded metadata.

(mmcif-input)=
### mmCIF and large structures

Every calculation command that accepts PDB also accepts `.cif` and `.mmcif`. Use mmCIF for chain IDs longer than one character, residue numbers beyond four digits, atom serial numbers beyond five digits, or structures with 10,000 or more residues.

`pdb2reaction` reads the first coordinate model and keeps one alternate location (altLoc) per residue, the one with the highest mean occupancy. During the calculation the atoms carry temporary chain IDs and residue numbers; output CIF files restore the original chain IDs, residue numbers, and insertion codes. Large or non-standard PDB files are handled the same way, for example files with 10,000 or more residues, 99,999 or more atoms, hybrid-36 numbering, or numbers that overflow their columns.

Residue and atom selectors use the original chain IDs and residue numbers. For the `.cif` files written next to each output, see [Output Directory Layout](output-layout.md).

(trajectory-one-frame)=
### Extract one frame from a trajectory

A `_trj.xyz` file is a plain multi-frame XYZ file, so frame k (counted from 1) can be extracted with:

```bash
N=$(head -1 scan_trj.xyz); k=12
sed -n "$(( (k-1)*(N+2)+1 )),$(( k*(N+2) ))p" scan_trj.xyz > frame_12.xyz
```

To continue with the PDB topology, pass the original PDB to `--ref-pdb` of the next command; the coordinates come from the frame.

(exit-codes)=

## Exit codes

| Code | Meaning |
|---|---|
| `0` | Success or usable partial results |
| `1` | Non-convergence, no usable result, runtime exception, or output failure |
| `2` | Invalid input, CLI arguments, or configuration |
| `130` | User interruption (SIGINT) |

Exit codes do not depend on JSON output. Exit code `0` covers both `success` and `partial`; tell them apart by `scientific_status`. `all` and `path-search` write it to `summary.log`, and `all` also prints `Scientific status:` on the console. The other commands record it in `result.json` when run with `--out-json` (see [Execution and requested-stage completion](json-output.md#execution-and-requested-stage-completion)). Without `--out-json`, read the console lines that the command's page lists, for example [Judging the IRC](irc.md#judging-the-irc) for `irc`.

(opt-mode-semantics)=

## `--opt-mode` (subcommand-dependent)

`--opt-mode` picks the optimizer. L-BFGS (limited-memory BFGS) and RFO (rational function optimization) find minima; Dimer, RS-P-RFO (restricted-step partitioned RFO), RS-I-RFO (restricted-step image RFO), and TRIM (trust-region image minimization) search for a TS (transition state).

| Subcommand | `grad` alias selects | `hess` alias selects | Default |
|---|---|---|---|
| `opt` | L-BFGS (`lbfgs`) | RFO (`rfo`) | `grad` (L-BFGS) |
| `tsopt` | Dimer (`dimer`) | RS-P-RFO (`rsprfo`) | `hess` (RS-P-RFO) |
| `path-opt` (endpoint preopt) | L-BFGS | RFO | `grad` |
| `path-search` (single-structure optimization of HEI±1 and kink nodes; HEI = highest-energy image) | L-BFGS | RFO | `grad` |
| `scan` / `scan2d` / `scan3d` (per-grid relaxation) | L-BFGS | RFO | `grad` |
| `all` (pre-optimization, `--opt-mode`) | L-BFGS | RFO | `grad` |
| `all` (TS optimization, `--opt-mode-post`) | Dimer | RS-P-RFO | `hess` |
| `all` (endpoint optimization after IRC, `--opt-mode-post`) | L-BFGS | RFO | `hess` |

The same `--opt-mode` value selects a **different algorithm** on each subcommand, and the defaults differ, so check the table before copying a recipe. Algorithm names are accepted on `opt` (`lbfgs` / `rfo`) and `tsopt` (`dimer` / `rsirfo` / `trim` / `rsprfo`); all other subcommands accept only `grad` / `hess`. On `tsopt`, `--opt-mode grad` is therefore a **Dimer** TS search, not an L-BFGS minimization, and this Dimer periodically computes the Hessian to update its direction. Write `--opt-mode dimer` or `rsirfo` on `tsopt` and `--opt-mode lbfgs` or `rfo` on `opt` to make a recipe unambiguous.

## CLI ↔ YAML name mismatches

A few CLI flags use slightly different names than their YAML counterparts, and a few are renamed when wrapped in `all`. The main flags and their YAML keys are in {ref}`YAML Reference › Common CLI-to-YAML mapping <common-cli-to-yaml-mapping>`. The two most-asked cases:

(pressure-vs-pressure-atm)=
- **`--pressure` (CLI) vs `pressure_atm` (YAML)** — on `freq` the flag is `--pressure FLOAT`; in `all` it is exposed as `--freq-pressure`. YAML key: `thermo.pressure_atm`. Both carry **atm** values (converted to Pa internally).

- **`--step-size` (CLI) vs `step_length` (YAML)** — on `irc` the flag is `--step-size FLOAT` (bohr); in `all` it is `--irc-step-size`. YAML key: `irc.step_length`.

```bash
pdb2reaction irc -i ts.pdb -q 0 --step-size 0.05
pdb2reaction all -i r.pdb p.pdb -c SAM -l 'SAM:1' --tsopt --irc-step-size 0.05
```

## YAML configuration

```bash
pdb2reaction all -i r.pdb p.pdb -q -1 --config my_settings.yaml --out-dir result/
```

(configuration-precedence)=

```
built-in defaults  <  --config (YAML)  <  CLI options
```

`pdb2reaction <subcmd> --help-advanced` and the [Command Reference](reference/commands/index.md) show the built-in default of each option (`[default: …]`). Only *explicitly supplied* CLI values override YAML; options left at their CLI default do not mask YAML values. This order holds for every command that takes `--config`. Full schema: [YAML Reference](yaml-reference.md).

## Output directory

`-o/--out-dir ./my_results/` sets the output directory of a calculation command; each command has its own default, listed in [Output Directory Layout](output-layout.md). `extract` instead takes one or more file paths with `-o/--output` and writes to the current directory by default.

## Notes

* **Flags with a value**: a flag followed by a value (`true` / `false`) is also accepted; write paired flags in commands and scripts.
* **Put the chain before a residue name and number.** `TYR:44` is read as chain `TYR`, residue 44, and stops with a "not found" error; write `A:TYR:44`.
* **Use one kind of residue selector per list.** The forms in the table fall into three kinds: names only (`SAM`), chain + name with or without a number (`A:SAM`, `A:TYR:44`), and numbers with or without a chain (`A:123`, `123`). Kinds cannot be mixed in one list: `A:TYR:44,A:SAM` works, while `A:SAM,SAM`, `A:44,A:SAM`, and `SAM,TYR:44` stop with an error.
* **Atom selectors on PDB files with an empty chain column** take three fields, such as `'SER:11:HG'` or `'SER 11 HG'`. `_` does not stand for an empty chain, so `'_:SER:11:HG'` matches no atom and stops with an error.
* **mmCIF and large structures**: up to 619,938 residues can be handled (62 one-character chain IDs × 9,999 residue numbers in the internal PDB used during the calculation). `fix-altloc` and `add-elem-info` read PDB only. For mmCIF, the altLoc is chosen and element symbols are taken from `_atom_site.type_symbol` when the file is read. A row without it stops with an error.
* **`--flatten` and YAML**: `--flatten` (`opt`, `tsopt`, and `all`; off by default) displaces the structure along extra imaginary modes and optimizes again. When neither `--flatten` nor `--no-flatten` is given, `tsopt` keeps a `hessian_dimer.flatten_max_iter` value set in YAML, and `opt` always uses up to 50 rounds with `--flatten`; `--no-flatten` forces 0 (see {ref}`When --flatten is on <flatten-precedence-caveat>`).

## See Also

- [Installation](installation.md) — setup and dependencies
- [Getting Started](getting-started.md) — the shortest run and which page to read next
- [Output Directory Layout](output-layout.md) — file names and default output directories
- [Troubleshooting](troubleshooting.md) — common errors and fixes
- [YAML Reference](yaml-reference.md) — all configuration options
