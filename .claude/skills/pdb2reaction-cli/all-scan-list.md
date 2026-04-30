# `pdb2reaction all` — scan-list mode

## When to use

You have **only the reactant** (no product structure) and you can
articulate the chemistry as a sequence of staged distance scans —
e.g. "first push the methyl from S of SAM to C7 of GPP, then snap H11
to OE2 of GLU 186". `pdb2reaction all` runs each stage in order, ties
the resulting trajectories into an MEP, and the recursive bond-change
segmentation slots in any intermediates it finds.

## Synopsis

```bash
pdb2reaction all -i 1.R.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --scan-lists \
        '[("CS1 SAM 320","C7 GPP 321",1.60)]' \
        '[("H11 GPP 321","OE2 GLU 186",0.90)]' \
    --tsopt --thermo \
    -o result_scan
```

Each `--scan-lists '...'` argument is **one stage**. Stages run
sequentially; the final geometry of stage *k* is the input geometry of
stage *k+1*.

## `--scan-lists` syntax

Each argument is a Python literal-eval expression: a list of bond
tuples, where each tuple is `(atom_a, atom_b, target_distance_Å)`.

```
[ ("<atom-spec>", "<atom-spec>", <float>) , ... ]
```

`<atom-spec>` is a string of three tokens (atom name, residue name, residue index) in **any order**, separated by whitespace, comma, slash, backtick, or backslash. Common conventions:

| Form | Example |
|---|---|
| `"NAME RESNAME RESID"` (whitespace) | `"CS1 SAM 320"` |
| `"NAME,RESNAME,RESID"` (comma) | `"CS1,SAM,320"` |
| `"RESNAME/RESID/NAME"` (slash) | `"SAM/320/CS1"` |

The parser (`utils.resolve_atom_spec_index`) auto-detects the role of each token by type (integer = resid, known residue name = resname, otherwise atom name); chain IDs are not part of the spec.

Multiple bonds per stage drive simultaneously. If you want them done
**sequentially**, split them into separate `--scan-lists` arguments.

Examples:

```bash
# One stage, two bonds driven together (concerted SN2):
--scan-lists '[("CS1 SAM 320","C7 GPP 321",1.60),("C7 GPP 321","S SAM 320",3.0)]'

# Two stages, one bond each (stepwise mechanism):
--scan-lists '[("CS1 SAM 320","C7 GPP 321",1.60)]' \
             '[("H11 GPP 321","OE2 GLU 186",0.90)]'
```

## Mode-specific flags

| Flag | Default | Meaning |
|---|---|---|
| `--scan-lists` | required | One or more stages of distance-restraint scans |
| `--mep-mode` | `gsm` | After scans complete, MEP refinement uses GSM unless `dmf` |

Unlike endpoint-MEP mode, `-i` is **a single PDB** (the reactant). The
toolkit synthesizes intermediate / product geometries from the scan
trajectories.

## Output

Same overall tree as in `all.md`, plus per-stage scan output:

```
result_scan/
├── scan/                           # at out_dir top level (NOT under path_search/)
│   ├── stage_01/  scan_*.xyz       # raw distance-restraint scan trajectory
│   ├── stage_02/  scan_*.xyz
│   └── ...
├── path_search/                    # (or path_opt/ when --refine-path False)
│   ├── mep_seg_NN.{pdb,xyz}        # per-segment MEP strings
│   └── post_seg_NN/                # per-segment refinements + energy diagrams
├── seg_NN/                         # canonical R/TS/P per segment (top-level)
└── summary.json
```

`summary.json["scan"]` carries the stage-by-stage record:

```python
import json
d = json.load(open("result_scan/summary.json"))
for stage in d["path_search"]["scan"]["stages"]:
    print(stage["index"], stage["bonds"], stage["status"], stage["final_distance_A"])
```

## Distinctive failure modes

| Symptom | Cause | Fix |
|---|---|---|
| Stage k goes to a different geometry than expected | Distance restraint not strong enough; SCF found a side product | Tighten the target distance, or split a complex stage into two simpler ones |
| `--scan-lists` triggers a Python literal-eval error | Quoting mistake | Wrap each stage in single quotes outside, double quotes inside; backticks survive bash without escaping |
| Path search reports more segments than expected | Bond-change detector found a "free" intermediate | This is usually correct; check the IM geometry in `seg_01/product.pdb` (= `seg_03/reactant.pdb`) |

## Caveats

- The atom specs must match the **exact** atom names in the input PDB
  (case sensitive). PyMOL/Maestro sometimes rename `CB` ↔ `CB1`.
- `--scan-lists` is incompatible with multiple `-i` inputs (the latter
  triggers `all-endpoint-mep.md`).
- Each stage can take longer than path-search itself; budget walltime
  accordingly.

## See also

- `all.md` — base orientation.
- `scan.md`, `scan2d.md`, `scan3d.md` — standalone distance scan
  subcommands (without the surrounding pipeline).
- `path-search.md` — what happens after all scans complete.
- Defaults: `import pdb2reaction.defaults as d; print(d.SEARCH_KW, d.STOPT_KW)`.