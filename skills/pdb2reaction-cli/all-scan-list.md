# `pdb2reaction all` — scan-list mode

## When to use

You have **only the reactant** (no product structure) and you can
articulate the chemistry as a sequence of staged distance scans —
e.g. "first push the methyl from S of SAM to C7 of GPP, then snap H11
to OE2 of GLU 186". `pdb2reaction all` runs each stage in order and ties
the resulting trajectories into an MEP. By default the MEP stage is
single-pass `path-opt`; pass `--refine-path` to run the recursive
bond-change segmentation that slots in any intermediates it finds.

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

Use exactly one `--scan-lists` flag. Each space-separated literal following
that flag is **one stage**. Stages run
sequentially; the final geometry of stage *k* is the input geometry of
stage *k+1*.

## `--scan-lists` syntax

Each argument is a Python literal-eval expression: a list of bond
tuples, where each tuple is `(atom_a, atom_b, target_distance_Å)`.

```
[ ("<atom-spec>", "<atom-spec>", <float>) , ... ]
```

`<atom-spec>` is either three tokens (atom name, residue name, residue index)
in **any order**, or the positional four-field form
`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`. The latter is required when repeated residue
numbers or names would otherwise be ambiguous. Three-field tokens may be
separated by whitespace, comma, slash, backtick, or backslash.

| Form | Example |
|---|---|
| `"NAME RESNAME RESID"` (whitespace) | `"CS1 SAM 320"` |
| `"NAME,RESNAME,RESID"` (comma) | `"CS1,SAM,320"` |
| `"RESNAME/RESID/NAME"` (slash) | `"SAM/320/CS1"` |
| `"CHAIN:RESNAME:RESSEQ[ICODE]:ATOM"` | `"A:SAM:320:CS1"` |

The parser (`utils.resolve_atom_spec_index`) auto-detects roles in the
three-field form. The four-field chain-qualified form is positional and must
use the order shown above; this keeps numeric or repeated chain IDs
unambiguous.

Multiple bonds in one stage are driven simultaneously. If you want them done
**sequentially**, split them into separate literal values after the same
`--scan-lists` occurrence. Repeating the flag is rejected.

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
| `--scan-lists` | required once | One or more following literal values, one per distance-restraint stage |
| `--mep-mode` | `gsm` | After scans complete, MEP refinement uses GSM unless `dmf` |

Unlike endpoint-MEP mode, `-i` is **one reactant structure** in
PDB/mmCIF/XYZ/GJF format. The toolkit synthesizes intermediate/product
geometries from the scan trajectories. XYZ/GJF inputs cannot use residue
selectors or `-c` extraction; use numeric atom selectors and provide a verified
explicit charge (or a valid GJF header).

## Output

Same overall layout as `all.md`, plus per-stage scan output:

| Path | When | Content |
|---|---|---|
| `<out_dir>/_work/scan/stage_NN/{scan_trj.xyz,result.{xyz,pdb,cif,gjf}}` | corresponding scan stage reaches output; companions depend on topology/template | raw distance-restraint scan trajectory + per-stage final geometry (scratch) |
| `<out_dir>/_work/path_opt/mep_seg_NN_trj.xyz`, `mep_seg_NN.{pdb,cif,gjf}` | corresponding MEP segment succeeds (`_work/path_search/` with `--refine-path`); companions depend on topology/template | per-segment MEP strings (scratch) |
| `<out_dir>/segments/seg_NN/` | a candidate segment enters post-processing | per-segment deliverables; may be partial after failure |
| `<out_dir>/segments/seg_NN/{reactant,ts,product}.{pdb,cif,xyz,gjf}` | successful `--tsopt` + IRC/endpoint processing; companions depend on topology/template | canonical R/TS/P per processed segment; CIF restores bridge-input IDs |
| `<out_dir>/summary.json` | pipeline summary/error handling reaches output | machine-readable result; check top-level and per-stage status |

For per-stage scan diagnostics (target distances, convergence, energies),
run `pdb2reaction scan` standalone with `--out-json` and parse
`result_scan/result.json` `["stages"]` (see [`scan.md`](scan.md)).
`pdb2reaction all` does **not** propagate the per-stage scan record into
its `summary.json`.

## Distinctive failure modes

| Symptom | Cause | Fix |
|---|---|---|
| Stage k goes to a different geometry than expected | The restrained MLIP optimization relaxed into another basin, or the chosen coordinate under-specifies the mechanism | Inspect the trajectory, revise/add a chemically meaningful coordinate, or split a complex stage into simpler stages; do not merely assume the side product is valid |
| `--scan-lists` triggers a Python literal-eval error | Quoting mistake | Wrap each stage in outer single quotes and use double quotes inside. Prefer whitespace/comma/slash atom selectors; a backtick is safe only while it remains inside those outer single quotes. |
| Path search reports more segments than expected (`--refine-path` only) | Bond-change segmentation proposed an additional candidate intermediate | Inspect and validate the IM and adjacent TS/IRC results; extra segmentation is not proof that the intermediate is chemically real. The default single-pass `path-opt` does not add segments. |

## Caveats

- The atom specs must match the **exact** atom names in the input PDB/mmCIF
  (case sensitive). PyMOL/Maestro sometimes rename `CB` ↔ `CB1`.
- If a three-field spec matches more than one atom, add the auth chain ID as
  `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`; do not guess from the first match.
- `--scan-lists` is incompatible with multiple `-i` inputs (the latter
  triggers `all-endpoint-mep.md`).
- Every stage adds a restrained optimization before the MEP calculation;
  benchmark a pilot stage and budget from measured timings rather than assuming
  a fixed scan-to-path-search cost ratio.

## See also

- `all.md` — base orientation.
- `scan.md`, `scan2d.md`, `scan3d.md` — standalone distance scan
  subcommands (without the surrounding pipeline).
- `path-search.md` — what happens after all scans complete.
- Defaults: `import pdb2reaction.core.defaults as d; print(d.SEARCH_KW, d.STOPT_KW)`.
