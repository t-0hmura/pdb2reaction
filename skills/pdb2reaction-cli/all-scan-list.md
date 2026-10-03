# `pdb2reaction all`: Single structure + scan

Give one reactant and the coordinates to drive with `-s`; `all` runs the scan
stages in order, searches the MEP through the stage ends, and, with `--tsopt`,
optimizes each TS candidate and runs IRC. It succeeded when the console prints
`Scientific status: success` under the last `====== Pipeline summary ======`
(and `[tsopt] Converged (n_imag=1).` for each TS).

## When to use

You have only the reactant and can write the chemistry as a sequence of scans
of distances, angles, or dihedrals, for example "first move the methyl from S
of SAM to C7 of GPP, then move H11 of GPP onto OE2 of Glu186". The start and
the stage ends become the inputs of the MEP search: single-pass `path-opt` by
default, or the recursive `path-search` with `--refine-path`, which can add
intermediates it finds.

## Minimal run

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
that flag is one stage. Stages run in order; the final geometry of stage k is
the input geometry of stage k+1.

## Writing --scan-lists

Each literal is a list of target tuples: distance `(i, j, target_Å)`, angle
`(i, j, k, target_deg)`, or dihedral `(i, j, k, l, target_deg)`.

An atom is either a 1-based atom number or a selector in double quotes. A
three-field selector gives the atom name, residue name, and residue number in
any order, separated by spaces, commas, colons, slashes, backticks, or
backslashes. The four-field form
`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` must use this order; use it when repeated
residue numbers or names would otherwise be ambiguous. In a PDB with an empty
chain column (the bundled examples), use three fields; `_` does not mean an
empty chain.

| Form | Example |
|---|---|
| spaces | `"CS1 SAM 320"` |
| commas | `"SAM,320,CS1"` |
| slashes | `"SAM/320/CS1"` |
| chain, four fields | `"A:SAM:320:CS1"` |

Several tuples in one literal move together. If you want them done
sequentially, split them into separate literal values after the same
`--scan-lists` occurrence. Repeating the flag is rejected.

```bash
# One stage, two bonds driven together (concerted SN2):
--scan-lists '[("CS1 SAM 320","C7 GPP 321",1.60),("CS1 SAM 320","SD SAM 320",3.0)]'

# Two stages, one bond each (stepwise mechanism):
--scan-lists '[("CS1 SAM 320","C7 GPP 321",1.60)]' \
             '[("H11 GPP 321","OE2 GLU 186",0.90)]'
```

To choose the coordinates for your own reaction, see
[Staged vs concerted scans](../pdb2reaction-overview/ts-strategy.md#staged-vs-concerted-scans).

## Judge success

Read the console, `summary.json`, and the endpoints as in
[all.md](all.md#judge-success). For the scan itself:

- **Stages**: open `_work/scan/stage_NN/scan_trj.xyz` and check that the coordinates change as intended; `_work/scan/stage_NN/result.*` is the restrained end of each stage. `summary.json["scan"]` holds the scan status and `stages`; the full record is `_work/scan/result.json`.
- **MEP and R/TS/P**: the MEP of each segment is `_work/path_opt/mep_seg_NN_trj.xyz` (`_work/path_search/` with `--refine-path`). With `--tsopt`, `segments/seg_NN/{reactant,ts,product}.*` are the R/TS/P of each processed segment.

## Pitfalls and recovery

- **A stage reaches an unexpected geometry.** The restrained optimization relaxed into another basin, or the coordinate under-specifies the mechanism. Inspect the trajectory, revise or add a chemically meaningful coordinate, or split a complex stage; do not assume the side product is valid.
- **Python literal error.** Wrap each stage in single quotes and use double quotes inside. Prefer space, comma, or slash selectors; a backtick is safe only inside the outer single quotes.
- **Atom not found or matched twice.** Names must match the input exactly, case included; editing tools sometimes rename atoms (`CB` to `CB1`). If a three-field selector matches more than one atom, add the chain with `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`; do not guess from the first match.
- **XYZ or GJF input.** Residue selectors and `-c` extraction are unavailable; use atom numbers and give `-q` (or a valid GJF header).
- **Several `-i` inputs.** `-s` takes exactly one structure; with two or more, the run stops with an error.
- **More segments than expected** (`--refine-path` only). Bond-change splitting proposed another candidate intermediate; validate it and the neighbouring TS/IRC. The default `path-opt` adds no segments.
- **Cost.** Each stage adds a restrained optimization before the MEP; time one pilot stage and budget from it.

## Next step

- [all.md](all.md): mode choice, success criteria, resume, output tree.
- [scan.md](scan.md): `scan`, `scan2d`, and `scan3d` on their own.
- [path.md](path.md): the MEP search after the scans.
- Defaults: `python -c "import pdb2reaction.core.defaults as d; print(d.SEARCH_KW, d.STOPT_KW)"`.
