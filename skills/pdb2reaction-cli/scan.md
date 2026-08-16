# `pdb2reaction scan`

## Purpose

1D bond-length-driven scan with staged harmonic restraints and inter-stage
relaxation. Drive one or more bonds toward target distances and inspect the
trajectory or select frames as later path endpoints. Use
`pdb2reaction all --scan-lists` when the integrated downstream MEP/TS/IRC
pipeline is wanted; use standalone `scan` when the scan itself is the task.

## Synopsis

```bash
pdb2reaction scan -i input.pdb \
    -s '[(idx_a, idx_b, target_A), ...]' \
    [-l 'RES:Q,...'] [-q / -m] \
    [-b uma|orb|mace|aimnet2] [-o ./result_scan/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Reactant `.pdb` / `.cif` / `.mmcif` / `.xyz` / `.gjf` |
| `-s, --scan-lists` | str | required | Inline Python literal `'[(a,b,target),...]'`, or YAML/JSON spec path. **Pass multiple stages as space-separated literals after a single `-s`** — repeating `-s` is rejected. |
| `-q` / `-l` / `-m` | — | — | Charge / spin |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_scan/` | Output directory |
| `--ref-pdb` | path | none | Residue context for XYZ/GJF inputs |
| `--config` / `--dry-run` / `--help-advanced` | — | — | Standard (`scan` has no `--show-config`) |

The tuple grammar in `-s` accepts atom-index ints (`(1, 5, 1.4)`) or atom
specs (`("CS1 SAM 320", "C7 GPP 321", 1.60)`). When residue numbering or
names repeat, use `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`, for example
`("A:SAM:320:CS1", "B:GPP:321:C7", 1.60)`. Multiple stages chain
sequentially as space-separated literals after a **single** `-s`; each
stage starts from the previous stage's final geometry.

## Examples

### Single stage by atom name

```bash
pdb2reaction scan -i 1.R.pdb -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.60)]' \
    -b uma -o result_scan
```

### Two sequential stages

```bash
pdb2reaction scan -i 1.R.pdb -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.60)]' \
       '[("H11 GPP 321","OE2 GLU 186",0.90)]' \
    -b uma -o result_scan_staged
```

Each space-separated literal after a single `-s` is one stage; do **not** repeat `-s` (rejected with `repeated flags are not accepted`).

## Output

| Path | When | Content |
|---|---|---|
| `<out_dir>/result.json` | `--out-json` | machine-readable result |
| `<out_dir>/preopt/result.{xyz,pdb,gjf}` | `--preopt` | pre-optimized starting geometry |
| `<out_dir>/stage_NN/result.xyz` | stage reaches its output-writing step | final attempted geometry; check the stage `converged` value |
| `<out_dir>/stage_NN/scan_trj.xyz` | stage is attempted | per-stage scan trajectory; it may be empty when the target already equals the starting distance |
| `<out_dir>/stage_NN/scan_*.xyz` | `--dump` | intermediate optimizer steps; run-scoped YAML `opt.dump` is ignored |
| `<out_dir>/scan_trj.xyz` | at least one scan step produces a frame | stitched scan trajectory across all stages |
| `<out_dir>/scan.pdb` | stitched XYZ exists, `--convert-files`, and PDB/mmCIF topology/reference available | normalized PDB companion used between pipeline stages |
| `<out_dir>/scan.cif` | stitched XYZ exists, `--convert-files`, and input/reference required the mmCIF or oversized-PDB bridge | public trajectory with original IDs |

`result.json` uses top-level `status: "completed"` when the runner returned;
that is not proof every optimizer converged. Inspect each
`stages[i]["converged"]`, target distance, final energy, and trajectory. Plot
the stitched trajectory with `trj2fig.md` when it exists.

## Caveats

- `-s` is Python literal-eval. Quote with single quotes outside,
  double quotes inside. Atom-name strings use `"NAME RESNAME RESID"`
  with single spaces.
- A chain-qualified selector is positional:
  `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`. The legacy three-field form remains
  order-flexible.
- Stage *k+1* starts from stage *k*'s final geometry; a diverged
  stage poisons downstream stages.
- For coupled multi-bond drives in one stage, put multiple tuples in
  one `-s` argument: `'[(a,b,1.6),(c,d,3.0)]'`.

## See also

- `scan2d.md`, `scan3d.md` — higher-dim analogs.
- `all-scan-list.md` — wraps `scan` inside the full pipeline.
- Defaults: `import pdb2reaction.core.defaults as d; print(d.BIAS_KW, d.BOND_KW, d.OUT_DIR_SCAN)`
