# `pdb2reaction path-search`

## Purpose

Recursive minimum-energy-path (MEP) search across two or more
endpoints. Detects bond changes along the candidate MEP and
**recursively re-segments** the path until each segment crosses
exactly one transition state. Output: one `seg_NN/` per elementary
step, plus a stitched `mep.pdb` and energy diagrams.

This is the engine behind `pdb2reaction all` in endpoint-MEP mode when
`--refine-path True` is set; by default `all` uses single-pass `path-opt`
for its MEP stage instead.

## Synopsis

```bash
pdb2reaction path-search -i 1.R.pdb 3.P.pdb [-i 1.R.pdb 2.IM.pdb 3.P.pdb] \
    [--mep-mode gsm|dmf] [--refine-mode peak|minima] \
    [--max-nodes 20] [-l 'RES:Q,...'] [-b uma|orb|mace|aimnet2] \
    [-o ./result_path_search/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path(s) | required (≥ 2) | Two or more endpoints in reaction order |
| `--mep-mode` | str | `gsm` | `gsm` (Growing String) or `dmf` (Direct Max Flux) |
| `--refine-mode` | str | mode-dep | `peak` (HEI±1) or `minima` (nearest local minima) |
| `--max-nodes` | int | 20 | Max internal nodes per segment string |
| `-q, --charge` / `-l` / `-m` | — | — | Charge / multiplicity (see common conventions) |
| `-b, --backend` | str | `uma` | MLIP backend |
| `--solvent` | str | none | xTB-ALPB solvent name |
| `-o, --out-dir` | path | `./result_path_search/` | Output directory |
| `--config` / `--show-config` / `--dry-run` | — | — | YAML config + preview |

## Examples

### 2-endpoint MEP, GSM, default refinement

```bash
pdb2reaction path-search -i 1.R.pdb 3.P.pdb \
    -l 'SAM:1,GPP:-3' -b uma \
    -o result_path_search
```

### 3-endpoint with explicit intermediate

```bash
pdb2reaction path-search -i 1.R.pdb 2.IM.pdb 3.P.pdb \
    -l 'SAM:1,GPP:-3' -b uma --max-nodes 30 \
    -o result_path_search
```

### DMF mode (sometimes better for ill-conditioned strings)

```bash
pdb2reaction path-search -i 1.R.pdb 3.P.pdb \
    --mep-mode dmf --refine-mode minima \
    -l 'SAM:1,GPP:-3' -b uma -o result_path_search
```

## Output

| Path | When | Content |
|---|---|---|
| `<out_dir>/summary.json`, `summary.log` | always | machine + human result |
| `<out_dir>/mep_trj.xyz` | always | stitched MEP across all segments |
| `<out_dir>/mep.pdb` | reference PDB available | PDB companion |
| `<out_dir>/mep.gjf` | input is `.gjf` | GJF companion |
| `<out_dir>/mep_w_ref.pdb`, `mep_w_ref_seg_NN.pdb` | reference PDB available | full / per-segment MEP merged with reference |
| `<out_dir>/mep_plot.png` | always | MEP energy plot |
| `<out_dir>/seg_NNN_<tag>/{final_geometries_trj.xyz, mep_plot.png, hei.{xyz,pdb,gjf}}` | always | per-string scratch (3-digit, `_mep` / `_maxdepth`) |
| `<out_dir>/mep_seg_NN_trj.xyz`, `mep_seg_NN.{pdb,gjf}` | always | canonical per-segment MEP frames |
| `<out_dir>/hei_seg_NN.{xyz,pdb,gjf}` | always | canonical HEI per segment (TS seed) |
| `<out_dir>/energy_diagram_MEP.png` | always | bare MEP energies |

Standalone `path-search` does **not** run post-processing —
tsopt / irc / freq / dft are `all`'s job (it writes them under
`segments/seg_NN/`).

`summary.json["segments"]` lists each elementary step with:

```python
{
  "index": 1,
  "tag": "seg_01",
  "kind": "seg",
  "barrier_kcal": 21.5,
  "delta_kcal": -0.7,
  "bond_changes": [
    {"Bond formed (1)": ["C508-C567 : 3.166 Å --> 1.675 Å"]},
    {"Bond broken (1)": ["S507-C508 : 1.798 Å --> 3.459 Å"]}
  ]
}
```

Top-level `summary.json["status"]` is `"success"` (energy diagram
written) or `"partial"` (no diagram — usually means MEP did not
complete cleanly).

## Caveats

- All `-i` inputs must have identical atom counts and ordering.
- Recursive segmentation can produce **more** segments than `len(-i) - 1`
  — that's the whole point: it finds intermediates the user didn't
  supply.
- `--max-nodes` bigger than 30 rarely helps; if a segment doesn't
  converge with 20 nodes, the chemistry is usually the problem.
- Output **does not** include refined TSs; those are
  `segments/seg_NN/ts/` produced by the `all` pipeline.

## See also

- `path-opt.md` — single-segment MEP optimization (the building block).
- `tsopt.md` — runs after path-search on each `hei_seg_NN.xyz` (TS seed).
- `bond-summary.md` — same bond-change algorithm used here, standalone.
- Defaults: `import pdb2reaction.core.defaults as d; print(d.SEARCH_KW, d.GS_KW, d.DMF_KW)`