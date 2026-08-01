# `pdb2reaction path-search`

## Purpose

Recursive minimum-energy-path (MEP) search across two or more
endpoints. Detects bond changes along the candidate MEP and
**recursively re-segments** the path using bond-change evidence to isolate
candidate elementary steps. It does not prove that each segment has exactly
one transition state; TSOPT, an independent frequency calculation, and IRC
provide that validation. Output includes per-segment HEI candidates and a
stitched MEP.

This is the engine behind `pdb2reaction all` in endpoint-MEP mode when
`--refine-path` is set; by default `all` uses single-pass `path-opt`
for its MEP stage.

## Synopsis

```bash
pdb2reaction path-search -i 1.R.pdb 3.P.pdb \
    [--mep-mode gsm|dmf] [--refine-mode peak|minima] \
    [--max-nodes 20] [-l 'RES:Q,...'] [-b uma|orb|mace|aimnet2] \
    [-o ./result_path_search/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path(s) | required (≥ 2) | Two or more endpoints in reaction order, with identical atom identities/elements/order |
| `--mep-mode` | str | `gsm` | `gsm` (Growing String) or `dmf` (Direct Max Flux) |
| `--refine-mode` | str | engine-dependent | `peak` (HEI±1; GSM default) or `minima` (nearest local minima; DMF default) |
| `--max-nodes` | int | 20 | Max internal nodes per segment string |
| `--thresh` | str | `gau` | Single-structure optimization convergence preset |
| `--thresh-gsm` | str | `gau_loose` | GSM string-optimizer convergence preset |
| `--thresh-dmf` | str/float | `tight` | DMF IPOPT dual-infeasibility tolerance: `tight`, `middle`, `loose`, or a positive float |
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

### Optional DMF mode

```bash
pdb2reaction path-search -i 1.R.pdb 3.P.pdb \
    --mep-mode dmf --refine-mode minima \
    -l 'SAM:1,GPP:-3' -b uma -o result_path_search
```

## Output

| Path | When | Content |
|---|---|---|
| `<out_dir>/summary.json`, `summary.log` | successful/partial path construction | machine + human result; handled runtime failures write the standard error envelope, while early CLI/input validation can fail before output exists |
| `<out_dir>/mep_trj.xyz` | successful/partial path construction | stitched candidate MEP across all path segments |
| `<out_dir>/mep.pdb` | `--convert-files` and PDB/mmCIF reference available | normalized PDB companion used between pipeline stages |
| `<out_dir>/mep.cif` | `--convert-files` and input/reference required the mmCIF or oversized-PDB bridge | public trajectory with original IDs |
| `<out_dir>/mep.gjf` | `--convert-files` and input is `.gjf` | GJF companion |
| `<out_dir>/mep_w_ref.{pdb,cif}`, `mep_w_ref_seg_NN.{pdb,cif}` | `--ref-full-pdb` merge requested and succeeds; CIF needs a bridged reference | active-site MEP merged into the full-system reference template(s) |
| `<out_dir>/mep_plot.png` | plotting succeeds | MEP energy plot |
| `<out_dir>/seg_NNN_<tag>/...` | corresponding recursive string is attempted | per-string scratch (3-digit; tags include `_mep`, `_maxdepth`, and `_bridge`; individual artifacts depend on engine progress and conversion inputs) |
| `<out_dir>/mep_seg_NN_trj.xyz`, `mep_seg_NN.{pdb,cif,gjf}` | bond-change segments; companions require `--convert-files` and a reference topology/template | canonical per-segment MEP frames |
| `<out_dir>/hei_seg_NN.{xyz,pdb,cif,gjf}` | bond-change segments; companions require `--convert-files` and a reference/template | canonical HEI per segment (TS seed) |
| `<out_dir>/energy_diagram_MEP.png` | diagram construction and PNG export succeed | bare MEP energies; PNG export requires Kaleido |

Standalone `path-search` does **not** run post-processing —
tsopt / irc / freq / dft are `all`'s job (it writes them under
`segments/seg_NN/`).

`summary.json["segments"]` lists path-segmentation candidates in MEP order,
including non-reactive bridge/kink segments. A `kind: "seg"` entry with bond
changes is a TS candidate, not proof of an elementary step; validate it with
TSOPT, an independent frequency calculation, and IRC. A representative entry is:

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

Top-level `summary.json["status"]` is `"success"` when energy-diagram
metadata was constructed, or `"partial"` when it was not. `success` does not
guarantee that PNG export succeeded or that the MEP converged; inspect the
expected files and segment diagnostics separately.

## Caveats

- All `-i` inputs must have identical atom identities, element sequence, and
  ordering; equal counts alone are insufficient.
- Recursive segmentation can produce **more** segments than `len(-i) - 1`
  by proposing intermediates the user did not supply. Validate every proposed
  intermediate and adjacent saddle; segmentation is geometry evidence, not a
  mechanism verdict.
- Raising `--max-nodes` increases calculator calls and may help an
  under-resolved path, but it does not repair inconsistent endpoints or a
  wrong mechanism. Inspect and pre-relax endpoints before increasing it.
- `--mep-mode dmf` cannot be combined with `--solvent`, because the DMF ASE
  path has no xTB correction while the rest of the run would use the corrected
  surface. Use GSM or omit the solvent correction; pdb2reaction rejects the
  mixed-PES combination.
- Output **does not** include refined TSs; those are
  `segments/seg_NN/ts/` produced by the `all` pipeline.

## See also

- `path-opt.md` — single-segment MEP optimization (the building block).
- `tsopt.md` — runs after path-search on each `hei_seg_NN.xyz` (TS seed).
- `bond-summary.md` — same bond-change algorithm used here, standalone.
- Defaults: `import pdb2reaction.core.defaults as d; print(d.SEARCH_KW, d.GS_KW, d.DMF_KW)`
