# `path-search` (recursive MEP through two or more structures)

## Overview

`path-search` builds one continuous minimum-energy path (MEP) through **two or more** structures given in reaction order (R → … → P). It refines the path recursively, only in the regions where covalent bonds change, and builds each piece with GSM (growing string method, the default) or DMF (direct max flux).

### What it is for

* **Splitting R → P into reactive segments**: when you do not know whether the reaction has one step or several, find the regions where bonds change.
* **A multistep path through intermediates**: give known intermediates between R and P and get one stitched path.
* **TS candidates per segment**: each reactive segment gets its own HEI (highest-energy image), `hei_seg_NN.xyz`, to optimize with [`tsopt`](tsopt.md).

For exactly two endpoints without recursive refinement, [`path-opt`](path-opt.md) is simpler.

---

## Examples

### 1. Two endpoints

Give the reactant and the product after one `-i`, with the charge and the spin multiplicity.

```bash
pdb2reaction path-search -i reactant.pdb product.pdb -q 0 -m 1 --out-dir ./result_path_search
```

When the run finishes, open the `[2] Segment-level MEP summary` section of `summary.log`, or read `summary.json`. There, `scientific_status` is `success` when the pre-optimizations and every path run converged, otherwise `partial` or `failed`. `segments` lists each segment's `index`, `tag`, `kind`, `converged`, and `barrier_kcal`; `kind` is `seg` (reactive), `kink` (conformation only), or `bridge` (short connecting path).

### 2. Add intermediates for a multistep path

List the structures in reaction order after one `-i`; each adjacent pair is searched and the pieces are stitched into one path.

```bash
pdb2reaction path-search -i R.pdb IM1.pdb IM2.pdb P.pdb -q -1 -m 1 \
  --out-dir ./result_path_search_multi
```

### 3. DMF with minima refinement

Build the paths with DMF and refine around each HEI from the nearest local minima.

```bash
pdb2reaction path-search -i reactant.pdb product.pdb -q 0 -m 1 \
  --mep-mode dmf --refine-mode minima --out-dir ./result_path_search_dmf
```

---

## How it works

Before the search, each input is pre-optimized (`--preopt`) and aligned to the one before it (`--align`), with frozen atoms matched step by step while the other atoms relax.

1. **A coarse MEP for each pair**:
Between each pair of adjacent inputs (A → B), GSM or DMF builds a coarse MEP and finds its HEI.
2. **Relaxing around the HEI**:
`--refine-mode peak` optimizes the images on either side of the HEI (HEI ± 1); `minima` searches outward from the HEI for the nearest local minimum on each side. The result is two nearby minima, End1 and End2. When `--refine-mode` is omitted, GSM uses `peak` and DMF uses `minima`.
3. **Kink or reactive segment**:
If no covalent bond changes between End1 and End2, the region is a *kink*: `path-search` inserts a few linear nodes and optimizes each one. Otherwise the region is a *reactive segment*, and a new GSM or DMF path between End1 and End2 sharpens its barrier.
4. **Recursing where bonds still change**:
The parts A → End1 and End2 → B are checked for bond changes, and only parts that still have them are searched again, down to `--max-depth` levels.
5. **Stitching**:
The pieces are joined into one path. Duplicate endpoints are dropped; where the ends of two neighboring pieces still differ in bonds, that gap is searched as a new segment, and any other gap is filled with a short connecting path.

Bond changes are judged with the thresholds in the YAML `bond` section, by the same rules as in {ref}`scan <section-bond>`.

---

## Reading the segments

| What you see | Meaning | Next step |
| --- | --- | --- |
| A segment with bond changes, with its `hei_seg_NN.xyz` | A TS candidate for that step | Optimize it with [`tsopt`](tsopt.md), check for one imaginary mode, then run [`irc`](irc.md) |
| A segment whose `tag` is `seg_NNN_maxdepth` | Splitting stopped there, at the depth limit or after repeated kinks | It may hold more than one step; check it as above, raise `--max-depth`, or give intermediates |
| Only `kink` segments, or the warning `HEI is at an endpoint` | No bond change was found, or the path has no peak between its ends | Check the inputs, or give intermediates (example 2) |

A successful TS optimization gives one imaginary mode along the reaction coordinate; confirm every HEI with `tsopt` (n_imag = 1) and IRC before you read it as a step of the mechanism.

---

## Output files

`path-search` writes these files to `--out-dir`:

```text
result_path_search/
├─ mep_trj.xyz               # The whole stitched MEP, energies on the comment lines
├─ mep_trj.pdb               # Same path as PDB (PDB/mmCIF input or --ref-pdb)
├─ mep_plot.png              # ΔE profile along the path (kcal/mol, relative to the reactant)
├─ energy_diagram_MEP.png    # State-energy diagram of the MEP (relative to the reactant)
├─ summary.json              # Barrier and classification summary for every segment
├─ summary.log               # The same summary as text
├─ mep_seg_NN_trj.xyz        # Path of reactive segment NN
├─ hei_seg_NN.xyz            # HEI of reactive segment NN (TS candidate)
├─ hei_mode_seg_NN.*         # Reaction-direction guess at that HEI; all passes it to tsopt --ref-mode (a standalone tsopt does not need it)
├─ mep_w_ref*.pdb, hei_w_ref_seg_NN.pdb  # Paths and HEIs placed in the full system (--write-ref-merge)
├─ align_refine/             # Alignment and relaxation files of the inputs (--align)
└─ seg_NNN_*/                # Working files of each GSM/DMF run
```

`summary.json` differs from the `result.json` of the other commands; see {ref}`summary.json for path-search and all <summary-json-path-search-all>`. Only segments with bond changes get `mep_seg_NN_*` and `hei_seg_NN.*` files. NN is the segment's `index` in `summary.json` (counted from 01 along the final path), while NNN in a `seg_NNN` tag or directory counts the GSM/DMF runs from 000, so the two numbers differ. The working files of a reactive segment are in `<tag>_mep/`, named after the segment's `tag` in `summary.json` (for example `seg_000_refine_mep/`).

For PDB, mmCIF, or `.gjf` input, the outputs are also written in that format under the same name (the whole path as `mep.gjf`); {ref}`mmCIF input <mmcif-input>`, and PDB input too large for the PDB columns, also get `.cif` files that keep the original identifiers.

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | paths | (required) | Two or more structures in reaction order, after one `-i` (`-i` may also be repeated for each file) |
| `-q, --charge` | integer | `None` | Total charge. Required unless `-l` is given or the input is `.gjf` |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1); a `.gjf` input supplies its own |
| `-l, --ligand-charge` | text | `None` | Total ligand charge (for example `-1`) or a charge per residue name (for example `'GPP:-3,SAM:1'`), used when `-q` is omitted (PDB/mmCIF input only) |
| `-b, --backend` | text | `uma` | Backend (`uma`, `orb`, `mace`, `aimnet2`, `dft`) |
| `-o, --out-dir` | path | `./result_path_search/` | Output directory |
| `--mep-mode` | `gsm` / `dmf` | `gsm` | Path method: growing string method / direct max flux |
| `--dmf-backend` | `gpu` / `cpu` | `gpu` | DMF compute backend (`--mep-mode dmf` only): PyTorch on CUDA / NumPy |
| `--refine-mode` | `peak` / `minima` | `peak` for GSM, `minima` for DMF | How the region around each HEI is relaxed: HEI ± 1 / nearest local minima |
| `--max-depth` | integer | `10` | Maximum levels of recursive subdivision; `0` turns subdivision off |
| `--max-nodes` | integer | `20` | Movable images per segment; a segment has `max_nodes + 2` images |
| `--preopt/--no-preopt` | flag | `True` | Pre-optimize each input before the search |
| `--align/--no-align` | flag | `True` | Align each input to the one before it before the search |
| `--write-ref-merge/--no-write-ref-merge` | flag | `False` | Place the paths and HEIs in a full-system template (`mep_w_ref*`, `hei_w_ref*`); needs `--align` and `--ref-full-pdb` |
| `--ref-full-pdb` | path | `None` | Full-system PDB/mmCIF template for `--write-ref-merge`; use the one that matches the first input |
| `--ref-pdb` | paths | `None` | Cluster-model PDB/mmCIF files for `.xyz` / `.gjf` inputs, one per input in the same order; used to write the PDB outputs and for `--write-ref-merge`, not for `-l` or `--freeze-links` |
| `--freeze-links/--no-freeze-links` | flag | `True` | Freeze the parent atoms of cap hydrogens (PDB/mmCIF input only) |
| `--climb/--no-climb` | flag | `True` | Run the GSM climbing-image search on the reactive segments; connecting paths never climb |

See the [generated CLI reference](reference/commands/path_search.md) for every option.

> **Note:** In YAML (`--config`), `search.max_depth` sets the depth limit when `--max-depth` is not given, `search.kink_max_nodes` (default `3`) sets the number of nodes inserted in a kink, and `bond.bond_factor` (default `1.20`) scales the covalent radii used to decide whether a bond has changed.

---

## Notes

* **Inputs**: give at least two structures, all with the same atoms in the same order; fewer than two stops with an error.
* **`--write-ref-merge` without its partners**: without `--align` or without `--ref-full-pdb`, `path-search` prints a warning and writes no `*_w_ref*` files.
* **`--ref-pdb` count**: when `--write-ref-merge` is in effect, the number of `--ref-pdb` files must match the number of inputs, or the run stops with an error.
* **Inputs are protected**: if a fixed output name (`mep_trj.*`, `mep_plot.png`, `energy_diagram_MEP.png`, `summary.json`, `summary.log`) would replace an input file, `path-search` stops before writing anything.
* **DMF and YAML limits**: the DMF and YAML notes in [path-opt → Notes](path-opt.md#notes) also apply here, when you use `--mep-mode dmf` or `--config`.
* **Segment boundaries**: the segmentation is a guide based on bond-distance criteria; one segment is not guaranteed to be one elementary step or to contain exactly one TS.

---

## See also

* [path-opt](path-opt.md) — single-pass MEP between two structures
* [scan](scan.md) — drive a bond step by step to make a path or a TS candidate
* [tsopt](tsopt.md) — optimize each segment HEI into a TS
* [extract](extract.md) — make the cluster-model PDBs used as inputs
* [all](all.md) — the full workflow; `all --refine-path` runs `path-search` for its MEP step
* [YAML Reference](yaml-reference.md) — every `search`, `bond`, `gs`, and `dmf` setting
* [Glossary](glossary.md) — MEP, GSM, DMF, HEI, kink, and other terms
* [Troubleshooting](troubleshooting.md) — when a run fails
* {ref}`Exit codes <exit-codes>` — what each exit status means
