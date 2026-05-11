# Quickstart: `pdb2reaction all`

## Goal

Run the end-to-end workflow once from two full PDB structures.

## Prerequisites

- pdb2reaction installed (see [Installation](installation.md))
- Two PDB files (reactant R and product P) with **hydrogen atoms** already added
- The same atoms in the same order across all input PDB files

> **About the example filenames:** `1.R.pdb` and `3.P.pdb` mirror the numbered reactant/product files shipped in the geranyl pyrophosphate (GPP) C6-methyltransferase BezA example directory ([`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples) — `1.R.pdb` = reactant state, `3.P.pdb` = product state, with intermediate `2.*.pdb` files available for runs that include additional reactant/product/intermediate structures). Replace them with the two (or more) full-system PDBs for your own reaction.

## Minimal command

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --out-dir ./result_all
```

### (Optional) Add post-processing in the same run

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft --out-dir ./result_all
```

> **VRAM warning:** `--dft` launches GPU4PySCF single-point jobs on the extracted cluster model and can easily OOM on GPUs with < 24 GB VRAM for clusters above ~200 atoms. If you hit `CUDA out of memory`, either drop `--dft` and run `pdb2reaction dft` separately with a smaller basis / trimmed cluster, or move the DFT step to a larger-VRAM node. The `[dft]` extra must also be installed (see [Installation](installation.md) Step 7).

## Expected output

A successful run produces a directory like:

```text
result_all/
├── summary.log                    # Human-readable summary
├── summary.json                   # Machine-readable results
├── path_search/
│   ├── mep.pdb                    # Merged MEP trajectory
│   ├── energy_diagram_UMA_all.png # Energy profile
│   ├── summary.json               # Path-search results
│   └── post_seg_01/               # Post-processing (if --tsopt)
│       ├── ts/final_geometry.pdb
│       ├── irc/finished_irc_trj.xyz
│       └── freq/
└── seg_01/                        # IRC-optimized R/TS/P structures
    ├── reactant.pdb
    ├── ts.pdb
    └── product.pdb
```

**What to check:**

1. `summary.json` — check the `status` field (`"success"`, `"partial"`, or `"failed"`) and the per-segment `barrier_kcal` values; `summary.log` mirrors the same information in human-readable form
2. `seg_01/*.pdb` — open in PyMOL to verify the R/TS/P structures make chemical sense
3. `energy_diagram_*.png` — the energy profile should show a clear barrier

**Sample terminal output (successful run):**

```
[all] Elapsed for Whole Pipeline: HH:MM:SS.sss
```

(Wall-clock varies with system size, GPU, and selected stages.)

If `--tsopt` is enabled, you should also see:

```
[Imaginary modes] n=1  ([-425.9])
```

A first-order saddle point shows exactly one imaginary mode along the reaction coordinate; IRC validation (run automatically as part of `--tsopt`) confirms it connects the expected reactant and product.

## Tips

- `pdb2reaction all --help` shows core options; `pdb2reaction all --help-advanced` shows the full list.

## Next step

- Single-structure staged scan route: [Quickstart: `pdb2reaction all --scan-lists`](quickstart-scan.md)
- TS candidate validation: [Quickstart: TS-only mode](quickstart-tsopt-freq.md)
- Full option reference: [all](all.md)
