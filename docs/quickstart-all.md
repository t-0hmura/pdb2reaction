# Quickstart: `pdb2reaction all` (Endpoint mode)

## Goal

Run the end-to-end workflow once from two full PDB structures.

## Prerequisites

- pdb2reaction installed (see [Installation](installation.md))
- Two PDB/mmCIF files (reactant R and product P) with **hydrogen atoms** already added
- The same atom identities in the same order across all reaction-ordered input files

> **About the example filenames:** `1.R.pdb` and `3.P.pdb` mirror the numbered reactant/product files shipped in the geranyl pyrophosphate (GPP) C6-methyltransferase BezA example directory ([`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples) — `1.R.pdb` = reactant state, `3.P.pdb` = product state, with intermediate `2.*.pdb` files for multi-step runs). Replace them with the two (or more) full-system PDBs for your own reaction. To run the commands below verbatim, first fetch the bundled example: `git clone https://github.com/t-0hmura/pdb2reaction && cd pdb2reaction/examples`.

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

> **VRAM warning:** `--dft` launches GPU4PySCF single-point jobs on the
> extracted cluster. Memory use depends on the structure, basis, functional,
> precision, and software stack; pilot a representative state and monitor peak
> memory on the target node. On OOM, run `pdb2reaction dft` separately with a
> smaller basis / trimmed cluster or use a larger node. The `[dft]` extra must
> also be installed (see [Installation](installation.md) Step 7).

## Expected output

A successful run produces a directory like:

```text
result_all/
├── summary.log                    # Human-readable summary
├── summary.json                   # Machine-readable results
├── mep.pdb                        # Concatenated MEP path (promoted to the root)
├── energy_diagram_MEP.png         # All-segment MEP energy profile
└── _work/                         # Pipeline scratch (safe to delete)
    └── path_opt/                  # Raw MEP-engine output (path_search/ with --refine-path)
        ├── hei_seg_01.{xyz,pdb}   # Highest-energy MEP image
        └── summary.json           # MEP engine results
```

The minimal command stops after the MEP stage and therefore does **not** create
`segments/`. With `--tsopt`, a successfully validated reactive segment adds
`segments/seg_01/{reactant.pdb,ts.pdb,product.pdb}`, `ts/`, and `irc/`; adding
`--thermo` also adds `freq/`.

### Output validation

1. `summary.json` — use `scientific_status` and `scientific_status_reasons` for usability. In path mode, `segments[].barrier_kcal` is the raw MEP electronic barrier; requested post-processing results are reported under `rate_limiting_step` and `post_segments`. `status` is retained for compatibility.
2. `_work/path_opt/hei_seg_01.pdb` — inspect the highest-energy image; with `--tsopt`, also inspect the canonical `segments/seg_01/*.pdb` R/TS/P structures
3. `energy_diagram_*.png` — the energy profile should show a clear barrier

**Sample terminal output (successful run):**

```
[time] Elapsed Time for Whole Pipeline: HH:MM:SS.sss
```

(Wall-clock varies with system size, GPU, and selected stages.)

If `--tsopt` is enabled, you should also see:

```
[Imaginary modes] n=1 ([-425.9])
```

A first-order saddle point shows exactly one imaginary mode along the reaction coordinate. IRC validation (run automatically as part of `--tsopt`) confirms it connects the expected reactant and product.

## Tips

- `pdb2reaction all --help` shows core options; `pdb2reaction all --help-advanced` shows the full list.

## Next step

- Scan-defined single-structure route: [Quickstart: `pdb2reaction all --scan-lists`](quickstart-scan.md)
- TS candidate validation: [Quickstart: TS-only mode](quickstart-tsopt-freq.md)
- Full option reference: [all](all.md)
