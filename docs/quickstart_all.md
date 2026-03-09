# Quickstart: `pdb2reaction all`

## Goal

Run the end-to-end workflow once from two full PDB structures.

## Prerequisites

- pdb2reaction installed (see [Installation](installation.md))
- Two PDB files (reactant R and product P) with **hydrogen atoms** already added
- The same atoms in the same order across all input PDB files

## Minimal command

```bash
pdb2reaction all -i R.pdb -i P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --out-dir ./result_all
```

If you want post-processing in the same run:

```bash
pdb2reaction all -i R.pdb -i P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft --out-dir ./result_all
```

## What to check

- `result_all/summary.log`
- `result_all/summary.yaml`
- `result_all/path_search/mep.pdb` (or segment outputs under `result_all/path_search/seg_*/`)

## Tips

- `pdb2reaction all --help` shows core options; `pdb2reaction all --help-advanced` shows the full list.

## Next step

- Single-structure staged scan route: [Quickstart: `pdb2reaction scan`](quickstart_scan.md)
- TS optimization and validation: [Quickstart: `pdb2reaction tsopt`](quickstart_tsopt_freq.md)
- Full option reference: [all](all.md)
