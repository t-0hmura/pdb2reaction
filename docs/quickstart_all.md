# Quickstart: `pdb2reaction all`

## Goal

Run the end-to-end workflow once from two full PDB structures.

## Minimal command

```bash
pdb2reaction all -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --out-dir ./result_all
```

If you want post-processing in the same run:

```bash
pdb2reaction all -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft --out-dir ./result_all
```

## What to check

- `result_all/summary.log`
- `result_all/summary.yaml`
- `result_all/path_search/mep.pdb` (or segment outputs under `result_all/path_search/seg_*/`)

## Tips

- Use `--dry-run` first to validate parsing and execution plan without running heavy stages.
- `pdb2reaction all --help` shows core options; `pdb2reaction all --help-advanced` shows the full list.

## Next step

- Single-structure scan route: [Quickstart: `pdb2reaction scan` with `--spec`](quickstart_scan_spec.md)
- TS validation route: [Quickstart: `pdb2reaction tsopt` -> `pdb2reaction freq`](quickstart_tsopt_freq.md)
- Full option reference: [all](all.md)
