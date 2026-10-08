# Small-molecule reactions

`pdb2reaction` also analyzes reaction paths of small molecules, not only enzyme cluster models. Give the reactant (R) and product (P) as XYZ/GJF/PDB/CIF structures and run `all` without `-c`; with no `-c`, the cluster-model extraction is skipped and the input structures are analyzed as they are.

## First example

[`examples/aromatic_claisen/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples/aromatic_claisen) holds the reactant and product of an aromatic Claisen rearrangement (allyl phenyl ether → 6-allylcyclohexa-2,4-dien-1-one).

```bash
pdb2reaction all -i examples/aromatic_claisen/reactant.xyz examples/aromatic_claisen/product.xyz -q 0 --tsopt --thermo
```

The run searches the minimum energy path (MEP) and continues to transition-state (TS) optimization, the intrinsic reaction coordinate (IRC), and vibrational analysis with thermochemistry. `Scientific status: success` under `====== Pipeline summary ======` near the end of the terminal output means success; `scientific_status` in `summary.json` has the same value.

## Input notes

* **Formats**: XYZ, GJF, PDB, and CIF are accepted. Several structures must have the same atoms in the same order.
* **Charge and spin**: set the charge with `-q` and the spin multiplicity with `-m` (default 1).
* **Other input modes**: the [Scan-list mode](quickstart-scan.md), which builds the path from one structure by a scan, and the [TS-only mode](quickstart-tsopt.md), which starts from a TS candidate, also work without `-c`.
* **Reactions in solution**: `--solvent` adds the xTB solvation correction ([xTB solvent correction](backends.md#xtb-solvent-correction)).

## See also

* [`all`](all.md) — options
* [Quickstart: `pdb2reaction all`](quickstart-all.md) — starting from the structures before and after the reaction
