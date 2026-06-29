# `bond-summary`

Detect and report covalent bond changes between consecutive molecular structures (R → TS → P or multi-intermediate chains) via element-specific covalent-radius perception. For *N* input files it produces *N − 1* comparison blocks (A→B, B→C, …) and prints them to stdout; no file is written. Use it to audit which covalent bonds form or break between sequential structures along a reaction path — e.g. validating an IRC endpoint pair, screening multistep mechanisms, or sanity-checking `all` post-processing manually. Supported input formats are **XYZ**, **PDB**, and **GJF** (auto-detected by extension); distances are reported in Ångström.

## Examples

Two-structure comparison (R → P):

```bash
pdb2reaction bond-summary -i 1.R.xyz 3.P.xyz
```

Multi-structure chain — produces three comparison blocks (R→IM1, IM1→IM2, IM2→P):

```bash
pdb2reaction bond-summary -i 1.R.xyz 3.IM1.xyz 5.IM2.xyz 7.P.xyz
```

## Outputs

`bond-summary` writes no files. It prints a text report to stdout (one comparison block per consecutive pair), or machine-readable JSON to stdout with `--json`. Each text block lists the formed and broken bonds with their before/after distances in Å:

```
============================================================
  1.R.xyz  →  3.P.xyz
============================================================
Bond formed (2):
  - O14-H106 : 1.502 Å --> 1.011 Å
  - P95-O107 : 3.477 Å --> 1.523 Å
Bond broken (2):
  - P95-O97 : 1.585 Å --> 3.270 Å
  - H106-O107 : 1.034 Å --> 1.673 Å
```

## Python API

The bond change detection functions can also be used programmatically:

```python
from pdb2reaction.domain.bond_changes import compare_structures, has_bond_change, summarize_changes
```

| Function | Description |
|----------|-------------|
| `compare_structures(geom1, geom2, device="cuda", bond_factor=1.20)` | Detect covalent bonds formed or broken between two pysisyphus geometries. Returns a `BondChangeResult` with `formed_covalent` and `broken_covalent` (sets of 0-based index pairs). |
| `has_bond_change(geom_start, geom_end, bond_cfg)` | Convenience wrapper: returns `(has_changes: bool, summary_text: str)`. `bond_cfg` accepts keys: `device`, `bond_factor`, `margin_fraction`, `delta_fraction`. |
| `summarize_changes(geom, result, one_based=True)` | Format bond changes as a text report with distances in Angstrom. |

### Example

```python
from pysisyphus.helpers import geom_loader
from pdb2reaction.domain.bond_changes import has_bond_change

geom_r = geom_loader("R.xyz")
geom_p = geom_loader("P.xyz")

changed, summary = has_bond_change(geom_r, geom_p, {"device": "cpu", "bond_factor": 1.20})
if changed:
    print(summary)
```

## CLI options

| Option | Description | Default |
|--------|-------------|---------|
| `-i, --input FILE` | Input structure file in XYZ, PDB, or GJF format (auto-detected by extension; repeat for each file, ≥ 2 required) | — |
| `--device TEXT` | Compute device (`cpu`, `cuda`) | `cpu` |
| `--bond-factor FLOAT` | Scaling factor for covalent radii sum | `1.20` |
| `--one-based / --zero-based` | Atom index convention in output | `--one-based` |
| `--json / --no-json` | **Emit machine-readable JSON to stdout** instead of the text report. See [JSON Output Schema → bond-summary](json-output.md#bond-summary) for the schema and stdout/persistence behavior. | `False` |

The full flag list is in the generated [command reference](reference/commands/index.md).

## Notes

- All input structures must have **identical atom counts and element ordering**.
- Bond detection uses the same algorithm the `all` workflow applies for IRC endpoint validation.
- To adjust sensitivity to borderline bonds (e.g., metal coordination at 2.0–2.4 Å), increase `--bond-factor` (e.g., `1.30`).

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [irc](irc.md) -- IRC trajectories whose endpoints are validated by bond detection
- [all](all.md) -- End-to-end workflow that uses bond-change validation internally
- [trj2fig](trj2fig.md) -- Visualize energy profiles from trajectories
