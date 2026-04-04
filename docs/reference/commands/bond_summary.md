# `pdb2reaction bond-summary`

```text

Usage: pdb2reaction bond-summary [OPTIONS] [EXTRA_INPUTS]...

  Detect bond changes between consecutive structures.

Options:
  --help-advanced             Show all options (including advanced settings) and
                              exit.
  -i, --input TEXT            Input structure files (XYZ/PDB/GJF). Repeat -i for
                              each file.
  --device TEXT               Compute device for distance calculations.
                              [default: cpu]
  --bond-factor FLOAT         Scaling factor for covalent radii sum.  [default:
                              1.2]
  --one-based / --zero-based  Use 1-based atom indices in output.  [default:
                              one-based]
  --out-json / --no-out-json  Print machine-readable JSON to stdout instead of
                              text.  [default: no-out-json]
  -h, --help                  Show this message and exit.
```
