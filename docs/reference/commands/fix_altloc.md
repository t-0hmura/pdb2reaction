# pdb2reaction fix-altloc

```
pdb2reaction ver. 0.2.1.dev116+gff86b942a

Usage: cli fix-altloc [OPTIONS]

  Blank PDB altLoc column (col 17) without shifting, and keep one altLoc per
  atom by default rule: highest occupancy, then earliest appearance.

Options:
  --help-advanced   Show all options (including advanced settings) and exit.
  -i, --input PATH  Input PDB file or directory.  [required]
  -o, --out PATH    Output file (if input is a file) or output directory (if
                    input is a directory).
  -h, --help        Show this message and exit.
```
