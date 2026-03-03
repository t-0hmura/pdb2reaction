# `pdb2reaction add-elem-info`

```text
pdb2reaction ver. 0.2.1.dev39+gd7f50241d

Usage: pdb2reaction add-elem-info [OPTIONS]

  Add/repair element columns (77–78) in a PDB using Biopython.

Options:
  --help-advanced               Show all options (including advanced settings)
                                and exit.
  -i, --input FILE              Input PDB filepath.  [required]
  -o, --out FILE                Output PDB filepath (default: replace ".pdb"
                                with "_add_elem.pdb"; when provided, --overwrite
                                is ignored).
  --overwrite / --no-overwrite  Overwrite the input file in-place when -o/--out
                                is omitted.  [default: no-overwrite]
  -h, --help                    Show this message and exit.
```

---

See also: [add-elem-info workflow](../../add_elem_info.md) | [extract](extract.md)
