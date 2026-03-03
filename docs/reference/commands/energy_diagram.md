# `pdb2reaction energy-diagram`

```text
pdb2reaction ver. 0.2.1.dev39+gd7f50241d

Usage: pdb2reaction energy-diagram [OPTIONS]

  Plot an energy diagram from numeric inputs only. Use -i for values (multiple
  numbers or list-like string).

Options:
  --help-advanced    Show all options (including advanced settings) and exit.
  -i, --input TEXT   Numeric sequence. Accepts: -i 0 12.5 4.3  or  -i "[0, 12.5,
                     4.3]"  or repeated -i.  [required]
  -o, --output FILE  Output image path.  [default: energy_diagram.png]
  --label-x TEXT     State labels on x-axis. Accepts: --label-x R TS P  or
                     --label-x "['R','TS','P']".
  --label-y TEXT     Y-axis label.  [default: ΔE (kcal/mol)]
  -h, --help         Show this message and exit.
```

---

See also: [energy-diagram workflow](../../energy_diagram.md) | [dft](dft.md) | [trj2fig](trj2fig.md)
