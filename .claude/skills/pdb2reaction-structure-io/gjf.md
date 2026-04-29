# Gaussian gjf format (gjf.md)

`gjf` (or `.com`) is the Gaussian input format. From `pdb2reaction`'s
perspective it is **XYZ with a 5-line header** that adds charge / spin /
route line / title; everything from the coordinate block onward is the
same as XYZ.

## XYZ → GJF in one diff

```
                               %nproc=8                ← optional %-block
                               %mem=8GB                ← (skip if unused)
                               %chk=run.chk
                                                        ← (blank)
                               # <route line>          ← method/basis/keywords
                                                        ← (blank)
                               Title (free text)       ← any string
<n_atoms>                                              ← (blank)
<comment line>                 0 1                      ← <charge> <multiplicity>
<element> x y z                <element> x y z          ← coordinates (Å)
<element> x y z                <element> x y z
...                            ...
                                                        ← (terminating blank)
```

Left = XYZ (2-line header). Right = GJF (~5-line header, same coordinates).

## What `pdb2reaction` reads from a `.gjf`

| gjf field | Used as |
|---|---|
| Element + x/y/z | Geometry |
| Charge | `-q` (CLI overrides) |
| Multiplicity | `-m` (CLI overrides) |
| Route line | **Ignored** — backend / `--func-basis` decides the QM method |
| Frozen flag (`-1` in column 2) | **Preserved on round-trip** but does **not** drive the runtime freeze list. Use `--freeze-atoms` / `--freeze-links` / YAML `geom.freeze_atoms`. |
| Connectivity / ECP / custom basis (after coords) | Discarded |

## Minimal CLI usage

```bash
# Charge / spin already in the gjf header
pdb2reaction tsopt -i ts.gjf -b uma -o result_tsopt

# CLI override
pdb2reaction tsopt -i ts.gjf -q -1 -m 2 -b uma -o result_tsopt

# DFT — note --func-basis takes precedence over the route line
pdb2reaction dft -i ts.gjf --func-basis 'wb97m-v/def2-tzvpd' --engine gpu
```

## Generating gjf from pdb2reaction outputs

`pdb2reaction extract` itself writes only PDB. Use one of:

- `--convert-files` (default on) on any geometry-touching subcommand
  reuses the user-supplied gjf template's header to emit `.gjf`
  companions to every output `.xyz`.
- Programmatic conversion:
  ```python
  from pathlib import Path
  from pdb2reaction.utils import parse_gjf_template, convert_xyz_to_gjf
  template = parse_gjf_template(Path("template.gjf"))
  convert_xyz_to_gjf(Path("ts.xyz"), template, Path("ts.gjf"))
  ```

If no template is on hand, hand-write the header once and reuse it; the
appropriate functional / basis depends on the downstream calculation.

## See also

- `pdb.md`, `xyz.md` — the other two input formats.
- `charge-multiplicity.md` — deriving `-q` / `-m` when the gjf header is
  missing or wrong.
- [`pdb2reaction-cli/dft.md`](../pdb2reaction-cli/dft.md) — gjf is the
  most natural input for `pdb2reaction dft` because it carries charge
  and spin already.