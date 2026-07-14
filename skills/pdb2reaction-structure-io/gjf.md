# Gaussian gjf format (gjf.md)

`gjf` is the Gaussian input format (only the `.gjf` extension is recognized for header parsing; rename `.com` → `.gjf` first). From `pdb2reaction`'s
perspective it provides Cartesian coordinates plus a Gaussian-style header.
The parser needs the standard blank-separated layout, charge/multiplicity
line, and coordinate block; it does not use the route line to configure the
calculation. Unlike XYZ,
the coordinate block has no atom-count/comment header and may be followed by
connectivity, ECP, or custom-basis sections.

## XYZ → GJF block-by-block

| Block | XYZ | GJF |
|---|---|---|
| Link0 (`%nproc`, `%mem`, `%chk`) | — | optional, top of file |
| Route line (`# <method> <basis> <keywords>`) | — | required by a standard Gaussian input; ignored by `pdb2reaction` |
| Title | — | required by a standard Gaussian input (free text) |
| Charge / multiplicity | — | required by Gaussian and by `pdb2reaction` (`<q> <m>`) |
| Atom count | line 1 (`<n_atoms>`) | — |
| Comment line | line 2 (free text) | — |
| Coordinates (`<element> x y z`) | lines 3 … | after charge/mult, terminated by a blank line |
| Connectivity / ECP (optional) | — | after coordinates, blank-separated |

## What `pdb2reaction` reads from a `.gjf`

| gjf field | Used as |
|---|---|
| Element + x/y/z | Geometry |
| Charge | `-q` (CLI overrides) |
| Multiplicity | `-m` (CLI overrides) |
| Route line | **Ignored** — backend / `--func-basis` decides the QM method |
| Frozen flag (`-1` in column 2) | **Preserved on round-trip** but does **not** drive the runtime freeze list. Use `--freeze-atoms` / `--freeze-links` / YAML `geom.freeze_atoms`. |
| Connectivity / ECP / custom basis (after coords) | **Ignored by the calculation, but preserved verbatim** — the whole post-coordinate block is stored on the gjf template (`core/utils.py`) and re-emitted at the bottom of every `.gjf` companion written by `--convert-files`. It is copied unchanged from the input, so a connectivity block always describes the *input* bonding pattern; strip it from the template when the output geometry (TS / IRC frames) has different bonding, so a downstream Gaussian job reads connectivity that matches its own coordinates. |

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

`pdb2reaction extract` does not generate GJF. It writes its internal PDB model
and, for mmCIF/oversized-PDB input, a CIF companion. To generate GJF, use one
of:

- `--convert-files` (default on) on a geometry-touching subcommand reuses the
  user-supplied gjf template's header to emit `.gjf` companions. A single-frame
  output is a normal template-based Gaussian input. For a multi-frame
  trajectory, the current converter writes one header followed by
  blank-separated coordinate blocks and one copied suffix; treat that file as
  a coordinate archive, **not** as a directly executable Gaussian QST/Link1
  input. Split the desired frame and generate a single-frame job first.
- Programmatic conversion:
  ```python
  from pathlib import Path
  from pdb2reaction.core.utils import parse_gjf_template, convert_xyz_to_gjf
  template = parse_gjf_template(Path("template.gjf"))
  convert_xyz_to_gjf(Path("ts.xyz"), template, Path("ts.gjf"))
  ```

If no template is on hand, hand-write the header once and reuse it; the
appropriate functional / basis depends on the downstream calculation.

## See also

- `pdb.md`, `cif.md`, `xyz.md` — the other input formats.
- `charge-multiplicity.md` — deriving `-q` / `-m` when the gjf header is
  missing or wrong.
- [`pdb2reaction-cli/dft.md`](../pdb2reaction-cli/dft.md) — gjf is the
  most natural input for `pdb2reaction dft` because it carries charge
  and spin already.
