# Gaussian gjf format (gjf.md)

`gjf` (or `.com`) is the input format for Gaussian. `pdb2reaction`
reads it as an alternative to PDB / XYZ when you have a Gaussian-style
input that already encodes charge, spin, and a route specification.

## Structure

```
%nproc=8
%mem=8GB
%chk=run.chk

# <route line: functional/basis + keywords; user-supplied>

  Title (one line, free text)

0 1
  C   0.000   0.000   0.000
  H   0.000   1.090   0.000
  ...

[blank line]
[optional: connectivity table, ECP, basis set, ...]
```

Layout:

| Block | Lines | Content |
|---|---|---|
| `%`-section | 0+ | Link0 commands: `%nproc`, `%mem`, `%chk`, … |
| Route line | 1 | Starts with `#`, declares method, basis, jobtype |
| (blank) | 1 | Required separator |
| Title | 1 | Free text |
| (blank) | 1 | Required separator |
| Charge / spin | 1 | Two integers separated by space: `<charge> <multiplicity>` |
| Coordinates | n | `<element>  <x>  <y>  <z>` (Å), or with frozen flag `<element> -1 <x> <y> <z>` |
| (blank) | 1 | Terminator |
| (optional) | 0+ | Connectivity, ECP block, custom basis, etc. |

`pdb2reaction` reads the geometry block plus the charge/spin line. It
**ignores** the route line (the MLIP backend / DFT engine selected on
the CLI overrides whatever's in the gjf).

## What `pdb2reaction` parses

| gjf field | Where it goes |
|---|---|
| Element + x/y/z | Geometry, fed to the calculator |
| Charge | `-q` (only if you don't override on CLI) |
| Multiplicity | `-m` (only if you don't override on CLI) |
| Frozen flag (`-1` in column 2) | Preserved verbatim by the gjf I/O layer (captured in the line prefix and re-emitted by `convert_xyz_to_gjf`). It does **not** drive the runtime freeze list — use `--freeze-links` / `--freeze-atoms` / YAML `geom.freeze_atoms` to declare frozen atoms. |

## CLI usage

```bash
pdb2reaction tsopt -i ts.gjf -b uma -o result_tsopt
```

No need to pass `-q`/`-m` when they're in the header. To override:

```bash
pdb2reaction tsopt -i ts.gjf -q -1 -m 2 -b uma -o result_tsopt
```

`pdb2reaction dft` reads gjf the same way, but its `--func-basis` flag
takes precedence over the route-line method:

```bash
pdb2reaction dft -i ts.gjf --func-basis 'wb97m-v/def2-tzvpd' --engine gpu
```

## Frozen atoms

A `-1` in the second column of a coordinate line marks the atom as
frozen during optimization. Gaussian convention:

```
  C  -1   1.234   5.678   9.012
```

**`pdb2reaction` does not act on this flag** — i.e. it does not add
the `-1`-marked atoms to the runtime freeze list. The flag itself is
preserved on round-trip: `_GJF_COORD_RE` puts everything before the
x/y/z floats into the line `prefix` capture group, and
`convert_xyz_to_gjf` re-emits the prefix verbatim.

To declare frozen atoms in pdb2reaction, use either:

- `--freeze-links` (CLI; default-on for `extract`-produced clusters
  with link-H caps),
- `--freeze-atoms 'i,j,k,...'` (CLI, 1-based), or
- a YAML `--config` with `geom.freeze_atoms: [<index>, ...]`.

If you have a Gaussian gjf with frozen markers you care about, extract
those indices in a preprocessing step and feed them via `--freeze-atoms`
or `geom.freeze_atoms`.

## Generating gjf from `pdb2reaction` output

`pdb2reaction extract` itself writes only PDB. To get a `.gjf` from a
`pdb2reaction` stationary point, the simplest path is to let
`pdb2reaction path-search` (or any geometry-touching subcommand) emit
GJF companions via `--convert-files/--no-convert-files` (default on);
that pipeline reuses the user-supplied gjf template's route line and
header, so you don't have to specify the QM method here.

To convert a single XYZ programmatically, use the in-tree helper
`pdb2reaction.utils.convert_xyz_to_gjf(xyz_path, template, out_path)`,
which writes the new geometry under the existing template's header
(charge / spin / route / link0 / connectivity sections all preserved):

```python
from pathlib import Path
from pdb2reaction.utils import parse_gjf_template, convert_xyz_to_gjf
template = parse_gjf_template(Path("template.gjf"))   # any prior gjf
convert_xyz_to_gjf(Path("ts.xyz"), template, Path("ts.gjf"))
```

If no template is available, hand-write the route line in a stub gjf
once and reuse it; do not embed a specific QM method in this skill,
since the appropriate functional/basis depends on the user's downstream
calculation.

## Validation

```bash
# Are charge and spin present?
awk '/^[ ]*-?[0-9]+ +[1-9][0-9]*$/{print "charge,spin:", $0; exit}' ts.gjf

# Atom count
awk 'BEGIN{at=0} /^[ ]*[A-Z][a-z]?[ ]+(-?[0-9]+|[0-9.]+)/{at++} END{print "atoms:", at}' ts.gjf

# Route line
grep -m1 '^#' ts.gjf
```

## Common gotchas

| Symptom | Cause | Fix |
|---|---|---|
| `pdb2reaction` says "unknown element symbol Mn"  | Two-letter elements written without correct casing (`MN`) | Use `Mn`, not `MN`, in the element column |
| Charge is read wrong | A blank line missing between title and charge/spin block | Insert the missing blank line |
| Frozen flag ignored | `-1` is in the wrong column or there's a tab instead of space | Use space-separated columns; `-1` must be the 2nd whitespace-separated token |
| Connectivity block confuses the parser | Trailing blocks (after coords) are not parsed by `pdb2reaction` and are discarded | Strip them or just ignore the warning |

## See also

- `pdb.md`, `xyz.md` — the other two formats.
- `charge-multiplicity.md` — when the gjf header is missing or wrong.
- [`pdb2reaction-cli/dft.md`](../pdb2reaction-cli/dft.md) — gjf is the most natural input for
  `pdb2reaction dft` because it carries charge / spin already.