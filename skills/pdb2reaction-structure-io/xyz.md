# XYZ format (xyz.md)

XYZ is a simple plain-text dump of element + Cartesian coordinate, with
no residue or charge information embedded. `pdb2reaction` writes XYZ for trajectories,
optimized stationary points, and IRC paths; it accepts XYZ as input
when `-q` and `-m` are supplied (or `--ref-pdb` plus `-m`).

## Layout

```
<n_atoms>
<comment line — free text, often used by ASE for energies and properties>
<element>  <x>  <y>  <z>
<element>  <x>  <y>  <z>
...
```

- Line 1: integer atom count (no decimals).
- Line 2: arbitrary string. ASE encodes per-frame energies / properties
  here in a `key=value` syntax (see below).
- Lines 3 ... 2 + n_atoms: one atom each. Element is the periodic-table
  symbol (`H`, `C`, `N`, `Mg`, `Zn`, …); coordinates are Å.
- Multiple frames (a trajectory) are concatenated: each frame restarts
  with the atom-count line.

## Comment line (line 2)

`pdb2reaction` writes `*_trj.xyz` files with a bare floating-point
energy as the comment line:

```
56
-11148.201817745587
H   0.123  4.567  8.901
...
```

The comment is parsed by `read_xyz_energies` (a regex that pulls the
first numeric token), and emitted on write by `xyz_string_with_energy`
(per-frame) / `write_xyz_trj_with_energy` (trajectory) in
`pdb2reaction.core.utils`. ASE-style extxyz tags (`Lattice="..."`,
`Properties=species:S:1:pos:R:3`, `pbc="F F F"`) are **not** produced
by pdb2reaction's writer; they may appear when the file was written by
ASE directly.

## Reading and writing in Python

```python
from ase.io import read, write
atoms = read("ts.xyz")          # single frame
trj   = read("mep.xyz", ":")    # all frames as a list
write("out.xyz", trj)           # round-trip
```

Pure-text reading without ASE:

```python
def read_xyz(path):
    with open(path) as f:
        n = int(f.readline())
        comment = f.readline().rstrip()
        coords = [f.readline().split() for _ in range(n)]
    return n, comment, coords
```

## Charge / multiplicity for XYZ inputs

XYZ has **no header for charge or spin**. When passing XYZ to
`pdb2reaction`, supply both:

```bash
pdb2reaction tsopt -i ts.xyz -q 0 -m 1 -b uma -o result_tsopt
pdb2reaction dft   -i ts.xyz -q -1 -m 1 --func-basis 'wb97m-v/def2-svp'
```

Or, if a corresponding PDB exists, point at it with `--ref-pdb` so the
residue-based charge mapping (`-l`) still resolves:

```bash
pdb2reaction tsopt -i ts.xyz --ref-pdb cluster.pdb -l 'SAM:1,GPP:-3' -m 1
```

If unsure about `-q` / `-m`, see `charge-multiplicity.md` before
running.

## Common edits

### Strip a frame from a trajectory

```python
from ase.io import read, write
trj = read("mep.xyz", ":")
write("ts.xyz", trj[5])         # frame index 5
```

### Convert XYZ → PDB so that residue-aware tools can consume it

```python
from ase.io import read, write
atoms = read("ts.xyz")
# ASE's PDB writer assigns residue 'MOL' to all atoms — overwrite if needed.
write("ts.pdb", atoms)
```

### Combine multiple stationary points into one trajectory

```bash
cat reactant.xyz ts.xyz product.xyz > rts.xyz
```

This works because each frame already starts with its own atom-count
line; `pdb2reaction trj2fig` consumes the result for energy plots.

## Validation hooks

```bash
# atom count consistent with line 1?
awk 'NR==1{n=$1; expected=n+2} END{if(NR!=expected) print "BAD: line count " NR " expected " expected}' ts.xyz

# any non-element symbols?
awk 'NR>2 && !/^[A-Z][a-z]?[ ]/{print "weird element: " $0}' ts.xyz

# multiple frames? frame count = lines / (natoms+2)
awk 'NR==1{n=$1; per=n+2} END{print "frames:", NR/per}' mep.xyz
```

## See also

- `pdb.md` — the format you'll typically convert XYZ back to when
  you need residue / chain context.
- `charge-multiplicity.md` — how to fill `-q` and `-m` when the XYZ
  doesn't say.
- `pdb2reaction-cli/{trj2fig,energy-diagram}.md` — plotting an XYZ
  trajectory directly.
- `pdb2reaction-workflows-output/SKILL.md` — where the toolkit's
  output XYZ files live.