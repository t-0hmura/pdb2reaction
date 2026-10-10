# Structure formats

Per-format layout for pdb2reaction inputs and outputs. Format choice, selectors,
and the charge rules are in [SKILL.md](SKILL.md).

## Fields at a glance

| Format | Coordinates | Identity | Charge and spin |
|---|---|---|---|
| PDB | Columns 31–54, Å | Atom name 13–16, residue name 18–20, chain 22, residue number 23–26, element 77–78 | Not read from the file; from residue names and `-l`, or `-q` / `-m` |
| mmCIF | `_atom_site` rows | Original chain IDs, residue numbers, insertion codes, residue and atom names, `type_symbol` | As for PDB |
| XYZ | Lines 3 onward: `<element> <x> <y> <z>` | Element only | None; `-q` / `-m`, or `--ref-pdb` with `-l` |
| GJF | After the charge and multiplicity line | Element only | Header line `<charge> <multiplicity>`, multiplicity = 2S+1 |

## PDB

| Record | Use in pdb2reaction |
|---|---|
| `ATOM` | Polymer atoms by PDB convention |
| `HETATM` | Ligands, metals, waters, cofactors, and cap hydrogens |
| `TER` | Not parsed; `extract` finds chain breaks from the C–N peptide distance |
| `END` | Informational only |
| `MODEL` / `ENDMDL` | Only the first model of an input is read, with a warning; `extract` writes multi-MODEL output for several inputs, and `fix-altloc` processes each MODEL block |
| `CRYST1` | Not used and omitted from cluster PDBs |

`ANISOU`, `LINK`, `SSBOND`, and other connectivity records are not used. Keep
the original file; if a parser needs a cleanup, write a separate derived PDB
and check atom identity and order afterwards.

Columns are 1-based and inclusive:

| Cols | Field | Format | Example |
|---|---|---|---|
| 1–6 | Record name | Left-justified | `ATOM  ` |
| 7–11 | Atom serial | Right-justified integer | `   42` |
| 13–16 | Atom name | Fixed width; alignment depends on element and name, so keep a parser's formatting | ` CB ` |
| 17 | altLoc | Character | ` ` or `A` |
| 18–20 | Residue name | Upper-case three letters | `SAM` |
| 22 | Chain ID | Character | `A` |
| 23–26 | Residue number | Right-justified integer | `  44` |
| 27 | Insertion code | Character | ` ` |
| 31–38, 39–46, 47–54 | X, Y, Z | Right-justified, 3 decimals | `   4.050` |
| 55–60 | Occupancy | Right-justified, 2 decimals | `  1.00` |
| 61–66 | B-factor | Right-justified, 2 decimals | `  0.00` |
| 77–78 | Element | Right-justified symbol | ` C`, `Fe` |
| 79–80 | Formal charge | Right-justified | `2+` |

The record type does not set the charge: pdb2reaction uses the residue name,
its internal amino-acid, ion, and water tables, and `-l` for unknown residues.
`pdb2reaction add-elem-info` fills blank element columns, which some structure
editors leave; run it only when elements are missing and check the symbols
before extraction. For 10,000 or more residues, atom serials past 99,999,
longer residue numbers, or multi-character chain IDs, use mmCIF.

### Per-residue charge (-l)

```bash
pdb2reaction extract -i complex.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -o cluster.pdb
```

Name only unknown/non-standard residues in `-l`. `MG` is already +2 from the
ion table; adding `MG:2` is accepted without effect, and a different value such
as `MG:3` is ignored with a warning. The total is the sum over all residues
kept after extraction.

Amber terminal residue names (`NPRO`, `CGLU`, …, PDB columns 18–21) are read as
the standard residue, and the terminal charge comes from the atoms present: +1
only when H1, H2, and H3 remain (two of them for the Pro or Hyp ring N), −1
only when OXT remains. Check the terminal H atoms and OXT before trusting the
summed charge, and select such a residue by the name written in the file
(`-c NPRO`; `-c PRO` does not match it).

### Cap hydrogens

When `extract` cuts an amino-acid CB–CA, CA–N, or CA–C bond (only CA–C for
proline and hydroxyproline), it adds a hydrogen 1.09 Å from the kept carbon toward the removed
partner. It does not cap other cuts. Each cap is a `HETATM` with atom name
`HL`, residue name `LKH`, and chain `L`, and carries no formal charge. Geometry commands freeze the
carbon parent of each cap by default (`--freeze-links`), finding it again from
the geometry at run time; the PDB holds no freeze list. Frozen-atom control is in
[extract](../pdb2reaction-cli/extract.md#freeze-atoms-at-the-cluster-boundary);
where to cut is in
[model-setup](SKILL.md#check-the-boundary-and-the-charge).

### Common edits

Rename residue 44 of chain A from CYS to CSS by column, which is safer than a
regex on fixed-width lines:

```bash
awk 'BEGIN{OFS=""} \
  ($1=="ATOM" || $1=="HETATM") && substr($0,18,3)=="CYS" \
    && substr($0,22,1)=="A" && substr($0,23,4)+0==44 \
    { $0 = substr($0,1,17) "CSS" substr($0,21) } \
  { print }' my.pdb > my_renamed.pdb
```

Fill the element column, or keep one altLoc per residue:

```bash
pdb2reaction add-elem-info -i my.pdb -o my_with_elem.pdb
pdb2reaction fix-altloc -i my.pdb -o my_clean.pdb
```

Four-point water models (OPC, TIP4P) carry a massless virtual site on each
water; `add-elem-info` labels it `EP` and `extract` keeps it, but it has no
nucleus and does not belong in a cluster model. Delete the virtual sites from
the full MD structure before `extract`, and keep a map between the new atom
order and the MD atom order.

Use Biopython for larger edits; it handles altLoc and `ANISOU`:

```python
from Bio.PDB import PDBParser, PDBIO
p = PDBParser(QUIET=True).get_structure("x", "my.pdb")
for atom in p.get_atoms():
    if atom.get_name() == "OD1" and atom.get_parent().get_resname() == "ASP":
        atom.set_bfactor(20.0)        # cosmetic only
io = PDBIO()
io.set_structure(p)
io.save("my_edited.pdb")
```

The B-factor column is never a freeze flag. `extract` and `add-elem-info` pass
it through unchanged. Freeze atoms with `--freeze-atoms`, `--freeze-links`, or
`geom.freeze_atoms` in YAML
([extract](../pdb2reaction-cli/extract.md#freeze-atoms-at-the-cluster-boundary)).

### Validation checks

```bash
# atom count, and unique residues (chain, number, insertion code)
grep -c '^ATOM\|^HETATM' my.pdb
awk '/^ATOM|^HETATM/{print substr($0,22,6)}' my.pdb | sort -u | wc -l

# missing element columns
awk '/^ATOM|^HETATM/{e=substr($0,77,2); if(e ~ /^[[:space:]]*$/) print NR, $0}' my.pdb

# duplicate atom names within one residue
awk '/^ATOM|^HETATM/{key=substr($0,22,6)"-"substr($0,13,4); print key}' my.pdb \
    | sort | uniq -c | awk '$1>1'
```

## mmCIF and very large structures

Use `.cif` or `.mmcif` when chain IDs are longer than one character, residue
numbers need more than four columns, atom serials need more than five, or the
model has 10,000 or more residues. Do not widen PDB fields instead; fixed-column
readers then misread atom names, residue numbers, or coordinates. A large or
non-standard PDB, with 10,000 or more residues, 99,999 or more atoms, hybrid-36
numbering, or overflowing numbers, is handled the same way as mmCIF.

How the input is read:

1. Only the first `_atom_site.pdbx_PDB_model_num` model is used.
2. Each residue keeps one altLoc, the one with the highest mean occupancy, and
   blank or shared atoms are kept.
3. Atoms get temporary one-character chain IDs and residue numbers 1–9999,
   over as many chains as needed. This covers up to 619,938 residues
   (62 × 9,999).
4. The calculation runs on that internal PDB. Its atom serials can repeat past
   99,999; atom order, not the serial, identifies atoms.
5. With `--convert-files` enabled (the default), coordinate outputs keep the
   internal PDB used between stages and also get a `.cif`; the `.cif` keeps
   the original chain IDs and residue numbers, and also the insertion codes,
   residue names, and atom names. Its `_atom_site.id` runs sequentially
   without the PDB width limit.

Selectors use the original identifiers:

```bash
pdb2reaction extract -i complex.cif -c 'SAM' -o model.pdb                    # every SAM
pdb2reaction extract -i complex.cif -c 'LONG_CHAIN:SAM' -o model.pdb         # every SAM in chain LONG_CHAIN
pdb2reaction extract -i complex.cif -c 'LONG_CHAIN:SAM:10001' -o model.pdb   # one SAM past the PDB limit
pdb2reaction extract -i complex.cif -c 'LONG_CHAIN:10001' -o model.pdb
```

For scan atoms, `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`, such as
`LONG_CHAIN:SAM:10001:C1` or `A:SAM:12B:C1`, separates repeated numbering.
Chain IDs are case-sensitive; `A` and `a` can be different chains.

A geometry command on `complex.cif` writes both files:

```text
final_geometry.pdb  # internal topology used between stages
final_geometry.cif  # coordinates with the original identifiers
```

The chain IDs and residue numbers in that PDB are temporary; report residues
from the `.cif`.

Limits:

- Only the first model is read. Multi-frame outputs are written as several CIF
  models.
- Categories other than `_atom_site` are not copied. Identity, occupancy,
  B-factor, formal charge, and atom labels are kept; refinement metadata is not.
- Inputs along one path still need identical atoms in identical order; format
  conversion does not fix atom mapping.
- `fix-altloc` and `add-elem-info` are PDB tools. For mmCIF, altLoc is chosen
  and elements are taken from `_atom_site.type_symbol` on reading; a row without
  `type_symbol` stops with an error.

## XYZ

XYZ holds elements and Cartesian coordinates only. pdb2reaction writes it for
trajectories, optimized stationary points, and IRC paths.

```text
<n_atoms>
<comment line>
<element>  <x>  <y>  <z>
...
```

Line 1 is the integer atom count, line 2 a free comment, and each following
line one atom with its element symbol and coordinates in Å. A trajectory
repeats this block per frame.

pdb2reaction writes `*_trj.xyz` with a bare energy as the comment line:

```text
56
-11148.201817745587
H   0.123  4.567  8.901
```

On reading, an explicit `E=` or `Energy:` token wins; otherwise a single number
is taken as the energy, and several bare numbers are rejected as ambiguous.
pdb2reaction does not write extended-XYZ tags such as `Lattice=`,
`Properties=`, or `pbc=`; they appear only in files written by ASE itself.

XYZ has no charge or spin, so give both, or point `--ref-pdb` at the matching
PDB so `-l` works:

```bash
pdb2reaction tsopt -i ts.xyz -q 0 -m 1 -b uma -o result_tsopt
pdb2reaction dft -i ts.xyz -q -1 -m 1 --func-basis 'wb97m-v/def2-svp'
pdb2reaction tsopt -i ts.xyz --ref-pdb cluster.pdb -l 'SAM:1,GPP:-3' -m 1
```

Common edits with ASE:

```python
from ase.io import read, write
trj = read("mep.xyz", ":")       # all frames
write("ts.xyz", trj[5])          # one frame, index 5
write("ts.pdb", read("ts.xyz"))  # every atom gets residue 'MOL'; fix names if needed
```

Without Python, frame k (counted from 1) is
`N=$(head -1 scan_trj.xyz); sed -n "$(( (k-1)*(N+2)+1 )),$(( k*(N+2) ))p" scan_trj.xyz`.
To keep the PDB topology for the next command, pass the original PDB to
`--ref-pdb`. `cat reactant.xyz ts.xyz product.xyz > rts.xyz` makes a valid
trajectory for `trj2fig`
([utilities](../pdb2reaction-cli/utilities.md#trj2fig)).

```bash
# atom count consistent with line 1
awk 'NR==1{n=$1; expected=n+2} END{if(NR!=expected) print "BAD: line count " NR " expected " expected}' ts.xyz
# non-element symbols
awk 'NR>2 && !/^[A-Z][a-z]?[ ]/{print "weird element: " $0}' ts.xyz
# frame count
awk 'NR==1{n=$1; per=n+2} END{print "frames:", NR/per}' mep.xyz
```

## GJF

Only the `.gjf` extension is parsed as Gaussian input; rename `.com` files
first. The file needs the standard blank-line-separated layout, the charge and
multiplicity line, and the coordinate block. Unlike XYZ, it has no atom-count
or comment line.

| Block | Required | Use in pdb2reaction |
|---|---|---|
| Link0 (`%nproc`, `%mem`, `%chk`) | Optional | Ignored |
| Route line (`# method basis keywords`) | Yes, for Gaussian | Ignored; the backend or `--func-basis` sets the method |
| Title | Yes, for Gaussian | Ignored |
| `<charge> <multiplicity>` | Yes | Total charge and multiplicity; `-q` / `-m` override |
| Coordinates (`<element> x y z`), ending with a blank line | Yes | Geometry |
| Frozen flag (`-1` in column 2) | Optional | Kept in written files but not used for freezing; use `--freeze-atoms`, `--freeze-links`, or `geom.freeze_atoms` |
| Connectivity, ECP, or custom basis after the coordinates | Optional | Not used, but copied unchanged into every `.gjf` written with `--convert-files` |

Because that trailing block is copied from the input, its connectivity
describes the input bonding. Remove it from the template when outputs such as
TS or IRC frames bond differently, so a later Gaussian job does not read
connectivity that contradicts its coordinates.

```bash
pdb2reaction tsopt -i ts.gjf -b uma -o result_tsopt                  # charge and spin from the header
pdb2reaction tsopt -i ts.gjf -q -1 -m 2 -b uma -o result_tsopt       # override
pdb2reaction dft -i ts.gjf --func-basis 'wb97m-v/def2-tzvpd' --dft-engine gpu
```

`extract` does not write GJF. With a `.gjf` input and `--convert-files`
enabled, geometry commands write `.gjf` outputs from the input header. A
single-frame output is a normal Gaussian input. A multi-frame output has one
header and blank-separated coordinate blocks; it is a coordinate archive, not
an executable QST or Link1 job, so take one frame out first. Without a
template, write the header by hand once and reuse it; the functional and basis
depend on the downstream calculation.

## Ligand, ion, and metal charges

Common ligand charges, not defaults; confirm the modeled protonation state for
the mechanism.

| Ligand | Residue name | Charge at pH 7 |
|---|---|---|
| S-Adenosylmethionine | `SAM` | +1 |
| S-Adenosylhomocysteine | `SAH` | 0 |
| Geranyl pyrophosphate | `GPP` | −3 |
| ATP / ADP / GTP | `ATP` / `ADP` / `GTP` | −4 / −3 / −4 |
| NADH / NAD⁺ | `NAI` or `NDH` / `NAD` | −2 / −1 |
| FAD | `FAD` | −2 |
| Pyridoxal phosphate | `PLP` | −2 |
| Heme | `HEM` | From the oxidation state, axial ligands, and propionate protonation |
| Free phosphate | `PO4` | −2 to −3 |

Monatomic ions follow the rules in
[Charge and multiplicity](SKILL.md#charge-and-multiplicity). The values follow the PDB CCD residue
names, so the oxidation state lives in the name:

| Ion | Residue name | Charge |
|---|---|---|
| Mg²⁺, Zn²⁺, Mn²⁺ | `MG`, `ZN`, `MN` | +2 |
| Fe³⁺ / Fe²⁺ | `FE` / `FE2` | +3 / +2 |
| Cu²⁺ / Cu⁺ | `CU` / `CU1` | +2 / +1 |
| Na⁺, K⁺ | `NA`, `K` | +1 |
| Cl⁻ | `CL` | −1 |

`FE3` is not in the table; to check, dump the table with the one-liner in
[SKILL.md](SKILL.md#charge-and-multiplicity). When the deposited name disagrees
with the oxidation state the mechanism needs, such as an `FE` that is really
Fe(II), set the cluster total with `-q`. For `all -c`, also fix the residue name
or mapping so the mismatch warning does not hide a model-building error.

| Metal | Common S | Multiplicity (2S+1) |
|---|---|---|
| Mn²⁺ (d⁵), high spin | 5/2 | 6 |
| Fe²⁺ (d⁶), high spin, tetrahedral or weak field | 2 | 5 |
| Fe²⁺ (d⁶), low spin, octahedral strong field | 0 | 1 |
| Fe³⁺ (d⁵), high spin | 5/2 | 6 |
| Co²⁺ (d⁷), high spin | 3/2 | 4 |
| Cu²⁺ (d⁹) | 1/2 | 2 |
| Zn²⁺ (d¹⁰) | 0 | 1 |

## See also

- [SKILL.md](SKILL.md): format choice, selectors, and the charge rules.
- [extract](../pdb2reaction-cli/extract.md): `-c`, `-l`, and frozen atoms.
- [utilities](../pdb2reaction-cli/utilities.md): `fix-altloc`, `add-elem-info`, `trj2fig`, and `energy-diagram`.
- [dft](../pdb2reaction-cli/dft.md): GJF is a natural `dft` input, since it carries charge and spin.
- [outputs](../pdb2reaction-overview/outputs.md): where output XYZ, PDB, and CIF files are written.
