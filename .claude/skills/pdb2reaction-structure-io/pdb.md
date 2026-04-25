# PDB format (pdb.md)

The Protein Data Bank format is the primary input to `pdb2reaction`. It
is **column-based** — each field has a fixed character range — so a one-
character drift breaks everything downstream. Treat editing as a
column-aware operation, not as plain-text find-and-replace.

## Record types `pdb2reaction` cares about

| Record | Used for |
|---|---|
| `ATOM` | Standard amino-acid / nucleic-acid atoms (residue ≤ 3 letters, in `pdb2reaction`'s AMINO_ACIDS table) |
| `HETATM` | Ligand, metal, water, cofactor, link-H atoms |
| `TER` | Chain terminator — used by `extract` to identify chain breaks |
| `END`, `ENDMDL` | File terminator — informational only |
| `CRYST1` | Unit cell — read but not written by `pdb2reaction` (cluster model only) |

`pdb2reaction` ignores `MODEL`, `ANISOU`, `LINK`, `SSBOND`, etc. Strip
them with `pdb2reaction add-elem-info` or `fix-altloc` if a downstream
step complains.

## Column layout (`ATOM` / `HETATM`)

Columns are 1-based, inclusive on both ends.

| Cols | Field | Width | Type | Example |
|---|---|---|---|---|
| 1–6 | Record name | 6 | left-just `ATOM  ` / `HETATM` | `ATOM  ` |
| 7–11 | Atom serial | 5 | right-just int | `   42` |
| 13–16 | Atom name | 4 | left-just (with 1-char element prefix) | ` CB ` |
| 17 | Alt-loc indicator | 1 | char | ` ` or `A`/`B` |
| 18–20 | Residue name | 3 | upper-case 3-letter code | `SAM` |
| 22 | Chain ID | 1 | char | `A` |
| 23–26 | Residue sequence | 4 | right-just int | `  44` |
| 27 | Insertion code | 1 | char | ` ` |
| 31–38 | X | 8 | right-just float, 3 decimals | `   4.050` |
| 39–46 | Y | 8 | right-just float, 3 decimals | `  -8.106` |
| 47–54 | Z | 8 | right-just float, 3 decimals | `   6.935` |
| 55–60 | Occupancy | 6 | right-just float, 2 decimals | `  1.00` |
| 61–66 | Temperature factor | 6 | right-just float, 2 decimals | `  0.00` |
| 77–78 | Element symbol | 2 | right-just upper-case | ` C` |
| 79–80 | Formal charge | 2 | right-just (e.g. `2+`, `1-`) | `  ` |

`pdb2reaction.add-elem-info` repairs columns 77-78 when they are blank,
which they often are after PyMOL/Maestro export. Always run it before
`extract` if the elements are missing.

## Residue selectors (the `-c / --center` flag)

`pdb2reaction extract` (and any subcommand that also extracts internally)
uses three forms:

```bash
# Form 1 — comma-separated residue names. Matches every residue with
# resName == SAM, GPP, or MG.
pdb2reaction extract -i complex.pdb -c 'SAM,GPP,MG' -o cluster.pdb

# Form 2 — a separate PDB containing only the substrate residues.
pdb2reaction extract -i complex.pdb -c substrate.pdb -o cluster.pdb

# Form 3 — chain-specific: `chainID:resSeq` or `chainID:resName`
pdb2reaction extract -i complex.pdb -c 'A:44,A:SAM,A:GPP' -o cluster.pdb
```

The pocket radius around the centers is set by `-r <Å>` (default 3.0 Å).
All residues with at least one heavy atom inside the radius are kept.

## Per-residue charge mapping (`-l / --ligand-charge`)

```bash
pdb2reaction extract -i complex.pdb \
    -c 'SAM,GPP,MG' \
    -l 'SAM:1,GPP:-3,MG:2' \
    -o cluster.pdb
```

Standard amino acids are looked up from `pdb2reaction`'s internal
`AMINO_ACIDS` table — you only need to provide ligand / metal /
non-standard residues in `-l`. The total cluster charge is the sum of
all residues kept (post extraction).

If you don't know a ligand's formal charge, see
`charge-multiplicity.md` for the lookup workflow.

## Link-hydrogen capping

When `extract` cuts a covalent bond between an in-cluster atom (`A`)
and an out-of-cluster atom (`B`), it places a hydrogen `H_link` along
the `A→B` direction at 1.09 Å (standard C-H length). The link hydrogen
is written as a `HETATM` with residue name `LKH` (or similar marker
in your version — check by reading `pdb2reaction.extract.AMINO_ACIDS`).

Link hydrogens carry **no formal charge**; they do not enter the
charge sum.

The atoms that **donate** hydrogens (i.e. atom `A`, the cluster-side
parent of each cap) are added to a `freeze_atoms` list so they do not
move during optimization. `pdb2reaction` writes the indices to the
extracted PDB's B-factor field for downstream consumption.

## Common edits

### Rename residue 44 in chain A from CYS to CSS

PDB columns: resName 18-20, chainID 22, resSeq 23-26. Match all three
explicitly via `awk substr` (the natural column-aware tool — sed regex
on PDB lines is error-prone because of the dot-counts):

```bash
awk 'BEGIN{OFS=""} \
  ($1=="ATOM" || $1=="HETATM") && substr($0,18,3)=="CYS" \
    && substr($0,22,1)=="A" && substr($0,23,4)+0==44 \
    { $0 = substr($0,1,17) "CSS" substr($0,21) } \
  { print }' my.pdb > my_renamed.pdb
```

For non-trivial edits use Biopython (handles altloc, ANISOU, …).

### Add element column with `add-elem-info`

```bash
pdb2reaction add-elem-info -i my.pdb -o my_with_elem.pdb
```

### Resolve alt-loc

```bash
pdb2reaction fix-altloc -i my.pdb -o my_clean.pdb
```

### Programmatic edits via Biopython

```python
from Bio.PDB import PDBParser, PDBIO
p = PDBParser(QUIET=True).get_structure("x", "my.pdb")
for atom in p.get_atoms():
    if atom.get_name() == "OD1" and atom.get_parent().get_resname() == "ASP":
        atom.set_bfactor(20.0)        # mark frozen, for example
io = PDBIO()
io.set_structure(p)
io.save("my_edited.pdb")
```

## Validation hooks

```bash
# atom count + residue count
grep -c '^ATOM\|^HETATM' my.pdb
awk '/^ATOM|^HETATM/{print substr($0,18,3)}' my.pdb | sort -u

# any missing element columns?
awk '/^ATOM|^HETATM/{e=substr($0,77,2); if(e=="  ") print NR, $0}' my.pdb

# any duplicate atom names within one residue (often breaks AMBER / parm7)?
awk '/^ATOM|^HETATM/{key=substr($0,22,5)"-"substr($0,13,4); print key}' my.pdb \
    | sort | uniq -c | awk '$1>1'
```

## See also

- `xyz.md`, `gjf.md` — alternative formats.
- `charge-multiplicity.md` — figuring out per-ligand charges.
- `pdb2reaction-cli/extract.md` — full `extract` flag set and examples.
- `pdb2reaction-cli/{add-elem-info,fix-altloc}.md` — utility subcommands.