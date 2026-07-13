# PDB format (pdb.md)

PDB is `pdb2reaction`'s primary input. Column-based, fixed-width fields.

## Record types `pdb2reaction` cares about

| Record | Used for |
|---|---|
| `ATOM` | Standard amino-acid / nucleic-acid atoms (residue ≤ 3 letters, in `pdb2reaction`'s AMINO_ACIDS table) |
| `HETATM` | Ligand, metal, water, cofactor, cap-H atoms |
| `TER` | Chain terminator; `extract` infers chain breaks from C–N peptide-adjacency distance, not by parsing TER directly |
| `END` | File terminator — informational only |
| `MODEL` / `ENDMDL` | `extract` writes multi-MODEL output for multi-structure input, and `fix-altloc` processes each MODEL block. Multi-MODEL **input** is not read — `extract` uses only the first model and warns. |
| `CRYST1` | Unit cell — read but not written by `pdb2reaction` (cluster model only) |

`pdb2reaction` ignores `ANISOU`, `LINK`, `SSBOND`, etc. If a
downstream step complains, strip them with a one-line `awk` / `grep`
filter.

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

# Form 3 — chain-aware: `chainID:resSeq` (numeric resSeq is honored)
pdb2reaction extract -i complex.pdb -c 'A:44' -o cluster.pdb
```

> **Caveat**: only `chainID:resSeq` (numeric) is parsed as chain-aware
> in `extract.py`. Tokens like `'A:SAM'` (chain:resName) fail the
> `_parse_res_tokens` regex (which requires a numeric resSeq), and
> `resolve_substrate_residues` then silently falls through to the
> resname splitter (`[,\s]+`); since the literal token `'A:SAM'` is
> not a known resname, the call raises
> `ValueError("Residue name 'A:SAM' not found in complex.")`. To
> restrict by chain you must supply numeric resSeq. To select all SAM
> in chain A, run `extract` with `-c 'SAM'` first then trim chains by
> hand, or use the substrate-PDB form (`-c <substrate.pdb>`).

The pocket radius around the centers is set by `-r <Å>` (default 2.6 Å).
All residues with at least one atom (any element, including H) inside the radius are kept.

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

## Cap-hydrogen capping

When `extract` cuts a covalent bond between an in-cluster atom (`A`)
and an out-of-cluster atom (`B`), it places a hydrogen `H_link` along
the `A→B` direction at 1.09 Å (standard C-H length). The cap hydrogen
is written as a `HETATM` with atom name `HL` and residue name `LKH`
(hard-coded in `_format_linkH_block`).

Cap hydrogens carry **no formal charge**; they do not enter the
charge sum.

The atoms that **donate** hydrogens (i.e. atom `A`, the cluster-side
parent of each cap) are frozen during optimization. The PDB itself
carries no encoded freeze list — `pdb2reaction.core.utils.detect_freeze_links`
re-derives the parent atom of each `LKH/HL` record at runtime via a
nearest-neighbor search in Cartesian space. The B-factor column on
the LKH atoms is hard-coded to 0.00 by `_format_linkH_block`.

## Manual cluster-boundary checklist

- Make a retained backbone fragment terminate consistently at alpha carbons
  (`CA`) at both main-chain ends, then cap the open valences.
- For side-chain/ligand/cofactor truncation, choose an aliphatic C–C single
  bond whenever possible (`CA–CB` or farther from the reactive center). Do not
  cut peptide C–N, polar C–N/C–O, aromatic/conjugated, disulfide, or
  metal-coordination bonds; move the boundary or include the bonded partner.
- One intended cap per severed valence; visually inspect bond lengths and
  valence. Keep cap parents frozen with `--freeze-links` in production.
- Apply exactly the same atom selection, order, and cap topology to R/IM/P.
  Never re-extract compared states independently.
- Recompute net charge/multiplicity after truncation. A geometrically neat
  cluster with the wrong electron count is still an invalid MLIP input.

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
        atom.set_bfactor(20.0)        # cosmetic edit only; see note below
io = PDBIO()
io.set_structure(p)
io.save("my_edited.pdb")
```

**Note:** `pdb2reaction` does **not** treat the PDB B-factor column as
a freeze flag. To freeze atoms during optimization, use the CLI
`--freeze-atoms` / `--freeze-links` options or the YAML
`geom.freeze_atoms` key (see `pdb2reaction-cli/freeze-atoms.md`).
B-factors edited above are passed through verbatim by `extract` /
`add-elem-info` but have no effect on geometry constraints.

## Validation hooks

```bash
# atom count + residue count
grep -c '^ATOM\|^HETATM' my.pdb
awk '/^ATOM|^HETATM/{print substr($0,18,3)}' my.pdb | sort -u

# any missing element columns?
awk '/^ATOM|^HETATM/{e=substr($0,77,2); if(e=="  ") print NR, $0}' my.pdb

# any duplicate atom names within one residue (often breaks downstream parsers)?
awk '/^ATOM|^HETATM/{key=substr($0,22,5)"-"substr($0,13,4); print key}' my.pdb \
    | sort | uniq -c | awk '$1>1'
```

## See also

- `xyz.md`, `gjf.md` — alternative formats.
- `charge-multiplicity.md` — figuring out per-ligand charges.
- [`pdb2reaction-cli/extract.md`](../pdb2reaction-cli/extract.md) — full `extract` flag set and examples.
- `pdb2reaction-cli/{add-elem-info,fix-altloc}.md` — utility subcommands.
