# PDB format (pdb.md)

PDB is `pdb2reaction`'s primary input. Column-based, fixed-width fields.
For 10,000 or more residues, atom serials beyond 99,999, residue IDs beyond
four columns, or multi-character chain IDs, use mmCIF; see `cif.md`.

## Record types `pdb2reaction` cares about

| Record | Used for |
|---|---|
| `ATOM` | Polymer atoms by PDB convention (normally amino acids or nucleic acids) |
| `HETATM` | Non-polymer atoms by PDB convention (normally ligand, metal, water, cofactor, or generated cap H) |
| `TER` | Chain terminator; `extract` infers chain breaks from C–N peptide-adjacency distance, not by parsing TER directly |
| `END` | File terminator — informational only |
| `MODEL` / `ENDMDL` | `extract` writes multi-MODEL output for multi-structure input, and `fix-altloc` processes each MODEL block. Multi-MODEL **input** is not read — `extract` uses only the first model and warns. |
| `CRYST1` | Unit-cell metadata is not used by the cluster calculations and is omitted from generated cluster PDBs |

`pdb2reaction` does not use `ANISOU`, `LINK`, `SSBOND`, and related
connectivity annotations for its calculations. Preserve the original file;
if a parser-specific cleanup is needed, write a separate derived PDB and
validate atom identity/order afterward.

## Column layout (`ATOM` / `HETATM`)

Columns are 1-based, inclusive on both ends.

| Cols | Field | Width | Type | Example |
|---|---|---|---|---|
| 1–6 | Record name | 6 | left-just `ATOM  ` / `HETATM` | `ATOM  ` |
| 7–11 | Atom serial | 5 | right-just int | `   42` |
| 13–16 | Atom name | 4 | fixed-width; PDB alignment depends on element/name, so preserve a parser's formatting | ` CB ` |
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
| 77–78 | Element symbol | 2 | right-justified canonical symbol (` C`, `Fe`, `Mg`) | ` C` |
| 79–80 | Formal charge | 2 | right-just (e.g. `2+`, `1-`) | `  ` |

The record type does not decide the charge: `pdb2reaction` uses the residue
name and its internal amino-acid/ion/water tables, plus `-l` for unknown
residues. `pdb2reaction add-elem-info` repairs columns 77-78 when they are
blank, which can happen after structure-editor export. Run it only when
elements are missing, then validate the resulting symbols before extraction.

## Residue selectors (the `-c / --center` flag)

`pdb2reaction extract` (and any subcommand that also extracts internally)
uses five forms:

```bash
# Form 1 — comma-separated residue names. Matches every residue with
# resName == SAM, GPP, or MG.
pdb2reaction extract -i complex.pdb -c 'SAM,GPP,MG' -o cluster.pdb

# Form 2 — a separate PDB containing only the substrate residues.
pdb2reaction extract -i complex.pdb -c substrate.pdb -o cluster.pdb

# Form 3 — chain-aware: `chainID:resSeq` (numeric resSeq is honored)
pdb2reaction extract -i complex.pdb -c 'A:44' -o cluster.pdb

# Form 4 — every residue with a name in one chain
pdb2reaction extract -i complex.pdb -c 'A:SAM' -o cluster.pdb

# Form 5 — chain + residue name + sequence number (one intended residue)
pdb2reaction extract -i complex.pdb -c 'A:SAM:44' -o cluster.pdb
```

`chainID:resName` selects every matching residue in that chain and warns when
multiple residues match. Add the numeric `:resSeq` component when exactly one
is intended. This is especially useful when the same ligand/cofactor name is
repeated across chains. mmCIF chain IDs may contain multiple characters.

The pocket radius around the centers is set by `-r <Å>` (default 2.6 Å).
All residues with at least one atom (any element, including H) inside the radius are kept.

## Per-residue charge mapping (`-l / --ligand-charge`)

```bash
pdb2reaction extract -i complex.pdb \
    -c 'SAM,GPP,MG' \
    -l 'SAM:1,GPP:-3' \
    -o cluster.pdb
```

Standard amino acids are looked up from `pdb2reaction`'s internal
`AMINO_ACIDS` table and recognized monatomic ions from `ION`; provide only
unknown/non-standard ligand residues in `-l`. For example, `MG` is already
`+2`; adding `MG:2` to `-l` does not override it and emits an unmatched-entry
warning. The total cluster charge is the sum of all retained residues after
extraction.

If you don't know a ligand's formal charge, see
`charge-multiplicity.md` for the lookup workflow.

## Cap-hydrogen capping

When `extract` truncates one of its supported carbon-boundary bonds, it
places a hydrogen `H_link` from the retained carbon (`A`) toward the removed
partner (`B`) at 1.09 Å (standard C-H length). It does not generically cap
arbitrary C–X/non-carbon cuts. The cap hydrogen
is written as a `HETATM` with atom name `HL` and residue name `LKH`
(hard-coded in `_format_linkH_block`).

Cap hydrogens carry **no formal charge**; they do not enter the
charge sum.

For geometry commands, cap-parent atoms (atom `A`, the cluster-side parent of
each cap) are detected and frozen by default through `--freeze-links`; users
can explicitly disable that behavior. The PDB itself carries no encoded
freeze list — `pdb2reaction.core.utils.detect_freeze_links` re-derives each
parent of an `LKH/HL` record at runtime via a nearest-neighbor search in
Cartesian space. The B-factor column on LKH atoms is hard-coded to 0.00 by
`_format_linkH_block` and is not a freeze flag.

## Manual cluster-boundary checklist

- Make a retained backbone fragment terminate consistently at alpha carbons
  (`CA`) at both main-chain ends, then cap the open valences.
- For side-chain/ligand/cofactor truncation, choose an aliphatic C–C single
  bond whenever possible (`CA–CB` or farther from the reactive center). Do not
  cut peptide C–N, polar C–N/C–O, aromatic/conjugated, disulfide, or
  metal-coordination bonds; move the boundary or include the bonded partner.
- One intended cap per severed valence; visually inspect bond lengths and
  valence. Keep cap parents frozen with `--freeze-links` in production.
- Within one R/IM/P reaction path, require exactly the same atom identities,
  order, cluster boundary, and cap topology. Prefer one multi-input extraction;
  if states must be prepared separately, explicitly harmonize and validate all
  four properties before path calculations.
- WT/mutant or other cross-variant models may necessarily differ in atom
  identity. Keep the modeling protocol and boundary definition comparable,
  document the differences, and compare reaction observables such as barriers
  rather than trying to combine unlike structures into one path.
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
# atom count + unique residue-identity count (chain, residue number, insertion code)
grep -c '^ATOM\|^HETATM' my.pdb
awk '/^ATOM|^HETATM/{print substr($0,22,6)}' my.pdb | sort -u | wc -l

# any missing element columns?
awk '/^ATOM|^HETATM/{e=substr($0,77,2); if(e ~ /^[[:space:]]*$/) print NR, $0}' my.pdb

# any duplicate atom names within one residue (often breaks downstream parsers)?
awk '/^ATOM|^HETATM/{key=substr($0,22,6)"-"substr($0,13,4); print key}' my.pdb \
    | sort | uniq -c | awk '$1>1'
```

## See also

- `xyz.md`, `gjf.md` — alternative formats.
- `charge-multiplicity.md` — figuring out per-ligand charges.
- [`pdb2reaction-cli/extract.md`](../pdb2reaction-cli/extract.md) — full `extract` flag set and examples.
- `pdb2reaction-cli/{add-elem-info,fix-altloc}.md` — utility subcommands.
