# `pdb2reaction extract`

## Purpose

Cuts an active-site cluster from a PDB around a substrate selection.
Severed covalent bonds are capped with hydrogens (link-H), residue
charges are summed (`-l 'RES:Q'` for non-standard residues), and the
extracted cluster is written as a PDB ready for any other subcommand.

## Synopsis

```bash
pdb2reaction extract -i complex.pdb -c <substrate-spec> [-l 'RES:Q,...'] \
    [-r <radius_A>] [-o cluster.pdb] [--out-json]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path(s) | required | Protein–substrate complex PDB(s); multi-input requires identical atom counts |
| `-c, --center` | str | required | Substrate selector: residue-name list `'GPP,SAM'`, residue-ID list `'A:44,B:321'`, or a PDB path. Chain-qualified residue *names* (`'B:SAM'`) are not supported — use the residue ID instead. |
| `-r, --radius` | float | 2.6 | Pocket radius (Å) around `-c` atoms |
| `--radius-het2het` | float | (live default) | Separate radius for HET-to-HET inclusion |
| `-l, --ligand-charge` | str | none | Per-residue charges (amino acids derived from internal table) |
| `-o, --output` | path | (auto: `<input>_extracted.pdb`) | Output PDB path; multiple inputs produce per-file outputs |
| `--include-h2o / --no-include-h2o` | flag | (live default) | Include water residues found within radius |
| `--exclude-backbone / --no-exclude-backbone` | flag | (live default) | Trim backbone atoms outside the active site |
| `--add-linkh / --no-add-linkh` | flag | (live default) | Cap severed bonds with link hydrogens |
| `--selected-resn` | str | none | Force-include extra residue IDs (`'A:123,B:456'`); IDs only — passing residue names raises `ValueError` |
| `--modified-residue` | str | none | Comma-separated residue names (with optional charge) to **treat as amino acids** for backbone truncation and charge assignment. Examples: `'HD1,HD2,HD3'` (charge defaults to 0) or `'HD1:0,SEP:-2'`. |
| `-v, --verbose / --no-verbose` | flag | (live default) | Echo per-residue inclusion + charge sums |
| `--out-json / --no-out-json` | flag | off | Write a JSON summary alongside the PDB |

`extract` does **not** accept `-q`, `-m`, `--multi-model`,
`--freeze-links`, `--convert-files`, `--show-config`, or `--dry-run`;
those flags live on the downstream subcommands (`all`, `opt`, `tsopt`,
`dft`, `path-opt`, `path-search`, `freq`, `irc`).

## Examples

### Minimal — extract around two residues

```bash
pdb2reaction extract -i 1abc.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' -r 4.0 -o cluster.pdb
```

### Multiple structures, identical ordering

```bash
pdb2reaction extract -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -o "1.R_clu.pdb" "3.P_clu.pdb"
```

### Substrate-only PDB driving extraction

```bash
pdb2reaction extract -i complex.pdb -c substrate.pdb -r 3.5 -o cluster.pdb
```

## Output

```
cluster.pdb                  # extracted cluster, link-H caps applied
result.json (if --out-json)  # extraction stats: charges, atom counts, link-H count, file list
```

```python
import json
d = json.load(open("result.json"))
print(d["status"], d["total_charge"])
print(d["n_atoms_raw"], "->", d["n_atoms_extracted"])
print(d["n_link_hydrogens"])
print(d["files"])              # {basename: full_path} for each written PDB
print(d["protein_charge"], d["ligand_total_charge"], d["ion_total_charge"])
```

Frozen-atom indices (link-H parents) are surfaced by the *downstream*
subcommands via `freeze_atoms` in their summary.json — `extract` itself
just writes the cluster PDB.

## Caveats

- Atom names must match exactly (case-sensitive) when using `'A:44:CA'`-style
  selectors. `add-elem-info` and `fix-altloc` should run before extract
  if the PDB came out of PyMOL / Maestro.
- `-r` < 2 Å usually leaves the cluster missing essential coordinating
  atoms; 3.0–4.5 Å is typical.
- Ligand charges come **only** from `-l`; the internal table covers
  standard amino acids and a small list of common cofactors. For
  uncommon ligands, see `pdb2reaction-structure-io/charge-multiplicity.md`.

## See also

- `pdb2reaction-structure-io/pdb.md` — PDB column layout, residue selectors.
- `add-elem-info.md`, `fix-altloc.md` — pre-clean a raw PDB.
- Defaults: `import pdb2reaction.defaults as d; print(d.GEOM_KW_DEFAULT)`