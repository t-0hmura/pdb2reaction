# `pdb2reaction extract`

## Purpose

Cuts an active-site cluster from a PDB/mmCIF around a substrate selection.
Supported carbon-boundary truncations are capped with hydrogens (cap-H); the
extractor does not generically cap every severed bond type. Residue
charges are summed (`-l 'RES:Q'` for non-standard residues). The internal
cluster is written as PDB; mmCIF/oversized-PDB inputs also emit
a CIF companion with original identifiers.

## Synopsis

```bash
pdb2reaction extract -i complex.pdb -c <substrate-spec> [-l 'RES:Q,...'] \
    [-r <radius_A>] [-o cluster.pdb] [--out-json]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path(s) | required | Protein–substrate complex PDB/mmCIF file(s); multi-input requires identical atom identities and order (and therefore identical counts) |
| `-c, --center` | str | required | Residue names (`'GPP,SAM'`), IDs (`'A:44,B:321'`), chain-qualified names (`'B:SAM'`), exact named ID (`'B:SAM:321'`), or a PDB/mmCIF path |
| `-r, --radius` | float | 2.6 | Pocket radius (Å) around `-c` atoms |
| `--radius-het2het` | float | `0` | Separate radius for HET-to-HET inclusion (`0` disables) |
| `-l, --ligand-charge` | str | none | Charge for unknown ligands: total (`-3`) or per-resname mapping (`'GPP:-3,SAM:1'`); amino acids and ions use internal tables |
| `-o, --output` | path | `model.pdb` (single input); `model_<filename>.pdb` (multi) | Internal/output PDB path; bridge inputs also emit the same stem as `.cif` |
| `--include-h2o / --no-include-h2o` | flag | `--include-h2o` | Include water residues found within radius |
| `--exclude-backbone / --no-exclude-backbone` | flag | `--no-exclude-backbone` | Trim backbone atoms outside the active site |
| `--add-linkh / --no-add-linkh` | flag | `--add-linkh` | Cap severed bonds with cap hydrogens |
| `--selected-resn` | str | none | Force-include extra residue IDs or names (`'A:123,B:456'`, `'A:SAM'`) |
| `--modified-residue` | str | none | Comma-separated residue names (with optional charge) to **treat as amino acids** for backbone truncation and charge assignment. Examples: `'HD1,HD2,HD3'` (charge defaults to 0) or `'HD1:0,SEP:-2'`. |
| `--out-json / --no-out-json` | flag | off | Write a JSON summary alongside the PDB |

`extract` does **not** accept `-q`, `-m`, `--freeze-links`, `--config`,
`--convert-files`, `--show-config`, or `--dry-run`;
those flags live on the downstream subcommands (`all`, `opt`, `tsopt`,
`dft`, `path-opt`, `path-search`, `freq`, `irc`).

## Examples

### Minimal — extract around two residues

```bash
pdb2reaction extract -i 1abc.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' -r 4.0 --out-json -o cluster.pdb
```

`--out-json` enables the `result.json` shown below; omit it for `cluster.pdb` only.

### Multiple structures, identical ordering

```bash
pdb2reaction extract -i 1.R.pdb -i 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -o 1.R_clu.pdb -o 3.P_clu.pdb
```

### Substrate-only PDB driving extraction

```bash
pdb2reaction extract -i complex.pdb -c substrate.pdb -r 3.5 -o cluster.pdb
```

### Repeated ligand name in a large mmCIF

```bash
pdb2reaction extract -i complex.cif -c 'LONG_CHAIN:SAM:10001' \
    -l 'SAM:1' -o cluster.pdb
# writes cluster.pdb for internal chaining and cluster.cif with original IDs
```

## Output

```
cluster.pdb                  # extracted cluster, cap-H caps applied
cluster.cif                  # bridge inputs only; original chain/residue IDs restored
result.json (if --out-json)  # extraction stats: charges, atom counts, cap-H count, file list
```

```python
import json
d = json.load(open("result.json"))
print(d["status"], d["total_charge"])
print(d["n_atoms_raw"], "->", d["n_atoms_extracted"])
print(d["n_link_hydrogens"])
print(d["files"])              # {basename: full_path} for each written PDB/CIF
print(d["protein_charge"], d["ligand_total_charge"], d["ion_total_charge"])
```

Frozen-atom indices (cap-H parents) are reported by the *downstream*
subcommands via `freeze_atoms` in their summary.json — `extract` itself
just writes the cluster PDB.

## Caveats

- When `-c` is a PDB/mmCIF path, atom names must match exactly
  (case-sensitive) between the substrate and complex for the
  `is_exact_match` coordinate check. Run `add-elem-info` for missing/wrong PDB
  element fields. The common bridge chooses one coherent altloc per residue;
  use standalone `fix-altloc` only when a cleaned PDB itself is required.
- There is no universally safe radius. Inspect retained hydrogen-bonding,
  catalytic, metal-coordination, and charge-compensating residues, then run a
  sensitivity check when the boundary could affect the mechanism. The 2.6 Å
  software default is a starting value, not a chemically validated cutoff.
- Unknown ligand charges come **only** from `-l`; the internal tables cover
  standard/recognized modified amino acids and ions, not a general set of
  cofactors. A mapping name that matches no unknown selected residue produces
  a warning and is ignored. For ligand and ion details, see
  [`pdb2reaction-structure-io/charge-multiplicity.md`](../pdb2reaction-structure-io/charge-multiplicity.md).

## See also

- [`pdb2reaction-structure-io/pdb.md`](../pdb2reaction-structure-io/pdb.md) — PDB column layout, residue selectors.
- `add-elem-info.md`, `fix-altloc.md` — pre-clean a raw PDB.
- Defaults: `import pdb2reaction.core.defaults as d; print(d.GEOM_KW_DEFAULT)`
