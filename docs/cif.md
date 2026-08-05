# mmCIF and large structures

`pdb2reaction` accepts `.cif` and `.mmcif` anywhere a coordinate workflow
accepts PDB. Use mmCIF for multi-character chain IDs, residue numbers beyond
the four-column PDB field, atom counts beyond the five-column serial field, or
models with 10,000 or more residues.

## Internal conversion

The calculation code remains PDB-based. At input, `pdb2reaction` therefore:

1. reads the first mmCIF coordinate model;
2. resolves altLoc at the residue level by mean occupancy;
3. assigns safe temporary one-character chains and residue numbers 1–9999;
4. retains the original auth chain, residue number, insertion code, residue
   name, atom name, occupancy, B-factor, element, and formal charge; and
5. sends the temporary PDB to the ordinary optimizer/path code.

The same bridge is activated for a PDB with at least 10,000 residues, at least
99,999 atoms, a hybrid-36 field, or a recognized over-width decimal field.
The temporary PDB IDs are implementation details.

## Output

With `--convert-files` enabled (the default), coordinate-producing workflows
retain the PDB needed between pipeline stages and add a CIF companion with the
original identifiers:

`extract` is an exception: it has no conversion toggle and automatically writes
the retained-template CIF companion for a bridged input.

```text
final_geometry.pdb  # normalized representation used between pipeline stages
final_geometry.cif  # original chain/residue identity restored
```

Multi-frame trajectories use `_atom_site.pdbx_PDB_model_num`. Calculated CIF
files contain the atom-site coordinate table; unrelated crystallographic
refinement categories from the source file are not copied.

## Residue and atom selectors

```bash
# all SAM residues in auth chain LONG_CHAIN
pdb2reaction extract -i complex.cif -c 'LONG_CHAIN:SAM' -o model.pdb

# exactly one SAM, even when its residue number exceeds the PDB field
pdb2reaction extract -i complex.cif -c 'LONG_CHAIN:SAM:10001' -o model.pdb

# numeric form remains available
pdb2reaction extract -i complex.cif -c 'LONG_CHAIN:10001' -o model.pdb
```

`CHAIN:RESNAME` deliberately selects all matches in that chain and warns when
there is more than one. Add `:RESSEQ` to select one. A chain-qualified scan atom
selector has four fields: `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`; append the
insertion code to the residue number when needed (for example,
`A:SAM:12B:C1`). Chain IDs are case-sensitive (`A` and `a` may be distinct).

## Limits

- Up to 619,938 residues can be represented by the internal bridge.
- All reaction-ordered inputs must still have identical atom identity and
  order. Conversion does not infer atom mapping.
- Only the first input coordinate model is calculated.
- `fix-altloc` and `add-elem-info` remain PDB-only utilities; mmCIF altLoc and
  `type_symbol` are handled during conversion.

See also [extract](extract.md) and [CLI conventions](cli-conventions.md).
