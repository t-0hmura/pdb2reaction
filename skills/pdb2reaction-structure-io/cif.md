# mmCIF input and large-structure bridge

Use `.cif` or `.mmcif` when chain IDs exceed one character, residue identifiers
exceed four PDB columns, atom serials exceed five PDB columns, or the model has
10,000 or more residues. Do not make a nominal “PDB” by widening fixed fields;
fixed-column readers will interpret shifted atom names, residue IDs, or
coordinates incorrectly.

## Runtime contract

1. `pdb2reaction` reads the first `_atom_site.pdbx_PDB_model_num` model.
2. It chooses one coherent nonblank altLoc label per residue using mean
   occupancy, while retaining blank/shared atoms.
3. It assigns safe one-character internal chains and residue numbers 1–9999
   across as many internal chains as necessary.
4. Calculation stages consume that temporary PDB.
5. With `--convert-files` enabled (the default), coordinate outputs retain the
   internal PDB needed by composite workflows and also write a `.cif`
   companion. The CIF restores `auth_asym_id`,
   `auth_seq_id`, insertion code, residue name, and atom name from the input.

The bridge supports up to 62 × 9999 = 619,938 residues. Atom serials in the
internal PDB may repeat after 99,999 because decimal PDB has no larger field;
atom order, not serial, is authoritative internally. Output CIF `_atom_site.id`
is sequential and unbounded by the PDB field width.

## Residue selection

All forms below work for PDB and mmCIF:

```bash
# Every SAM in every chain
pdb2reaction extract -i complex.cif -c 'SAM' -o model.pdb

# Every SAM in auth chain LONG_CHAIN
pdb2reaction extract -i complex.cif -c 'LONG_CHAIN:SAM' -o model.pdb

# One SAM, including a residue number beyond the PDB limit
pdb2reaction extract -i complex.cif -c 'LONG_CHAIN:SAM:10001' -o model.pdb

# Numeric ID form remains valid
pdb2reaction extract -i complex.cif -c 'LONG_CHAIN:10001' -o model.pdb
```

`CHAIN:RESNAME` intentionally selects all matching residues in that chain and
warns when there is more than one. Add `:RESSEQ` when one residue is intended.
For scan atom selectors, use
`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`, for example
`LONG_CHAIN:SAM:10001:C1`, to disambiguate repeated numbering. Append an
insertion code to `RESSEQ` when needed, as in `A:SAM:12B:C1`.
Treat chain IDs as case-sensitive: `A` and `a` can identify different chains.

## Outputs

Given `complex.cif`, a coordinate-producing geometry workflow with conversion
enabled may contain both:

```text
final_geometry.pdb  # normalized topology used between pipeline stages
final_geometry.cif  # public coordinate output with original identifiers
```

The PDB's synthetic identifiers are implementation details. Use the CIF for
downstream analysis and provenance when the input required the bridge.

## Limitations and checks

- Only the first coordinate model is used as an input geometry. Multi-frame
  calculated trajectories are written as multiple CIF model numbers.
- Non-`_atom_site` categories from the source CIF are not copied to calculated
  outputs. Coordinate identity, occupancy, B-factor, formal charge, and atom
  labels are retained; crystallographic refinement metadata is not.
- All R/IM/P inputs still require identical atom identity and order. Format
  conversion does not repair atom mapping.
- `fix-altloc` and `add-elem-info` remain PDB utilities. CIF altLoc and
  `type_symbol` are handled by the bridge.
