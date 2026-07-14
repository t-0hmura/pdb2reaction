# `extract`

`pdb2reaction extract` carves an active-site cluster from a protein–ligand
PDB/mmCIF. Select by residue name (`GPP,SAM`), ID (`A:123A`), chain-qualified
name (`A:SAM`), exact named ID (`A:SAM:123`), or a PDB/mmCIF file. mmCIF and
PDBs at the fixed-column size boundary are safely reindexed internally and
also produce CIF output with the original identifiers.

## Examples

Command form:

```bash
pdb2reaction extract -i COMPLEX.pdb [COMPLEX2.pdb ...] \
    -c SUBSTRATE_SPEC \
    [-o MODEL.pdb [MODEL2.pdb ...]] \
    [--radius Å] [--radius-het2het Å] \
    [--include-h2o / --no-include-h2o] \
    [--exclude-backbone / --no-exclude-backbone] \
    [--add-linkh / --no-add-linkh] \
    [--selected-resn LIST] [--modified-residue LIST] \
    [-l, --ligand-charge MAP_OR_NUMBER] \
    [--out-json / --no-out-json] \
    [-v LEVEL]
```

Minimal (ID-based substrate) with explicit total ligand charge:

```bash
# Minimal (ID-based substrate) with explicit total ligand charge
pdb2reaction extract -i complex.pdb -c '123' -o model.pdb -l -3
```

Substrate provided as a PDB; per-resname charge mapping (others remain 0):

```bash
# Substrate provided as a PDB; per-resname charge mapping (others remain 0)
pdb2reaction extract -i complex.pdb -c substrate.pdb -o model.pdb -l 'GPP:-3,SAM:1'
```

Name-based substrate selection — all matches included (WARNING logged):

```bash
# Name-based substrate selection — all matches included (WARNING logged)
pdb2reaction extract -i complex.pdb -c 'GPP,SAM' -o model.pdb -l 'GPP:-3,SAM:1'
```

Multi-structure → single multi-MODEL output with hetero-hetero proximity:

```bash
# Multi-structure → single multi-MODEL output with hetero-hetero proximity
pdb2reaction extract -i complex1.pdb -i complex2.pdb -c 'GPP,SAM' \
    -o model_multi.pdb --radius-het2het 2.6 -l 'GPP:-3,SAM:1'
# Multi-structure → multiple outputs: pass -o model1.pdb -o model2.pdb instead
```

## Workflow

### Residue inclusion

- Always include the substrate residues from `-c/--center`.
- **Standard cutoff (`--radius`, default 2.6 Å)**: with `--no-exclude-backbone` (default), any atom within the cutoff makes its residue qualify (i.e. includes the residue). With `--exclude-backbone`, amino-acid residues must contact the substrate with a **non-backbone** atom (not N / H* / CA / HA* / C / O / OXT). Non-amino acids always use any atom.
- **Independent hetero–hetero cutoff (`--radius-het2het`)**: adds residues when a substrate hetero atom (non C/H) lies within the specified Å of a protein hetero atom. With backbone exclusion enabled, the protein atom must be non-backbone.
- **Water handling**: HOH / WAT / H2O / DOD / TIP / TIP3 / SOL are included by default (`--include-h2o`).
- **Forced inclusion**: `--selected-resn` accepts the same ID/name/chain-qualified selectors as `--center`; see {ref}`selected-resn-takes-ids`.
- **Chain/name disambiguation**: `A:SAM` selects every SAM in chain A and
  warns on multiple matches; `A:SAM:123` selects one intended residue.
- **Neighbor safeguards**:
  - When backbone exclusion is off and a residue contacts the substrate with a backbone atom, the peptide-adjacent N / C neighbors (C–N ≤ 1.9 Å) are auto-included; termini keep caps (N/H* or C/O/OXT).
  - Disulfide bonds (SG–SG ≤ 2.5 Å) bring both cysteines.
  - Non-terminal PRO residues always pull in the N-side amino acid; CA is preserved even when backbone atoms are removed, and under `--exclude-backbone` the neighbor's C / O / OXT remain to maintain the peptide bond.

### Truncation and capping

- Isolated residues retain only side-chain atoms; amino-acid backbone atoms (N, CA, C, O, OXT plus N/CA hydrogens) are removed except for PRO / HYP safeguards.
- Continuous peptide stretches keep internal backbone atoms; only terminal caps (N/H* or C/O/OXT) are removed. TER awareness prevents capping across chain breaks.
- With `--exclude-backbone`, main-chain atoms on all **non-substrate** amino acids are stripped (subject to PRO / HYP safeguards and PRO neighbor retention).
- Non-amino-acid residues never lose atoms named like backbone (N / CA / HA / H / H1 / H2 / H3).

### Cap hydrogens (`--add-linkh`)

- Cap hydrogens are placed at 1.09 Å along severed bond vectors at carbon boundaries only (CB–CA, CA–N, CA–C; PRO / HYP use CA–C only); non-carbon boundaries are not capped.
- Inserted after a `TER` as contiguous `HETATM` records named `HL` in residue `LKH` (chain `L`). Serial numbers continue from the main block.
- In multi-structure mode the same bonds are capped across all models; coordinates remain model-specific.

### Building or auditing a cluster model manually

Treat every severed bond as a chemical modeling decision, not just a geometric
radius cutoff.

- For a retained protein-backbone fragment, choose the residue span so its two
  main-chain ends terminate consistently at alpha carbons (PDB atom name
  `CA`), then satisfy the terminal valences with caps.
- At side-chain, ligand, or cofactor boundaries, cut a non-polar **C–C single
  bond** whenever possible (commonly `CA–CB` or a more distant aliphatic C–C
  bond). Do not cut a peptide C–N bond, C–O/C–N polar bond, aromatic/conjugated
  bond, disulfide, or metal-coordination bond merely because it crosses the
  distance cutoff; include the bonded partner or move the boundary.
- Verify every cluster-side boundary atom has the intended valence and exactly
  one cap. Keep cap parents frozen (`--freeze-links`, the default) in
  production optimizations.
- R/IM/P structures must contain the same atoms in the same order. Apply the
  selection/capping pattern from one structure to all states rather than
  re-extracting them independently.
- Recalculate the total charge and multiplicity after truncation and inspect
  the boundary visually before starting an MEP/Hessian calculation.

The automatic extractor covers common protein cases, but unusual covalent
cofactors, modified residues, metal sites, or a boundary that violates these
rules should be corrected manually.

### Charge summary (`--ligand-charge/-l`)

Amino acids and common ions draw charges from internal dictionaries; waters are zero. Unknown residues default to 0 unless `--ligand-charge/-l` supplies either a total charge (distributed across unknown substrate residues, or across all unknown residues when there is no unknown substrate residue) or a per-resname mapping like `GPP:-3,SAM:1`.

### Multi-structure ensembles

`extract` accepts multiple PDB/mmCIF inputs. Full per-atom identity and order
are validated; each structure is processed independently and the **union** of
selected residues is applied to every model.

| Output policy | Layout |
|---|---|
| No `-o`, multiple inputs | per-file `model_<original_basename>.pdb` |
| One `-o` path | single multi-MODEL PDB |
| N outputs matching N inputs | N individual PDBs |

Diagnostics report raw vs. kept atom counts per model along with residue IDs.

## Outputs

```text
<output>.pdb        # Active-site model PDB(s) with optional cap hydrogens after a TER record
<output>.cif        # mmCIF/oversized-PDB input: original auth chain/residue IDs restored
                    # Single input → model.pdb by default
                    # Multiple inputs without -o → model_<original_basename>.pdb per structure
                    # One -o path with multiple inputs → single multi-MODEL PDB
                    # Parent output directories are created automatically
```

Charge summary (protein / ligand / ion / total) is logged for model #1 when verbose mode is enabled. Programmatic use (`extract_api`) returns `{"outputs": [...], "counts": [...], "charge_summary": {...}, "n_link_hydrogens": N}`.

## CLI options

Defaults shown are used when the option is not specified. The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | One or more protein–ligand PDB/mmCIF files (identical atom identity/order required). | Required |
| `-c, --center SPEC` | PDB/mmCIF path, residue IDs/names, `CHAIN:RESNAME`, or `CHAIN:RESNAME:RESSEQ`. | Required |
| `-o, --output PATH...` | Active-site model PDB output(s); see [Outputs](#outputs) for naming/layout. | Auto (`model.pdb` or `model_<input>.pdb`) |
| `-r, --radius FLOAT` | Atom–atom distance cutoff (Å) for inclusion (internally `0.001 Å` when zero). | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero–hetero cutoff (Å, non C/H). | `0.0` (internally `0.001 Å` when zero) |
| `--include-h2o / --no-include-h2o` | Include HOH / WAT / H2O / DOD / TIP / TIP3 / SOL waters. | `True` |
| `--exclude-backbone / --no-exclude-backbone` | Remove backbone atoms on non-substrate amino acids (PRO / HYP safeguards). | `False` |
| `--add-linkh / --no-add-linkh` | Add cap hydrogens at 1.09 Å along severed bonds at carbon boundaries only (non-carbon boundaries are not capped). | `True` |
| `--selected-resn TEXT` | Force-include by the same ID/name/chain-qualified selectors as `--center`. | `""` |
| `--modified-residue TEXT` | Comma-separated residue names (with optional per-residue charge) to treat as amino acids for backbone truncation and charge assignment (e.g. `HD1,HD2,HD3` or `HD1:0,SEP:-2`). A residue given without a trailing `:charge` defaults to charge 0. | `""` |
| `-l, --ligand-charge TEXT` | Total charge or per-resname mapping (e.g. `GPP:-3,SAM:1`). | _None_ |
| `--out-json / --no-out-json` | Write a machine-readable `result.json` alongside the extracted PDB(s). Schema: [JSON Output Schema](json-output.md). | `False` |

### Substrate specification (`-c/--center`)

- **PDB/mmCIF path**: coordinates must match the first input exactly (tolerance 1e-3 Å); residue IDs propagate to other structures.
- **Residue IDs**: `'123,124'`, `'A:123,B:456'`, `'123A'`, `'A:123A'` (insertion codes supported).
- **Residue names**: comma-separated, case-insensitive. If multiple residues share a name, **all** matches are included and a warning is logged.
- **Chain + name**: `'A:SAM'` selects all matching residues in chain A;
  `'A:SAM:123'` narrows it to one residue.

## Notes

- Validate extraction-radius convergence and chemical completeness for the
  target system. Increasing `-r` includes more environment but does not
  guarantee monotonic accuracy and increases computational cost.
- INFO logs summarize residue selection, truncation counts, and charge breakdowns.

## Systems with non-standard residues (MCPB, etc.)

When metal-coordinating amino-acid parameters are generated by tools such as Amber's `MCPB.py` (Metal Center Parameter Builder), the coordinating residues are assigned non-standard names (e.g. `HD1`, `HE1`, `CM1`, `AP1`). These are not in `extract`'s internal `AMINO_ACIDS` dictionary, so **backbone truncation and cap-hydrogen placement will not be applied correctly**, and a warning is emitted:

```text
[extract] WARNING: Residue HD1 83 may be an amino acid (has N, CA, C, O)
but is not recognized as a standard residue name.
Backbone truncation was not applied.
Consider preparing the active site model manually.
```

Use `--modified-residue` to register non-standard residue names as amino acids so backbone truncation and charge assignment apply automatically — useful for phosphoserine, methylated residues, D-amino acids with unusual names, and MCPB-renamed metal-coordinating residues:

```bash
# Treat HD1, HD2, HD3 as amino acids (charge defaults to 0)
pdb2reaction extract -i complex.pdb -c 'SUB' -o model.pdb \
    --modified-residue 'HD1,HD2,HD3'

# Specify explicit charges per modified residue
pdb2reaction extract -i complex.pdb -c 'SUB' -o model.pdb \
    --modified-residue 'HD1:0,SEP:-2'
```

```{important}
If `--modified-residue` does not cover your use case, **manual active-site model construction is recommended**:

1. Select residues around the active site and determine truncation points.
2. Add a cap hydrogen on the parent atom (the atom that remains) of each severed covalent bond.
3. Use residue name `LKH` (chain `L`) and atom name `HL` for the cap hydrogen.
4. Place it at **1.09 Å** along the original bond direction.
```

## Appendix: PDB naming requirements and reference lists

This appendix exists mainly for debugging cases where `extract` misclassifies residues due to **non-standard residue or atom naming**. If your inputs follow standard PDB conventions, you can usually skip it.

```{important}
For `extract` to work correctly, **residue and atom names in PDB/mmCIF must follow the expected PDB chemical-component naming conventions**. The tool relies on internal dictionaries to recognize amino acids, ions, waters, and backbone atoms. Non-standard naming can misclassify residues or assign incorrect charges.
```

### `AMINO_ACIDS`

A dictionary mapping residue names to their nominal integer charges. Membership determines whether a residue is treated as an amino acid for backbone handling, truncation, and charge calculation.

**Standard 20** (charges reflect physiological pH):

- Neutral: `ALA`, `ASN`, `CYS`, `GLN`, `GLY`, `HIS`, `ILE`, `LEU`, `MET`, `PHE`, `PRO`, `SER`, `THR`, `TRP`, `TYR`, `VAL`
- Positive (+1): `ARG`, `LYS`
- Negative (−1): `ASP`, `GLU`

**Canonical extras:** `SEC` (selenocysteine, 0), `PYL` (pyrrolysine, +1).

**Protonation / tautomer variants** (Amber / CHARMM style): `HIP` (+1, fully protonated His), `HID` (0, Nδ-protonated His), `HIE` (0, Nε-protonated His), `ASH` (0, neutral Asp), `GLH` (0, neutral Glu), `LYN` (0, neutral Lys), `ARN` (0, neutral Arg), `TYM` (−1, deprotonated Tyr phenolate).

**Phosphorylated:** dianionic (−2) `SEP`, `TPO`, `PTR`; monoanionic (−1) `S1P`, `T1P`, `Y1P`; phospho-His (phosaa19SB) `H1D` (0), `H2D` (−1), `H1E` (0), `H2E` (−1).

**Cysteine variants:** `CYX` (0, disulfide), `CSO` (0, sulfenic acid), `CSD` (−1, sulfinic acid), `CSX` (0, generic), `OCS` (−1, cysteic acid), `CYM` (−1, deprotonated Cys).

**Lysine variants / carboxylation:** `MLY` (+1), `LLP` (+1), `KCX` (−1, Nz-carboxylic acid).

**D-amino acids** (19): `DAL`, `DAR`, `DSG`, `DAS`, `DCY`, `DGN`, `DGL`, `DHI`, `DIL`, `DLE`, `DLY`, `MED`, `DPN`, `DPR`, `DSN`, `DTH`, `DTR`, `DTY`, `DVA`.

**Other modified:** `CGU` (−2, γ-carboxy-glutamate), `CGA` (−1), `PCA` (0, pyroglutamate), `MSE` (0, selenomethionine), `OMT` (0, methionine sulfone), `HYP` (0, hydroxyproline); also `ASA`, `CIR`, `FOR`, `MVA`, `IIL`, `AIB`, `HTN`, `SAR`, `NMC`, `PFF`, `NFA`, `ALY`, `AZF`, `CNX`, `CYF`.

**N-terminal variants** (`N` prefix): `NALA` (+1), `NARG` (+2), `NASP` (0), `NGLU` (0), `NLYS` (+2), … plus `ACE` (0), `NTER` (+1, generic).
**C-terminal variants** (`C` prefix): `CALA` (−1), `CARG` (0), `CASP` (−2), `CGLU` (−2), `CLYS` (0), … plus `NHE` (0), `NME` (0), `CTER` (−1, generic).

### `BACKBONE_ATOMS`

Atom names treated as backbone for amino acids; under `--exclude-backbone` these are removed from non-substrate residues:

```
N, C, O, CA, OXT, H, H1, H2, H3, HN, HA, HA2, HA3
```

### `ION`

Recognized ion residue names with formal charges:

| Charge | Residue names |
|---|---|
| +1 | `LI`, `NA`, `K`, `RB`, `CS`, `TL`, `AG`, `CU1`, `K+`, `NA+`, `NH4`, `H3O+`, `HE+`, `HZ+` |
| +2 | `MG`, `CA`, `SR`, `BA`, `MN`, `FE2`, `CO`, `NI`, `CU`, `ZN`, `CD`, `HG`, `PB`, `BE`, `PD`, `PT`, `SN`, `RA`, `YB2`, `V2+` |
| +3 | `FE`, `AU3`, `AL`, `GA`, `IN`, `CE`, `CR`, `DY`, `EU`, `EU3`, `ER`, `GD3`, `LA`, `LU`, `ND`, `PR`, `SM`, `TB`, `TM`, `Y`, `PU` |
| +4 | `U4+`, `TH`, `HF`, `ZR` |
| −1 | `F`, `CL`, `BR`, `I`, `CL-`, `IOD` |

### `WATER_RES`

Recognized water residue names (included by default with `--include-h2o`, assigned zero charge):

```
HOH, WAT, H2O, DOD, TIP, TIP3, SOL
```

(link-hydrogen-and-frozen-atoms)=

## Cap hydrogen and frozen atoms

When `pdb2reaction` extracts an active-site model from a larger structure, severed bonds are capped with **cap hydrogens**. By default (`--freeze-links`), the parent atoms of these cap hydrogens are frozen during optimization and path searches to prevent unphysical rearrangement at the boundary.

- **Forces** — frozen atoms receive zeroed forces.
- **Hessian** — frozen degrees of freedom are either removed (`return_partial_hessian: true`) or zeroed in the full matrix.
- **Vibrational analysis** — when frozen atoms are present, `freq` automatically performs Partial Hessian Vibrational Analysis (PHVA), diagonalizing only the active block of the Hessian.

Frozen atoms can also be set manually via the `geom.freeze_atoms` YAML key (1-based indices). CLI-detected cap atoms are merged with YAML-specified atoms.

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [all](all.md) — End-to-end workflow that calls extract internally via `-c/--center`
- [path-search](path-search.md) — MEP search on extracted active-site models
- [scan](scan.md) — Staged scan on extracted active site models
- [add-elem-info](add-elem-info.md) — Fix missing PDB element columns before extraction
- [Troubleshooting](troubleshooting.md) — Common extraction errors
- [Glossary](glossary.md) — Definitions of Active Site Model, Cluster Model, Cap Hydrogen
