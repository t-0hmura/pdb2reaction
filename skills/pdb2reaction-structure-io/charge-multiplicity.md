# Charge and multiplicity (charge-multiplicity.md)

Every run needs a total charge and a multiplicity, and a wrong value silently produces
wrong-chemistry trajectories rather than an error.

**For PDB/mmCIF input, give the charge with `-l 'RES:Q'` — the per-residue mapping — and let
`pdb2reaction` derive the total.** Name only unknown/non-standard ligand residues; standard
amino acids and recognized ions come from internal tables,
and waters and cap hydrogens are neutral. This route makes the reported total
follow the residues retained by `extract`. It is still only chemically correct
when residue names, protonation/oxidation states, ligand mappings, and the
chosen cluster boundary are correct; inspect the charge breakdown whenever the
model changes.

Use `-q INTEGER` when there are no residues to sum (notably `.xyz` without
`--ref-pdb`), or to deliberately override the derivation in ordinary geometry
commands. A valid `.gjf` already carries charge and multiplicity in its header;
CLI `-q` / `-m` override it. The `all` workflow is intentionally stricter when
`-c/--center` triggers extraction: there `-q` is an assertion and must match the
extract-derived total. Set multiplicity explicitly whenever it is not the known
singlet default.

## Multiplicity (`-m`)

| Default | Use case |
|---|---|
| **1 (singlet, closed shell)** | Use only when the modeled electron count and electronic state are known to be closed-shell. Do not infer singlet merely because the structure is biological or metal-bound. |
| 2 (doublet) | Radical species, unpaired-electron transition states (e.g. radical SAM enzymes, Fe(III) low-spin) |
| 3 (triplet) | O₂, some carbenes, Ni(II) (d⁸) high-spin tetrahedral / weak-field octahedral |
| 4 (quartet) | Co²⁺ (d⁷) high-spin, Cr³⁺ / V²⁺ (d³) |
| 5 (quintet) | Mn(III), Fe(II) high-spin |
| 6 (sextet) | Mn(II) high-spin, S=5/2 ferric |

> The entries are examples, not a spin-state calculator. For metals, radicals,
> antiferromagnetically coupled centers, or uncertain protonation/oxidation
> states, derive total charge and multiplicity for the **whole retained
> cluster** from the mechanism and primary literature.

## Charge — derive it with `-l 'RES:Q'`

1. **Per-residue mapping (the route to use for PDB)** — pass `-l 'RES1:Q1,RES2:Q2,...'`
   and let `pdb2reaction` sum amino-acid + ion + ligand charges into the total.
   The mapping is honored whether or not extraction runs: with `-c` it feeds the
   extractor's charge summary; with `-c` omitted (a pre-carved model passed as-is)
   the same mapping is applied to the full input PDB/mmCIF structure.
2. **Direct total / override** — pass `-q INTEGER` for an input with no residues to
   sum (`.xyz` without `--ref-pdb`), to replace a GJF header deliberately, or
   to override the `-l` derivation in ordinary geometry commands. In those
   commands `-q` wins over `-l` when both are given. Exception: `all` with
   `-c/--center` compares `-q` with the extracted total and aborts on a mismatch;
   `-q` does not silently replace extraction's charge result.

The amino-acid table is internal:

```bash
python -c "from pdb2reaction.workflows.extract import AMINO_ACIDS, ION; print(dict(AMINO_ACIDS)); print(dict(ION))"
```

(`AMINO_ACIDS` and `ION` are `Dict[str, int]` mapping residue/ion
name to formal charge.)

For unknown/non-standard ligand residues, supply `-l`. Recognized monatomic
ions use the internal `ION` table and must not be repeated in `-l`. A mapping
does not override a standard amino-acid or recognized-ion entry. For a
non-default protonation/oxidation state represented by such a resname, either
use an appropriate distinct residue name during model construction or provide
the verified total with `-q` as a deliberate override.

## Lookup workflow for an unfamiliar substrate

When you don't know a ligand's formal charge:

- **Lookup**: primary mechanism paper → PubChem / ChEBI `Formal Charge` → RCSB ligand summary (e.g. `https://www.rcsb.org/ligand/SAM`).
- **Derive from SMILES** if needed: `sum(a.GetFormalCharge() for a in Chem.MolFromSmiles(smi).GetAtoms())`.
- **Sanity-check before a long job**: use `--dry-run` on calculation
  subcommands that expose it; `extract` itself has no dry-run mode. It exits
  before MLIP/DFT stages. `all -c/--center ... --dry-run` deliberately performs
  extraction in a temporary directory so it can print and validate the
  extracted model's charge and electron parity; the temporary directory is
  removed afterward. `--show-config` prints the resolved configuration but then
  proceeds with the full run, so it is not a preview. `pdb2reaction extract ...
  --verbose 2` prints the per-residue charge sum used for `cluster.pdb`.

### Protonation state at physiological pH

Many ligands have multiple protonation states. Common rule of thumb:

| Group | At pH 7 | Typical formal charge contribution |
|---|---|---|
| Carboxylate (`-COO⁻`) | deprotonated | −1 each |
| Phosphate, monoester | mostly `-OPO₃²⁻` | −2 |
| Phosphate, diester | mostly `-OPO₂⁻` | −1 |
| Triphosphate (e.g. ATP) | fully deprotonated | −4 |
| Sulfonium (e.g. SAM cofactor) | quaternary | +1 |
| Lysine / Arginine side chain | protonated | +1 |
| Aspartate / Glutamate side chain | deprotonated | −1 |
| Histidine | mostly neutral, possibly +1 | 0 or +1 |

Check the literature for the cluster you are modeling — biological
mechanisms sometimes invoke an unusual protonation state.

## Source policy for an unknown value

For an unknown charge or multiplicity, use authoritative sources in this order:
the mechanism's primary paper or deposited structure documentation, then
PubChem/ChEBI/RCSB CCD. Cite the source and state the modeled protonation and
oxidation state. If the sources do not determine one unambiguous state, stop
and ask the user; never silently default a metal/radical cluster to
`-q 0 -m 1`.

## Quick-reference ligand charges (commonly seen)

These are common examples, not defaults. Confirm the exact deposited/modeled
protonation state against the relevant mechanism before using a value.

| Ligand | Resname (PDB) | Charge at pH 7 |
|---|---|---|
| Methionine sulfonium (SAM) | `SAM` | +1 |
| Adenosylhomocysteine | `SAH` | 0 |
| Geranyl pyrophosphate | `GPP` | −3 |
| ATP | `ATP` | −4 |
| ADP | `ADP` | −3 |
| GTP | `GTP` | −4 |
| NADH | `NAI`/`NDH` | −2 |
| NAD⁺ | `NAD` | −1 |
| FAD | `FAD` | −2 |
| Pyridoxal phosphate (PLP) | `PLP` | −2 |
| Heme | `HEM` | Derive from the modeled oxidation state, axial ligands, and propionate protonation |
| Phosphate ion (free) | `PO4` | −2 to −3 |

### Monatomic ions are summed from the internal `ION` table

`-l` is applied only to residues that are in **none** of `AMINO_ACIDS`,
`ION` or the water set: an ion resname takes its charge from `ION`. A token
such as `-l 'MG:2'` or `-l 'FE:2'` therefore matches no unknown residue,
produces a warning, and is ignored. List only true ligands / non-standard
residues in `-l` — the ions are already counted in the total that `extract`
reports.

The built-in values follow the PDB CCD resnames, so the oxidation state
lives in the resname:

| Ion | Resname (PDB) | `ION` value |
|---|---|---|
| Mg²⁺ | `MG` | +2 |
| Zn²⁺ | `ZN` | +2 |
| Mn²⁺ | `MN` | +2 |
| Fe³⁺ | `FE` | +3 |
| Fe²⁺ | `FE2` | +2 |
| Cu²⁺ / Cu⁺ | `CU` / `CU1` | +2 / +1 |
| Na⁺ / K⁺ | `NA` / `K` | +1 |
| Cl⁻ | `CL` | −1 |

Dump the whole table with the `python -c` one-liner above (`FE3` is not a
key in it). When the deposited resname disagrees with the oxidation state
the mechanism requires — an `FE` record that is really Fe(II), say — set
the cluster total explicitly with `-q` in an ordinary geometry command. For
`all -c/--center`, correct the residue naming/mapping instead because `-q` is a
consistency assertion, not an override.

## Multiplicity for metals (look-up shortcuts)

| Metal | Common high-spin S | Multiplicity (2S+1) |
|---|---|---|
| Mn²⁺ (d⁵) | 5/2 | 6 |
| Fe²⁺ (d⁶), high-spin tetrahedral / weak field | 2 | 5 |
| Fe²⁺ (d⁶), low-spin octahedral / strong field | 0 | 1 |
| Fe³⁺ (d⁵) high-spin | 5/2 | 6 |
| Co²⁺ (d⁷) high-spin | 3/2 | 4 |
| Cu²⁺ (d⁹) | 1/2 | 2 |
| Zn²⁺ (d¹⁰) | 0 | 1 |

## See also

- `pdb.md` — `-l 'RES:Q'` syntax and where it parses from.
- `xyz.md` — XYZ has no header, so `-q`/`-m` must be on the CLI.
- `gjf.md` — gjf encodes charge / spin in the header.
- [`pdb2reaction-cli/extract.md`](../pdb2reaction-cli/extract.md) — the subcommand that consumes
  `-l` and `-q` first.
