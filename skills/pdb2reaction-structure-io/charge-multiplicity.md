# Charge and multiplicity (charge-multiplicity.md)

`pdb2reaction` needs total charge (`-q`) and multiplicity (`-m`) as
integers. Wrong values silently produce wrong-chemistry trajectories.

## Multiplicity (`-m`)

| Default | Use case |
|---|---|
| **1 (singlet, closed shell)** | The default for almost every organic / biological / metal-coordination system whose textbook description is closed-shell. |
| 2 (doublet) | Radical species, unpaired-electron transition states (e.g. radical SAM enzymes, Fe(III) low-spin) |
| 3 (triplet) | O₂, some carbenes, Ni(II) (d⁸) high-spin tetrahedral / weak-field octahedral |
| 4 (quartet) | Co²⁺ (d⁷) high-spin, Cr³⁺ / V²⁺ (d³) |
| 5 (quintet) | Mn(III), Fe(II) high-spin |
| 6 (sextet) | Mn(II) high-spin, S=5/2 ferric |

> If the system contains a known paramagnetic metal, look up the
> oxidation state and use the high-spin/low-spin assignment from the
> primary literature for that enzyme.

## Charge (`-q`, or summed via `-l 'RES:Q'`)

**For PDB inputs, prefer `-l`**: give only the non-standard-residue charges
and let the total be auto-derived (standard AAs from the internal table +
ions + your ligand charges; waters / cap-H are neutral). It matches
`extract`'s reported total and stays correct when the cluster (radius / rep)
changes — so you never hand-enter a per-pocket total. Reserve `-q` for
`.xyz` / `.gjf` inputs (no residues to sum) or to deliberately override.

1. **Per-residue mapping (recommended for PDB)** — pass `-l 'RES1:Q1,RES2:Q2,...'`
   and let `pdb2reaction` sum amino-acid + ligand charges.
2. **Direct total / override** — pass `-q INTEGER` if you already know the
   total charge of the cluster (or to override the `-l` derivation).

The amino-acid table is internal:

```bash
python -c "from pdb2reaction.workflows.extract import AMINO_ACIDS, ION; print(dict(AMINO_ACIDS)); print(dict(ION))"
```

(`AMINO_ACIDS` and `ION` are `Dict[str, int]` mapping residue/ion
name to formal charge.)

For non-standard residues / ligands / metals, you must supply `-l`.

## Lookup workflow for an unfamiliar substrate

When you don't know a ligand's formal charge:

- **Lookup**: primary mechanism paper → PubChem / ChEBI `Formal Charge` → RCSB ligand summary (e.g. `https://www.rcsb.org/ligand/SAM`).
- **Derive from SMILES** if needed: `sum(a.GetFormalCharge() for a in Chem.MolFromSmiles(smi).GetAtoms())`.
- **Sanity-check**: `pdb2reaction extract … --verbose` echoes the per-residue charge sum used for `cluster.pdb`; `--show-config` / `--dry-run` on `all` / `opt` / `tsopt` / `dft` / `path-opt` / `path-search` / `freq` / `irc` (not `extract`) print the resolved charge before running.

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

## Permission to web-search

When the agent does not know a charge / multiplicity:

- **Confirm with the user before running a web search.** Many users
  prefer to point to the relevant paper themselves.
- If web search is allowed, prefer authoritative sources in this order:
  primary research paper → PubChem / ChEBI → general databases →
  general web. Cite the source in the agent's output.
- If neither user input nor a clean web source is available, **stop
  and ask** — do not silently default to `-q 0 -m 1` for a metal
  cluster.

## Quick-reference ligand charges (commonly seen)

Always confirm against the relevant mechanism.

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
| Mg²⁺ | `MG` | +2 |
| Zn²⁺ | `ZN` | +2 |
| Mn²⁺ / Mn³⁺ | `MN` | +2 / +3 |
| Fe²⁺ / Fe³⁺ | `FE` / `FE2` / `FE3` | +2 / +3 |
| Heme (Fe(III) protoporphyrin) | `HEM` | +1 (with Fe³⁺ + porphyrin²⁻) |
| Phosphate ion (free) | `PO4` | −2 to −3 |

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