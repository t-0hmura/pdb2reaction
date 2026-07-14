---
name: pdb2reaction-structure-io
description: PDB / mmCIF / XYZ / GJF input-file reference for pdb2reaction, including large-structure identifiers and the charge / multiplicity decision workflow. TRIGGER on editing or inspecting a structure file, deciding `-q` / `-l` / `-m`, or interpreting residue / charge / spin in an input. SKIP for subcommand syntax, output parsing, install, or HPC questions.
---

# pdb2reaction Structure I/O

## Purpose

`pdb2reaction` accepts four input formats; each carries different
information and is preferred for different stages of the workflow:

| Format | Carries | Preferred for |
|---|---|---|
| **PDB** | atom name, residue name, chain, occupancy, B-factor, element | Initial input from PDB Bank, residue-aware extraction (`-c`, `-l`) |
| **mmCIF** (`.cif` / `.mmcif`) | PDB metadata without one-character chain, four-digit residue, or five-digit atom-serial limits | Large structures, multi-character chain IDs, structures distributed as mmCIF |
| **XYZ** | element + Cartesian coordinates only | Trajectories, post-IRC outputs, when residue info is unnecessary |
| **GJF** | element + coords + charge / spin / route line | Re-running a Gaussian-style input through the MLIP pipeline |

All four formats use Å for coordinates and the conventional periodic-table
element symbols. Per-format details are in:

| File | Topic |
|---|---|
| `pdb.md` | PDB column-by-column layout, residue selectors, cap-H placement |
| `cif.md` | mmCIF conversion contract, identifier preservation, large-structure limits |
| `xyz.md` | XYZ format, ASE extension comment line |
| `gjf.md` | Gaussian gjf header (`%link0 → route → charge multiplicity → coords`) |
| `charge-multiplicity.md` | Deciding `-q` and `-m` for an unfamiliar substrate (literature lookup workflow) |

## Decision tree: which format to feed `pdb2reaction`

| Input situation | Format | How to set `-q` / `-m` |
|---|---|---|
| Fresh extraction, standard-size model | **PDB** | `-l 'RES:Q,...'` for unknown/non-standard ligand charges; standard amino acids and recognized ions use internal tables |
| PDB Bank mmCIF, multi-character chains, or ≥10,000 residues | **mmCIF** | Same residue-aware charge rules as PDB; prefer mmCIF instead of inventing over-width PDB columns |
| Single-segment optimized geometry (TS candidate, IRC endpoint) | **XYZ** | pass `-q TOTAL_CHARGE` and `-m MULT` explicitly; or use `--ref-pdb` pointing back to the original PDB so `-l` still works |
| Gaussian gjf with route line, charge, spin in header | **GJF** | `pdb2reaction` parses the header automatically; `-q` / `-m` inferred unless you override |

## Editing approach (agent-side)

When an agent must edit a structure file, the basic approach is:

1. **Read the file first** to understand current layout (residues, atom
   counts, charge / multiplicity if present).
2. **Identify the change** and confirm it does not violate format
   conventions (e.g. PDB column widths, XYZ first-line atom count).
3. For unknown charge / multiplicity values, **confirm with the user
   or do a literature lookup** before guessing — see
   `charge-multiplicity.md` for the workflow.
4. Make the smallest possible edit (single residue rename, single
   charge change). Avoid wholesale rewrites.

## Subcommand × format compatibility

| Subcommand | PDB | mmCIF | XYZ | GJF |
|---|---|---|---|---|
| `extract` | ✓ | ✓ | — | — |
| `path-search` | ✓ | ✓ | ✓ | ✓ |
| `path-opt` | ✓ | ✓ | ✓ | ✓ |
| `opt` | ✓ | ✓ | ✓ | ✓ |
| `tsopt` | ✓ | ✓ | ✓ | ✓ |
| `freq` | ✓ | ✓ | ✓ | ✓ |
| `irc` | ✓ | ✓ | ✓ | ✓ |
| `sp` | ✓ | ✓ | ✓ | ✓ |
| `dft` | ✓ | ✓ | ✓ | ✓ |
| `scan`, `scan2d`, `scan3d` | ✓ | ✓ | ✓ | ✓ |
| `all` | ✓ | ✓ | ✓ (no residue-aware extraction) | ✓ |
| `bond-summary` | ✓ | ✓ | ✓ | ✓ |

PDB-utility subcommands (`fix-altloc`, `add-elem-info`) take PDB only;
`trj2fig` takes trajectory XYZ; `energy-diagram` takes no structure.

mmCIF and a PDB at/above the fixed-column size boundary are converted to a
temporary, safely reindexed PDB before calculation. The original chain,
residue, insertion-code, residue-name, and atom-name metadata remain attached
to the topology. With `--convert-files` enabled, coordinate outputs include a
`.cif` companion with the original identifiers restored. Internal `.pdb` files
can remain as workflow intermediates; do not use their synthetic chain/residue
IDs for scientific reporting.

If you pass an XYZ to a subcommand that needs residue context (e.g.
`-l 'GLU:-1'`), supply `--ref-pdb <path>` so the residue mapping can be
recovered. `sp` is the exception: it reads PDB / mmCIF / XYZ / GJF but exposes no
`--ref-pdb`, so feed it a PDB/mmCIF directly when you need residue context for
`-l`, or pass explicit `-q` / `-m` with an XYZ.

## Quick reference: which fields where

PDB ATOM / HETATM record (column-positions, 1-indexed):

| Cols | Field |
|---|---|
| 13–16 | atom name |
| 17 | altLoc |
| 18–20 | resName |
| 22 | chainID |
| 23–26 | resSeq |
| 31–38 / 39–46 / 47–54 | X / Y / Z |
| 55–60 | occupancy |
| 61–66 | B-factor |
| 77–78 | element |

XYZ:

| Line | Content |
|---|---|
| 1 | `<natoms>` |
| 2 | comment, optional ASE `Properties=...` |
| 3+ | `<element>  <x>  <y>  <z>` |

GJF (top-to-bottom block order):

| Block | Content |
|---|---|
| Link0 | `%nproc=...`, `%mem=...` |
| Route | `# <functional/basis  options>` |
| Title | `<title>` |
| Charge/Multiplicity | `<charge> <multiplicity>` (multiplicity = 2S+1) |
| Coords | `<element>  <x>  <y>  <z>` … |
| Optional | connectivity / ECP blocks |

Full detail: see `pdb.md`, `cif.md`, `xyz.md`, and `gjf.md`.

## Charge / multiplicity defaults

- `-m 1` is the software default, but it is chemically valid only for a
  closed-shell system. Confirm metal oxidation/spin states and radicals.
- Use `-m 2` for a verified doublet radical. For a high-spin state, pass the
  verified integer multiplicity at least 3 (for example `-m 3` or `-m 5`);
  `-m` accepts an integer, not the literal text `3+`.
- `-q` (total charge) should be explicitly given for XYZ inputs unless a
  matching PDB is supplied with `--ref-pdb`; then `-l 'RES:Q'` can derive it
  from PDB/mmCIF residue metadata. The reference must have the same atom count.
- A valid GJF supplies charge and multiplicity in its header; CLI `-q` / `-m`
  override those values when the modeled state deliberately differs.

If you're not sure about charge or spin, do **not** guess silently —
follow `charge-multiplicity.md`.

## See also

- [`pdb2reaction-cli/extract.md`](../pdb2reaction-cli/extract.md) — residue selectors and cap-H caps.
- `pdb2reaction-cli/SKILL.md` — common flag conventions across
  subcommands.
- `pdb2reaction-workflows-output/SKILL.md` — what comes out of the
  pipeline (XYZ / PDB / mmCIF).
