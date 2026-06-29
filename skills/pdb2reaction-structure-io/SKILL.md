---
name: pdb2reaction-structure-io
description: PDB / XYZ / GJF input-file reference for pdb2reaction, plus the charge / multiplicity decision workflow for arbitrary substrates. TRIGGER on editing or inspecting a structure file, deciding `-q` / `-l` / `-m`, or interpreting residue / charge / spin in an input. SKIP for subcommand syntax, output parsing, install, or HPC questions.
---

# pdb2reaction Structure I/O

## Purpose

`pdb2reaction` accepts three input formats; each carries different
information and is preferred for different stages of the workflow:

| Format | Carries | Preferred for |
|---|---|---|
| **PDB** | atom name, residue name, chain, occupancy, B-factor, element | Initial input from PDB Bank, residue-aware extraction (`-c`, `-l`) |
| **XYZ** | element + Cartesian coordinates only | Trajectories, post-IRC outputs, when residue info is unnecessary |
| **GJF** | element + coords + charge / spin / route line | Re-running a Gaussian-style input through the MLIP pipeline |

All three formats use Å for coordinates and the conventional periodic-table
element symbols. Per-format details are in:

| File | Topic |
|---|---|
| `pdb.md` | PDB column-by-column layout, residue selectors, cap-H placement |
| `xyz.md` | XYZ format, ASE extension comment line |
| `gjf.md` | Gaussian gjf header (`%link0 → route → charge multiplicity → coords`) |
| `charge-multiplicity.md` | Deciding `-q` and `-m` for an unfamiliar substrate (literature lookup workflow) |

## Decision tree: which format to feed `pdb2reaction`

| Input situation | Format | How to set `-q` / `-m` |
|---|---|---|
| Fresh extraction from PDB Bank or model from PyMOL / Maestro | **PDB** | `-l 'RES:Q,...'` for per-residue ligand charges; `pdb2reaction` reads residue names directly |
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

| Subcommand | PDB | XYZ | GJF |
|---|---|---|---|
| `extract` | ✓ (input + output) | — | — |
| `path-search` | ✓ | ✓ | ✓ |
| `path-opt` | ✓ | ✓ | ✓ |
| `opt` | ✓ | ✓ | ✓ |
| `tsopt` | ✓ | ✓ | ✓ |
| `freq` | ✓ | ✓ | ✓ |
| `irc` | ✓ | ✓ | ✓ |
| `dft` | ✓ | ✓ | ✓ |
| `scan`, `scan2d`, `scan3d` | ✓ | ✓ | ✓ |
| `all` | ✓ | ✓ (single segment) | ✓ |
| `bond-summary` | ✓ | ✓ | ✓ |

PDB-utility subcommands (`fix-altloc`, `add-elem-info`) take PDB only;
`trj2fig` takes trajectory XYZ; `energy-diagram` takes no structure.

If you pass an XYZ to a subcommand that needs residue context (e.g.
`-l 'GLU:-1'`), supply `--ref-pdb <path>` so the residue mapping can be
recovered.

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

Full byte-by-byte / per-keyword detail: see `pdb.md`, `xyz.md`, `gjf.md`.

## Charge / multiplicity defaults

- `-m 1` (singlet, closed shell) is the default for almost every
  organic / biological cluster.
- Use `-m 2` for radicals, `-m 3+` for unusual high-spin metal centers.
- `-q` (total charge) must be explicitly given for XYZ inputs (XYZ has
  no header). `-l 'RES:Q'` derives `-q` for PDB input (or XYZ with
  `--ref-pdb`) from per-residue charges plus `pdb2reaction`'s amino-acid
  table.

If you're not sure about charge or spin, do **not** guess silently —
follow `charge-multiplicity.md`.

## See also

- [`pdb2reaction-cli/extract.md`](../pdb2reaction-cli/extract.md) — residue selectors and cap-H caps.
- `pdb2reaction-cli/SKILL.md` — common flag conventions across
  subcommands.
- `pdb2reaction-workflows-output/SKILL.md` — what comes out of the
  pipeline (also XYZ / PDB).