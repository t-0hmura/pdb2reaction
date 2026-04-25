---
name: pdb2reaction-structure-io
description: How to read, edit, and write the three molecular-structure formats pdb2reaction handles — PDB, XYZ, GJF — plus the charge / multiplicity decision workflow for arbitrary substrates.
---

# pdb2reaction Structure I/O

## Overview

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
| `pdb.md` | PDB column-by-column layout, residue selectors, link-H placement |
| `xyz.md` | XYZ format, ASE extension comment line |
| `gjf.md` | Gaussian gjf header (`%link0 → route → charge spin → coords`) |
| `charge-multiplicity.md` | Deciding `-q` and `-m` for an unfamiliar substrate (literature lookup workflow) |

## Decision tree: which format to feed `pdb2reaction`

```
Is the input a fresh extraction from the PDB Bank or a model from PyMOL/Maestro?
  └── PDB. pdb2reaction reads residue names directly; use -l 'RES:Q,...'
          to assign per-residue ligand charges.

Is the input a single-segment optimized structure (TS candidate, IRC endpoint)?
  └── XYZ is fine. Pass -q TOTAL_CHARGE and -m MULT explicitly.
      Or use --ref-pdb pointing back to the original PDB so -l still works.

Is the input a Gaussian gjf (with route line, charge, spin in header)?
  └── GJF. pdb2reaction parses the header automatically; -q and -m
      are inferred unless you override.
```

## Editing approach (agent-side)

When an agent must edit a structure file, the basic posture is:

1. **Read the file first** to understand current layout (residues, atom
   counts, charge / multiplicity if present).
2. **Identify the change** and confirm it does not violate format
   conventions (e.g. PDB column widths, XYZ first-line atom count).
3. For unknown charge / multiplicity values, **confirm with the user
   or do a literature lookup** before guessing — see
   `charge-multiplicity.md` for the workflow.
4. Make the smallest possible edit (single residue rename, single
   charge change). Avoid wholesale rewrites.

> Future expansion: this skill is intentionally a base layer.
> Subsequent rounds will add literature-database integration
> (PubChem / ChEBI / UniProt) and ML-based charge inference.

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
| `scan`, `scan2d`, `scan3d` | ✓ | — | — |
| `all` | ✓ | ✓ (single segment) | ✓ |
| `bond-summary` | ✓ | ✓ | — |

If you pass an XYZ to a subcommand that needs residue context (e.g.
`-l 'GLU:-1'`), supply `--ref-pdb <path>` so the residue mapping can be
recovered.

## Quick reference: which fields where

```
PDB  ATOM/HETATM record
     ┌───────────────────────────────────────────────────────┐
     │ name(13-16) altloc(17) resName(18-20) chainID(22)     │
     │ resSeq(23-26)  X(31-38)  Y(39-46)  Z(47-54)           │
     │ occupancy(55-60) bfactor(61-66) element(77-78)        │
     └───────────────────────────────────────────────────────┘

XYZ  line 1: <natoms>
     line 2: <comment, optional ASE Properties=…>
     line 3+: <element>  <x>  <y>  <z>

GJF  %nproc=...  %mem=...
     # <route line:  functional/basis  options>

     <title>

     <charge> <spin>
     <element>  <x>  <y>  <z>
     ...

     [optional connectivity, ECP blocks, ...]
```

Full byte-by-byte / per-keyword detail: see `pdb.md`, `xyz.md`, `gjf.md`.

## Charge / multiplicity defaults

- `-m 1` (singlet, closed shell) is the default for almost every
  organic / biological cluster.
- Use `-m 2` for radicals, `-m 3+` for unusual high-spin metal centers.
- `-q` (total charge) must be explicitly given for XYZ inputs (XYZ has
  no header). `-l 'RES:Q'` derives `-q` for PDB / GJF inputs from
  per-residue charges plus `pdb2reaction`'s amino-acid table.

If you're not sure about charge or spin, do **not** guess silently —
follow `charge-multiplicity.md`.

## See also

- `pdb2reaction-cli/extract.md` — residue selectors and link-H caps.
- `pdb2reaction-cli/SKILL.md` — common flag conventions across
  subcommands.
- `pdb2reaction-workflows-output/SKILL.md` — what comes out of the
  pipeline (also XYZ / PDB).