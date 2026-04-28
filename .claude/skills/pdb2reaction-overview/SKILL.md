---
name: pdb2reaction-overview
description: What pdb2reaction is, when to use it, and the design choices that distinguish it from generic QM/MLIP path-search tools (PDB-native input, GPU-accelerated pysisyphus fork, recursive bond-change-driven path search).
---

# pdb2reaction Overview

## What it is

`pdb2reaction` is a command-line toolkit that takes a protein–ligand PDB and
runs the entire MLIP-driven reaction-path workflow — active-site extraction,
minimum-energy-path (MEP) search, transition-state (TS) optimization, IRC
validation, vibrational analysis, and an optional DFT single-point — through
a single `pdb2reaction all` invocation.

Three things make it different from gluing together generic tools:

1. **PDB-native automation.** A residue-aware extractor cuts an active-site
   cluster, sums residue/ligand formal charges, and places link hydrogens
   along severed covalent bonds without manual atom mapping.
2. **GPU-accelerated pysisyphus fork (bundled).** Geometry optimizers, TS
   searches (RS-I-RFO default, Dimer alternative), and IRC integrators keep the heavy tensor
   work on the same device as the MLIP — no CPU round-trip per step.
3. **Recursive bond-change-driven path search.** When the reactant and
   product differ by more than one elementary step, the path search
   detects bond changes along the MEP and recursively re-segments until
   each segment crosses exactly one TS.

## When to use it

| Goal | Fit |
|---|---|
| Cluster-model enzyme reaction mechanism (single or multi-step) | Primary use case |
| Validate a TS candidate with IRC + thermochemistry on MLIP | `pdb2reaction tsopt → irc → freq` |
| DFT//MLIP barrier refinement | `pdb2reaction dft -i <ts.pdb>` after IRC |
| Single-point energies on an arbitrary geometry (MLIP or DFT) | `pdb2reaction opt` / `pdb2reaction dft` |

## When *not* to use it

- Pure QM (DFT-only) without MLIP: stick with ORCA/Gaussian/Q-Chem directly.
- Explicit-solvent QM/MM: see ML/MM-style frameworks (e.g. `mlmm_toolkit`),
  not `pdb2reaction` (which is cluster-model only).
- Free-energy simulations (umbrella sampling, metadynamics): out of scope.

## Quick check

```bash
pdb2reaction --version           # confirm install
pdb2reaction --help              # list subcommands
pdb2reaction all --help          # end-to-end flag list
```

If `pdb2reaction` is not on PATH, see the `pdb2reaction-install-backends`
skill (`SKILL.md` plus `core.md`) before doing anything else.

## Pipeline at a glance

```
PDB(s)
  │
  ▼
[extract]      active-site cluster + link-H caps + total charge
  │
  ▼
[path-search]  MEP (GSM or DMF), recursive bond-change segmentation
  │            → seg_01, seg_02, ... (one per elementary step)
  ▼
[tsopt]        TS refinement per segment (RS-I-RFO default; Dimer alternative)
  │
  ▼
[irc]          forward/backward IRC, endpoint LBFGS optimization
  │
  ▼
[freq]         Hessian, vibrational frequencies, QRRHO thermochemistry
  │
  ▼
[dft]          (optional) ωB97M-V/def2-TZVPD single point on R, TS, P
```

`pdb2reaction all` chains all of these. Each stage is also available as
its own subcommand (`pdb2reaction tsopt -i ts.xyz`) if you only want one
piece.

## Backend choices (MLIP)

`pdb2reaction` ships with four MLIP backends; pick with `-b <name>`:

| `-b` | Model family | Strength |
|---|---|---|
| `uma` (default) | UMA-s-1.1, UMA-s-1.2, UMA-m-1.1 (Meta FAIR) | Broadest coverage, default for most workflows |
| `mace` | MACE-OMOL-0 | Strong on organic + 1st-row metals; needs separate env (e3nn conflict) |
| `orb` | Orb-v3-omol | Fast screening; lower TS accuracy than UMA/MACE |
| `aimnet2` | AIMNet2 | Element coverage limited; small organics |

Backend-specific install notes live in
`pdb2reaction-install-backends/{uma,mace,orb,aimnet2}.md`. The default-value
dictionaries are in `pdb2reaction.defaults` (read live, not transcribed):

```bash
python -c "import pdb2reaction.defaults as d; print(sorted(n for n in dir(d) if not n.startswith('_')))"
```

## Sibling project: mlmm-toolkit

`mlmm_toolkit` shares the same `pysisyphus` and `thermoanalysis` core but
targets ML/MM ONIOM workflows (3-layer ML/movable-MM/frozen, AmberTools
parm7 topology, microiteration). Choose:

| Use case | Toolkit |
|---|---|
| Cluster model, no MM environment | `pdb2reaction` |
| Solvated enzyme with MM force field around ML region | `mlmm_toolkit` |
| Need automatic Amber parameterization | `mlmm_toolkit` (`mm-parm` subcommand) |
| Need recursive multi-step path-search | `pdb2reaction` |

The two toolkits **cannot share a Python environment** (incompatible
pysisyphus forks and `e3nn` versions). Keep separate `conda env`s.

## Where the code lives

| File | What's there |
|---|---|
| `pdb2reaction/cli.py` | Click entry point, subcommand registry |
| `pdb2reaction/defaults.py` | All default kwarg dicts (UMA_CALC_KW, RSIRFO_KW, IRC_KW, …) — single source of truth |
| `pdb2reaction/backends/__init__.py` | `BACKEND_REGISTRY`, `create_calculator(...)` factory |
| `pdb2reaction/all.py` | End-to-end orchestration for `pdb2reaction all` |
| `pdb2reaction/extract.py` | PDB → cluster, residue table, link-H placement |
| `pdb2reaction/path_search.py` | Recursive MEP search, bond-change segmentation |
| `pdb2reaction/tsopt.py` | RS-I-RFO (default) / Dimer (alternative) transition-state search |
| `pdb2reaction/irc.py` | EulerPC IRC + endpoint optimization |
| `pdb2reaction/freq.py` | Hessian, frequencies, QRRHO thermochemistry |
| `pdb2reaction/dft.py` | PySCF / GPU4PySCF single-point driver |
| bundled `pysisyphus/` | GPU-tensor pysisyphus fork |
| bundled `thermoanalysis/` | QRRHO thermochemistry |

## Navigation map of the skill set

| You want to … | Read |
|---|---|
| Pick a subcommand and run it | `pdb2reaction-cli/SKILL.md` then the per-subcommand md |
| Read or edit a `.pdb` / `.xyz` / `.gjf` input | `pdb2reaction-structure-io/{SKILL,pdb,xyz,gjf}.md` |
| Decide charge / multiplicity for a substrate | [`pdb2reaction-structure-io/charge-multiplicity.md`](../pdb2reaction-structure-io/charge-multiplicity.md) |
| Install the toolkit or a specific backend | `pdb2reaction-install-backends/SKILL.md` + the relevant backend md |
| Build a recipe (multi-step / scan-list / endpoint MEP) | `pdb2reaction-workflows-output/SKILL.md` |
| Submit on a PBS or SLURM cluster | `pdb2reaction-hpc/SKILL.md` |
| Detect what cluster / GPU / scheduler you are on | `pdb2reaction-env-detect/SKILL.md` |
