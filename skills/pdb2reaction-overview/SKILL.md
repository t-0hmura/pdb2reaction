---
name: pdb2reaction-overview
description: Orientation for pdb2reaction — what it is, when to use it, and how it differs from generic QM/MLIP path-search tools (PDB-native input, GPU-accelerated pysisyphus fork, recursive bond-change-driven path search). TRIGGER on first-touch / "what is pdb2reaction" / "should I use it" questions. SKIP when the user already named a subcommand, an install issue, an output file, a structure format, or a cluster — sibling skills cover those.
---

# pdb2reaction Overview

## Purpose

`pdb2reaction` is a command-line toolkit that drives an MLIP reaction-path
workflow from a single `pdb2reaction all` invocation. With `-c/--center`, `all`
first extracts an active-site cluster from a protein–ligand PDB; without it,
the supplied PDB/mmCIF/XYZ/GJF model is used as-is. It then pre-optimizes and runs a
minimum-energy-path (MEP) search, stopping at the MEP's highest-energy image
(a TS *candidate*). The post-processing stages are opt-in flags: `--tsopt`
adds transition-state (TS) optimization + IRC validation, `--thermo` adds
vibrational analysis + QRRHO thermochemistry, `--dft` adds the DFT
single-point — e.g. `pdb2reaction all -i R.pdb P.pdb -q -1 --tsopt --thermo`.

Three things make it different from gluing together generic tools:

1. **PDB-native automation.** A residue-aware extractor cuts an active-site
   cluster, sums residue/ligand formal charges, and places cap hydrogens at
   supported carbon truncation boundaries without manual atom mapping.
2. **GPU-aware pysisyphus fork (bundled).** Geometry optimizers, TS searches
   (RS-P-RFO default, Dimer alternative), and IRC integrators include locally
   maintained tensor/GPU paths designed for the MLIP backends. Exact device
   placement is backend and operation dependent; verify it by profiling rather
   than assuming every operation stays on the GPU.
3. **Recursive bond-change-driven path search.** When the reactant and
   product differ by more than one elementary step, the path search
   detects bond changes along the MEP and recursively proposes narrower
   reactive segments. This is a geometry-based segmentation heuristic, not
   proof that a segment contains exactly one TS; validate every exported HEI
   with TS optimization, one-imaginary-mode analysis, and IRC connectivity.

## When to use it

| Goal | Fit |
|---|---|
| Cluster-model enzyme reaction mechanism (single or multi-step) | Primary use case |
| Validate a TS candidate with IRC + thermochemistry on MLIP | `pdb2reaction tsopt → irc → freq` |
| DFT//MLIP barrier evaluation | Run `pdb2reaction dft` consistently on the IRC-refined R, TS, and P geometries; one TS single point alone is not a barrier |
| Single-point energies on an arbitrary geometry (MLIP or DFT) | `pdb2reaction sp` (MLIP energy + forces, optional Hessian via `--hess`) / `pdb2reaction dft` |

## When *not* to use it

- Pure QM (DFT-only) without MLIP: stick with ORCA/Gaussian/Q-Chem directly.
- Explicit-solvent QM/MM with full force-field embedding: out of scope
  (`pdb2reaction` is cluster-model only).
- Free-energy simulations (umbrella sampling, metadynamics): out of scope.

## Quick check

```bash
pdb2reaction --version           # confirm install
pdb2reaction --help              # list subcommands
pdb2reaction all --help          # end-to-end primary flags
pdb2reaction all --help-advanced # every flag (--mep-mode, --opt-mode, --precision, ...)
```

If `pdb2reaction` is not on PATH, see the `pdb2reaction-install-backends`
skill (`SKILL.md` plus `core.md`) before doing anything else.

## Pipeline at a glance

| Stage | Role |
|---|---|
| `extract` | active-site cluster + cap-H atoms + total charge |
| `path-opt` / `path-search` | MEP (GSM or DMF): single-pass `path-opt` by default; `--refine-path` runs recursive `path-search` with bond-change-based candidate segmentation → `seg_01`, `seg_02`, … |
| `tsopt` | TS refinement per segment (RS-P-RFO default; Dimer alternative) |
| `irc` | forward / backward EulerPC IRC (caches endpoint Hessians) |
| `freq` | Hessian, vibrational frequencies, QRRHO thermochemistry |
| `dft` | (optional) ωB97M-V/def2-TZVPD single point on R, TS, P |

`pdb2reaction all` orchestrates the stages selected by its input mode and
flags; extraction (`-c`), TS/IRC (`--tsopt`), thermo (`--thermo`), and DFT
(`--dft`) are conditional. Each stage is also available as its own subcommand.

## Backend choices (MLIP)

`pdb2reaction` ships with four MLIP backends; pick with `-b <name>`:

| `-b` | Model family | Strength |
|---|---|---|
| `uma` (default) | UMA-s-1.1, UMA-s-1.2, UMA-m-1.1 (Meta FAIR) | Default integration; gated checkpoint access is required |
| `mace` | MACE-OMOL-0 | OMol model integration; currently needs a separate environment because its supported package stack conflicts with `fairchem-core` |
| `orb` | `orb_v3_conservative_omol` | Conservative model integration; pdb2reaction defaults it to fp64 because explicit fp32/TF32 can make finite-difference Hessians noisy. Validate `n_imag=1` and IRC like every backend. |
| `aimnet2` | AIMNet2 family | Model-specific element/state domains: default 14-element organic model, with separately selected radical, Pd, and reactive models |

Backend-specific install notes live in
`pdb2reaction-install-backends/{uma,mace,orb,aimnet2}.md`. The default-value
dictionaries are in `pdb2reaction.core.defaults` (read live, not transcribed):

```bash
python -c "import pdb2reaction.core.defaults as d; print(sorted(n for n in dir(d) if not n.startswith('_')))"
```

## Where the code lives

| File | What's there |
|---|---|
| `pdb2reaction/cli/app.py` | Click entry point, subcommand registry |
| `pdb2reaction/core/defaults.py` | Primary source for shared calculation-default dictionaries (UMA_CALC_KW, RSIRFO_KW, IRC_KW, …); some Click presentation/workflow-local defaults remain inline, so verify live help too |
| `pdb2reaction/backends/__init__.py` | `BACKEND_REGISTRY`, `create_calculator(...)` factory |
| `pdb2reaction/workflows/all.py` | End-to-end orchestration for `pdb2reaction all` |
| `pdb2reaction/workflows/extract.py` | PDB → cluster, residue table, cap-H placement |
| `pdb2reaction/workflows/path_opt.py` | Single-pass pairwise MEP (GSM or DMF) between adjacent endpoints — the default MEP route |
| `pdb2reaction/workflows/path_search.py` | Recursive MEP search with bond-change segmentation — the `--refine-path` route |
| `pdb2reaction/workflows/tsopt.py` | RS-P-RFO (default) / Dimer (alternative) transition-state search |
| `pdb2reaction/workflows/irc.py` | EulerPC IRC (caches endpoint Hessians) |
| `pdb2reaction/workflows/freq.py` | Hessian, frequencies, QRRHO thermochemistry |
| `pdb2reaction/workflows/dft.py` | PySCF / GPU4PySCF single-point driver |
| bundled `pysisyphus/` | GPU-tensor pysisyphus fork |
| bundled `thermoanalysis/` | QRRHO thermochemistry |

## Navigation map of the skill set

| You want to … | Read |
|---|---|
| Pick a subcommand and run it | `pdb2reaction-cli/SKILL.md` then the per-subcommand md |
| Read or edit a `.pdb` / `.cif` / `.mmcif` / `.xyz` / `.gjf` input | `pdb2reaction-structure-io/{SKILL,pdb,cif,xyz,gjf}.md` |
| Decide charge / multiplicity for a substrate | [`pdb2reaction-structure-io/charge-multiplicity.md`](../pdb2reaction-structure-io/charge-multiplicity.md) — for PDB/mmCIF, name unknown ligand charges with `-l 'RES:Q'`; standard amino acids and recognized ions use internal tables. Explicit `-q` sets the total, including in `all -c`; disagreement with the extraction-derived value produces a warning. |
| Install the toolkit or a specific backend | `pdb2reaction-install-backends/SKILL.md` + the relevant backend md |
| Build a recipe (multi-step / scan-list / endpoint MEP) | `pdb2reaction-workflows-output/SKILL.md` |
| Choose a TS-search strategy, or fix a bad imaginary-mode count | `pdb2reaction-ts-strategy/SKILL.md` |
| Submit on a PBS or SLURM cluster | `pdb2reaction-hpc/SKILL.md` |
| Detect what cluster / GPU / scheduler you are on | `pdb2reaction-env-detect/SKILL.md` |
| Drive the tool from an MCP client (18 MCP tools, `SubcmdResult` schema) | `pdb2reaction-mcp/SKILL.md` |
| Find which package layer to grep before touching code | `pdb2reaction-architecture/SKILL.md` |
