# pdb2reaction Documentation Refactoring Plan

This document outlines a comprehensive plan to improve the readability, clarity, and consistency of the pdb2reaction documentation.

---

## Current State Assessment

### Scoring (out of 10)

| Metric | Score | Notes |
|--------|-------|-------|
| Accuracy | 8.9 | Documentation matches code behavior well |
| Clarity | 6.9 | Entry points are scattered; role overlap causes confusion |
| Readability | 6.5 | Inconsistent structure; unnatural English in places |

### Key Problems Identified

1. **TL;DR + Overview duplication**: 9 of 21 pages have TL;DR blocks that largely repeat the Overview paragraph
2. **Landing page overlap**: README.md, docs/index.md, and getting-started.md contain duplicate content (examples, subcommand tables, requirements)
3. **Unnatural English**: Terms like "mental model" and "orchestrator" feel overly formal for user documentation
4. **Inconsistent templates**: Subcommand pages vary in structure (some have TL;DR, some don't; section order varies)
5. **Heavy content in wrong places**: HPC scripts in uma_pysis.md make the page harder to navigate

---

## Phase 1: TL;DR Removal and Opening Paragraph Integration

**Goal**: Replace TL;DR blocks with integrated opening paragraphs
**Impact**: High
**Effort**: Low

### Target Files (9 files)

| File | Current TL;DR | Integrated Opening (proposed) |
|------|---------------|-------------------------------|
| `extract.md` | "Extract a cluster model (active-site pocket)..." | `pdb2reaction extract` creates an active-site pocket from a protein–ligand PDB. Specify substrates with `-c` (residue names, IDs, or a PDB path), and use `--ligand-charge` for non-standard residue charges. Link hydrogens cap severed bonds automatically. |
| `path_search.md` | "Build a continuous MEP across 2+ structures..." | `pdb2reaction path-search` builds a continuous minimum-energy path (MEP) across two or more structures using GSM (default) or DMF. It automatically refines regions with bond changes and identifies the highest-energy image as a TS candidate. |
| `path_opt.md` | "Find the MEP between exactly two structures..." | `pdb2reaction path-opt` finds the minimum-energy path between exactly two structures using GSM (default) or DMF. It outputs the path trajectory and identifies the highest-energy image (HEI). For multi-structure workflows with automatic refinement, use `path-search` instead. |
| `tsopt.md` | "Optimize transition states using Dimer..." | `pdb2reaction tsopt` optimizes a transition state using Dimer (`--opt-mode light`) or RS-I-RFO (`--opt-mode heavy`). When VRAM permits, use `--hessian-calc-mode Analytical` for faster convergence. The optimized TS has exactly one imaginary frequency. |
| `opt.md` | "Optimize a single structure to a local minimum..." | `pdb2reaction opt` optimizes a single structure to a local minimum using L-BFGS (`--opt-mode light`, default) or RFO (`--opt-mode heavy`). For PDB inputs, link-hydrogen parents are automatically frozen. |
| `freq.md` | "Compute vibrational frequencies and thermochemistry..." | `pdb2reaction freq` computes vibrational frequencies and thermochemistry (ZPE, Gibbs energy, etc.) using UMA. Use `--hessian-calc-mode Analytical` for faster Hessian evaluation when VRAM permits. Imaginary frequencies appear as negative values. |
| `irc.md` | "Trace the reaction path from a TS..." | `pdb2reaction irc` traces the intrinsic reaction coordinate from a transition state toward reactant and product. Both forward and backward branches run by default. Use `--hessian-calc-mode Analytical` for faster Hessian evaluation when VRAM permits. |
| `scan.md` | "Drive a reaction coordinate by scanning bond distances..." | `pdb2reaction scan` drives a reaction coordinate by scanning bond distances with harmonic restraints. Use `--scan-lists` to specify target distances. Multiple stages run sequentially, each starting from the previous result. |
| `dft.md` | "Run single-point DFT calculations..." | `pdb2reaction dft` runs single-point DFT calculations using GPU4PySCF (or CPU PySCF as fallback). The default functional/basis is wB97M-V/def2-TZVPD. Results include energy and population analysis (Mulliken, meta-Lowdin, IAO charges). |

### Edit Pattern

```diff
  ## Overview

- > **TL;DR:** <summary content>
-
- <longer explanation that often repeats TL;DR content>
+ `pdb2reaction <cmd>` <does X for Y>. <key inputs/outputs>. <when to use key option>.
+
+ <remaining unique details from the original longer explanation>
```

---

## Phase 2: Unnatural English Replacement

**Goal**: Replace formal/technical jargon with natural, readable English
**Impact**: Medium
**Effort**: Low

### Replacement Table

| Current Expression | Replacement | Files Affected |
|--------------------|-------------|----------------|
| "mental model" | "overview" or "how the pieces fit together" | concepts.md, docs/index.md |
| "orchestrator" / "orchestrates" | "runs" / "ties together" / "coordinates" | getting-started.md, concepts.md, yaml-reference.md |
| "format-aware conversions mirror" | "also writes" / "generates companion" | path_search.md, all.md |
| "End-to-end ensemble" | "Multi-structure workflow" | all.md |
| "acts as an orchestrator" | "runs the full pipeline" | getting-started.md |

### Specific Edits

#### concepts.md:3
```diff
- This page gives a **mental model** of how `pdb2reaction` is structured: what "pockets", "templates", "segments", and "images" mean, and how the top-level `all` workflow orchestrates subcommands.
+ This page explains the key terms in pdb2reaction (pockets, templates, segments, images) and how the `all` command ties together the subcommands.
```

#### getting-started.md:191
```diff
- The `all` workflow acts as an **orchestrator**: it chains cluster extraction, MEP search, TS optimization, vibrational analysis, and optional single-point DFT calculations into a single command.
+ The `all` command runs the full pipeline—cluster extraction, MEP search, TS optimization, vibrational analysis, and optional DFT—in a single invocation.
```

#### docs/index.md:57
```diff
- [**Concepts & Workflow**](concepts.md) - Mental model of pockets, templates, segments, and stages
+ [**Concepts & Workflow**](concepts.md) - Key terms: pockets, templates, segments, and stages
```

---

## Phase 3: Landing Page Deduplication

**Goal**: Establish clear roles for README, index, and getting-started
**Impact**: High
**Effort**: Medium

### Role Definitions

| Page | Purpose | Content Scope |
|------|---------|---------------|
| `README.md` | GitHub landing page | 1-paragraph overview, 1 minimal example, link to docs |
| `docs/index.md` | Documentation home | Quick Start by Goal table, navigation to sections |
| `getting-started.md` | Tutorial | Installation, first run, workflow explanation |

### Content Migration

| Content | Current Location(s) | Target Location |
|---------|---------------------|-----------------|
| System Requirements | README, index, getting-started | getting-started.md only |
| CLI Subcommands table | README, index, getting-started | index.md only (abbreviated in README) |
| Quick Examples (4) | README, index | README (2 examples), remove from index |
| Output Structure tree | README, index, all.md | all.md only |
| Installation steps | README, getting-started | getting-started.md only (brief in README) |

### README.md Target Structure

```markdown
# pdb2reaction

One-paragraph overview.

## Quick Start

```bash
pip install ...
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

## Key Commands

| Command | Purpose |
|---------|---------|
| `all` | End-to-end workflow |
| `extract` | Pocket extraction |
| `path-search` | MEP search |
| `tsopt` | TS optimization |
| `freq` | Vibrational analysis |

[Full documentation](docs/index.md) | [Getting Started](docs/getting-started.md)

## License / Citation
```

---

## Phase 4: Subcommand Template Unification

**Goal**: All subcommand pages follow the same structure
**Impact**: Medium
**Effort**: Medium

### Standard Template

```markdown
# `<command>`

## Overview
<1-3 sentences: what it does, key inputs/outputs, when to use>

## When to use
- Use this when...
- Prefer `<alternative>` when...

## Usage
```bash
pdb2reaction <command> -i INPUT [options]
```

## Examples
### Minimal
```bash
<shortest working example>
```

### With common options
```bash
<practical example with frequently-used options>
```

## CLI options
<Table of most-used options only; full table in collapsed section or at end>

## Outputs
```text
<output directory tree>
```

## Notes
<Important caveats, tips>

## YAML configuration
<If applicable: key YAML sections, minimal example>

## See Also
- [Related command](related.md)
```

### Files to Update (13 subcommand pages)

- all.md
- extract.md
- opt.md
- tsopt.md
- path_opt.md
- path_search.md
- scan.md
- scan2d.md
- scan3d.md
- irc.md
- freq.md
- dft.md
- trj2fig.md

---

## Phase 5: CLI Conventions Consolidation

**Goal**: Single source of truth for CLI conventions
**Impact**: Medium
**Effort**: Low

### New File: `docs/cli-conventions.md`

Content to consolidate from:
- getting-started.md "CLI conventions" section
- concepts.md "A few CLI conventions worth knowing"
- Scattered explanations in subcommand pages

### Proposed Content

```markdown
# CLI Conventions

## Boolean Options
All boolean options require explicit `True` or `False`:
```bash
--tsopt True --thermo True --dft False
```

## Residue Selectors
...

## Charge Mapping
...

## Atom Selectors
...

## Input Requirements
...
```

---

## Phase 6: Heavy Content Separation (Future)

**Goal**: Move HPC/advanced content to dedicated pages
**Impact**: Low
**Effort**: High

### Proposed Changes

| Current | Proposed |
|---------|----------|
| uma_pysis.md with HPC scripts | uma_pysis.md (basics) + advanced/hpc-scaling.md |
| yaml-reference.md (all in one) | Keep as-is (reference pages can be long) |

---

## Implementation Schedule

| Week | Phase | Files | Estimated Changes |
|------|-------|-------|-------------------|
| 1 | Phase 1 | 9 files | ~90 line edits |
| 1 | Phase 2 | 5 files | ~10 line edits |
| 2 | Phase 3 | 3 files | ~200 line edits |
| 2 | Phase 5 | 1 new file | ~50 lines new |
| 3 | Phase 4 | 13 files | ~300 line edits |
| 4 | Japanese sync | 13 files | Mirror English changes |

---

## Success Metrics

After implementation:

| Metric | Current | Target |
|--------|---------|--------|
| Clarity Score | 6.9 | 8.0+ |
| Readability Score | 6.5 | 7.5+ |
| TL;DR blocks | 9 | 0 |
| "mental model" occurrences | 2 | 0 |
| "orchestrator" occurrences | 3 | 0 |
| Duplicate content sections | ~8 | ~2 |

---

## Changelog

| Date | Phase | Status |
|------|-------|--------|
| 2026-01-31 | Plan created | Complete |
| 2026-01-31 | Phase 1: TL;DR removal | Complete |
| 2026-01-31 | Phase 2: Unnatural English replacement | Complete |
| 2026-01-31 | Phase 3: Landing page deduplication | Complete |
| 2026-01-31 | Phase 4: Subcommand template review | Complete |
| 2026-01-31 | Phase 5: cli-conventions.md created | Complete |

