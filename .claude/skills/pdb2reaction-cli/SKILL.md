---
name: pdb2reaction-cli
description: Index and recipes for the 17 pdb2reaction subcommands. Each subcommand has its own md (extract.md / tsopt.md / …); this SKILL.md is the orientation, the cross-cutting flag conventions, and a small library of canonical recipes.
---

# pdb2reaction CLI

## Subcommand index

Each row points to the full per-subcommand md in this skill directory.

| md | subcommand | role (2 lines) |
|---|---|---|
| `all.md` | `all` | End-to-end pipeline: extract → MEP → TS → IRC → freq → (DFT) in one invocation.<br>Delegates to a base orientation; specific modes are in `all-{endpoint-mep,scan-list,ts-only}.md`. |
| `all-endpoint-mep.md` | `all` (mode 1) | Drives the pipeline from N reaction-ordered structures (R, optionally IM₁ … IMₙ, P).<br>Path search runs GSM/DMF between adjacent endpoints; recursion handles multi-step mechanisms. |
| `all-scan-list.md` | `all` (mode 2) | Drives the pipeline from a single reactant + a list of staged distance scans.<br>The scan list seeds the MEP; recursion handles intermediate states like in mode 1. |
| `all-ts-only.md` | `all` (mode 3) | Skips path search and starts from a TS candidate; runs `tsopt → irc → freq → dft`.<br>Use when you already have a transition-state guess (from a different code or a prior run). |
| `extract.md` | `extract` | Cuts an active-site cluster from a PDB around the substrate residues.<br>Handles residue selection, per-residue charge mapping (`-l`), and link-H placement in one pass. |
| `path-search.md` | `path-search` | Recursive MEP search (GSM or DMF) across N endpoints with bond-change segmentation.<br>Splits multi-step paths into one-TS-per-segment automatically. |
| `path-opt.md` | `path-opt` | MEP optimization for a **single** segment between two endpoints.<br>Building block of `path-search`; also useful for refining one segment without re-running the whole search. |
| `opt.md` | `opt` | Single-structure geometry optimization with LBFGS or RFO.<br>`--opt-mode grad` (LBFGS, default) is fast; `--opt-mode hess` (RFO) is robust on tricky surfaces. |
| `tsopt.md` | `tsopt` | TS optimization: RS-I-RFO (default) or Hessian-Guided Dimer.<br>`--opt-mode hess/rsirfo` for full-Hessian RS-I-RFO (default); `--opt-mode grad/dimer` for cheaper Dimer. |
| `freq.md` | `freq` | Vibrational analysis: Hessian, frequencies, normal-mode visualization, QRRHO thermochemistry.<br>Default temperature/pressure 298.15 K / 1 atm; partial-Hessian variant when `freeze_atoms` is non-empty. |
| `irc.md` | `irc` | IRC integration with EulerPC in mass-weighted Cartesians.<br>Forward + backward from a TS, plus LBFGS optimization of each endpoint. |
| `dft.md` | `dft` | Single-point DFT through PySCF (CPU) or GPU4PySCF (CUDA, x86_64).<br>`--engine gpu` is default when available; falls back to CPU on aarch64. |
| `scan.md` | `scan` | 1D distance scan with harmonic restraints to seed a path search.<br>Useful when neither endpoint nor TS guess is available — drives the bond manually. |
| `scan2d.md` | `scan2d` | 2D analog of `scan` with two restrained distances.<br>Generates a grid; pdb2reaction interpolates the MEP through the grid minima. |
| `scan3d.md` | `scan3d` | 3D analog with three restrained distances.<br>Rare but supported; output volume grows quickly, plan resources. |
| `trj2fig.md` | `trj2fig` | Plot an energy profile from an XYZ trajectory.<br>Reads ASE-style energies in the comment line and writes a static PNG/HTML. |
| `energy-diagram.md` | `energy-diagram` | Build an ad-hoc energy diagram from a list of state names + energies.<br>For composing diagrams that combine multiple `pdb2reaction` runs. |
| `add-elem-info.md` | `add-elem-info` | Repair / add the element column (PDB cols 77-78).<br>Run before `extract` if your PDB came out of PyMOL or Maestro and elements are missing. |
| `fix-altloc.md` | `fix-altloc` | Resolve PDB alternate locations (`altloc` field).<br>Pick a single conformation per residue; needed before `extract` on raw RCSB downloads. |
| `bond-summary.md` | `bond-summary` | Detect bond changes between two structures (e.g. R vs P).<br>Uses the same bond-change algorithm `path-search` invokes for segmentation. |

## Pipeline at a glance

```
PDB(s) ──► extract ──► path-search ──► tsopt ──► irc ──► freq ──► (dft)
                          │   │            │       │       │
                          ▼   └─ path-opt  ▼       ▼       ▼
                       seg_NN/                  result_freq/
                                                  result_irc/
                                                  result_tsopt/
```

`pdb2reaction all` chains the whole pipeline; each box is also available
as its own subcommand.

## Common flag conventions

These flags appear on most subcommands (canonical list:
`pdb2reaction <subcommand> --help`):

| Flag | Meaning |
|---|---|
| `-i, --input` | Input file(s); accepts `.pdb`, `.xyz`, `.gjf` |
| `-q, --charge` | Total charge (integer) |
| `-l, --ligand-charge` | `'RES1:Q1,RES2:Q2'` per-residue mapping (PDB inputs) |
| `-m, --multiplicity` | Spin multiplicity (2S+1), default 1 |
| `-b, --backend` | MLIP backend: `uma` / `orb` / `mace` / `aimnet2` |
| `-o, --out-dir` | Output directory, subcommand-specific default |
| `--config` | YAML configuration file applied before CLI flags |
| `--show-config` / `--dry-run` | Print resolved config without running |
| `--help-advanced` | Reveal hidden / advanced flags |
| `--ref-pdb` | Reference PDB used to derive residue context for XYZ inputs |
| `--solvent` | xTB-ALPB solvent name (e.g. `water`); `none` to disable |

Charge precedence: explicit `-q` > `-l 'RES:Q'` derivation > `--config` YAML > `defaults.py`.

## Canonical recipes

### Multi-input MEP for a 1-step reaction

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    --out-dir result_mep
```

### Single-input scan-list (when only the reactant is available)

```bash
pdb2reaction all -i 1.R.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --scan-lists '[("CS1 SAM 320","GPP 321 C7",1.60)]' \
                 '[("GPP`321/H11","GLU`186/OE2",0.90)]' \
    --tsopt --thermo \
    --out-dir result_scan
```

### Validate a TS candidate without re-running path search

```bash
pdb2reaction tsopt -i ts_guess.xyz -q -1 -m 1 -b uma -o result_tsopt
pdb2reaction freq  -i result_tsopt/final_geometry.xyz -q -1 -m 1 -b uma -o result_freq
pdb2reaction irc   -i result_tsopt/final_geometry.xyz -q -1 -m 1 -b uma -o result_irc
```

### DFT//MLIP single point on the rate-limiting TS

```bash
pdb2reaction dft -i seg_01/ts.pdb \
    -l 'SAM:1,GPP:-3' \
    --func-basis 'wb97m-v/def2-tzvpd' \
    --engine gpu
```

### Bond-change report between R and P

```bash
pdb2reaction bond-summary -i reactant.pdb product.pdb
```

## Cross-cutting caveats

| Pitfall | Fix |
|---|---|
| `--scan-lists` syntax error | The list is a Python literal-eval expression. Quote with single-quotes outside, double-quotes inside, and watch ` ` vs ``\``. |
| Wrong charge silently | Always run `--show-config` once before a long job; it prints the resolved charge. |
| Forgetting to pin `-b` for production | The default is `-b uma`; specify `-b uma` / `-b orb` / `-b mace` / `-b aimnet2` explicitly so a future default change cannot silently re-route the run. |
| `--config` YAML ignored | YAML is read **after** built-in defaults but **before** explicit CLI flags. Anything also given on CLI overrides YAML. |
| `--help-advanced` flags differ between versions | They are subject to change; if a flag isn't in `--help`, check `--help-advanced` and version-pin if the workflow is shared. |
| OOM on the Hessian step | Reduce `hessian_calc_mode` to `'FiniteDifference'`, set `return_partial_hessian=True`, or downgrade backend (UMA-m → UMA-s). |

## Defaults

Every default value is exported from `pdb2reaction.defaults` (read with
`import` — the skill does not transcribe values that change between
releases):

```bash
python -c "import pdb2reaction.defaults as d; print([n for n in dir(d) if n.endswith('_KW') or n.startswith('OUT_DIR')])"

# Examples:
python -c "import pdb2reaction.defaults as d; print(d.LBFGS_KW)"
python -c "import pdb2reaction.defaults as d; print(d.RSIRFO_KW)"
python -c "import pdb2reaction.defaults as d; print(d.IRC_KW)"
python -c "import pdb2reaction.defaults as d; print(d.UMA_CALC_KW)"
```

Each per-subcommand md points at the relevant `_KW` dict in the
"See also" section.

## See also

- `pdb2reaction-overview/SKILL.md` — what `pdb2reaction` is and when to
  use it.
- `pdb2reaction-structure-io/` — input file formats and charge / spin.
- `pdb2reaction-install-backends/` — `<tool>` / backend installation.
- `pdb2reaction-workflows-output/SKILL.md` — what comes out of each
  invocation, summary.json schema, R/TS/P canonical paths.
- `pdb2reaction-hpc/SKILL.md` — running these recipes on PBS / SLURM.
