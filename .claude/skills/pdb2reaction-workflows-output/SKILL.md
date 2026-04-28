---
name: pdb2reaction-workflows-output
description: Canonical pdb2reaction workflows (cluster + 1-step / multistep / scan-list / endpoint-MEP / TS-only / DFT//MLIP), plus the summary.json schema, R/TS/P/IM canonical paths, bond-change interpretation, and energy-diagram conventions used to extract numerical results.
---

# pdb2reaction Workflows and Output Parsing

## Six canonical workflows

### 1. Cluster + 1-step reaction (multi-input MEP)

You have R and P PDBs (from a published QM study). One step.

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -o result_mep
```

Result: `result_mep/seg_NN/{reactant,ts,product}.pdb`,
`summary.json["segments"][0]["barrier_kcal"]`.

### 2. Multi-step recursive (multi-input MEP, recursive segmentation)

You have R and P. The mechanism is multi-step, but you don't have
intermediates handy. Path-search recursively segments by detecting
bond changes.

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c '...' -l '...' \
    --tsopt --thermo \
    -o result_mep
```

The output `summary.json["n_segments"]` may be > 1 — that's the
recursion finding intermediates the inputs didn't contain.

### 3. Single-input scan-list (when only R is available)

You have just the reactant. Articulate the reaction as a sequence of
distance-restraint scans.

```bash
pdb2reaction all -i 1.R.pdb \
    -c '...' -l '...' \
    --scan-lists '[("CS1 SAM 320","C7 GPP 321",1.60)]' \
                 '[("H11 GPP 321","OE2 GLU 186",0.90)]' \
    --tsopt --thermo \
    -o result_scan
```

Each `--scan-lists` argument is one stage. See
[`pdb2reaction-cli/all-scan-list.md`](../pdb2reaction-cli/all-scan-list.md) for syntax details.

### 4. Endpoint-MEP with explicit intermediates

You have R, IM₁, IM₂, P from the literature. Pass them all in order:

```bash
pdb2reaction all -i 1.R.pdb 2.IM1.pdb 3.IM2.pdb 4.P.pdb \
    -c '...' -l '...' \
    --tsopt --thermo \
    -o result_mep_4pt
```

Recursive segmentation still runs *between* adjacent endpoints, so
you don't have to provide every elementary step.

### 5. TS-only validation (existing TS candidate)

You have a TS guess from another code or a prior run. Skip
extract / path-search:

```bash
pdb2reaction tsopt -i ts.xyz -q -1 -m 1 -b uma -o result_tsopt
pdb2reaction freq  -i result_tsopt/final_geometry.xyz -q -1 -m 1 -b uma -o result_freq
pdb2reaction irc   -i result_tsopt/final_geometry.xyz -q -1 -m 1 -b uma -o result_irc
```

Or use `pdb2reaction all` with a single `-i` (collapses to TS-only
mode automatically; see [`pdb2reaction-cli/all-ts-only.md`](../pdb2reaction-cli/all-ts-only.md)).

### 6. DFT//MLIP refinement

After any of the above, refine R / TS / P energies at DFT level:

```bash
pdb2reaction dft -i seg_01/reactant.pdb \
    -l 'SAM:1,GPP:-3' \
    --func-basis 'wb97m-v/def2-tzvpd' \
    --engine gpu \
    -o dft_R
pdb2reaction dft -i seg_01/ts.pdb       -l '...' --func-basis '...' -o dft_TS
pdb2reaction dft -i seg_01/product.pdb  -l '...' --func-basis '...' -o dft_P
```

Composite the energies with `energy-diagram` (see below).

## summary.json schema (`pdb2reaction all` output)

Top-level keys:

| Key | Description |
|---|---|
| `command` | The subcommand (`"all"`, `"tsopt"`, …) |
| `pdb2reaction_version` | Toolkit version that produced this output |
| `status` | `"success"`, `"partial"`, or `"failed"` |
| `charge` / `spin` | Resolved cluster charge / multiplicity |
| `environment` | `{device, gpu_name, gpu_vram_gb, cuda_version, cpu, n_cpus, ram_gb}` |
| `config` | Full effective config after CLI + YAML + defaults merge |
| `freeze_atoms` | Indices held fixed during optimization (link-H parents). 0-based internal indices; `summary.log` echoes them as 1-based. |
| `n_images` | MEP image count in path-search / path-opt; IRC trajectory frame count in tsopt-only mode |
| `n_segments` | Total MEP segment count (reactive + bridge) |
| `n_segments_reactive` | Reactive segment count (kind != "bridge") |
| `rate_limiting_step` | Dict `{segment, barrier_kcal, method}` for the highest-barrier segment |
| `overall_reaction_energy_kcal` | R → P total energy difference |
| `segments` | Lightweight per-segment summary (see below) |
| `post_segments` | Per-segment post-processing details (tsopt / IRC / freq / DFT outputs) |
| `key_output_files` | Map of role → path (str) for top-level files; `seg_NN` entries are `{description, files}` dicts |
| `pipeline_mode` | One of `"path-search"`, `"path-opt"`, `"tsopt-only"` |
| `mlip_backend` | Which backend produced the energies |
| `energy_diagrams` | List of energy-diagram entries (PNG / HTML paths and metadata) |
| `out_dir` | Output directory absolute path |

Per-segment keys (`summary.json["segments"][i]`) — lightweight, MEP-level:

| Key | Description |
|---|---|
| `index` | Segment index (1-based int; written zero-padded as `seg_01/`, `seg_02/`, …) |
| `tag` | Segment tag (`"reactive"` / `"non-reactive"`) |
| `kind` | Segment kind (`"seg"`, `"bridge"`, or `"tsopt"`) |
| `barrier_kcal` | TS – R energy (kcal/mol) — the rate constant input |
| `delta_kcal` | P – R energy (kcal/mol) |
| `bond_changes` | List of single-key dicts (one per detected change): `[{"Bond formed (k)": ["A-B : 3.17 Å --> 1.68 Å", ...]}, {"Bond broken (k)": [...]}]`. Cutoff 1.20× covalent radii with margin 0.05 — see [`pdb2reaction-cli/bond-summary.md`](../pdb2reaction-cli/bond-summary.md). |

Per-segment post-processing keys (`summary.json["post_segments"][i]`) — when `--tsopt`, `--thermo`, or `--dft` was passed:

| Key | Description |
|---|---|
| `index` / `tag` / `kind` / `bond_changes` | Mirror of the corresponding `segments[i]` row |
| `mep_barrier_kcal` / `mep_delta_kcal` | Refined barriers from the post-IRC / tsopt re-evaluation |
| `post_dir` | Subdirectory holding tsopt / freq / IRC outputs for this segment |
| `irc_plot` / `irc_traj` | Paths to the IRC trace PNG and trajectory XYZ |
| `uma` | Per-stage MLIP energy block (or whichever backend was used) |
| `ts_imag` | Dict `{n_imag, nu_imag_max_cm, min_abs_imag_cm, min_freq_cm}` describing the TS spectrum |
| `ts_imag_freq_cm` | Peak imaginary frequency (cm⁻¹); same as `ts_imag.nu_imag_max_cm` |
| `gibbs_uma` | QRRHO Gibbs energies (when `--thermo`) |

## R/TS/P canonical paths

Two locations get written for each elementary step:

```
result_all/
├── seg_NN/                                 # CANONICAL — top-level, post-IRC re-optimized (RFO by default; LBFGS via --opt-mode-post grad)
│   ├── reactant.{xyz,pdb}                  # IRC backward endpoint, re-optimized
│   ├── ts.{xyz,pdb}                        # tsopt'd transition state
│   └── product.{xyz,pdb}                   # IRC forward endpoint, re-optimized
└── path_search/post_seg_NN/structures/
    ├── reactant.{xyz,pdb}                  # same as above (canonical) — nested copy
    ├── reactant_irc.{xyz,pdb}              # raw IRC backward end (pre-RFO/LBFGS re-optimization)
    ├── ts.{xyz,pdb}                        # same as above
    ├── product.{xyz,pdb}                   # same as above (canonical)
    └── product_irc.{xyz,pdb}               # raw IRC forward end (pre-RFO/LBFGS re-optimization)
```

**Rule of thumb**: read from `seg_NN/` for downstream stages. Use
`reactant_irc.xyz` / `product_irc.xyz` only when debugging
IRC vs. post-IRC re-optimization divergence.

`bond_changes` are computed from `reactant.xyz` / `product.xyz`
(post-IRC re-optimization), not from the raw IRC endpoints.

## Programmatic key extraction

```python
import json

# Per-segment barriers
d = json.load(open("result_all/summary.json"))
for seg in d["segments"]:
    print(f"seg_{seg['index']:02d}: ΔE‡ = {seg['barrier_kcal']:.1f} kcal/mol, "
          f"ΔE = {seg['delta_kcal']:.1f} kcal/mol")

# Rate-limiting barrier
rls = d["rate_limiting_step"]
print(f"rate-limiting: seg_{rls['segment']:02d}, barrier = "
      f"{d['segments'][rls-1]['barrier_kcal']:.1f} kcal/mol")

# Imaginary-mode check on every TS (post-processing data)
for ps in d.get("post_segments", []):
    ts = ps.get("ts_imag") or {}
    n = ts.get("n_imag")
    if n is not None and n != 1:
        print(f"WARNING: seg_{ps['index']:02d} has {n} imaginary modes "
              f"(freq {ps.get('ts_imag_freq_cm')!r})")

# Backend-level energies / Gibbs (when --tsopt and/or --thermo was passed)
for ps in d.get("post_segments", []):
    print(ps["index"], ps.get("uma"), ps.get("gibbs_uma"))
```

## Bond-change interpretation

`segments[i]["bond_changes"]` is a list of single-key dicts, one per
detected change, with the change kind encoded in the key name:

```json
[
  {"Bond formed (1)": ["C508-C567 : 3.166 Å --> 1.675 Å"]},
  {"Bond broken (1)": ["S507-C508 : 1.798 Å --> 3.459 Å"]}
]
```

The `(k)` integer is the consecutive-frame-pair index for multi-step
IRC traces. Each string carries `<atom>-<atom> : <d_R> Å --> <d_P> Å`.

Reading rules:

- "Bond formed (k)" entries list bonds that exist in P but not R
  (covalent-radius cutoff 1.20×).
- "Bond broken (k)" entries list bonds that exist in R but not P.
- For a single reactive segment you usually expect 1–4 entries combined.
- If a single segment shows > 8 bond changes, the recursive
  segmentation may have failed — inspect the geometries before trusting
  the barrier.

## Failed-run output

When `summary.json["status"] != "success"`, look at:

1. `summary.log` — human-readable, prints the failure point first.
2. `path_search/post_seg_NN/<stage>/result.json` (or `path_opt/...`
   when `--refine-path False`) — per-stage status (which step
   crashed). Note: top-level `post_seg_NN/` is a re-export; the
   primary copy lives under `path_search/` or `path_opt/`.
3. The per-stage log produced by the underlying tool (e.g. PySCF /
   pysisyphus) is silenced by default (see
   `default_group._silence_pysisyphus_loggers`); rerun with `--dump`
   to keep optimizer trajectories.

Even on failed runs, partial outputs are kept:

- `path_search/seg_NN/` exists for any segment that completed
  path-search (even if downstream stages failed).
- `seg_NN/` (top-level) is **only populated** for fully-successful
  segments.

## Energy diagrams

`pdb2reaction all` writes `path_search/energy_diagram_*.png` (or
`path_opt/energy_diagram_*.png` when `--refine-path False`):

- `energy_diagram_MEP.png` — bare MEP energies from the path-search
  string (MLIP, no thermochemistry).
- `energy_diagram_UMA_all.png` (etc.) — per-segment energies for
  whichever backend was used.
- `energy_diagram_G_UMA_all.png` — Gibbs free-energy diagram with
  QRRHO thermochemistry (when `--thermo`).

To compose a custom diagram from energies of multiple runs, use
[`pdb2reaction-cli/energy-diagram.md`](../pdb2reaction-cli/energy-diagram.md):

```bash
pdb2reaction energy-diagram \
    -i 0.0 21.5 -0.7 2.2 -18.2 \
    --label-x R TS1 IM TS2 P \
    -o my_diagram.png
```

## Cross-references

- [`pdb2reaction-cli/all.md`](../pdb2reaction-cli/all.md) and the three `all-*.md` mode files.
- `pdb2reaction-cli/{tsopt,freq,irc,dft}.md` — per-stage
  `result.json` schemas.
- [`pdb2reaction-cli/bond-summary.md`](../pdb2reaction-cli/bond-summary.md) — same bond-change algorithm,
  standalone.
- `pdb2reaction-structure-io/SKILL.md` — input file formats that feed
  these workflows.
