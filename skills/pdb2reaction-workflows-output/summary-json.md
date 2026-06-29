# `summary.json` — schema and parsing

`pdb2reaction all` writes `summary.json` at the top of `--out-dir`. It
is the canonical machine-readable artifact for downstream analysis
(barriers, ΔE, Gibbs, TS imaginary modes, file paths).

## Top-level keys

| Key | Description |
|---|---|
| `command` | The subcommand (`"all"`, `"tsopt"`, …) |
| `pdb2reaction_version` | Toolkit version that produced this output |
| `status` | `"success"`, `"partial"`, or `"failed"` |
| `charge` / `spin` | Resolved cluster charge / multiplicity |
| `environment` | `{device, gpu_name, gpu_vram_gb, cuda_version, cpu, n_cpus, ram_gb}` |
| `config` | Full effective config after CLI + YAML + defaults merge |
| `freeze_atoms` | Indices held fixed during optimization (cap-H parents). 0-based internal indices; `summary.log` echoes them as 1-based. |
| `n_images` | MEP image count in path-search / path-opt; IRC trajectory frame count in tsopt-only mode |
| `n_segments` | Total MEP segment count (reactive + bridge) |
| `n_segments_reactive` | Reactive segment count (kind != "bridge") |
| `rate_limiting_step` | Dict `{segment, barrier_kcal, method}` for the highest-barrier segment |
| `overall_reaction_energy_kcal` | R → P total energy difference |
| `segments` | Lightweight per-segment summary (see below) |
| `post_segments` | Per-segment post-processing details (tsopt / IRC / freq / DFT outputs) |
| `key_output_files` | Top-level entries map filename → description (str); `seg_NN` entries are `{description, files}` dicts |
| `pipeline_mode` | One of `"path-search"`, `"path-opt"`, `"tsopt-only"` |
| `mlip_backend` | Which backend produced the energies |
| `energy_diagrams` | List of energy-diagram entries (PNG / HTML paths and metadata) |
| `out_dir` | Output directory absolute path |

## Per-segment keys (`summary.json["segments"][i]`)

Lightweight, MEP-level:

| Key | Description |
|---|---|
| `index` | Segment index (1-based int; written zero-padded as `seg_01/`, `seg_02/`, …) |
| `tag` | Segment-name string (`seg_NN` / `*_bridge`). To select reactive segments, filter `kind != "bridge"` |
| `kind` | Segment kind (`"seg"`, `"bridge"`, or `"tsopt"`) |
| `barrier_kcal` | TS – R energy (kcal/mol) — the rate constant input |
| `delta_kcal` | P – R energy (kcal/mol) |
| `bond_changes` | List of single-key dicts: `[{"Bond formed (k)": ["A-B : 3.17 Å --> 1.68 Å", ...]}, {"Bond broken (k)": [...]}]`. **Standalone `irc result.json["bond_changes"]` uses a different shape**: `{"formed": [str], "broken": [str]}` (flat dict). The list-of-dicts form here is the `path_search` / `all` summary.json shape. See "Bond-change interpretation" below for cutoff / algorithm. |

## Per-segment post-processing keys (`summary.json["post_segments"][i]`)

Present when `--tsopt`, `--thermo`, or `--dft` was passed:

| Key | Description |
|---|---|
| `index` / `tag` / `kind` / `bond_changes` | Mirror of the corresponding `segments[i]` row |
| `mep_barrier_kcal` / `mep_delta_kcal` | MEP-level (un-refined) barrier / ΔE, mirroring `segments[i]` (refined post-IRC/tsopt energies live in the per-segment `uma` block) |
| `post_dir` | Subdirectory holding tsopt / freq / IRC outputs for this segment |
| `irc_plot` / `irc_traj` | Paths to the IRC trace PNG and trajectory XYZ |
| `uma` | Per-stage MLIP energy block (or whichever backend was used) |
| `ts_imag` | Dict `{n_imag, nu_imag_max_cm, min_abs_imag_cm, min_freq_cm}` describing the TS spectrum |
| `ts_imag_freq_cm` | Peak imaginary frequency (cm⁻¹); same as `ts_imag.nu_imag_max_cm` |
| `gibbs_uma` | QRRHO Gibbs energies (when `--thermo`) |

## R/TS/P canonical paths

Two locations get written for each elementary step:

| Path | Content |
|---|---|
| `<out_dir>/segments/seg_NN/{reactant,ts,product}.{xyz,pdb}` | **CANONICAL** — post-IRC re-optimized (RFO default; LBFGS via `--opt-mode-post grad`). `reactant` = backward endpoint, `ts` = tsopt'd, `product` = forward endpoint, all re-optimized |
| `<out_dir>/segments/seg_NN/structures/{reactant,ts,product}.{xyz,pdb}` | nested copy of canonical |
| `<out_dir>/segments/seg_NN/structures/{reactant,product}_irc.{xyz,pdb}` | raw IRC endpoints (pre-RFO/LBFGS re-optimization); for debugging IRC vs. post-IRC re-optimization divergence |

**Rule of thumb**: read from `segments/seg_NN/` for downstream stages. Use
`reactant_irc.xyz` / `product_irc.xyz` only when debugging
IRC vs. post-IRC re-optimization divergence.

`bond_changes` are computed from the MEP segment endpoints (first/last MEP
frames) in the default path-opt and path-search modes; only in tsopt-only
mode are they computed from the post-IRC `reactant.xyz` / `product.xyz`.

## Programmatic key extraction

```python
import json

# Per-segment barriers
d = json.load(open("result_all/summary.json"))
for seg in d["segments"]:
    print(f"seg_{seg['index']:02d}: ΔE‡ = {seg['barrier_kcal']:.1f} kcal/mol, "
          f"ΔE = {seg['delta_kcal']:.1f} kcal/mol")

# Rate-limiting barrier
rls = d["rate_limiting_step"]                # {"segment", "barrier_kcal", "method"}
idx = rls["segment"] - 1                     # segments[] is 0-indexed; rls["segment"] is 1-based
print(f"rate-limiting: seg_{rls['segment']:02d}, barrier = "
      f"{d['segments'][idx]['barrier_kcal']:.1f} kcal/mol")

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

Reading rules (cutoff: 1.20× covalent radii, margin 0.05; algorithm in [`pdb2reaction-cli/bond-summary.md`](../pdb2reaction-cli/bond-summary.md)):

- "Bond formed (k)" entries list bonds that exist in P but not R.
- "Bond broken (k)" entries list bonds that exist in R but not P.
- For a single reactive segment you usually expect 1–4 entries combined.
- If a single segment shows > 8 bond changes, the recursive
  segmentation may have failed — inspect the geometries before trusting
  the barrier.

## Failed-run output

When `summary.json["status"] != "success"`, look at:

1. `summary.log` — human-readable, prints the failure point first.
2. `segments/seg_NN/<stage>/result.json` — per-stage status (which step
   crashed) for the post-processing stages (ts / irc / freq / dft).
3. The per-stage log produced by the underlying tool (e.g. PySCF /
   pysisyphus) is silenced by default (see
   `default_group._silence_pysisyphus_loggers`); rerun with `--dump`
   to keep optimizer trajectories.

Even on failed runs, partial outputs are kept:

- `_work/path_opt/seg_NNN_<tag>/` (3-digit scratch; `_work/path_search/`
  with `--refine-path True`) exists for any segment that completed the MEP
  stage (even if downstream stages failed).
- `segments/seg_NN/` (deliverable) is **only populated** for
  fully-successful segments.

## See also

- `SKILL.md` — six canonical workflows and energy-diagram outputs.
- [`pdb2reaction-cli/all.md`](../pdb2reaction-cli/all.md) and the three `all-*.md` mode files — the subcommand that writes this `summary.json`.
- `pdb2reaction-cli/{tsopt,freq,irc,dft}.md` — per-stage `result.json` schemas (different shape from `summary.json`).
- [`pdb2reaction-cli/bond-summary.md`](../pdb2reaction-cli/bond-summary.md) — standalone bond-change tool with the same algorithm.