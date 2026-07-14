# `summary.json` — schema and parsing

This page describes the aggregate file written by **`pdb2reaction all`** at
the top of `--out-dir`. It
is the canonical machine-readable artifact for downstream analysis
(barriers, ΔE, Gibbs, TS imaginary modes, file paths).

Do not apply this table to a standalone leaf command's `result.json` /
`summary.json`; those smaller stage and error envelopes are documented on the
corresponding command page and in [`pdb2reaction-cli`](../pdb2reaction-cli/SKILL.md).

## Top-level keys

| Key | Description |
|---|---|
| `schema_version` | Aggregate summary schema version; consumers should validate this before assuming keys or nested shapes |
| `command` | Full recorded invocation string for this `all` run (for example, `pdb2reaction all -i R.pdb -i P.pdb ...`) |
| `pdb2reaction_version` | Toolkit version that produced this output |
| `status` | `"success"`, `"partial"`, or `"failed"` |
| `charge` / `spin` | Resolved cluster charge / multiplicity |
| `environment` | `{device, gpu_name, gpu_vram_gb, cuda_version, cpu, n_cpus, ram_gb}` |
| `config` | Full effective config after CLI + YAML + defaults merge |
| `freeze_atoms` | Resolved union of YAML `geom.freeze_atoms` and auto-detected cap-H parents. These are 0-based internal indices; `summary.log` echoes them as 1-based. |
| `n_images` | MEP image count in path-search / path-opt; IRC trajectory frame count in tsopt-only mode |
| `n_segments` | Total MEP segment count (reactive + bridge) |
| `n_segments_reactive` | Non-bridge segment count used by the schema; validation is still required before calling each one an elementary reaction |
| `rate_limiting_step` | When available, `{segment, barrier_kcal, method, ...}` for the highest barrier at the highest method covering every reactive segment: DFT//MLIP Gibbs > DFT > MLIP Gibbs > MLIP > MEP |
| `overall_reaction_energy_kcal` | R → P total energy difference from the best available all-segment diagram |
| `overall_reaction_energy_method` | Method for `overall_reaction_energy_kcal` (`MEP`, `MLIP`, `MLIP_Gibbs`, `DFT`, or `DFT//MLIP_Gibbs`) |
| `segments` | Lightweight per-segment summary (see below) |
| `post_segments` | Per-segment post-processing details (tsopt / IRC / freq / DFT outputs) |
| `key_output_files` | Top-level entries map filename → description (str); `seg_NN` entries are `{description, files}` dicts |
| `pipeline_mode` | One of `"path-search"`, `"path-opt"`, `"tsopt-only"` |
| `mlip_backend` | Which backend produced the energies |
| `mlip_model` | Exact model/checkpoint identifier, kept separate from `mlip_backend` |
| `energy_diagrams` | Energy data plus optional image metadata. Check `image_written` and that `image` is non-null before treating a PNG as an artifact; energies remain available after renderer failure. |
| `out_dir` | Output directory absolute path |

## Per-segment keys (`summary.json["segments"][i]`)

Lightweight, MEP-level:

| Key | Description |
|---|---|
| `index` | Segment index (1-based int; written zero-padded as `seg_01/`, `seg_02/`, …) |
| `tag` | Segment-name string (`seg_NN` / `*_bridge`). To select reactive segments, filter `kind != "bridge"` |
| `kind` | Segment kind (`"seg"`, `"bridge"`, or `"tsopt"`) |
| `barrier_kcal` | MEP-level band: highest MEP image − first MEP image (kcal/mol). In `tsopt-only` mode it is instead refined TS − assigned R. Whenever `--tsopt` ran, match a `post_segments` entry by its `index` and report its `mlip.barrier_kcal` (or `gibbs_mlip.barrier_kcal` with thermo); list positions need not match when bridge/skipped segments exist. |
| `delta_kcal` | MEP-level last − first energy. Match `post_segments` by `index` for the refined `mlip.delta_kcal`. |
| `bond_changes` | List of single-key dicts: `[{"Bond formed (k)": ["A-B : 3.17 Å --> 1.68 Å", ...]}, {"Bond broken (k)": [...]}]`. **Standalone `irc result.json["bond_changes"]` uses a different shape**: `{"formed": [str], "broken": [str]}` (flat dict). The list-of-dicts form here is the `path_search` / `all` summary.json shape. See "Bond-change interpretation" below for cutoff / algorithm. |

## Per-segment post-processing keys (`summary.json["post_segments"][i]`)

Present when `--tsopt`, `--thermo`, or `--dft` was passed:

| Key | Description |
|---|---|
| `index` / `tag` / `kind` / `bond_changes` | Identify the matching `segments` row; match by `index`, not list position |
| `mep_barrier_kcal` / `mep_delta_kcal` | MEP-level (un-refined) barrier / ΔE, mirroring `segments[i]` (refined post-IRC/tsopt energies live in the per-segment `mlip` block) |
| `post_dir` | Subdirectory holding tsopt / freq / IRC outputs for this segment |
| `irc_plot` / `irc_traj` | Paths to the IRC trace PNG and trajectory XYZ |
| `mlip` | Selected MLIP backend's electronic-energy block. Read top-level `mlip_backend` / `mlip_model` for exact provenance. |
| `ts_imag` | Dict `{n_imag, nu_imag_max_cm, min_abs_imag_cm, min_freq_cm}` describing the TS spectrum |
| `ts_imag_freq_cm` | Peak imaginary frequency (cm⁻¹); same as `ts_imag.nu_imag_max_cm` |
| `gibbs_mlip` | MLIP electronic energy + QRRHO thermal correction (when `--thermo`) |
| `dft` | DFT//MLIP single-point energies (when `--dft`). Same shape as `mlip`: `{labels, energies_au, energies_kcal, barrier_kcal, delta_kcal, diagram, structures}`. If DFT failed for any of R/TS/P, the block is `{"status": "failed", "failed_states": [...]}` instead, and no DFT diagram is written. |
| `gibbs_dft_mlip` | DFT electronic energy + MLIP QRRHO thermal correction (when `--dft` **and** `--thermo`, and all three DFT single-points succeeded). Same shape as `mlip`; read `barrier_kcal` here for the DFT//MLIP ΔG‡. |

`mlip`, `gibbs_mlip`, and `gibbs_dft_mlip` are the only emitted identifiers;
no compatibility aliases are written. Energy-diagram filenames use `MLIP` as
the backend-neutral layer label, while top-level `mlip_backend` / `mlip_model`
record the exact provenance.

## R/TS/P canonical paths

Two locations are written for each successfully post-processed candidate
segment:

| Path | Content |
|---|---|
| `<out_dir>/segments/seg_NN/{reactant,ts,product}.{xyz,pdb,cif}` | **CANONICAL** — TSOPT TS plus post-IRC endpoint optimization (RFO default; LBFGS via `--opt-mode-post grad`). CIF is present only for bridged mmCIF/oversized-PDB topology and restores original IDs. In MEP modes, IRC ends are oriented to the MEP left/right states by bond pattern then RMSD. In TS-only mode, lacking endpoint references, the higher-energy end is named reactant and the lower product; inspect chemical identity. |
| `<out_dir>/segments/seg_NN/structures/{reactant,ts,product}.{xyz,pdb,cif}` | nested copy of canonical; CIF follows the same bridge condition |
| `<out_dir>/segments/seg_NN/structures/{reactant,product}_irc.{xyz,pdb,cif}` | raw IRC endpoints (pre-RFO/LBFGS re-optimization); for debugging IRC vs. post-IRC re-optimization divergence |

**Rule of thumb**: read from `segments/seg_NN/` for downstream stages. Use
`reactant_irc.xyz` / `product_irc.xyz` only when debugging
IRC vs. post-IRC re-optimization divergence.

`bond_changes` are computed from the MEP segment endpoints (first/last MEP
frames) in the default path-opt and path-search modes; only in tsopt-only
mode are they computed from the post-IRC `reactant.xyz` / `product.xyz`.

## Programmatic key extraction

```python
import json

# Reportable per-segment barriers (TSOPT + IRC refined; present when --tsopt ran)
d = json.load(open("result_all/summary.json"))
for ps in d.get("post_segments", []):
    mlip = ps.get("mlip") or {}
    gibbs = ps.get("gibbs_mlip") or {}
    if mlip.get("barrier_kcal") is None:
        continue
    line = (f"seg_{ps['index']:02d}: ΔE‡ = {mlip['barrier_kcal']:.1f} kcal/mol, "
            f"ΔE = {mlip['delta_kcal']:.1f} kcal/mol")
    if gibbs.get("barrier_kcal") is not None:
        line += f", ΔG‡ = {gibbs['barrier_kcal']:.1f} kcal/mol"
    print(line)

# MEP-level band (un-refined) — keep the "MEP" label when printing it
for seg in d["segments"]:
    print(f"seg_{seg['index']:02d}: MEP band = {seg['barrier_kcal']:.1f} kcal/mol (un-refined)")

# Rate-limiting barrier — highest method complete across all reactive segments;
# may fall back to MLIP or the MEP band when a higher layer is incomplete
rls = d["rate_limiting_step"]   # {"segment", "barrier_kcal", "method", "mep_barrier_kcal"}
print(f"rate-limiting: seg_{rls['segment']:02d}, barrier = "
      f"{rls['barrier_kcal']:.1f} kcal/mol ({rls['method']})")
# Always branch/report by rls["method"]; "MEP" is explicitly un-refined.

# Imaginary-mode check on every TS (post-processing data)
for ps in d.get("post_segments", []):
    ts = ps.get("ts_imag") or {}
    n = ts.get("n_imag")
    if n is not None and n != 1:
        print(f"WARNING: seg_{ps['index']:02d} has {n} imaginary modes "
              f"(freq {ps.get('ts_imag_freq_cm')!r})")

# Backend-level energies / Gibbs (when --tsopt and/or --thermo was passed)
for ps in d.get("post_segments", []):
    print(ps["index"], ps.get("mlip"), ps.get("gibbs_mlip"))
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

The `(k)` integer is the number of bonds listed in that section, so
`"Bond formed (2)"` means two bonds were formed and its list holds two
strings. A section with no change is emitted without `(k)`, as
`{"Bond formed": ["None"]}`. Each string carries
`<atom>-<atom> : <d_R> Å --> <d_P> Å`.

Reading rules (cutoff: 1.20× covalent radii, margin 0.05; algorithm in [`pdb2reaction-cli/bond-summary.md`](../pdb2reaction-cli/bond-summary.md)):

- "Bond formed (k)" entries list bonds that exist in P but not R.
- "Bond broken (k)" entries list bonds that exist in R but not P.
- An unexpected or large set of changes is not a numeric validator. Inspect
  the endpoint geometries and covalent-radius cutoff, then establish segment
  validity with TS optimization, an independent frequency calculation, and
  IRC connectivity.

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
  with `--refine-path`) exists for any segment that completed the MEP
  stage (even if downstream stages failed).
- `segments/seg_NN/` is created when post-processing begins and may contain
  partial files after a failed stage. Directory presence is not a success
  signal; use the aggregate and per-stage JSON statuses.

## See also

- `SKILL.md` — six canonical workflows and energy-diagram outputs.
- [`pdb2reaction-cli/all.md`](../pdb2reaction-cli/all.md) and the three `all-*.md` mode files — the subcommand that writes this `summary.json`.
- `pdb2reaction-cli/{tsopt,freq,irc,dft}.md` — per-stage `result.json` schemas (different shape from `summary.json`).
- [`pdb2reaction-cli/bond-summary.md`](../pdb2reaction-cli/bond-summary.md) — standalone bond-change tool with the same algorithm.
