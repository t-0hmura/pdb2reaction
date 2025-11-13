# `scan` subcommand

## Purpose
Performs staged bond-length scans with harmonic restraints and single-structure relaxations after each step, using UMA plus pysisyphus optimizers.

## Usage
```bash
pdb2reaction scan -i INPUT -q CHARGE --scan-lists "[(i,j,target), ...]" [...]
                  [--one-based/--zero-based] [--max-step-size ΔÅ] [--bias-k k]
                  [--relax-max-cycles N] [--opt-mode light|lbfgs|heavy|rfo]
                  [--freeze-links/--no-freeze-links] [--dump/--no-dump]
                  [--out-dir DIR] [--preopt/--no-preopt] [--endopt/--no-endopt]
                  [--args-yaml FILE]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge. | Required |
| `-s, --spin INT` | Spin multiplicity (2S+1). | `1` |
| `--scan-lists TEXT` | Repeatable Python-like list literal of `(i, j, targetÅ)` triples defining each scan stage. | Required |
| `--one-based / --zero-based` | Interpret `(i, j)` indices as 1-based (default) or 0-based. | `--one-based` |
| `--max-step-size FLOAT` | Maximum change in any scanned bond per integration step (Å). | `0.20` |
| `--bias-k FLOAT` | Harmonic bias strength `k` (eV·Å⁻²). Overrides `bias.k`. | `100` |
| `--relax-max-cycles INT` | Maximum optimizer cycles per step (overrides `opt.max_cycles`). | `10000` |
| `--opt-mode TEXT` | Relaxation optimizer (`light|lbfgs` or `heavy|rfo`). | `light` |
| `--freeze-links / --no-freeze-links` | Freeze link-hydrogen parents for PDB inputs. | `--freeze-links` |
| `--dump / --no-dump` | Dump stage trajectories. | `--no-dump` |
| `--out-dir TEXT` | Output directory. | `./result_scan/` |
| `--args-yaml FILE` | YAML overrides (see below). | _None_ |
| `--preopt / --no-preopt` | Pre-optimize the initial structure before scanning. | `--preopt` |
| `--endopt / --no-endopt` | Unbiased relaxation after each stage. | `--endopt` |

## YAML configuration (`--args-yaml`)
The YAML root must be a mapping. CLI parameters override YAML. Shared sections reuse the definitions documented for [`opt`](opt.md#yaml-configuration-args-yaml).

### Shared sections
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: see [`opt`](opt.md#yaml-configuration-args-yaml) for all keys and defaults. `opt.dump` is internally forced to `False`; dumping is controlled by `--dump`.

### Section `bias`
Controls the harmonic restraint applied during each scan stage.

- `k` (`100`): Harmonic strength in eV·Å⁻².

### Section `bond`
Parameters for UMA-based bond-change detection (mirrors `path_search`).

- `device` (`"cuda"`): UMA device for bond analysis.
- `bond_factor` (`1.20`): Covalent-radius scale factor for bond cutoff.
- `margin_fraction` (`0.05`): Fractional tolerance when comparing bond lengths.
- `delta_fraction` (`0.05`): Minimum fractional change to classify a bond as forming/breaking.

## Outputs
- `<out-dir>/stage_XX/scan.trj` (and `.pdb` when the input was PDB and dumping is enabled).
- `<out-dir>/stage_XX/scan_summary.yaml` with per-step metadata.
- Final unbiased geometries after `--endopt` stored per stage.
- Console summaries of resolved `geom`, `calc`, `opt`, `bias`, `bond`, and optimizer blocks.

## Notes
- `--scan-lists` accepts multiple literals; each defines one stage, and tuple indices are normalized to 0-based internally.
- When `--one-based` is used (default), indices follow PDB conventions; invalid indices or non-positive target distances raise `click.BadParameter`.
- Pre- and end-of-stage optimizations share UMA calculator instances for efficiency.
- Stage dumping writes one trajectory per stage; `--dump` also triggers `.pdb` conversion for PDB inputs.
