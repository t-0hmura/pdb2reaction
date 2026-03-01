# `init`

## Overview

`pdb2reaction init` writes a starter YAML configuration file for `pdb2reaction all`.

Use this when you want a reproducible config-first workflow:

```bash
pdb2reaction init --out pdb2reaction_all.config.yaml
pdb2reaction all --config pdb2reaction_all.config.yaml --dry-run
```

## Usage

```bash
pdb2reaction init [--out PATH] [--force]
```

## Options

| Option | Description | Default |
| --- | --- | --- |
| `-o, --out PATH` | Destination YAML path. | `pdb2reaction_all.config.yaml` |
| `--force/--no-force` | Overwrite existing file. | `False` |

## Generated template

The generated file contains the following starter configuration:

```yaml
# Starter config for `pdb2reaction all`

extract:
  radius: 2.6
  radius_het2het: 0.0

calc:
  workers: 1
  workers_per_node: 1

path_search:
  mep_mode: gsm
  max_nodes: 10
  max_cycles: 300

scan:
  max_step_size: 0.2
  bias_k: 300.0
  relax_max_cycles: 10000

tsopt:
  max_cycles: 10000

freq:
  max_write: 10
  amplitude_ang: 0.8
  n_frames: 20
  sort: value
  temperature: 298.15
  pressure: 1.0

dft:
  func_basis: wb97m-v/def2-tzvpd
  max_cycle: 100
  conv_tol: 1.0e-9
  grid_level: 3
```

## Next steps

1. Edit the generated YAML to match your system (charge, multiplicity, substrate, etc.).
2. Validate with a dry run: `pdb2reaction all --config pdb2reaction_all.config.yaml --dry-run -i R.pdb P.pdb -c "SUBSTRATE"`.
3. Once satisfied, remove `--dry-run` to execute the full workflow.
4. For advanced tuning, see [YAML Reference](yaml_reference.md) for the complete set of configurable keys.

## Notes

- The generated file is a starting point, not an exhaustive schema.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [YAML Reference](yaml_reference.md) -- Full configurable keys
