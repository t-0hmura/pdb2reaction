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

## Notes

- The generated file is a starting point, not an exhaustive schema.
- Runtime precedence is: `defaults < --config < CLI < --override-yaml`.

---

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [all](all.md) -- Run with `--config` / `--override-yaml`
- [YAML Reference](yaml-reference.md) -- Full configurable keys
