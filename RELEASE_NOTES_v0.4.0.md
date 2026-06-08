# pdb2reaction v0.4.0

`pdb2reaction` is a Python CLI for elucidating **enzymatic reaction pathways**
from **PDB structures** using machine-learning interatomic potentials (MLIPs).
Given (i) two or more PDB files (R → ... → P), (ii) one PDB with `--scan-lists`,
or (iii) one TS candidate with `--tsopt`, it extracts an **active-site cluster
model**, runs an **MEP search**, and optionally chains **TS optimization → IRC →
frequencies → DFT single-point**. Each stage is also exposed as an individual
subcommand.

An initial reaction path is one command:

```bash
# Multi-PDB mode (R + P endpoints → MEP, with TS optimisation + thermo)
pdb2reaction all -i R.pdb P.pdb -q 0 --tsopt --thermo
```

## What's new in v0.4.0

Major release: the package is refactored into a clean layered architecture
(cli / workflows / domain / backends / io / core). The CLI surface is unchanged,
so existing commands keep working.

### Breaking
- Removed the flat-top compatibility shim. Import paths are now layered
  (`pdb2reaction.workflows`, `.domain`, `.backends`, `.io`, `.core`, `.cli`);
  top-level imports of stage modules no longer exist.
- Removed the `--trust-band` / `--hessian-window` / `--weighted-trust` flags.

### Added
- `--deterministic` on every compute subcommand for bit-reproducible GPU runs.
- `--precision fp32|fp64` on every calculator-constructing subcommand
  (fp64 also forces the Hessian to fp64).
- `--irc-pos-def` IRC convergence guard.
- Structured `result.json` / `summary.json` envelope (`schema_version`,
  error-class chain), mirrored to a single `summary.json` per subcommand.
- New docs: `reproducibility.md`, `output-layout.md`.

### Changed
- `--help` groups subcommands into semantic sections.
- AIMNet2 rejects `--precision fp64` and `--deterministic` with clear errors.
- Cleaner CLI error rendering with a recovery hint.

See `CHANGELOG.md` for the complete list.
