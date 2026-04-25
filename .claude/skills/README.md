# Agent Skills for `pdb2reaction`

This folder contains a set of skills that common AI agent interfaces
will recognize and help speed up code development by providing
concise instructions on how to use the `pdb2reaction` API and CLI for
common tasks. Inspired by the [nvalchemi-toolkit][nvalchemi] skill
pattern.

[nvalchemi]: https://github.com/NVIDIA/nvalchemi-toolkit

- `pdb2reaction-overview`: what `pdb2reaction` is, when to use it,
  and how it differs from generic QM/MLIP path-search tools.
- `pdb2reaction-cli`: index of the 17 subcommands plus per-subcommand
  mds (each with synopsis, key flags, examples, output, caveats).
- `pdb2reaction-structure-io`: PDB / XYZ / GJF format references and
  the charge / multiplicity decision workflow.
- `pdb2reaction-install-backends`: install pdb2reaction itself, MLIP
  backends (UMA / Orb / MACE / AIMNet2), DFT (PySCF / GPU4PySCF), and
  xtb; CUDA + PyTorch pairing.
- `pdb2reaction-workflows-output`: canonical workflows (cluster /
  multistep / scan-list / endpoint-MEP / TS-only / DFT//MLIP), output
  schema, and R/TS/P canonical paths.
- `pdb2reaction-hpc`: PBS / SLURM preamble templates with placeholders,
  walltime guidance, monitoring, plus a flock+pbsdsh dynamic-dispatch
  recipe.
- `pdb2reaction-env-detect`: fallback for detecting scheduler / GPU /
  CUDA / conda env when the environment is unknown.

The skills are **self-contained**: copying this `.claude/skills/`
directory into another project (or `~/.claude/skills/`) gives an agent
everything it needs to work with `pdb2reaction` without consulting the
main documentation.
