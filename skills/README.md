# Agent Skills for `pdb2reaction`

This folder contains a set of skills that common AI agent interfaces
will recognize and help speed up code development by providing
concise instructions on how to use the `pdb2reaction` CLI for common
tasks.

- `pdb2reaction-overview`: what `pdb2reaction` is, when to use it,
  and how it differs from generic QM/MLIP path-search tools.
- `pdb2reaction-cli`: index of the 17 subcommands plus per-subcommand
  mds (each with synopsis, key flags, examples, output, caveats).
- `pdb2reaction-structure-io`: PDB / XYZ / GJF format references and
  the charge / multiplicity decision workflow.
- `pdb2reaction-install-backends`: install pdb2reaction itself, MLIP
  backends (UMA / Orb / MACE / AIMNet2), DFT (PySCF / GPU4PySCF), and
  xtb (ALPB implicit-solvent correction, not an MLIP backend);
  CUDA + PyTorch pairing.
- `pdb2reaction-workflows-output`: canonical workflows (cluster /
  multistep / scan-list / endpoint-MEP / TS-only / DFT//MLIP), output
  schema, and R/TS/P canonical paths.
- `pdb2reaction-hpc`: PBS / SLURM preamble templates with placeholders,
  walltime guidance, monitoring, plus a flock+pbsdsh dynamic-dispatch
  recipe.
- `pdb2reaction-env-detect`: fallback for detecting scheduler / GPU /
  CUDA / conda env when the environment is unknown.

Copying this `skills/` directory into another project (e.g.\ as
`.claude/skills/` for Claude Code, `~/.claude/skills/` for a
user-global install, or your agent's configured skill location)
gives an agent the orientation it needs to work with `pdb2reaction`. For version-sensitive defaults, advanced flags,
and backend model identifiers, agents should verify against the
installed CLI (`pdb2reaction <subcommand> --help-advanced`) and
`pdb2reaction.defaults`.
