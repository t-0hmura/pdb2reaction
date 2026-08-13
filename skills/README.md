# Agent Skills for `pdb2reaction`

This folder contains a set of skills that common AI agent interfaces
will recognize and help speed up code development by providing
concise instructions on how to use the `pdb2reaction` CLI for common
tasks.

- `pdb2reaction-overview`: what `pdb2reaction` is, when to use it,
  and how it differs from generic QM/MLIP path-search tools.
- `pdb2reaction-architecture`: 6-layer package map (cli / workflows /
  domain / backends / io / core) + bundled-fork policy; tells an agent
  which dir to grep for a given concern before touching code.
- `pdb2reaction-cli`: index of the 18 subcommands plus per-subcommand
  mds (each with synopsis, key flags, examples, output, caveats).
- `pdb2reaction-ts-strategy`: decision know-how for a reaction-barrier
  run — backend- and purpose-aware precision (fp32 vs fp64), the two TS-candidate
  routes (`path-search` MEP vs distance-restrained `scan`), fixing a bad
  imaginary-mode count (`--precision fp64` / `--coord-type dlc` /
  `--flatten`), reading a P-start barrier as the reverse direction,
  staged vs concerted scans, and the distinction between within-path atom-order
  requirements and controlled cross-variant comparisons.
- `pdb2reaction-mcp`: how to drive `pdb2reaction` from any MCP client
  (Claude Desktop / Claude Code / Cursor / custom SDK) via the bundled
  `pdb2reaction-mcp` server; lists the 18 MCP tools and the shared
  `SubcmdResult` return schema.
- `pdb2reaction-structure-io`: PDB / mmCIF / XYZ / GJF format references and
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
- `colab-local-gpu-runtime`: Windows setup and operation for running the Colab
  interface on a local NVIDIA GPU through WSL2 and Docker Desktop.

Copying this `skills/` directory into another project (e.g. as
`.claude/skills/` for Claude Code, `~/.claude/skills/` for a
user-global install, or your agent's configured skill location)
gives an agent the orientation it needs to work with `pdb2reaction`. For version-sensitive defaults, advanced flags,
and backend model identifiers, agents should verify against the
installed CLI (`pdb2reaction <subcommand> --help-advanced`) and
`pdb2reaction.core.defaults`.
