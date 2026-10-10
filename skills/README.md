# Agent Skills for `pdb2reaction`

These skills tell an AI agent which `pdb2reaction` command to run and when,
how to judge whether a run succeeded, which pitfalls to expect and how to
recover, and how to read the outputs.

- `pdb2reaction-overview` (start here): what `pdb2reaction` is, picking one of the three
  `all` modes (Endpoint mode, Scan-list mode, TS-only mode), step-by-step runs, a code map,
  TS strategy (`ts-strategy.md`), and reading outputs (`outputs.md`).
- `pdb2reaction-cli`: the 18 subcommands in 14 files (heavy commands one each,
  the rest grouped): when to use each, how to judge success, and pitfalls.
  Exact flags are left to `--help-advanced` and the command reference.
- `pdb2reaction-mcp`: how to drive `pdb2reaction` from any MCP client
  (Claude Desktop / Claude Code / Cursor / custom SDK) via the bundled
  `pdb2reaction-mcp` server; lists the 18 MCP tools and the result shape
  shared by every tool.
- `pdb2reaction-model-setup`: PDB / mmCIF / XYZ / GJF formats, the charge /
  multiplicity decision, and building the cluster model: which residues
  `extract` keeps, cutting and capping the boundary, trimming or enlarging it,
  and keeping the same atoms across states and variants
  ([full guide](../docs/model-setup.md)).
- `pdb2reaction-install`: install `pdb2reaction` itself, MLIP
  backends (UMA / ORB / MACE / AIMNet2), DFT (PySCF / GPU4PySCF), and xtb
  (ALPB implicit-solvent correction, not an MLIP backend); CUDA + PyTorch
  pairing; probing the scheduler, GPU, CUDA, and conda env when the
  environment is unknown.
- `pdb2reaction-hpc`: PBS / SLURM preamble templates with placeholders,
  walltime guidance, monitoring, plus a flock+pbsdsh dynamic-dispatch
  recipe.
- `colab-local-gpu-runtime`: Windows setup and operation for running the Colab
  interface on a local NVIDIA GPU through WSL2 and Docker Desktop.

## Install

Each folder here is one skill (`<name>/SKILL.md`). Copy the folders into your agent's skill directory:

- Claude Code: `.claude/skills/` in a project, or `~/.claude/skills/` for all projects.
- Codex: `.agents/skills/` in a repository, or `~/.agents/skills/` for all repositories.

For example, from the root of this repository:

```bash
mkdir -p ~/.claude/skills
cp -r skills/pdb2reaction-* skills/colab-local-gpu-runtime ~/.claude/skills/
```

Links from the skills to `../../docs/` open only inside a repository checkout; elsewhere, use the docs at <https://t-0hmura.github.io/pdb2reaction/>.

For exact flags and defaults, check the installed CLI (`pdb2reaction <subcommand> --help-advanced`) and the [command reference](../docs/reference/commands/index.md).
