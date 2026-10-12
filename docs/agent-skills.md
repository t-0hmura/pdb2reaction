# Skills for AI Agents

`pdb2reaction` ships instructions for AI agents (Claude Code, Codex, Cursor, and others) in the `skills/` directory. With the skills loaded, an agent can run and analyze calculations from plain-language instructions.

## First example

To install them, tell your AI agent:

> Import `https://github.com/t-0hmura/pdb2reaction/tree/main/skills` as skills, and install pdb2reaction by following `pdb2reaction-install`.

If you cloned the GitHub repository, you can give the local `skills/` path instead. Then you can ask, for example:

> Read *the paper*, build a model from the PDB structure *PDB ID*, and study the mechanism of *the reaction step* with the pdb2reaction skills.

## Skills

Each folder in `skills/` is one skill (`<name>/SKILL.md`).

| Skill | What it covers |
| --- | --- |
| `pdb2reaction-overview` (start here) | What pdb2reaction is, how to pick one of the three `all` modes, step-by-step runs, TS strategy, and how to read the outputs |
| `pdb2reaction-cli` | When to use each subcommand, how to judge success, and pitfalls |
| `pdb2reaction-model-setup` | PDB / mmCIF / XYZ / GJF formats, the charge and multiplicity, and building the cluster model |
| `pdb2reaction-install` | Installing pdb2reaction, the MLIP backends, DFT, and xTB, and probing an unknown environment |
| `pdb2reaction-hpc` | PBS / SLURM job scripts, walltime, monitoring, and dynamic dispatch |
| `colab-local-gpu-runtime` | Running the Colab interface on a local NVIDIA GPU through WSL2 and Docker Desktop on Windows |

## Install by hand

Copy the skill folders into your agent's skill directory:

* Claude Code: `.claude/skills/` in a project, or `~/.claude/skills/` for all projects.
* Codex: `.agents/skills/` in a repository, or `~/.agents/skills/` for all repositories.

For example, from the root of the repository:

```bash
mkdir -p ~/.claude/skills
cp -r skills/pdb2reaction-* skills/colab-local-gpu-runtime ~/.claude/skills/
```

The skills leave exact flags and defaults to the installed CLI (`pdb2reaction <subcommand> --help-advanced`) and the [command reference](reference/commands/index.md).

## See also

* [Getting Started](getting-started.md) — the shortest run, and which page to read next
* [Installation](installation.md) — setting up the environment
* [`skills/README.md`](https://github.com/t-0hmura/pdb2reaction/blob/main/skills/README.md) — the skill list on GitHub
