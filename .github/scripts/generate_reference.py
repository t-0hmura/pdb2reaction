#!/usr/bin/env python3
"""Generate auto-derived CLI reference pages for docs.

Only the per-subcommand pages under ``docs/reference/commands/`` are
auto-generated. The YAML reference is the curated narrative at
``docs/yaml-reference.md`` and is NOT regenerated here, because the
narrative carries design rationale, cross-references, and per-section
context that cannot be recovered from ``defaults.py`` introspection.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import click
from click.testing import CliRunner

REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS_ROOT = REPO_ROOT / "docs"
REF_ROOT = DOCS_ROOT / "reference"
COMMANDS_ROOT = REF_ROOT / "commands"
TOOL_NAME = "pdb2reaction"

sys.path.insert(0, str(REPO_ROOT))

from pdb2reaction.cli import cli as root_cli  # noqa: E402


@dataclass(frozen=True)
class RenderedFile:
    path: Path
    content: str


@dataclass(frozen=True)
class CommandDoc:
    name: str
    stem: str


def _collect_subcommands() -> list[str]:
    ctx = click.Context(root_cli)
    return sorted(root_cli.list_commands(ctx))


def _doc_stem(command_name: str) -> str:
    # Keep command spelling in headings/tables, but standardize filenames.
    return command_name.replace("-", "_")


def _collect_command_docs() -> list[CommandDoc]:
    docs: list[CommandDoc] = []
    used_stems: dict[str, str] = {}
    for name in _collect_subcommands():
        stem = _doc_stem(name)
        prev = used_stems.get(stem)
        if prev is not None and prev != name:
            raise RuntimeError(
                f"Command doc stem collision: '{prev}' and '{name}' both map to '{stem}'."
            )
        used_stems[stem] = name
        docs.append(CommandDoc(name=name, stem=stem))
    return docs


def _capture_help(command_name: str, *, advanced: bool) -> str:
    runner = CliRunner()
    args = [command_name, "--help-advanced"] if advanced else [command_name, "--help"]
    # The CLI start-header guard inspects sys.argv to detect a help/version
    # request. Under CliRunner the help flag lives in `args`, not sys.argv, so
    # without this the banner / [command] / [mode] lines leak into the captured
    # reference. Point sys.argv at the help invocation so the guard suppresses
    # them (real-invocation CLI behavior is unchanged).
    saved_argv = sys.argv
    sys.argv = [TOOL_NAME, *args]
    try:
        result = runner.invoke(
            root_cli,
            args,
            catch_exceptions=False,
            prog_name=TOOL_NAME,
        )
    finally:
        sys.argv = saved_argv
    if result.exit_code != 0:
        raise RuntimeError(
            f"Failed to collect help for '{TOOL_NAME} {command_name}' "
            f"(advanced={advanced}):\n{result.output}"
        )
    if "[Unavailable]" in result.output:
        raise RuntimeError(
            f"Cannot generate help for '{TOOL_NAME} {command_name}': the lazy "
            "subcommand could not be imported. Install the repository's "
            "development/runtime dependencies and retry."
        )
    # Strip pre-banner noise (version line, bundled-pysisyphus
    # rc-file warning) so docs stay stable across environments where
    # `~/.pysisyphusrc` may or may not exist.
    lines = result.output.splitlines(keepends=True)
    lines = [
        ln for ln in lines
        if not ln.startswith(f"{TOOL_NAME} ver. ")
        and not ln.startswith("Couldn't find configuration file. Expected it at ")
    ]
    return "".join(lines).rstrip() + "\n"


def _render_command_page(command_name: str, help_text: str) -> str:
    return (
        f"# `{TOOL_NAME} {command_name}`\n\n"
        "```text\n"
        f"{help_text.rstrip()}\n"
        "```\n"
    )


def _render_commands_index(command_docs: list[CommandDoc]) -> str:
    toctree_body = "\n".join(doc.stem for doc in command_docs)
    rows = "\n".join(
        f"| `{TOOL_NAME} {doc.name}` | [{doc.name}]({doc.stem}.md) |" for doc in command_docs
    )
    return (
        "# CLI Command Reference\n\n"
        "These pages are the generated `--help` reference for each command. For usage "
        "guidance and worked examples, see the command's narrative guide in the main "
        "documentation (e.g. [all](../../all.md), [opt](../../opt.md)).\n\n"
        "```{toctree}\n"
        ":maxdepth: 1\n"
        ":hidden:\n\n"
        f"{toctree_body}\n"
        "```\n\n"
        "| Command | Page |\n"
        "|---|---|\n"
        f"{rows}\n"
    )


def _render() -> list[RenderedFile]:
    command_docs = _collect_command_docs()
    rendered: list[RenderedFile] = []

    total = len(command_docs)
    for index, doc in enumerate(command_docs, start=1):
        print(
            f"[reference] Rendering {index}/{total}: {TOOL_NAME} {doc.name}",
            flush=True,
        )
        try:
            help_text = _capture_help(doc.name, advanced=True)
        except RuntimeError:
            help_text = _capture_help(doc.name, advanced=False)
        rendered.append(
            RenderedFile(
                path=COMMANDS_ROOT / f"{doc.stem}.md",
                content=_render_command_page(doc.name, help_text),
            )
        )

    rendered.append(
        RenderedFile(
            path=COMMANDS_ROOT / "index.md",
            content=_render_commands_index(command_docs),
        )
    )
    return rendered


def _check_or_write(rendered: list[RenderedFile], *, check: bool) -> int:
    expected_paths = {item.path.resolve() for item in rendered}
    existing_paths = set(COMMANDS_ROOT.glob("*.md"))
    stale_extra = sorted(p.resolve() for p in existing_paths if p.resolve() not in expected_paths)

    stale: list[Path] = []
    for item in rendered:
        current = item.path.read_text(encoding="utf-8") if item.path.exists() else None
        if current != item.content:
            stale.append(item.path.resolve())

    if check:
        if stale or stale_extra:
            print("Reference files are out of date. Run: python .github/scripts/generate_reference.py")
            for path in sorted(stale):
                print(f"  MISSING/DIFF: {path.relative_to(REPO_ROOT)}")
            for path in stale_extra:
                print(f"  EXTRA: {path.relative_to(REPO_ROOT)}")
            return 1
        print("Reference files are up to date.")
        return 0

    for item in rendered:
        item.path.parent.mkdir(parents=True, exist_ok=True)
        item.path.write_text(item.content, encoding="utf-8")

    for path in stale_extra:
        path.unlink()

    print(f"Generated {len(rendered)} files under {REF_ROOT.relative_to(REPO_ROOT)}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Fail when generated files are stale.")
    args = parser.parse_args()
    rendered = _render()
    return _check_or_write(rendered, check=args.check)


if __name__ == "__main__":
    raise SystemExit(main())
