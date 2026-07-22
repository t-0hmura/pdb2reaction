#!/usr/bin/env python3
"""Contract-check the advertised example shell scripts (M65).

Asserts the exact ``PUBLIC_SHELL_EXAMPLES`` paths exist, runs ``bash -n`` on
each, and validates their ``pdb2reaction`` invocations against the live CLI via
the shared command contract. Scientific execution stays in the explicit smoke
lanes; this lane only proves the advertised scripts parse and use real options.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import docs_command_contract as contract  # noqa: E402


def main() -> int:
    errors: list[str] = []
    total_cmds = 0
    for path in contract.PUBLIC_SHELL_EXAMPLES:
        rel = contract.rel_to_repo(path)
        if not path.exists():
            errors.append(f"advertised example script missing: {rel}")
            continue
        text = path.read_text(encoding="utf-8")
        lines = text.splitlines()
        if not lines or lines[0] != "#!/usr/bin/env bash":
            errors.append(f"{rel}: first line must be '#!/usr/bin/env bash'")
        if "set -euo pipefail" not in lines[:5]:
            errors.append(f"{rel}: missing 'set -euo pipefail' near the top")
        if not path.stat().st_mode & 0o111:
            errors.append(f"{rel}: script is not executable")
        if "BASH_SOURCE[0]" not in text:
            errors.append(f"{rel}: fixture paths are not anchored to the script directory")
        proc = subprocess.run(
            ["bash", "-n", str(path)], text=True, capture_output=True
        )
        if proc.returncode != 0:
            errors.append(f"bash -n failed for {rel}:\n{proc.stderr.strip()}")
            continue
        commands = contract.extract_shell_commands(path)
        if not commands:
            errors.append(f"no {contract.TOOL_NAME} invocations found in {rel}")
            continue
        total_cmds += len(commands)
        errors.extend(contract.validate_option_names(commands))

    if errors:
        print("[example-scripts] failed:")
        for err in errors:
            print(f"- {err}")
        return 1

    print(
        f"[example-scripts] validated {len(contract.PUBLIC_SHELL_EXAMPLES)} "
        f"advertised scripts, {total_cmds} CLI invocations."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
