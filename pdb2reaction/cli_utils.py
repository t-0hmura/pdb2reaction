from __future__ import annotations

from typing import List, Sequence


def collect_option_values(argv: Sequence[str], names: Sequence[str]) -> List[str]:
    """
    Collect values following any of the option names, including sloppy/eat-all forms.

    Example:
        collect_option_values(["--scan-lists", "stage1", "stage2"], ("--scan-list", "--scan-lists"))
        -> ["stage1", "stage2"]
    """
    vals: List[str] = []
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok in names:
            j = i + 1
            while j < len(argv) and not argv[j].startswith("-"):
                vals.append(argv[j])
                j += 1
            i = j
        else:
            i += 1
    return vals
