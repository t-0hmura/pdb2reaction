"""Public output declarations must stay inside the pipeline root on disk."""

from __future__ import annotations

from pathlib import Path

import pytest

from pdb2reaction.workflows._run_session import (
    InvocationManifest,
    declare_public_output,
)


@pytest.mark.parametrize("linked", ["segments", "summary.json"])
def test_public_declaration_rejects_symlink_resolving_outside_root(
    tmp_path: Path, linked: str,
) -> None:
    root = tmp_path / "out"
    external = tmp_path / "external"
    root.mkdir()
    external.mkdir()
    if linked == "segments":
        (root / "segments").symlink_to(external, target_is_directory=True)
        destination = root / "segments" / "seg_01" / "result.json"
    else:
        (root / "summary.json").symlink_to(external / "summary.json")
        destination = root / "summary.json"
    manifest = InvocationManifest()

    with pytest.raises(ValueError, match="resolves outside pipeline root"):
        declare_public_output(manifest, root, destination)
    assert manifest.expected == {}
