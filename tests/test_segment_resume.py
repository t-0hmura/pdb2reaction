from pathlib import Path

import click
import pytest

from pdb2reaction.workflows._run_session import ArtifactClaimError, InvocationManifest
from pdb2reaction.workflows._segment_resume import (
    build_resume_identity,
    invalidate_from_segment,
    retained_energy_diagrams,
    retained_post_segments,
    validate_resume_identity,
)


def _identity(tmp_path: Path):
    source = tmp_path / "r.xyz"
    source.write_text("1\nr\nH 0 0 0\n", encoding="utf-8")
    return build_resume_identity(
        inputs=[source], ref_pdb=None, pipeline_mode="path-opt",
        path_request={"radius": 2.6}, calculator_identity={"backend": "uma"},
        charge=0, spin=1,
    )


def test_resume_identity_requires_exact_path_state(tmp_path: Path) -> None:
    identity = _identity(tmp_path)
    validate_resume_identity(identity, identity, segment=2, available_segments=[1, 2])
    changed = {**identity, "charge": 1}
    with pytest.raises(click.UsageError, match="charge changed"):
        validate_resume_identity(identity, changed, segment=2, available_segments=[1, 2])
    with pytest.raises(click.BadParameter, match="not present"):
        validate_resume_identity(identity, identity, segment=3, available_segments=[1, 2])


def test_manifest_adopts_only_digest_verified_artifact(tmp_path: Path) -> None:
    path = tmp_path / "hei.xyz"
    path.write_text("current", encoding="utf-8")
    from pdb2reaction.workflows._segment_resume import file_sha256

    manifest = InvocationManifest()
    assert manifest.adopt_existing("path.hei.01", path, sha256_expected=file_sha256(path)) == path
    path.write_text("changed", encoding="utf-8")
    other = InvocationManifest()
    with pytest.raises(ArtifactClaimError, match="does not match"):
        other.adopt_existing("path.hei.01", path, sha256_expected="0" * 64)


def test_resume_invalidation_preserves_mep_and_earlier_segments(tmp_path: Path) -> None:
    for index in (1, 2, 3):
        segment = tmp_path / "segments" / f"seg_{index:02d}"
        segment.mkdir(parents=True)
        (segment / "result.txt").write_text(str(index), encoding="utf-8")
    (tmp_path / "mep_trj.xyz").write_text("mep", encoding="utf-8")
    (tmp_path / "energy_diagram_MLIP_all.png").write_bytes(b"old")
    invalidate_from_segment(tmp_path, 2)
    assert (tmp_path / "segments" / "seg_01").is_dir()
    assert not (tmp_path / "segments" / "seg_02").exists()
    assert not (tmp_path / "segments" / "seg_03").exists()
    assert (tmp_path / "mep_trj.xyz").is_file()
    assert not (tmp_path / "energy_diagram_MLIP_all.png").exists()


def test_resume_retains_only_completed_prefix_metadata(tmp_path: Path) -> None:
    summary = {
        "post_segments": [{"index": 1}, {"index": 2}],
        "energy_diagrams": [
            {"image": str(tmp_path / "segments/seg_01/a.png")},
            {"image": str(tmp_path / "segments/seg_02/a.png")},
            {"image": str(tmp_path / "energy_diagram_MLIP_all.png")},
        ],
    }
    assert retained_post_segments(summary, 2) == [{"index": 1}]
    assert retained_energy_diagrams(summary, 2) == [summary["energy_diagrams"][0]]
