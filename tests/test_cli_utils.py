# tests/test_cli_utils.py
"""Tests for pdb2reaction.cli_utils."""

from pathlib import Path

import click
import pytest
import yaml

from pdb2reaction.cli_utils import (
    resolve_yaml_sources,
    load_merged_yaml_cfg,
    link_or_copy_file,
)


class TestResolveYamlSources:
    def test_no_sources(self):
        config, override, legacy = resolve_yaml_sources(None, None, None)
        assert config is None
        assert override is None
        assert legacy is False

    def test_config_only(self):
        p = Path("config.yaml")
        config, override, legacy = resolve_yaml_sources(p, None, None)
        assert config is p
        assert override is None
        assert legacy is False

    def test_override_only(self):
        p = Path("override.yaml")
        config, override, legacy = resolve_yaml_sources(None, p, None)
        assert config is None
        assert override is p
        assert legacy is False

    def test_legacy_only(self):
        p = Path("legacy.yaml")
        config, override, legacy = resolve_yaml_sources(None, None, p)
        assert config is None
        assert override is p
        assert legacy is True

    def test_config_and_override(self):
        c = Path("config.yaml")
        o = Path("override.yaml")
        config, override, legacy = resolve_yaml_sources(c, o, None)
        assert config is c
        assert override is o
        assert legacy is False

    def test_conflict_raises(self):
        with pytest.raises(click.BadParameter):
            resolve_yaml_sources(None, Path("a.yaml"), Path("b.yaml"))


class TestLoadMergedYamlCfg:
    def test_no_files(self):
        merged, cfg, ovr = load_merged_yaml_cfg(None, None)
        assert merged == {}
        assert cfg == {}
        assert ovr == {}

    def test_single_config(self, tmp_path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(yaml.dump({"opt": {"max_cycles": 100}}))
        merged, cfg, ovr = load_merged_yaml_cfg(cfg_file, None)
        assert merged["opt"]["max_cycles"] == 100
        assert ovr == {}

    def test_merge_override(self, tmp_path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(yaml.dump({"opt": {"max_cycles": 100, "thresh": "gau"}}))
        ovr_file = tmp_path / "override.yaml"
        ovr_file.write_text(yaml.dump({"opt": {"max_cycles": 50}}))
        merged, cfg, ovr = load_merged_yaml_cfg(cfg_file, ovr_file)
        assert merged["opt"]["max_cycles"] == 50  # Override wins
        assert merged["opt"]["thresh"] == "gau"  # Config preserved


class TestLinkOrCopyFile:
    def test_creates_symlink(self, tmp_path):
        src = tmp_path / "src.txt"
        src.write_text("hello")
        dst = tmp_path / "dst.txt"
        result = link_or_copy_file(src, dst)
        assert result is True
        assert dst.exists()
        assert dst.read_text() == "hello"

    def test_overwrites_existing(self, tmp_path):
        src = tmp_path / "src.txt"
        src.write_text("new")
        dst = tmp_path / "dst.txt"
        dst.write_text("old")
        result = link_or_copy_file(src, dst)
        assert result is True
        assert dst.read_text() == "new"

    def test_directory_destination_returns_false(self, tmp_path):
        src = tmp_path / "src.txt"
        src.write_text("hello")
        dst = tmp_path / "dst_dir"
        dst.mkdir()
        result = link_or_copy_file(src, dst)
        assert result is False
