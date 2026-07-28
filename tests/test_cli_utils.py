"""Tests for pdb2reaction.cli.decorators."""

from pathlib import Path

import click
import pytest
import yaml

from pdb2reaction.cli.decorators import (
    resolve_yaml_sources,
    load_merged_yaml_cfg,
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

    def test_merge_keeps_source_layers_and_effective_tree_independent(self, tmp_path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(
            yaml.dump({"opt": {"max_cycles": 100, "nested": {"left": 1}}})
        )
        ovr_file = tmp_path / "override.yaml"
        ovr_file.write_text(
            yaml.dump({"opt": {"max_cycles": 50, "nested": {"right": 2}}})
        )

        merged, cfg, ovr = load_merged_yaml_cfg(cfg_file, ovr_file)

        assert cfg["opt"] == {"max_cycles": 100, "nested": {"left": 1}}
        assert ovr["opt"] == {"max_cycles": 50, "nested": {"right": 2}}
        assert merged["opt"] == {
            "max_cycles": 50,
            "nested": {"left": 1, "right": 2},
        }

        merged["opt"]["nested"]["left"] = 99
        assert cfg["opt"]["nested"]["left"] == 1
        assert "left" not in ovr["opt"]["nested"]
