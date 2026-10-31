"""Tests for pdb2reaction.cli.decorators."""

import yaml

from pdb2reaction.cli.decorators import load_merged_yaml_cfg


class TestLoadMergedYamlCfg:
    def test_no_files(self):
        merged, cfg = load_merged_yaml_cfg(None)
        assert merged == {}
        assert cfg == {}

    def test_single_config(self, tmp_path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(yaml.dump({"opt": {"max_cycles": 100}}))
        merged, cfg = load_merged_yaml_cfg(cfg_file)
        assert merged["opt"]["max_cycles"] == 100
        assert cfg == merged

    def test_effective_tree_is_independent_of_the_config_layer(self, tmp_path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(
            yaml.dump({"opt": {"max_cycles": 100, "nested": {"left": 1}}})
        )

        merged, cfg = load_merged_yaml_cfg(cfg_file)

        assert merged == cfg == {"opt": {"max_cycles": 100, "nested": {"left": 1}}}
        merged["opt"]["nested"]["left"] = 99
        assert cfg["opt"]["nested"]["left"] == 1

    def test_unknown_top_level_section_warns(self, tmp_path, capsys):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text("clac: {charge: 0}\n1: {}\nopt: {max_cycles: 10}\n")
        load_merged_yaml_cfg(cfg_file)
        err = capsys.readouterr().err
        assert "YAML section(s) 1, clac are not recognized" in err

    def test_shared_config_sections_do_not_warn(self, tmp_path, capsys):
        sections = (
            "geom", "calc", "opt", "lbfgs", "rfo", "rsirfo", "stopt", "gs", "dmf",
            "sp", "freq", "thermo", "hessian_dimer", "irc", "dft", "bond", "bias", "search",
        )
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(yaml.dump({name: {} for name in sections}))
        load_merged_yaml_cfg(cfg_file)
        assert capsys.readouterr().err == ""
