import json
from pathlib import Path

import click
import pytest
from click.testing import CliRunner

from pdb2reaction.cli.app import cli


def test_dft_rejects_redundant_backend_selector() -> None:
    result = CliRunner().invoke(cli, ["dft", "-b", "dft"])

    assert result.exit_code == 2
    assert "already selects the DFT evaluator" in result.output
    assert "sp -b dft" in result.output


def test_calculator_leaf_help_describes_backend_specific_solvent() -> None:
    result = CliRunner().invoke(cli, ["sp", "--help-advanced"])

    assert result.exit_code == 0, result.output
    normalized_help = " ".join(result.output.split())
    assert "xTB solvent delta" in normalized_help
    assert "native PySCF PCM/SMD" in normalized_help


def test_dft_resource_defaults_are_lowmem_and_explicit_values_normalize(
    monkeypatch,
) -> None:
    from pdb2reaction.core.dft_settings import resolve_dft_settings

    monkeypatch.setenv("OMP_NUM_THREADS", "6")
    settings = resolve_dft_settings({"backend": "dft"})
    explicit = resolve_dft_settings({
        "backend": "dft",
        "dft": {"nprocs": 3, "memory": "64GB"},
    })

    assert settings.func_basis == "wb97m-v/def2-svp"
    assert settings.lowmem is True
    assert settings.density_fit is False
    assert settings.nprocs <= 6
    assert explicit.nprocs == 3
    assert explicit.memory_mb == 64000


def test_standalone_dft_accepts_pyscf_object_defaults_without_false_conflict(
    tmp_path,
) -> None:
    config = tmp_path / "pyscf.yaml"
    config.write_text(
        "dft:\n"
        "  pyscf:\n"
        "    mf:\n"
        "      conv_tol: 2.0e-8\n"
        "    grids:\n"
        "      level: 1\n"
        "    mol:\n"
        "      max_memory: 4096\n",
        encoding="utf-8",
    )
    gjf = Path(__file__).resolve().parent / "smoke" / "h2.gjf"

    result = CliRunner().invoke(
        cli,
        [
            "dft", "-v", "3", "-i", str(gjf), "--engine", "cpu", "--config", str(config),
            "--show-config", "--dry-run", "--out-dir", str(tmp_path / "result"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "conv_tol: 2.0e-08" in result.output
    assert "grid_level: 1" in result.output
    assert "memory_mb: 4096" in result.output


def test_dft_resources_are_provenance_but_not_checkpoint_identity() -> None:
    from pdb2reaction.core.dft_settings import resolve_dft_settings
    from pdb2reaction.core.utils import calculator_provenance

    small = resolve_dft_settings({
        "backend": "dft", "dft": {"nprocs": 2, "memory": "4GB"}
    })
    large = resolve_dft_settings({
        "backend": "dft", "dft": {"nprocs": 16, "memory": "64GB"}
    })
    provenance = calculator_provenance({
        "backend": "dft", "dft_settings": large.to_dict()
    })

    assert small.scientific_identity() == large.scientific_identity()
    assert provenance["dft_resources"] == {
        "memory_mode": "gpu4pyscf_rks_lowmem",
        "nprocs": 16,
        "nprocs_source": "explicit",
        "memory_mb": 64000,
        "memory_source": "explicit",
    }


@pytest.mark.parametrize(
    ("dft", "spin", "expected"),
    [
        ({}, 1, "gpu4pyscf_rks_lowmem"),
        ({"engine": "cpu"}, 1, "direct_jk"),
        ({}, 2, "direct_jk"),
        ({"lowmem": False}, 1, "density_fit"),
    ],
)
def test_dft_memory_mode_tracks_the_effective_driver(dft, spin, expected) -> None:
    from pdb2reaction.core.dft_settings import resolve_dft_settings

    calc = {"backend": "dft", "spin": spin, "dft": dft}
    assert resolve_dft_settings(calc).memory_mode == expected


@pytest.mark.parametrize(
    "key", ["charge", "multiplicity", "embedcharge", "embedcharge_cutoff"]
)
def test_calculator_dft_mapping_rejects_top_level_owned_state(key) -> None:
    from pdb2reaction.core.dft_settings import resolve_dft_settings

    with pytest.raises(click.BadParameter, match=f"calc.dft.{key}"):
        resolve_dft_settings({"backend": "dft", "dft": {key: 1}})


def test_dft_yaml_nprocs_reports_a_click_validation_error() -> None:
    from pdb2reaction.core.dft_settings import resolve_dft_settings

    with pytest.raises(click.BadParameter, match="positive integer"):
        resolve_dft_settings({"backend": "dft", "dft": {"nprocs": "many"}})


def test_def2_ecp_is_only_auto_enabled_for_covered_elements() -> None:
    from pdb2reaction.workflows.dft import _def2_ecp_required

    atomic_numbers = {"H": 1, "C": 6, "O": 8, "Rb": 37, "I": 53}
    charge = atomic_numbers.__getitem__

    assert not _def2_ecp_required(
        [("H", (0, 0, 0)), ("C", (0, 0, 1)), ("O", (0, 1, 0))],
        charge,
    )
    assert _def2_ecp_required([("Rb", (0, 0, 0)), ("I", (0, 0, 3))], charge)


def test_detailed_pyscf_ecp_is_preserved_and_conflicts_fail_closed() -> None:
    from pdb2reaction.workflows.dft import _resolve_effective_ecp

    assert _resolve_effective_ecp(
        None,
        {},
        "sto-3g",
        [("H", (0, 0, 0))],
        lambda _symbol: 1,
    ) is None

    configured = {"I": "def2-tzvp"}
    effective = _resolve_effective_ecp(
        None,
        configured,
        "def2-svp",
        [("I", (0, 0, 0))],
        lambda _symbol: 53,
    )
    assert effective is configured

    with pytest.raises(click.BadParameter, match="conflicts"):
        _resolve_effective_ecp(
            "def2-svp",
            configured,
            "def2-svp",
            [("I", (0, 0, 0))],
            lambda _symbol: 53,
        )


@pytest.mark.parametrize(
    ("model", "detail", "expected"),
    [
        ("pcm", {"eps": 12.5}, ("pcm", 12.5)),
        ("smd", {"solvent": "acetonitrile"}, ("smd", "acetonitrile")),
    ],
)
def test_detailed_pyscf_solvent_attributes_reach_native_model(
    model, detail, expected
) -> None:
    from pdb2reaction.workflows.dft import _apply_implicit_solvent

    class Solvent:
        eps = None
        solvent = None

    class MF:
        def __init__(self):
            self.with_solvent = Solvent()
            self.applied = None

        def PCM(self):
            self.applied = "pcm"
            return self

        def SMD(self):
            self.applied = "smd"
            return self

    mf = _apply_implicit_solvent(
        MF(),
        {
            "solvent": "water",
            "solvent_model": model,
            "pyscf": {"with_solvent": detail},
        },
    )

    attribute = "eps" if model == "pcm" else "solvent"
    assert (mf.applied, getattr(mf.with_solvent, attribute)) == expected


def test_dft_rejects_unknown_yaml_engine_before_execution(tmp_path) -> None:
    xyz = tmp_path / "h2.xyz"
    xyz.write_text(
        "2\nH2\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n",
        encoding="utf-8",
    )
    config = tmp_path / "config.yaml"
    config.write_text("dft:\n  engine: quantum-potato\n", encoding="utf-8")

    result = CliRunner().invoke(
        cli,
        [
            "dft",
            "-i",
            str(xyz),
            "-q",
            "0",
            "--config",
            str(config),
            "--dry-run",
        ],
    )

    assert result.exit_code == 2
    assert "dft.engine must be either 'cpu' or 'gpu'" in result.output


def test_dft_nonconvergence_commits_json_before_exit(tmp_path) -> None:
    from pdb2reaction.workflows.dft import _finalize_dft_result

    with pytest.raises(SystemExit) as caught:
        _finalize_dft_result(
            out_json=True,
            out_dir=tmp_path,
            payload={"status": "not_converged", "converged": False},
            elapsed_seconds=1.0,
        )

    assert caught.value.code == 1
    payload = json.loads((tmp_path / "result.json").read_text())
    assert payload["scientific_status"] == "failed"
    assert payload["converged"] is False


def test_prepare_dft_output_dir_invalidates_prior_public_results(tmp_path) -> None:
    from pdb2reaction.workflows.dft import _prepare_dft_output_dir

    for name in ("result.yaml", "result.json", "summary.json"):
        (tmp_path / name).write_text("stale\n", encoding="utf-8")

    _prepare_dft_output_dir(tmp_path)

    assert all(
        not (tmp_path / name).exists()
        for name in ("result.yaml", "result.json", "summary.json")
    )


def test_prepare_dft_output_dir_rejects_input_alias_before_mutation(tmp_path) -> None:
    from pdb2reaction.workflows.dft import _prepare_dft_output_dir

    result = tmp_path / "result.yaml"
    result.write_text("input: retained\n", encoding="utf-8")
    alias = tmp_path / "input.yaml"
    alias.hardlink_to(result)

    with pytest.raises(click.UsageError, match="collides with reserved DFT output"):
        _prepare_dft_output_dir(tmp_path, protected_inputs=(alias,))

    assert result.read_text(encoding="utf-8") == "input: retained\n"


def test_dft_unexpected_config_failure_uses_yaml_effective_output(
    monkeypatch,
    tmp_path,
) -> None:
    from pdb2reaction.workflows import dft

    xyz = tmp_path / "h2.xyz"
    xyz.write_text(
        "2\nH2\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n",
        encoding="utf-8",
    )
    effective_out = tmp_path / "configured-output"
    config = tmp_path / "config.yaml"
    config.write_text(
        "dft:\n"
        f"  out_dir: {effective_out}\n"
        "  func_basis: wb97m-v/def2-tzvpd\n",
        encoding="utf-8",
    )

    def fail_after_output_resolution(_value):
        raise RuntimeError("config probe failed")

    monkeypatch.setattr(dft, "_parse_func_basis", fail_after_output_resolution)
    result = CliRunner().invoke(
        cli,
        [
            "dft",
            "-i",
            str(xyz),
            "-q",
            "0",
            "--config",
            str(config),
            "--dry-run",
        ],
    )

    assert result.exit_code == 1
    payload = json.loads((effective_out / "result.json").read_text())
    assert payload["execution_status"] == "failed"
    assert payload["command"] == "dft"
    assert payload["error"] == "config probe failed"


@pytest.mark.parametrize(
    ("label", "config_text", "expected"),
    [
        ("syntax", "calc: [1,\n", "invalid YAML"),
        ("root", "- 1\n- 2\n", "must be a mapping"),
        ("section", "geom: 5\n", "YAML section 'geom' must be a mapping"),
    ],
)
def test_malformed_config_yaml_is_an_input_error(
    tmp_path, label, config_text, expected
) -> None:
    """A malformed --config file exits 2 with a message, not a raw traceback.

    The YAML layer is loaded before each command's own exception rendering, so
    the loader and the section resolver own these diagnostics.
    """
    xyz = tmp_path / "h2.xyz"
    xyz.write_text("2\nH2\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n", encoding="utf-8")
    config = tmp_path / f"{label}.yaml"
    config.write_text(config_text, encoding="utf-8")

    result = CliRunner().invoke(
        cli,
        [
            "sp",
            "-i",
            str(xyz),
            "-q",
            "0",
            "--config",
            str(config),
            "--dry-run",
            "-o",
            str(tmp_path / f"out_{label}"),
        ],
    )

    assert result.exit_code == 2
    assert expected in result.output
    assert "Traceback" not in result.output


def test_leaf_dft_checkpoint_defaults_under_output_directory(tmp_path) -> None:
    from pdb2reaction.core.dft_settings import (
        DFT_CLI_META_KEY,
        finalize_dft_calculator_config,
    )

    command = click.Command("sp")
    ctx = click.Context(command, info_name="sp")
    ctx.params["out_dir"] = tmp_path / "result"
    ctx.meta[DFT_CLI_META_KEY] = {"save_scf_checkpoint": True}
    calc_cfg = {"backend": "dft", "charge": 0, "spin": 1}

    finalize_dft_calculator_config(ctx, calc_cfg)

    assert calc_cfg["dft_settings"]["checkpoint_path"] == str(
        tmp_path / "result" / "_work" / "dft_scf" / "state.chk"
    )


def test_leaf_dft_checkpoint_uses_yaml_effective_output_directory(tmp_path) -> None:
    from pdb2reaction.core.dft_settings import (
        DFT_CLI_META_KEY,
        finalize_dft_calculator_config,
    )

    ctx = click.Context(click.Command("sp"), info_name="sp")
    ctx.params["out_dir"] = tmp_path / "click-default"
    ctx.meta[DFT_CLI_META_KEY] = {"save_scf_checkpoint": True}
    calc_cfg = {"backend": "dft", "charge": 0, "spin": 1}
    effective_out = tmp_path / "yaml-output"

    finalize_dft_calculator_config(ctx, calc_cfg, output_dir=effective_out)

    assert calc_cfg["dft_settings"]["checkpoint_path"] == str(
        effective_out / "_work" / "dft_scf" / "state.chk"
    )


@pytest.mark.parametrize(
    "dft_config",
    [
        {"lowmem": "false"},
        {"density_fit": "false"},
        {"save_scf_checkpoint": "false"},
        {"scf_stepwise_grid": "false"},
        {"pyscf": {"density_fit": {"enabled": "false"}}},
    ],
)
def test_dft_yaml_booleans_require_yaml_boolean_type(dft_config) -> None:
    from pdb2reaction.core.dft_settings import resolve_dft_settings

    with pytest.raises(click.BadParameter, match="must be true or false"):
        resolve_dft_settings({"backend": "dft", "dft": dft_config})


@pytest.mark.parametrize("pbs_var", ["PBS_NP", "PBS_NUM_PPN"])
def test_dft_resources_honor_pbs_cpu_counts(monkeypatch, pbs_var) -> None:
    from pdb2reaction.core import dft_settings

    for name in (
        "OMP_NUM_THREADS",
        "SLURM_CPUS_PER_TASK",
        "NSLOTS",
        "PBS_NP",
        "PBS_NUM_PPN",
        "PBS_NODEFILE",
    ):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv(pbs_var, "3")
    monkeypatch.setattr(dft_settings, "_affinity_count", lambda: None)

    assert dft_settings.resolve_dft_settings({"backend": "dft"}).nprocs == 3


def test_dft_resources_honor_pbs_nodefile(monkeypatch, tmp_path) -> None:
    from pdb2reaction.core import dft_settings

    for name in (
        "OMP_NUM_THREADS",
        "SLURM_CPUS_PER_TASK",
        "NSLOTS",
        "PBS_NP",
        "PBS_NUM_PPN",
        "PBS_NODEFILE",
    ):
        monkeypatch.delenv(name, raising=False)
    nodefile = tmp_path / "pbs_nodes"
    nodefile.write_text("node02\nnode02\nnode03\n", encoding="utf-8")
    monkeypatch.setenv("PBS_NODEFILE", str(nodefile))
    monkeypatch.setattr(dft_settings, "_affinity_count", lambda: None)

    assert dft_settings.resolve_dft_settings({"backend": "dft"}).nprocs == 3


def test_dft_shorthand_conflict_with_pyscf_object_setting_is_rejected() -> None:
    from pdb2reaction.core.dft_settings import resolve_dft_settings

    with pytest.raises(click.BadParameter, match="with_solvent.solvent"):
        resolve_dft_settings(
            {
                "backend": "dft",
                "dft": {
                    "solvent": "water",
                    "pyscf": {"with_solvent": {"solvent": "methanol"}},
                },
            }
        )


def test_iterative_dft_does_not_inherit_legacy_calc_solvent_scope() -> None:
    from pdb2reaction.core.dft_settings import resolve_dft_settings

    settings = resolve_dft_settings(
        {"backend": "dft", "solvent": "water", "solvent_model": "alpb"}
    )

    assert settings.solvent == "none"
    assert settings.solvent_model == "none"


def test_standalone_solvent_config_cannot_replace_pyscf_callable() -> None:
    from types import SimpleNamespace

    from pdb2reaction.backends.base import BackendError
    from pdb2reaction.workflows.dft import _apply_implicit_solvent

    solvent = SimpleNamespace(eps=None, reset=lambda: None)
    mf = SimpleNamespace(PCM=lambda: mf, with_solvent=solvent)

    with pytest.raises(BackendError, match="callable attribute"):
        _apply_implicit_solvent(
            mf,
            {
                "solvent": "water",
                "solvent_model": "pcm",
                "pyscf": {"with_solvent": {"reset": "disabled"}},
            },
        )


def test_standalone_density_fit_forwards_canonical_auxbasis() -> None:
    from pdb2reaction.workflows.dft import _configure_scf_object

    class Mutable:
        pass

    class MF:
        grids = Mutable()
        with_df = Mutable()

        def density_fit(self, **kwargs):
            self.density_kwargs = kwargs
            return self

    mf = _configure_scf_object(
        MF(),
        {
            "conv_tol": 1.0e-9,
            "max_cycle": 20,
            "grid_level": 1,
            "auxbasis": "def2-universal-jkfit",
            "pyscf": {"density_fit": {"enabled": True}},
        },
        "pbe",
    )

    assert mf.density_kwargs == {"auxbasis": "def2-universal-jkfit"}


@pytest.mark.parametrize(
    "pyscf_config, field",
    [
        ({"mol": {"basis": "sto-3g"}}, "basis"),
        ({"mf": {"xc": "pbe"}}, "xc"),
    ],
)
def test_dft_rejects_resolver_owned_pyscf_method_fields(
    pyscf_config, field
) -> None:
    from pdb2reaction.core.dft_settings import resolve_dft_settings

    with pytest.raises(click.BadParameter, match=field):
        resolve_dft_settings(
            {"backend": "dft", "dft": {"pyscf": pyscf_config}}
        )


@pytest.mark.parametrize("conv_tol", [float("nan"), float("inf"), float("-inf")])
def test_dft_rejects_nonfinite_convergence_tolerance(conv_tol) -> None:
    from pdb2reaction.core.dft_settings import resolve_dft_settings

    with pytest.raises(click.BadParameter, match="conv_tol"):
        resolve_dft_settings(
            {"backend": "dft", "dft": {"conv_tol": conv_tol}}
        )


@pytest.mark.parametrize(
    "key",
    [
        "typo_setting",
        "charge",
        "multiplicity",
        "save_scf_checkpoint",
        "checkpoint_path",
        "nprocs_source",
        "memory_source",
        "embedcharge",
        "embedcharge_cutoff",
    ],
)
def test_standalone_dft_mapping_fails_closed_on_unowned_key(key) -> None:
    from pdb2reaction.core.dft_settings import standalone_dft_settings_mapping

    with pytest.raises(click.BadParameter, match=f"dft.{key}"):
        standalone_dft_settings_mapping({"out_dir": "result", key: 1})


@pytest.mark.parametrize(
    "calc_cfg, message",
    [
        ({"charge": 0.5}, "calc.charge"),
        ({"spin": 0}, "calc.spin"),
        ({"dft": {"embedcharge_cutoff": float("nan")}}, "embedcharge_cutoff"),
        ({"dft": {"nprocs": 1.5}}, "DFT nprocs"),
    ],
)
def test_dft_rejects_invalid_canonical_charge_spin_and_cutoff(
    calc_cfg, message
) -> None:
    from pdb2reaction.core.dft_settings import resolve_dft_settings

    with pytest.raises(click.BadParameter, match=message):
        resolve_dft_settings({"backend": "dft", **calc_cfg})


def test_dmf_solvent_guard_is_backend_capability_aware() -> None:
    from pdb2reaction.workflows.path_opt import _validate_dmf_solvent_compatibility

    _validate_dmf_solvent_compatibility(
        {"backend": "dft", "solvent": "water", "solvent_model": "pcm"}
    )
    with pytest.raises(click.ClickException, match="gas-phase ASE PES"):
        _validate_dmf_solvent_compatibility(
            {"backend": "uma", "solvent": "water", "solvent_model": "alpb"}
        )


def test_scf_stepwise_grid_is_on_by_default_for_both_dft_paths(tmp_path) -> None:
    from pdb2reaction.core.dft_settings import resolve_dft_settings

    assert resolve_dft_settings({"backend": "dft"}).scf_stepwise_grid is True
    assert resolve_dft_settings(
        {"backend": "dft", "dft": {"scf_stepwise_grid": False}}
    ).scf_stepwise_grid is False

    xyz = tmp_path / "h2o.xyz"
    xyz.write_text("3\n\nO 0 0 0\nH 0 0.76 0.59\nH 0 -0.76 0.59\n")
    out_dir = tmp_path / "dft"
    result = CliRunner().invoke(
        cli,
        [
            "dft", "-i", str(xyz), "-q", "0", "-m", "1",
            "--func-basis", "lda/sto-3g", "--dft-engine", "cpu",
            "--scf-stepwise-grid", "true", "--out-json", "-o", str(out_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads((out_dir / "result.json").read_text())
    assert payload["converged"] is True


def test_scf_stepwise_grid_shortens_the_final_scf_of_the_dft_command() -> None:
    from pyscf import dft, gto

    from pdb2reaction.workflows.dft import _run_scf_kernel

    mol = gto.M(atom="O 0 0 0; H 0 0.76 0.59; H 0 -0.76 0.59", basis="sto-3g", verbose=0)

    def build_mf():
        mf = dft.RKS(mol)
        mf.xc = "lda"
        return mf

    normal, e_normal = _run_scf_kernel(build_mf, False)
    staged, e_staged = _run_scf_kernel(build_mf, True)

    assert e_staged == pytest.approx(e_normal, abs=1.0e-8)
    assert staged.cycles < normal.cycles


def test_calculator_dft_cli_flags_reach_the_dft_settings(tmp_path) -> None:
    xyz = Path(__file__).parent / "smoke" / "r.xyz"
    result = CliRunner().invoke(
        cli,
        [
            "sp", "-i", str(xyz), "-q", "-1", "-m", "1", "-b", "dft",
            "--dft-engine", "cpu", "--no-dft-low-memory", "--scf-stepwise-grid",
            "--show-config", "--dry-run", "-o", str(tmp_path / "sp"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "    engine: cpu\n" in result.output
    assert "    lowmem: false\n" in result.output
    assert "    scf_stepwise_grid: true\n" in result.output
