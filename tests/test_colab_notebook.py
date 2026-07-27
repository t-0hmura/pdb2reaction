"""Static contracts for the release-matched Colab GUI notebook."""

from __future__ import annotations

import ast
import csv
import datetime
import glob
import hashlib
import html
import json
import os
import shlex
import subprocess
import sys
import types
import zipfile
from pathlib import Path

import pytest


NOTEBOOK = Path(__file__).parents[1] / "examples" / "pdb2reaction_colab.ipynb"


def _execute_app(monkeypatch, tmp_path: Path) -> tuple[dict, list]:
    """Execute the complete app cell with real widgets and captured HTML."""
    rendered: list = []
    import IPython.display as ipd
    monkeypatch.setattr(ipd, "display", lambda *args, **kwargs: rendered.extend(args))
    monkeypatch.setattr(ipd, "clear_output", lambda *args, **kwargs: None)
    monkeypatch.setattr(ipd, "HTML", lambda value: value)
    monkeypatch.setattr(ipd, "Image", lambda *args, **kwargs: (args, kwargs))
    monkeypatch.chdir(tmp_path)
    namespace = {"TOOL": "pdb2reaction", "BACKEND": "mace", "REPO_DIR": "unused"}
    source = _notebook()["cells"][2]["source"]
    exec(compile(source, str(NOTEBOOK), "exec"), namespace)
    return namespace, rendered


def _notebook() -> dict:
    return json.loads(NOTEBOOK.read_text(encoding="utf-8"))


def _root_normalized_subcommand_argv(app: dict, subcommand: str, argv: list[str]) -> list[str]:
    """Apply the root CLI's compatibility normalization before direct parsing."""
    root = app["PRODUCT_CLI"]
    root_context = app["click"].Context(root)
    bool_values, bool_toggles, negative_aliases, single_flags = (
        root._resolve_bool_options(root_context, subcommand)
    )
    normalized, _legacy = root._normalize_bool_argv(
        [subcommand, *argv],
        {subcommand: bool_values},
        {subcommand: bool_toggles},
        {subcommand: negative_aliases},
        {subcommand: single_flags},
    )
    return normalized[1:]


def _output_contract() -> dict:
    """Execute only the pure output-tracking helpers from the app cell."""
    source = _notebook()["cells"][2]["source"]
    wanted = {
        "_normalized_scope_argv", "_effective_out_dir",
        "_click_subcommand_params", "_click_output_default",
        "_trj2fig_output_targets", "_flag_enabled", "_force_dry_run",
        "_grouped_option_values", "_parsed_path",
        "_input_needs_cif_companion", "_exact_output_scope",
        "_output_scope", "_effective_result_root", "_matches_output_scope",
        "_snapshot_files", "_snapshot_output_scope",
        "_output_scope_collision", "_structured_current_paths",
    }
    tree = ast.parse(source)
    module = ast.Module(
        body=[node for node in tree.body
              if isinstance(node, ast.FunctionDef) and node.name in wanted],
        type_ignores=[],
    )
    namespace = {
        "glob": glob, "json": json, "os": os, "shlex": shlex,
        "Path": Path, "CLI": "pdb2reaction",
        "S": {"out_dir": "result", "subcmd": "all"},
    }
    exec(compile(module, str(NOTEBOOK), "exec"), namespace)
    assert wanted <= namespace.keys()
    return namespace


def _viewer_contract() -> dict:
    """Execute pure selection/highlight and trajectory-label helpers."""
    source = _notebook()["cells"][2]["source"]
    wanted = {
        "_pick_text",
        "_resolve_click_meta",
        "_residue_id_selector",
        "_rich_residue_selector",
        "_center_cli_selectors",
        "_input_count_error",
        "_artifact_kind",
        "_csv_preview_html",
        "_text_preview_html",
        "_atom_signatures",
        "_resolve_atom_query",
        "_trajectory_result_metadata",
        "_trajectory_segment_ranges",
        "_trajectory_semantics",
        "_trajectory_frame_context",
        "_stationary",
    }
    tree = ast.parse(source)
    module = ast.Module(
        body=[node for node in tree.body
              if isinstance(node, ast.FunctionDef) and node.name in wanted],
        type_ignores=[],
    )
    namespace = {
        "os": os, "S": {}, "_TRAJ": {}, "Path": Path, "csv": csv, "html": html,
        "json": json, "_TEXT_PREVIEW_LIMIT": 512 * 1024,
        "SPEC": {}, "_ARTIFACT_KINDS": {
            ".png": "image", ".jpg": "image", ".jpeg": "image", ".svg": "SVG",
            ".html": "interactive HTML", ".csv": "CSV table", ".pdf": "PDF",
            ".gjf": "text", ".com": "text", ".inp": "text", ".prms": "text",
            ".pdb": "structure", ".ent": "structure", ".cif": "structure",
            ".mmcif": "structure",
        },
    }
    exec(compile(module, str(NOTEBOOK), "exec"), namespace)
    assert wanted <= namespace.keys()
    return namespace


def _advanced_widget_values(app: dict, param) -> list:
    """Return browser-widget values that exercise each option serializer."""
    click = app["click"]
    if param.name == "verbose":
        return [0, 3]
    if param.is_bool_flag or isinstance(param.type, click.types.BoolParamType):
        return [True, False]
    if isinstance(param.type, click.Choice):
        choices = list(param.type.choices)
        assert choices, f"{param.name} has an empty Click choice"
        return [choices[0]]
    nargs = max(1, int(getattr(param, "nargs", 1)))
    groups = 2 if param.multiple else 1
    return [" ".join(str(index + 1) for index in range(nargs * groups))]


def _widget_descendants(widget):
    yield widget
    for child in getattr(widget, "children", ()):
        yield from _widget_descendants(child)


def _widget_with_description(widget, description: str, *, enabled: bool = True):
    matches = [
        child for child in _widget_descendants(widget)
        if getattr(child, "description", "") == description
        and (not enabled or not getattr(child, "disabled", False))
    ]
    assert matches, f"missing widget: {description}"
    return matches[0]


def _assert_advanced_argv_parses(app: dict, subcommand: str, argv: list[str]) -> None:
    """Use Click's option parser without invoking a scientific workflow."""
    click = app["click"]
    command = app["_advanced_command"](subcommand)
    assert command is not None
    parser = command.make_parser(click.Context(command, info_name=subcommand))
    _values, leftovers, _order = parser.parse_args(list(argv))
    assert leftovers == [], (subcommand, argv, leftovers)


def test_colab_notebook_has_valid_code_cells_and_gpu_metadata() -> None:
    notebook = _notebook()
    introduction = notebook["cells"][0]["source"]

    assert notebook["nbformat"] == 4
    assert notebook["metadata"]["accelerator"] == "GPU"
    assert len(notebook["cells"]) == 3
    assert "[GitHub](https://github.com/t-0hmura/pdb2reaction)" in introduction
    assert "PDB/mmCIF or small-molecule XYZ/GJF structures" in introduction
    assert "**1 Input → 2 Setup → 3 Options → 4 Results**" in introduction
    assert (
        "[ChemRxiv](https://chemrxiv.org/doi/full/"
        "10.26434/chemrxiv.15003538/v1)"
    ) in introduction
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            ast.parse(cell["source"])


def test_colab_setup_is_pinned_to_matching_release_and_one_backend() -> None:
    setup = _notebook()["cells"][1]["source"]

    assert 'pdb2reaction_version = "v0.4.12"' in setup
    assert "first run ~8 min" in setup
    assert "first run ~5-10 min" not in setup
    # The release notebook installs the pinned wheel from PyPI, the same way a
    # normal user does, so the version guard compares what pip actually resolved
    # against the requested tag.
    assert "pip('pdb2reaction==' + pdb2reaction_version.lstrip('v'))" in setup
    # The DFT extra installs with a visible log; a quiet pip looked stalled.
    assert "pip('pdb2reaction[dft]==' + pdb2reaction_version.lstrip('v'))" in setup
    assert "pip_logged" not in setup
    assert "install_dft is ticked" in setup
    # Gated UMA sign-in is the last step, so no install phase waits on a prompt.
    assert "Hugging Face sign-in runs at the end of Installation" in setup
    assert setup.index("notebook_login()") > setup.index("_phase_done('version verified')")
    assert "git clone" not in setup
    assert "v0.4.4" not in setup
    assert "Restart the Colab runtime first" in setup
    assert "mace-torch>=0.3.8" in setup
    assert "HF_TOKEN" in setup
    assert "installed_version != pdb2reaction_version[1:]" in setup
    assert "install_dft = True" in setup
    assert "INSTALL_DFT = install_dft" in setup
    assert "Installation activates the selected backend." in setup
    assert "Only the **selected backend** is installed" not in setup
    assert "DFT_SETUP_READY = True" in setup
    assert "_dft_packages = {'pyscf': 'pyscf', 'gpu4pyscf': 'gpu4pyscf-cuda12x'}" in setup


def test_colab_setup_dft_branch_installs_extra_and_checks_gpu(monkeypatch, capsys) -> None:
    import importlib.metadata

    setup = _notebook()["cells"][1]["source"]   # install_dft now defaults to True
    calls: list[list[str]] = []

    def fake_run(argv, **_kwargs):
        calls.append([str(value) for value in argv])
        return types.SimpleNamespace(stdout="GPU 0", stderr="", returncode=0)

    popen_calls: list[list[str]] = []

    class _FakePopen:
        """The [dft] extra installs through pip_logged, which streams pip output."""

        def __init__(self, argv, **_kwargs):
            popen_calls.append([str(value) for value in argv])
            self.stdout = iter([
                "Collecting pyscf\n",
                "Downloading pyscf-2.13.0-cp311.whl (30.2 MB)\n",
                "a quiet line that must stay unechoed\n",
                "Installing collected packages: pyscf\n",
                "Successfully installed pyscf-2.13.0\n",
            ])

        def wait(self):
            return 0

    real_import_module = importlib.import_module
    fake_cupy = types.SimpleNamespace(
        cuda=types.SimpleNamespace(
            runtime=types.SimpleNamespace(getDeviceCount=lambda: 1),
        ),
    )

    def fake_import_module(name):
        if name == "cupy":
            return fake_cupy
        if name in {"pyscf", "basis_set_exchange", "gpu4pyscf.dft"}:
            return types.SimpleNamespace()
        return real_import_module(name)

    versions = {
        "pyscf": "2.11.0",
        "gpu4pyscf-cuda12x": "1.5.2",
        "pdb2reaction": "0.4.12",
    }
    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(subprocess, "Popen", _FakePopen)
    monkeypatch.setattr(importlib.util, "find_spec", lambda _name: object())
    monkeypatch.setattr(importlib, "import_module", fake_import_module)
    monkeypatch.setattr(importlib.metadata, "version", lambda name: versions[name])
    namespace: dict = {}
    exec(compile(setup, str(NOTEBOOK), "exec"), namespace)

    installs = [argv for argv in calls if "install" in argv]
    assert any("pdb2reaction==0.4.12" in argv for argv in installs)
    assert any("pdb2reaction[dft]==0.4.12" in argv for argv in installs)
    assert popen_calls == []          # no streamed pip log, only the announcement
    logged = capsys.readouterr().out
    assert "install_dft is ticked" in logged
    assert "DFT support installed: PySCF 2.11.0" in logged
    assert "Collecting" not in logged and "Downloading" not in logged
    assert namespace["DFT_SETUP_READY"] is True
    assert namespace["INSTALL_DFT"] is True
    assert namespace["BACKEND"] == "mace"
    assert "importlib.util.find_spec(module)" in setup
    assert "_dft_imports = ('pyscf', 'basis_set_exchange', 'gpu4pyscf.dft')" in setup
    assert "for _module in _dft_imports: importlib.import_module(_module)" in setup
    assert "_cupy.cuda.runtime.getDeviceCount()" in setup
    assert "DFT packages installed but failed their import/GPU check" in setup
    assert "DFT support installed: PySCF %s · GPU4PySCF %s" in setup
    assert "anywidget==0.11.0" in setup
    assert "py3Dmol" not in setup
    assert "pip('ipywidgets','anywidget==0.11.0','matplotlib')" in setup
    assert "[%%d/5] %%s".replace("%%", "%") in setup
    assert "time.monotonic()" in setup
    assert ".pysisyphusrc" not in setup
    assert "nglview" not in setup


@pytest.mark.parametrize(
    ("backend", "use_token", "expected_login"),
    [
        ("orb", False, None),
        ("uma", True, "token"),
        ("uma", False, "notebook"),
    ],
)
def test_colab_setup_operates_orb_and_uma_branches(
    monkeypatch, backend, use_token, expected_login,
) -> None:
    import importlib.metadata

    setup = _notebook()["cells"][1]["source"].replace(
        'backend = "mace"', f'backend = "{backend}"', 1,
    ).replace("install_dft = True", "install_dft = False", 1)
    calls: list[list[str]] = []
    logins: list[tuple] = []

    def fake_run(argv, **_kwargs):
        calls.append([str(value) for value in argv])
        return types.SimpleNamespace(stdout="GPU 0", stderr="", returncode=0)

    fake_hf = types.ModuleType("huggingface_hub")
    fake_hf.login = lambda **kwargs: logins.append(("token", kwargs))
    fake_hf.notebook_login = lambda: logins.append(("notebook", {}))
    monkeypatch.setitem(sys.modules, "huggingface_hub", fake_hf)
    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(
        importlib.metadata,
        "version",
        lambda name: "0.4.12" if name == "pdb2reaction" else "test",
    )
    monkeypatch.delenv("HF_TOKEN", raising=False)
    if use_token:
        monkeypatch.setenv("HF_TOKEN", "test-token")

    namespace: dict = {}
    exec(compile(setup, str(NOTEBOOK), "exec"), namespace)

    joined = [" ".join(argv) for argv in calls]
    assert namespace["BACKEND"] == backend
    assert not any("mace-torch" in command for command in joined)
    assert not any("uninstall -y -q fairchem-core" in command for command in joined)
    if backend == "orb":
        assert any("orb-models" in command for command in joined)
        assert logins == []
    else:
        assert not any("orb-models" in command for command in joined)
        assert [entry[0] for entry in logins] == [expected_login]
        if use_token:
            assert logins[0][1]["token"] == "test-token"


def test_colab_gui_tracks_current_structure_and_execution_contracts() -> None:
    app = _notebook()["cells"][2]["source"]

    assert isinstance(app, str)
    assert ".pdb,.ent,.cif,.mmcif,.xyz,.gjf" in app
    assert "prepare_input_structure" in app
    assert "load_pdb_atom_metadata" in app
    assert "CHAIN:RESNAME:RESSEQ" in app
    assert "String(item.sourceIndex)" in app
    assert "callback_ns = 'pdb2reaction_gui'" in app
    assert "_co.register_callback('pdb2reaction_gui.on_click', on_click)" in app
    assert "return bridge.api.invokeFunction(cfg.callback+'.'+suffix,args,kwargs);" in app
    assert "demo.on_click" not in app
    assert "return shlex.split(line)" in app
    assert "dry_argv = _force_dry_run(list(a))" in app
    assert "if not _validate_command(a): return" in app
    assert "real_run = not _flag_enabled(effective, '--dry-run', '--no-dry-run')" in app
    assert "def _utility_autofill_complete(argv):" in app
    assert (
        app.index("input_records = (")
        < app.index("rc, transcript = _stream(a)")
        < app.index("'inputs': input_records")
    )
    assert "b_run = W.Button(description='Run'" in app
    assert "RUN exact command" not in app
    assert "def _preflight_output_scope(argv):" in app
    assert "Run not started: refusing to overwrite the existing output:" in app
    assert '"""Return result(1), result(2), ... without touching prior output."""' in app
    assert "uma-s-1p2" in app
    assert "MODELS = {'mace': ['MACE-OMOL-0']" in app
    assert "options=['auto', 'fp32', 'fp64']" in app
    assert "Appended to the command" not in app
    assert "Every remaining option" not in app
    assert "_MOLSTAR_VERSION = '5.6.1'" in app
    assert "molstar@%s/build/viewer/molstar.js" in app
    assert "layoutShowSequence:cfg.showSequence" in app
    assert "layoutShowControls:true" in app
    assert "collapseRightPanel:true" in app
    assert "rxworkspace" in app
    assert "rxviewer" in app
    assert "rxinspector" in app
    assert ":focus-visible" in app
    assert ".rxapp .widget-button .fa, .rxapp .widget-upload .fa { display:none !important; }" in app
    assert "overflow-wrap:anywhere" in app
    assert "max_width='100%'" in app
    assert "flex:0 0 auto; min-width:0" in app
    assert "layout=W.Layout(width='260px', max_width='100%')" in app
    assert "Mol* owns representation, colour, camera, measurement, and screenshot controls." in app
    assert "view_controls = W.HBox([" in app
    assert "cb_water," in app
    assert "last_pick_row = W.HBox([last_pick_html]," in app
    assert "btn_clear_pick" not in app
    assert "sel_lang" not in app
    assert "No structure loaded" in app
    # Colab renders ipywidgets' Tab as an empty block, so the tab strip is Buttons
    # over four persistently mounted panes. Mount persistence protects native
    # upload queues and the live WebGL canvas while navigating.
    assert "_TAB_PAGES = [('1 Input', input_box), ('2 Setup', viewer_box)," in app
    assert "('3 Options', options_box), ('4 Results', results_box)]" in app
    assert "_tab_body = W.VBox([page for _label, page in _TAB_PAGES])" in app
    assert "_tab_body.children = [_TAB_PAGES[i][1]]" not in app
    assert "_pane.layout.display = '' if _j == i else 'none'" in app
    assert "Options (optional)" not in app
    assert "def _tab_go(i):" in app
    assert "W.Tab(" not in app
    # Accordion has the same Colab rendering defect, so every collapsible uses
    # the same Button + VBox pattern.
    assert "def _collapsible(title, child, on_open=None):" in app
    assert "W.Accordion(" not in app
    assert "adv_acc = _collapsible(" in app
    assert "def _advanced_options(sub):" in app
    assert "adv_acc = _collapsible('Advanced flags', adv_box)" in app
    assert "every CLI option accounted for" not in app
    assert "freeze_acc = _collapsible(" in app
    assert "ngl_acc" not in app
    assert "nglview" not in app
    # "all workflow" mode and the depth switches only apply to `all`.
    assert "all_mode_box = all_mode" in app
    assert "('Scan 1', 'scan')" in app
    assert "style={'button_width': '74px', 'description_width': '0px'}" in app
    assert "depth_label = W.HTML('<b>Optional stages</b>')" in app
    assert "depth_box = W.VBox([depth_label, W.HBox([" in app
    assert "_flag_row(w_ts" in app and "_flag_row(w_th" in app
    assert "for _name in ('all_mode_box', 'depth_box'):" in app
    # -r/--radius is an extraction option; small-molecule `all` skips extraction.
    assert "radius_applies = (sub == 'extract' or" in app
    assert "(sub == 'all' and S.get('mode') in ('pdb', 'mmcif')))" in app
    assert "if radius_applies and r and r > 0: cmd += ['-r', str(r)]" in app
    assert "adv_radius.disabled = not radius_applies" in app
    assert "_set_flag_visible(adv_radius, radius_applies)" in app
    assert "_set_flag_visible(adv_dftfb, _dftfb_applicable)" in app
    assert "req='1 file + scan-lists" in app
    assert "No results yet." in app
    assert "frame_slider.disabled = last == 0" in app
    assert "trajectory_box.layout.display = 'none'" in app
    # `--tr-projection legacy-active` is deprecated (it warns and must not be used
    # for pass/HOSP transition-state certification): never offer it in the GUI.
    assert "legacy-active" not in app
    assert "--tr-projection" not in app
    # SPEC is the single source of truth for per-subcommand branching: input
    # requirement, Select panels, Options flag visibility and Results lookup all
    # read it, so a new subcommand is one table entry rather than N conditionals.
    assert "SPEC = {" in app
    assert "SUBREQ = {k: v['req'] for k, v in SPEC.items()}" in app
    assert "FLAG_SUBS = {" in app
    assert "'mep_mode': FLAG_SUBS['adv_mep']," in app
    assert "'threshold': FLAG_SUBS['adv_thresh']," in app
    assert "'adv_thresh':  {'all', 'opt', 'tsopt', 'scan', 'scan2d', 'scan3d', 'path-opt', 'path-search'}," in app
    assert "if 'freeze' in SPEC.get(sub, {}).get('panels', ()) and S['freeze_atoms']:" in app
    assert "_panels = SPEC.get(sub, {}).get('panels', ())" in app
    # Surfaced key flags + their emission.
    assert "key_opts_box = _collapsible('Key options', key_opts_content)" in app
    assert "cmd += ['--flatten']" in app
    assert "cmd += ['--refine-path']" in app
    assert "cmd += ['--max-cycles', str(int(mc))]" in app
    # Help is a browser-native disclosure. Opening it needs no Python round trip,
    # which keeps it responsive in Colab while the panel remains in normal flow.
    assert "def _info_markup(tip, revision=0, rich=False):" in app
    assert "def _info_control(tip, target=None, rich=False):" in app
    assert "def _set_info_text(control, tip):" in app
    assert "def _close_info_target(target):" in app
    assert "def _hdr(content, tip):" in app
    assert "def _flag_row(widget, tip, info_target=None, rich=False):" in app
    assert "'<summary aria-label=\"More information: %s\" title=\"%s\">&#9432;</summary>'" in app
    assert "'<details class=\"rxinfo-details\" data-revision=\"%d\">'" in app
    assert "'<div class=\"rxhelp-panel\" role=\"note\"><small>%s</small>" in app
    assert "_INFO_CONTROLS = weakref.WeakSet()" in app
    assert "description='Show information', icon='info-circle'" not in app
    assert "_rx_info_button" not in app
    assert "rxinfo-popover" not in app
    assert "📥 needs" not in app
    assert "w_show_run_log = W.Checkbox(value=True, description='Show run log', indent=False)" in app
    assert "logbox.layout.display = '' if change['new'] else 'none'" in app
    assert "cmdline_box = W.VBox([command_editor, command_footer, logbox])" in app
    assert "command_editor = W.VBox([" in app
    assert "W.HTML('<b>Command line</b>')" in app
    assert "command_editor = _collapsible('Command line'" not in app
    assert "viewer_toolbar = W.HBox(" in app
    assert "[view_input, pick_action, last_pick_info, view_controls, viewer_more]" in app
    assert "viewer_box = W.VBox([workflow_box, viewer_toolbar, view_input_note, selection_box])" in app
    assert "'<div class=\"rxmolstar-embed\">%s</div>' % _document_iframe(" in app
    assert "workflow_contract_row = W.HBox([subreq, outputs_html]" in app
    assert "workflow_controls.add_class('rxworkflow-controls')" in app
    assert "grid-template-columns:minmax(230px,270px) minmax(140px,1fr) 246px;" in app
    assert "height:auto; min-height:0; overflow-x:clip; overflow-y:visible;" in app
    assert "@media (min-width: 821px) and (max-height: 900px)" in app
    assert ".rxapp-main { flex:0 0 auto !important; min-height:0; overflow:visible; }" in app
    assert ".rxpages { flex:0 0 auto !important; min-height:0; overflow:visible; }" in app
    assert "height:auto; max-height:none; min-height:0; overflow:visible;" in app
    assert ".rxviewer-page { min-height:800px !important; }" in app
    assert "width:100% !important; min-width:0; max-width:none;" in app
    assert "max-width:clamp(600px,calc(250dvh - 1320px),1000px);" in app
    assert ".rxpath-panel svg, .rxpath-panel img, .rxpath-panel canvas {" in app
    assert "traj_out = W.HTML(layout={'width': '100%', 'min_width': '0'})" in app
    assert "plot_out = W.HTML(layout={'width': '100%', 'min_width': '0'})" in app
    assert "'flex': '1 1 440px'" not in app
    assert ".rxcommand-dock {" in app
    assert "rootbox = W.VBox([header, app, cmdline_box])" in app
    assert "rootbox = W.VBox([header, app, W.HTML('<hr" not in app
    assert "role=\"tooltip\"" not in app
    # Standalone opt/tsopt/path-opt users need a cluster model first.
    assert "b_extract = W.Button(description='Prepare cluster & use it'" in app
    assert "prep_radius = W.FloatText(" in app
    assert "keep_subcmd=True" in app
    # The wheel ships no examples, so Load example resolves them from the git
    # tag matching the installed release (a source checkout is used when present).
    assert "def _example_file(relpath):" in app
    assert "raw.githubusercontent.com/t-0hmura/pdb2reaction" in app
    assert "Aromatic Claisen rearrangement - small molecule" in app
    assert "aromatic_claisen/reactant.xyz" in app
    assert "aromatic_claisen/product.xyz" in app
    assert "HCN -> HNC" not in app
    assert "Run Setup first" not in app
    assert "cmd = [CLI, 'extract', '-i', *S['inputs'], '-o', *outs, '-c', ','.join(cen)]" in app
    assert "cen = _center_cli_selectors()" in app
    assert "command_footer = W.VBox([command_actions, w_show_run_log]" in app
    assert "cmdline_box.add_class('rxcommand-dock')" in app
    assert "for _label, _page in _TAB_PAGES:" in app
    assert "_page.add_class('rxpage')" in app
    assert "_tab_body.add_class('rxpages')" in app
    assert "app.add_class('rxapp-main')" in app
    assert "'MLIP_Gibbs': 'gibbs_mlip', 'MLIP': 'mlip'" in app
    assert "ps.get('uma')" not in app
    assert "ps.get('gibbs_uma')" not in app
    assert "def _effective_out_dir(argv=None):" in app
    assert "flags = ('-o', '--out', '--out-dir', '--output')" in app
    assert "def _trj2fig_output_targets(argv):" in app
    assert "command.make_parser(sub_ctx).parse_args" in app
    assert "def _snapshot_output_scope(scope):" in app
    assert "def _output_scope_collision(scope):" in app
    assert "scope = _output_scope(argv)" in app
    assert "target = scope['target']; out = scope['root']" in app
    assert "_results(out)" in app
    assert "'adv_mep':     {'all', 'path-opt', 'path-search'}," in app
    assert "'adv_dmf':     {'all', 'path-opt', 'path-search'}," in app
    assert "cmd += ['--dmf-backend', dmf_backend]" in app
    assert "adv_thresh.disabled = sub not in TOOL_CAPABILITIES['threshold']" in app
    assert "elif sub == 'dft':" in app
    assert "cmd += ['--func-basis', fb]" in app
    assert "adv_mep.disabled = sub not in TOOL_CAPABILITIES['mep_mode']" in app
    assert "_set_flag_visible(adv_dmf, sub in FLAG_SUBS['adv_dmf'] and adv_mep.value == 'dmf')" in app
    assert "d['all_mode'] = _wv('all_mode', 'mep')" in app
    assert "all_mode.value = d['all_mode']" in app
    assert "def _validate_and_normalize_session(payload):" in app
    assert "_SESSION_APPLY = {'active': False}" in app
    assert "bytes(c).decode('utf-8')" in app
    assert "def _ingest_saved_files(" in app
    assert "input_file_rows" in app
    assert "class _DropUpload(anywidget.AnyWidget):" in app
    assert "upl = _DropUpload(accept=_acc, formats=_drop_formats)" in app
    assert "upl = W.FileUpload(accept=_acc, multiple=True, description='Upload files'" in app
    assert "_drop = W.VBox(_drop_children," in app
    assert "align_items='center', justify_content='center'" in app
    assert "el.addEventListener('drop', onDrop)" in app
    assert "file.arrayBuffer()" in app
    assert "model.send({" in app
    assert "_DROP_STATE = {'generation': 0, 'seen': set()}" in app
    assert "def _claim_drop_batch(batch, generation):" in app
    assert "def _bump_drop_generation():" in app
    assert "queue.push({id: opaqueId(), files: selected, generation})" in app
    assert "if (busy || !queue.length || signal.aborted) return;" in app
    assert "message.batch !== active" in app
    assert "if clear:\n        _bump_drop_generation()" in app
    assert "pdb2reaction_gui.on_drop" not in app
    # Colab draws no AnyWidget and drops a FileUpload's binary buffers, so there
    # every upload ships its bytes through invokeFunction from its own zone.
    assert "_UPLOAD_MODE = ('colab' if IN_COLAB else" in app
    assert "_drop_children = ([upl] if _UPLOAD_MODE == 'anywidget'" in app
    assert "_cwm.register_callback('pdb2reaction_gui.upload_files', _on_colab_upload)" in app
    # callback_ns is local to _molstar_document; the module-level block needs a literal.
    assert "callback_ns" not in app.split("display(rootbox)")[1]
    assert "bridge.invokeFunction(CONFIG.callback, [spec.role, payload], {})" in app
    assert "reader.readAsDataURL(file)" in app
    assert "def _decode_colab_batch(files):" in app
    assert "base64.b64decode(entry['b64'], validate=True)" in app
    for _role in ("'input'", "'model'", "'session'"):
        assert "role=%s" % _role in app or "role == %s" % _role in app
    assert "_accept_upload_pairs(pairs, 'browser upload')" in app
    assert "_accept_upload_pairs(pairs, 'browser upload')" in app
    assert "def _delete_owned_uploads(paths):" in app
    assert "_HAS_DROP_WIDGET" in app
    assert "if _UPLOAD_MODE == 'anywidget': upl.on_msg(_on_drop_upload)" in app
    assert "description='Remove file'" in app
    assert "description='Move earlier'" in app
    assert "description='Move later'" in app
    assert "tooltip='Move earlier'" in app
    assert "tooltip='Move later'" in app
    assert "def _advanced_coverage(" in app
    assert "adv_extra" not in app


def test_colab_viewer_persists_exact_atom_and_residue_context() -> None:
    app = _notebook()["cells"][2]["source"]

    for marker in (
        "'_last_pick': None", "'_pick_history': []", "def _remember_pick",
        "_VIEWER_GENERATION", "_MOLSTAR_VERSION = '5.6.1'",
        "layoutShowSequence:cfg.showSequence", "layoutShowControls:true",
        "collapseRightPanel:true", "layoutShowLog:false",
        "layoutShowLeftPanel:false", "layoutShowRemoteState:false",
        "await molstar.Viewer.create", "await setCartoonAlpha()",
        "alpha:0.4", "structure-component-static-water",
        "await configureWater()", "tryCreateComponentStatic(",
        "'mmcif' if fmt in ('cif', 'mmcif')",
        "ignoreStartupEmpty", "loci.kind==='empty-loci'",
        "clickQueue=clickQueue.then(()=>handleClick(event))",
        "await invoke('clear_highlights',[cfg.generation])",
        "cfg.generation,exact",
        "exact atom", "set current pick", "view_input", "_view_mapping_ok",
        "last_pick_info", "Generated file preview", "Download current run (.zip)",
        "results_box.add_class('rxresults')", "overflow-x:auto",
        "colab_run.log", "energy unavailable", "Command was cancelled",
        "Command failed", "_frame_link = W.jslink", "linked structure + energy",
        "host.on('plotly_click'", "Plotly.restyle", "artifact_fold._rx_set_open",
        "message.type!=='rx-set-frame'", "update.to(model).update",
        "channel='trajectory', generation=generation",
        "''.join(_TRAJ['frames'])",
        "display(HTML(_molstar_iframe(source, fmt, show_sequence=(fmt != 'xyz'))))",
    ):
        assert marker in app
    assert "py3Dmol" not in app
    assert "v.addSphere(" not in app
    assert "v.addBox(" not in app
    assert "_pick_box_edges" not in app
    assert "'greenCarbon'" not in app
    assert "viewportFocusBehavior" not in app
    assert "tryCreateComponentFromExpression" not in app
    assert "lociSelects.select" not in app
    assert "canvas.setProps" not in app
    assert "managers.interactivity.setProps" not in app
    assert "_viewer_seed_picks" not in app
    assert "skipInitialClick" not in app
    assert "style_pick" not in app
    assert "('Measure (dist/angle/dihedral)'," not in app
    assert "measure_panel" not in app
    assert "freeze_acc = _collapsible('Freezing', freeze_panel)" in app
    assert "_PICK_HINT.get(pick_action.value, 'Choose a click action.')" in app
    assert "_PICK_HINT[pick_action.value]" not in app
    assert "pick_action.value = 'scanB'" in app
    assert "pick_action.value = 'freezeB'" in app
    assert "return bridge.api.invokeFunction(cfg.callback+'.'+suffix,args,kwargs);" in app
    assert "String(item.sourceIndex)" in app
    assert "cfg.generation,exact" in app
    assert "if exact in (False, 0, 'false', 'False')" in app
    assert "Mol* focused this residue" in app
    click_completion = app[
        app.index("def on_click("):
        app.index("    _co.register_callback('pdb2reaction_gui.on_click', on_click)")
    ]
    assert "if not live_marked:" in click_completion
    assert "selected in (False, 0, 'false', 'False')" not in click_completion
    assert "current_generation != _VIEWER_GENERATION['value']" in click_completion
    assert app.count("pdb2reaction_gui.clear_highlights") == 1
    assert "def _resolve_click_meta(" in app
    assert "'viewer_index': viewer_index" in app
    assert "def _load_small_view_structure(path):" in app
    assert "_view_format='xyz' if is_small else 'pdb'" in app
    assert "def _invalidate_last_run(" in app
    assert "def _artifact_kind(" in app

    contract = _viewer_contract()
    metadata = [
        {"index": 0, "chain": "A", "resname": "LIG", "resseq": 10,
         "icode": "", "name": "N", "xyz": (0, 0, 0)},
        {"index": 1, "chain": "A", "resname": "LIG", "resseq": 10,
         "icode": "", "name": "C1", "xyz": (1, 0, 0)},
        {"index": 2, "chain": "B", "resname": "LIG", "resseq": 10,
         "icode": "", "name": "N", "xyz": (5, 0, 0)},
        {"index": 3, "chain": "A", "resname": "LIG", "resseq": 10,
         "icode": "A", "name": "N", "xyz": (8, 0, 0)},
    ]
    serial_metadata = [dict(row, serial=100 + i) for i, row in enumerate(metadata)]
    contract["S"].update(_atom_meta=serial_metadata, _view_format="pdb")
    # Stable PDB serial/signature wins over a drifted browser-side source index.
    assert contract["_resolve_click_meta"](
        0, "101", "LIG", "10", "A", "C1", "",
    )["index"] == 1
    contract["S"]["_atom_meta"] = metadata
    assert contract["_resolve_atom_query"]("A:LIG:10:C1") == 1
    assert contract["_resolve_atom_query"]("2") == 1
    contract["S"].update(
        _view_format="xyz",
        _atom_meta=[{"index": 0, "serial": 1, "resname": "MOL",
                     "resseq": 1, "name": "C1", "xyz": (0, 0, 0)}],
    )
    resolved_xyz = contract["_resolve_click_meta"](
        "0", "0", "", "undefined", "", "C", "",
    )
    assert resolved_xyz["name"] == "C1"
    contract["S"].update(center=["LIG"], center_ids=["B:LIG:10"],
                         _primary_atom_meta=metadata)
    assert contract["_center_cli_selectors"]() == [
        "A:LIG:10", "B:LIG:10", "A:LIG:10A",
    ]
    contract["S"].update(center=[], center_ids=["A:LIG:10"])
    with pytest.raises(ValueError, match="insertion-code sibling"):
        contract["_center_cli_selectors"]()
    contract["SPEC"].update({
        "bond-summary": {"n_in": (2, None)},
        "trj2fig": {"n_in": (1, 1)},
    })
    assert contract["_input_count_error"]("bond-summary", 1)
    assert not contract["_input_count_error"]("bond-summary", 2)
    assert contract["_input_count_error"]("trj2fig", 2)
    assert contract["_artifact_kind"]("profile.html") == "interactive HTML"
    assert contract["_artifact_kind"]("plot.jpeg") == "image"
    assert contract["_artifact_kind"]("final.pdb") == "structure"
    assert contract["_artifact_kind"]("final.cif") == "structure"
    assert contract["_artifact_kind"]("final.xyz") == "structure"
    assert contract["_artifact_kind"]("job.gjf") == "text"
    assert contract["_artifact_kind"]("job.com") == "text"
    assert contract["_artifact_kind"]("job.inp") == "text"
    assert contract["_artifact_kind"]("optimization_trj.xyz") is None
    opt = contract["_trajectory_semantics"]("opt", "optimization_trj.xyz")
    path_semantics = contract["_trajectory_semantics"](
        "path-opt", "mep_trj.xyz",
    )
    vibration = contract["_trajectory_semantics"](
        "tsopt", "vib/imag_120i_trj.xyz",
    )
    assert vibration["title"] == "Vibrational-mode animation"
    assert vibration["x"] == "phase frame" and not vibration["extrema"]
    assert contract["_stationary"]([0.0, 2.0, 0.0], opt) == [
        (0, "initial"), (2, "optimized"),
    ]
    # A trajectory plot labels R and P only (user decision, 2026-07-27). An
    # energy-only extremum is neither a certified transition state nor a
    # certified intermediate; the energy diagram is where the profile is read.
    assert contract["_stationary"]([0.0, 2.0, 0.0], path_semantics) == [
        (0, "R"), (2, "P"),
    ]
    assert contract["_stationary"]([2.0, 0.0, 2.0], path_semantics) == [
        (0, "R"), (2, "P"),
    ]


def test_colab_results_bind_irc_truth_and_skip_bridge_extrema(tmp_path: Path) -> None:
    contract = _viewer_contract()
    irc_dir = tmp_path / "irc"
    irc_dir.mkdir()
    irc_path = irc_dir / "finished_irc_trj.xyz"
    irc_path.write_text("", encoding="utf-8")
    (irc_dir / "result.json").write_text(
        json.dumps({
            "n_frames_forward": 2,
            "n_frames_backward": 2,
            "n_frames_total": 5,
            "forward_converged": True,
            "backward_converged": False,
            "scientific_status": "partial",
            "scientific_status_reasons": ["backward IRC did not converge"],
        }),
        encoding="utf-8",
    )
    irc = contract["_trajectory_semantics"]("irc", str(irc_path), n_frames=5)
    assert irc["ts_index"] == 2
    assert "partial" in irc["trajectory_status"]
    assert "forward ✓" in irc["trajectory_status"]
    assert "backward ✗" in irc["trajectory_status"]
    mismatch = contract["_trajectory_semantics"](
        "irc", str(irc_path), n_frames=4,
    )
    assert "ts_index" not in mismatch
    assert mismatch["metadata_warning"] == "IRC frame metadata mismatch"

    mep_path = tmp_path / "mep_trj.xyz"
    mep_path.write_text("", encoding="utf-8")
    (tmp_path / "summary.json").write_text(
        json.dumps({
            "segments": [
                {"index": 1, "kind": "seg", "frame_ranges": [[0, 2]]},
                {"index": 2, "kind": "bridge", "frame_ranges": [[2, 5]]},
                {"index": 3, "kind": "seg", "frame_ranges": [[5, 7]]},
            ]
        }),
        encoding="utf-8",
    )
    path = contract["_trajectory_semantics"](
        "path-search", str(mep_path), n_frames=7,
    )
    assert path["bridge_ranges"] == [(2, 5)]
    assert contract["_trajectory_frame_context"](1, path)["label"] == "seg 01"
    assert contract["_trajectory_frame_context"](3, path)["label"] == "bridge"
    assert contract["_trajectory_frame_context"](6, path)["label"] == "seg 03"
    irc_context = contract["_trajectory_frame_context"](2, irc)
    assert "backward IRC did not converge" in " ".join(irc_context["details"])
    points = contract["_stationary"](
        [0.0, 0.5, 0.0, 3.0, 0.0, 2.0, 0.0], path,
    )
    assert points == [(0, "R"), (6, "P")]


def test_colab_drop_protocol_appends_and_rejects_incomplete_batches(
    tmp_path: Path, monkeypatch,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)

    class FakeDropWidget:
        def __init__(self, generation):
            self.generation = generation
            self.messages = []

        def send(self, message):
            self.messages.append(dict(message))

    generation = app["_DROP_STATE"]["generation"]
    widget = FakeDropWidget(generation)
    pdb = (
        b"HETATM    1  C1  LIG A  10       0.000   0.000   0.000"
        b"  1.00  0.00           C\nEND\n"
    )

    def upload(batch, name, payload=pdb, *, size=None, buffers=None):
        content = {
            "event": "upload",
            "batch": batch,
            "generation": widget.generation,
            "files": [{"name": name, "size": len(payload) if size is None else size}],
        }
        app["_on_drop_upload"](
            widget, content, [payload] if buffers is None else buffers,
        )

    upload("batch-1", "structure.pdb")
    upload("batch-2", "structure.pdb")
    assert [Path(path).name for path in app["S"]["inputs"]] == [
        "structure.pdb",
        "structure_2.pdb",
    ]
    assert [message["ok"] for message in widget.messages[-2:]] == [True, True]

    before = list(app["S"]["inputs"])
    upload("bad-size", "broken.pdb", size=len(pdb) + 1)
    upload("bad-count", "missing.pdb", buffers=[])
    assert app["S"]["inputs"] == before
    assert [message["ok"] for message in widget.messages[-2:]] == [False, False]
    assert "existing files were kept" in app["input_msg"].value

    old_generation = widget.generation
    app["_bump_drop_generation"]()
    widget.generation = app["_DROP_STATE"]["generation"]
    stale = {
        "event": "upload",
        "batch": "stale",
        "generation": old_generation,
        "files": [{"name": "late.pdb", "size": len(pdb)}],
    }
    app["_on_drop_upload"](widget, stale, [pdb])
    assert widget.messages[-1]["ok"] is False
    assert app["S"]["inputs"] == before


def test_colab_app_executes_atomic_view_and_result_transitions(
    tmp_path: Path, monkeypatch,
) -> None:
    app, calls = _execute_app(monkeypatch, tmp_path)
    drop_generation = app["_DROP_STATE"]["generation"]
    assert app["_claim_drop_batch"]("batch-a", drop_generation) == (True, "")
    assert app["_claim_drop_batch"]("batch-a", drop_generation) == (
        False, "duplicate batch",
    )
    app["b_clear_inputs"].click()
    assert app["_DROP_STATE"]["generation"] == drop_generation + 1
    assert app["_claim_drop_batch"]("late-batch", drop_generation) == (
        False, "stale generation",
    )
    primary = tmp_path / "primary.pdb"
    secondary = tmp_path / "secondary.pdb"
    primary_text = (
        "HETATM    1  C1  LIG A  10       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O1  LIG A  10       1.200   0.000   0.000  1.00  0.00           O\n"
        "HETATM    3  C1  COF B  20       4.000   0.000   0.000  1.00  0.00           C\nEND\n"
    )
    secondary_text = (
        "HETATM    1  C1  ALT A  10       0.100   0.000   0.000  1.00  0.00           C\nEND\n"
    )
    primary.write_text(primary_text, encoding="utf-8")
    secondary.write_text(secondary_text, encoding="utf-8")
    calc_file = tmp_path / "calculator.py"
    ref_mode = tmp_path / "reference-mode.npy"
    calc_file.write_text("calculator = 1\n", encoding="utf-8")
    ref_mode.write_bytes(b"reference mode")
    validation_argv = [
        "pdb2reaction", "tsopt", "-i", str(primary),
        "--calc-file", str(calc_file), "--ref-mode", str(ref_mode),
    ]
    app["cmd_box"].value = shlex.join(validation_argv)
    assert set(app["_command_input_files"](validation_argv)) == {
        str(primary), str(calc_file), str(ref_mode),
    }
    first_fingerprint = app["_validation_fingerprint"](validation_argv)
    calc_file.write_text("calculator = 2\n", encoding="utf-8")
    assert app["_validation_fingerprint"](validation_argv) != first_fingerprint
    extra_structure = tmp_path / "extra.pdb"
    extra_structure.write_text(secondary_text, encoding="utf-8")
    positional_argv = [
        "pdb2reaction", "bond-summary", str(primary), str(extra_structure),
    ]
    assert set(app["_command_input_files"](positional_argv)) == {
        str(primary), str(extra_structure),
    }
    scan_spec = tmp_path / "scan.yaml"
    scan_spec.write_text("scan: first\n", encoding="utf-8")
    scan_argv = [
        "pdb2reaction", "scan", "-i", str(primary), "-s", str(scan_spec),
    ]
    app["cmd_box"].value = shlex.join(scan_argv)
    assert set(app["_command_input_files"](scan_argv)) == {str(primary), str(scan_spec)}
    scan_fingerprint = app["_validation_fingerprint"](scan_argv)
    scan_spec.write_text("scan: second\n", encoding="utf-8")
    assert app["_validation_fingerprint"](scan_argv) != scan_fingerprint
    metadata = {
        str(primary): [
            {"chain": "A", "resname": "LIG", "resseq": 10, "icode": "", "name": "C1"},
            {"chain": "A", "resname": "LIG", "resseq": 10, "icode": "", "name": "O1"},
            {"chain": "B", "resname": "COF", "resseq": 20, "icode": "", "name": "C1"},
        ],
        str(secondary): [
            {"chain": "A", "resname": "ALT", "resseq": 10, "icode": "", "name": "C1"},
        ],
    }

    def load_view(path):
        path = str(path)
        return Path(path).read_text(encoding="utf-8"), [dict(row) for row in metadata[path]], path

    app["_load_view_structure"] = load_view
    app["load_pdb"]([str(primary), str(secondary)], center=["COF"], lcharge={"COF": 1})
    primary_widget = app["center_widget"]
    app["pick_action"].value = "center"
    app["on_click"]("0", "LIG", "10", "A", "C1")
    assert app["_center_cli_selectors"]() == ["A:LIG:10", "B:COF:20"]

    from Bio.PDB import PDBParser
    from pdb2reaction.workflows.extract import resolve_substrate_residues
    structure = PDBParser(QUIET=True).get_structure("primary", primary)
    resolved = resolve_substrate_residues(structure, "B:COF:20,A:LIG:10")
    assert {(residue.get_parent().id, residue.id[1], residue.resname) for residue in resolved} == {
        ("A", 10, "LIG"), ("B", 20, "COF"),
    }

    saved = (list(app["S"]["center"]), list(app["S"]["center_ids"]), dict(app["S"]["lcharge"]))
    app["S"]["scan_atoms"] = [
        {"index": 1, "chain": "A", "resn": "LIG", "resi": "10", "atom": "O1",
         "xyz": (1.2, 0.0, 0.0)},
        None,
    ]
    app["S"]["freeze_atoms"] = [2]
    app["S"]["measure_atoms"] = [
        {"index": 0, "chain": "A", "resn": "LIG", "resi": "10", "atom": "C1",
         "xyz": (0.0, 0.0, 0.0)},
    ]
    calls.clear()
    app["view_input"].value = 1
    assert app["S"]["_view_mapping_ok"] is False
    assert app["center_widget"] is primary_widget
    assert primary_widget.disabled is True
    assert (app["S"]["center"], app["S"]["center_ids"], app["S"]["lcharge"]) == saved
    assert any(
        isinstance(value, str) and '"type":"rx-load-structure"' in value
        and '"generation":%d' % app["_VIEWER_GENERATION"]["value"] in value
        for value in calls
    )
    app["dd_subcmd"].value = "scan"
    def _descendants(widget):
        yield widget
        for child in getattr(widget, "children", ()):
            yield from _descendants(child)
    clear_pair = next(
        widget for widget in _descendants(app["scan_panel"])
        if getattr(widget, "description", "") == "Clear current pair"
    )
    assert clear_pair.disabled
    before_scan = json.dumps(app["S"]["scan_atoms"], sort_keys=True)
    clear_pair.click()
    assert json.dumps(app["S"]["scan_atoms"], sort_keys=True) == before_scan
    app["view_input"].value = 0
    assert app["center_widget"].disabled is False

    old_text = app["S"]["_pdb_text"]
    app["_load_view_structure"] = lambda path: (_ for _ in ()).throw(ValueError("view failed")) \
        if str(path) == str(secondary) else load_view(path)
    app["view_input"].value = 1
    assert app["S"]["_view_input_index"] == 0
    assert app["view_input"].value == 0
    assert app["S"]["_pdb_text"] == old_text

    app["S"].update(_last_out_dir="old", _last_subcmd="opt", _last_argv=["old"],
                    _last_files=[str(primary)], _last_manifest={"status": "success"},
                    _last_log="old log")
    app["artifact_choice"].options = [("old", str(primary))]
    app["dl_btn"].disabled = False
    validate_out = tmp_path / "validate-only"
    app["cmd_box"].value = shlex.join([
        "pdb2reaction", "tsopt", "-i", str(primary), "--calc-file", str(calc_file),
        "--ref-mode", str(ref_mode), "-o", str(validate_out),
    ])
    assert app["S"]["_last_manifest"] == {} and app["S"]["_last_files"] == []
    assert app["artifact_choice"].disabled and app["dl_btn"].disabled
    app["S"]["_last_log"] = "old log"
    app["_stream"] = lambda argv: (0, "validation transcript")
    app["_do_validate"](None)
    assert app["S"]["_last_log"] == "old log"
    assert app["_RUN_STATE"]["validation_log"] == "validation transcript"
    app["_stream"] = lambda argv: (2, "invalid options")
    app["_do_validate"](None)
    assert app["w_show_run_log"].value is True
    assert app["logbox"].layout.display == ""

    # One path-position change drives both the Mol* structure and the energy
    # cursor/status from the same slider value.
    trajectory = tmp_path / "linked_path_trj.xyz"
    trajectory.write_text(
        "2\n-1.000000\n"
        "H 0.000 0.000 0.000\nH 0.700 0.000 0.000\n"
        "2\n-0.990000\n"
        "H 1.500 0.000 0.000\nH 2.200 0.000 0.000\n",
        encoding="utf-8",
    )
    calls.clear()
    app["S"]["_last_subcmd"] = "path-opt"
    app["_load_trajectory"](str(trajectory), str(tmp_path))
    assert app["frame_slider"].max == 1 and not app["frame_slider"].disabled
    assert app["trajectory_box"].layout.display == ""
    assert "Frame 1 of 2" in app["frame_state"].value
    assert "ΔE = 0.0 kcal/mol" in app["frame_state"].value
    assert not [
        value for value in calls
        if isinstance(value, str) and 'class="rxmolstar-frame"' in value
    ]
    trajectory_viewer = app["traj_out"].value
    assert trajectory_viewer.count('class="rxmolstar-frame"') == 1
    assert "H 0.000 0.000 0.000" in trajectory_viewer
    assert "H 1.500 0.000 0.000" in trajectory_viewer
    viewer_generation = app["_TRAJ"]["generation"]

    calls.clear()
    app["frame_slider"].value = 1
    assert "Frame 2 of 2" in app["frame_state"].value
    assert "ΔE = 6.3 kcal/mol" in app["frame_state"].value
    assert app["traj_out"].value.count('class="rxmolstar-frame"') == 1
    assert not [
        value for value in calls
        if isinstance(value, str) and 'class="rxmolstar-frame"' in value
    ]
    assert html.unescape(app["traj_signal_out"].value).count(
        '"type":"rx-set-frame"') == 1
    assert app["_TRAJ"]["generation"] == viewer_generation
    assert app["_TRAJ"]["frame_message"] == {
        "type": "rx-set-frame", "generation": viewer_generation, "index": 1,
    }
    assert not app["frame_prev"].disabled and app["frame_next"].disabled

    app["_invalidate_last_run"]("Input identity changed; run again.")
    assert app["S"]["_last_manifest"] == {} and app["S"]["_last_files"] == []
    assert app["artifact_choice"].disabled and app["dl_btn"].disabled
    assert "Input identity changed" in app["results_empty"].value


def test_colab_run_fails_closed_when_implicit_validation_fails(
    tmp_path: Path, monkeypatch,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)
    streamed = []
    app["_sync_run"] = lambda: None
    app["_argv"] = lambda: [
        "pdb2reaction", "opt", "-i", "missing.pdb", "-o", "result",
    ]
    app["_stream"] = lambda argv: (streamed.append(list(argv)) or (2, "invalid"))

    app["_do_run"](None)

    assert len(streamed) == 1
    assert "--dry-run" in streamed[0]
    assert app["_RUN_STATE"]["validated_fingerprint"] is None
    assert "validation failed" in app["run_status"].value


def test_colab_manifest_hashes_inplace_input_before_execution(
    tmp_path: Path, monkeypatch,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)
    structure = tmp_path / "structure.pdb"
    structure.write_text("ORIGINAL\n", encoding="utf-8")
    original_digest = app["_sha256"](str(structure))
    argv = [
        "pdb2reaction", "add-elem-info", "-i", str(structure), "--overwrite",
    ]
    app["_sync_run"] = lambda: None
    app["_argv"] = lambda: list(argv)
    app["_output_scope"] = lambda _argv: {
        "target": str(structure), "root": str(tmp_path), "shallow": True,
        "stdout_only": False, "direct_current": True,
    }
    app["_output_scope_collision"] = lambda _scope: False
    app["_snapshot_output_scope"] = lambda _scope: {}
    app["_begin_results_attempt"] = lambda *_args: None
    app["_results"] = lambda *_args: None

    def mutate_input(_argv):
        structure.write_text("MUTATED\n", encoding="utf-8")
        return 0, "updated in place"

    app["_stream"] = mutate_input
    app["_do_run"](None)

    assert app["S"]["_last_manifest"]["inputs"] == [{
        "path": str(structure), "sha256": original_digest,
    }]
    assert app["_sha256"](str(structure)) != original_digest


def test_colab_compact_selection_upload_viewer_and_advanced_contracts(
    tmp_path: Path, monkeypatch,
) -> None:
    app, calls = _execute_app(monkeypatch, tmp_path)
    primary = tmp_path / "primary.pdb"
    secondary = tmp_path / "secondary.pdb"
    primary.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  H   ALA A   1       0.000   1.000   0.000  1.00  0.00           H\n"
        "HETATM    3 MG    MG A   2       2.000   0.000   0.000  1.00  0.00          MG\n"
        "HETATM    4  O   WAT A  10       3.000   0.000   0.000  1.00  0.00           O\n"
        "HETATM    5  C1  LIG A   3       4.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    6  O1  LIG A   3       5.200   0.000   0.000  1.00  0.00           O\n"
        "END\n",
        encoding="utf-8",
    )
    secondary.write_text(
        "HETATM    1  C1  ALT B   4       8.000   0.000   0.000  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    metadata = {
        str(primary): [
            {"serial": 1, "chain": "A", "resname": "ALA", "resseq": 1, "icode": "", "name": "N"},
            {"serial": 2, "chain": "A", "resname": "ALA", "resseq": 1, "icode": "", "name": "H"},
            {"serial": 3, "chain": "A", "resname": "MG", "resseq": 2, "icode": "", "name": "MG"},
            {"serial": 4, "chain": "A", "resname": "WAT", "resseq": 10, "icode": "", "name": "O"},
            {"serial": 5, "chain": "A", "resname": "LIG", "resseq": 3, "icode": "", "name": "C1"},
            {"serial": 6, "chain": "A", "resname": "LIG", "resseq": 3, "icode": "", "name": "O1"},
        ],
        str(secondary): [
            {"serial": 1, "chain": "B", "resname": "ALT", "resseq": 4, "icode": "", "name": "C1"},
        ],
    }

    def load_view(path):
        path = str(path)
        return Path(path).read_text(encoding="utf-8"), [dict(row) for row in metadata[path]], path

    app["_load_view_structure"] = load_view
    assert app["_ingest_saved_files"]([str(primary)], "test drop")
    assert app["all_mode"].value == "scan"
    assert app["pick_action"].value == "scanA"
    assert app["_ingest_saved_files"]([str(secondary)], "test drop")
    assert app["all_mode"].value == "mep"
    assert app["S"]["inputs"] == [str(primary), str(secondary)]
    assert len(app["input_file_rows"].children) == 2
    assert app["prep_radius"].value == pytest.approx(2.6)
    assert app["adv_radius"].value == pytest.approx(2.6)
    assert app["dd_col"].value == "element"
    assert app["dd_rep"].value == "cartoon"
    assert app["dd_size"].value == 720
    assert app["S"]["viewer_width"] == 720
    assert app["S"]["viewer_height"] == 540

    # Information is a native details/summary disclosure: Colab opens it in the
    # browser without waiting for a Python widget callback.
    info = app["last_pick_info"]
    assert '<details class="rxinfo-details"' in info.value
    assert "<summary" in info.value and "&#9432;" in info.value
    assert 'role="note"' in info.value
    assert " open" not in info.value
    second_info = app["_info_control"]("second")
    assert "second" in second_info.value
    assert 'aria-label="More information: second"' in second_info.value
    assert 'title="second"' in second_info.value
    app["_set_info_text"](second_info, "updated")
    assert "updated" in second_info.value and "second" not in second_info.value
    before_close = second_info.value
    app["_close_info"](second_info)
    assert second_info.value != before_close
    assert "<details" in second_info.value and " open" not in second_info.value
    assert app["key_opts_box"].children[1].layout.display == "none"
    assert app["command_editor"].children[1] is app["cmd_box"]
    assert app["command_editor"].layout.display != "none"
    assert app["w_show_run_log"].value is True
    assert app["logbox"].layout.display == ""
    app["w_show_run_log"].value = False
    assert app["logbox"].layout.display == "none"
    app["w_show_run_log"].value = True
    assert app["_input_box_children"][1] is app["_drop"]
    assert app["_input_box_children"][2] is app["input_msg"]
    assert app["_input_box_children"][-1] is app["example_msg"]
    assert "Click sets:" in app["last_pick_html"].value
    assert "rxcard" in app["depth_box"]._dom_classes

    # The name picker contains only chemically useful hetero residues; waters
    # and canonical amino acids remain available through exact 3D clicks.
    assert set(app["_center_values"]()) == {"LIG", "MG"}
    mg_row = app["charge_rows"]["MG"]
    assert mg_row["auto"] and mg_row["use"].disabled and mg_row["val"].disabled
    assert mg_row["val"].value == 2
    assert "MG" not in app["S"]["lcharge"]

    app["center_widget"].value = ("MG",)
    with pytest.raises(ValueError, match="Confirm each ligand charge"):
        app["build_cmd"]()
    app["charge_rows"]["LIG"]["use"].value = True
    # New charge routing (author decision, 2026-07-27): a ligand-charge source
    # means the CLI derives the total charge, so the notebook emits no -q at
    # all. The checkbox became an explicit override, live only because -l is
    # present, and it keeps -q read-only until it is ticked.
    derived_command = app["build_cmd"]()
    assert "-q" not in derived_command
    assert app["w_charge_ok"].description == "overwrite system charge"
    assert not app["w_charge_ok"].disabled and app["w_q"].disabled
    app["w_q"].value = -1
    app["w_charge_ok"].value = True
    mg_row["val"].value = 3  # disabled widgets are still mutable from Python
    # Locked ion rows are content-bound to the built-in table, not to a
    # programmatic mutation of the disabled display widget.
    assert app["w_charge_ok"].value is True
    charged_command = app["build_cmd"]()
    assert charged_command[charged_command.index("-q") + 1] == "-1"
    ligand_values = dict(
        item.split(":", 1)
        for item in charged_command[charged_command.index("-l") + 1].split(",")
    )
    assert ligand_values == {"LIG": "0", "MG": "2"}
    confirmed_scope = app["S"]["_charge_scope"]
    app["all_mode"].value = "scan"
    assert app["w_charge_ok"].value is True
    assert app["S"]["_charge_scope"] == confirmed_scope
    app["all_mode"].value = "mep"
    app["_set_advanced_override"]("all", "workers", "2")
    assert app["w_charge_ok"].value is True
    app["_set_advanced_override"]("all", "include_h2o", False)
    assert app["w_charge_ok"].value is False
    mg_row["val"].value = 2
    app["charge_rows"]["LIG"]["use"].value = False
    app["center_widget"].value = ()
    app["w_charge_ok"].value = False

    calls.clear()
    app["cb_water"].value = True
    assert app["S"]["show_water"] is True
    updates = [
        value for value in calls
        if isinstance(value, str) and '"type":"rx-load-structure"' in value
    ]
    assert updates and '"showWater":true' in updates[-1]
    document = app["_molstar_document"](
        app["S"]["_pdb_text"], "pdb", show_water=True, interactive=True,
    )
    assert '"showWater":true' in document
    assert "layoutShowControls:true" in document
    assert "layoutShowSequence:cfg.showSequence" in document
    assert "structure-component-static-water" in document
    assert "WebGL is unavailable in this browser session." in document
    cif_document = app["_molstar_document"](
        "data_demo\n_entry.id demo\n", "cif", show_sequence=True,
    )
    assert '"format":"mmcif"' in cif_document
    assert '"showSequence":true' in cif_document

    # Serial/signature mapping wins over a drifted browser source index. The
    # keyboard fallback rebuilds; live Mol* clicks mutate browser state only.
    calls.clear()
    app["pick_action"].value = "center"
    app["on_click"]("1", "LIG", "3", "A", "C1", "5", "")
    assert app["S"]["_last_pick"]["index"] == 4
    assert app["S"]["_last_pick"]["viewer_index"] == 1
    assert any(
        isinstance(value, str) and '"type":"rx-load-structure"' in value
        for value in calls
    )
    assert [pick["index"] for pick in app["S"]["_pick_history"]] == [4]

    # A second Mol* click updates workflow state without publishing a
    # replacement iframe; Mol* retains ownership of its standard focus state.
    calls.clear()
    app["on_click"]("2", "MG", "2", "A", "MG", "3", "", live_marked=True)
    assert [pick["index"] for pick in app["S"]["_pick_history"]] == [4, 2]
    assert not any(
        isinstance(value, str)
        and ('class="rxmolstar-frame"' in value or '"type":"rx-load-structure"' in value)
        for value in calls
    )
    committed_centers = list(app["S"]["center_ids"])
    app["_clear_highlights_from_browser"](app["_VIEWER_GENERATION"]["value"])
    assert app["S"]["_pick_history"] == []
    assert app["S"]["center_ids"] == committed_centers

    # A later click is recorded in workflow state; visual focus remains Mol*'s
    # default and is intentionally not serialized by the notebook.
    app["on_click"]("2", "MG", "2", "A", "MG", "3", "", live_marked=True)
    assert [pick["index"] for pick in app["S"]["_pick_history"]] == [2]

    # A row-local × removes only that file and keeps primary-bound selections.
    close_secondary = app["input_file_rows"].children[1].children[1].children[-1]
    close_secondary.click()
    assert app["S"]["inputs"] == [str(primary)]
    assert app["S"]["center_ids"] == ["A:LIG:3", "A:MG:2"]
    assert secondary.exists()  # source/example paths are detached, never deleted
    owned = tmp_path / "owned.pdb"
    owned.write_text(primary.read_text(encoding="utf-8"), encoding="utf-8")
    app["_remember_uploaded_paths"]([str(owned)])
    app["S"]["inputs"].append(str(owned))
    app["_queue_change"](path=str(owned), remove=True)
    assert not owned.exists()

    # Every live Click option is classified, while editable rows and
    # serializers are generated from this repository's own command objects.
    all_statuses = {}
    root = app["PRODUCT_CLI"]
    for subcommand in root.list_commands(app["click"].Context(root)):
        options = app["_advanced_options"](subcommand)
        coverage = app["_advanced_coverage"](subcommand)
        assert set(coverage) == {param.name for param in options}
        assert set(coverage.values()) <= {"owned", "generated", "blocked", "rendered"}
        all_statuses.update(coverage)
    assert all_statuses["tr_projection"] == "blocked"
    assert all_statuses["verbose"] == "rendered"
    assert app["_advanced_coverage"]("dft")
    assert app["_advanced_coverage"]("sp")["hessian_calc_mode"] == "rendered"
    verbose_param = next(
        param for param in app["_advanced_options"]("all")
        if param.name == "verbose"
    )
    app["S"]["advanced_overrides"]["all"] = {"verbose": "7"}
    verbose_row = app["_advanced_widget"]("all", verbose_param)
    assert verbose_row.children[0].value is None
    assert "verbose" not in app["S"]["advanced_overrides"]["all"]

    app["dd_subcmd"].options = app["_sub_options"](app["SUBS"])
    app["dd_subcmd"].value = "add-elem-info"
    app["S"]["advanced_overrides"]["add-elem-info"] = {}
    safe_add_elem = app["build_cmd"]()
    assert "-o" in safe_add_elem and "--overwrite" not in safe_add_elem
    app["S"]["advanced_overrides"]["add-elem-info"] = {"overwrite": True}
    inplace_add_elem = app["build_cmd"]()
    assert "-o" not in inplace_add_elem and "--overwrite" in inplace_add_elem

    app["dd_subcmd"].value = "all"
    app["all_mode"].value = "mep"
    app["_render_advanced_rows"]()
    editable = [
        param for param in app["_advanced_options"]("all")
        if app["_advanced_status"]("all", param) == "rendered"
        and app["_advanced_semantic_applicable"]("all", param.name)
    ]
    assert len(app["adv_rows_box"].children) == len(editable)
    app["S"]["advanced_overrides"]["all"] = {
        "radius_het2het": "4.5", "include_h2o": False,
    }
    advanced_argv = app["_advanced_argv"]("all")
    assert advanced_argv[advanced_argv.index("--radius-het2het") + 1] == "4.5"
    assert advanced_argv[advanced_argv.index("--include-h2o") + 1] == "false"

    sp_hessian = next(param for param in app["_advanced_options"]("sp")
                      if param.name == "hessian_calc_mode")
    hessian_value = (list(sp_hessian.type.choices)[0]
                     if isinstance(sp_hessian.type, app["click"].Choice) else "analytical")
    app["S"]["advanced_overrides"]["sp"] = {"hessian_calc_mode": hessian_value}
    assert "--hessian-calc-mode" in app["_advanced_argv"]("sp")

    for utility, flag in (("trj2fig", "--out-json"),
                          ("energy-diagram", "--out-json"),
                          ("bond-summary", "--json")):
        out_json = next(param for param in app["_advanced_options"](utility)
                        if param.name == "out_json")
        assert app["_advanced_status"](utility, out_json) == "rendered"
        app["S"]["advanced_overrides"][utility] = {"out_json": True}
        assert flag in app["_advanced_argv"](utility)

    app["dd_subcmd"].value = "energy-diagram"
    app["S"]["advanced_overrides"]["energy-diagram"] = {}
    app["refresh"]()
    assert not app["_ACTION_STATE"]["auto_ready"]
    assert app["b_run"].disabled
    assert "Advanced flags" in app["subcmd_note"].value
    for invalid_values in (
        {"input_values": "0"},
        {"input_values": "0 oops"},
        {"input_values": "0 1", "label_x": "R"},
    ):
        app["S"]["advanced_overrides"]["energy-diagram"] = invalid_values
        invalid_command = app["build_cmd"]()
        assert not app["_utility_autofill_complete"](invalid_command)
        app["refresh"]()
        assert app["b_run"].disabled
    app["S"]["advanced_overrides"]["energy-diagram"] = {
        "input_values": "0 12.5 4.3",
        "output_path": str(tmp_path / "profile.png"),
    }
    energy_command = app["build_cmd"]()
    assert energy_command.count("--input") == 3
    assert app["_utility_autofill_complete"](energy_command)
    app["refresh"]()
    assert app["_ACTION_STATE"]["auto_ready"]
    assert not app["b_run"].disabled

    atom_a = {"index": 0, "chain": "A", "resn": "ALA", "resi": "1", "atom": "N"}
    atom_b = {"index": 1, "chain": "A", "resn": "ALA", "resi": "1", "atom": "H"}
    app["S"]["freeze_pairs"] = [{"a": atom_a, "b": atom_b, "t": None}]
    app["dd_subcmd"].value = "opt"
    opt_picks = {value for _label, value in app["pick_action"].options}
    assert {"freezeA", "freezeB", "freezeatom"} <= opt_picks
    assert "Distance restraints" in " ".join(
        getattr(child, "value", "") for child in app["freeze_panel"].children
    )
    app["dd_subcmd"].value = "tsopt"
    tsopt_picks = {value for _label, value in app["pick_action"].options}
    assert "freezeatom" in tsopt_picks
    assert not {"freezeA", "freezeB"} & tsopt_picks
    assert "Distance restraints" not in " ".join(
        getattr(child, "value", "") for child in app["freeze_panel"].children
    )
    assert len(app["S"]["freeze_pairs"]) == 1
    app["_render_summary"]()
    assert "freeze:" not in app["summary_html"].value
    app["center_widget"].value = ()
    app["S"]["center_ids"] = []
    app["charge_rows"]["LIG"]["use"].value = False
    app["w_q"].value = 0
    app["w_charge_ok"].value = True
    assert "--dist-freeze" not in app["build_cmd"]()
    app["dd_subcmd"].value = "opt"
    assert "Distance restraints" in " ".join(
        getattr(child, "value", "") for child in app["freeze_panel"].children
    )
    app["_render_summary"]()
    assert "freeze:" in app["summary_html"].value

    app["dd_subcmd"].value = "sp"
    assert app["key_opts_box"].layout.display == ""
    assert [value for _label, value in app["pick_action"].options] == [
        "center", "ligand",
    ]
    assert app["center_panel"].layout.display == ""
    assert app["charge_panel"].layout.display == ""
    assert app["extract_panel"].layout.display == ""
    assert app["adv_radius"]._rx_flag_row.layout.display == "none"
    assert app["adv_dftfb"]._rx_flag_row.layout.display == "none"
    app["dd_subcmd"].value = "freq"
    assert app["key_opts_box"].layout.display == ""
    app["dd_subcmd"].options = app["_sub_options"](app["SUBS"])
    app["dd_subcmd"].value = "extract"
    assert app["key_opts_box"].layout.display == ""
    assert app["adv_radius"]._rx_flag_row.layout.display == ""
    app["center_widget"].value = ()
    app["S"]["center_ids"] = []
    with pytest.raises(ValueError, match="needs a center residue"):
        app["build_cmd"]()
    app["center_widget"].value = ("LIG",)
    app["dd_subcmd"].value = "add-elem-info"
    assert app["center_panel"].layout.display == "none"
    assert app["charge_panel"].layout.display == "none"
    assert app["extract_panel"].layout.display == "none"
    app["dd_subcmd"].value = "all"
    app["all_mode"].value = "scan"
    assert {"scanA", "scanB"} <= {value for _label, value in app["pick_action"].options}
    app["all_mode"].value = "mep"
    assert not {"scanA", "scanB"} & {value for _label, value in app["pick_action"].options}

    # TS-only is a semantic mode, not only a visibility change: stale path
    # controls must neither remain visible nor leak into the built command.
    app["adv_refine"].value = True
    app["adv_mep"].value = "gsm"
    app["adv_thresh"].value = "gau"
    app["adv_maxcyc"].value = 9
    app["S"]["advanced_overrides"]["all"] = {
        "reject_uphill": False,
        "opt_mode": "grad",
        "preopt": False,
        "irc_step_size": "0.05",
        "opt_mode_post": "grad",
        "thresh_post": "baker",
        "hessian_calc_mode": "FiniteDifference",
    }
    app["w_th"].value = False
    app["adv_dft"].value = False
    app["w_ts"].value = False
    assert app["w_ts"].value is False and not app["w_ts"].disabled
    off_argv = app["_advanced_argv"]("all")
    for flag in (
        "--irc-step-size", "--opt-mode-post", "--thresh-post",
        "--hessian-calc-mode", "--no-reject-uphill",
    ):
        assert flag not in off_argv
    app["all_mode"].value = "tsonly"
    assert app["adv_refine"].layout.display == "none"
    assert app["w_ts"].value and app["w_ts"].disabled
    app["charge_rows"]["LIG"]["use"].value = True
    app["w_charge_ok"].value = True
    command = app["build_cmd"]()
    for flag in ("--refine-path", "--mep-mode", "--thresh", "--max-cycles"):
        assert flag not in command
    assert "--no-reject-uphill" in command
    assert command[command.index("--opt-mode") + 1] == "grad"
    assert command[command.index("--irc-step-size") + 1] == "0.05"
    assert command[command.index("--opt-mode-post") + 1] == "grad"
    assert command[command.index("--thresh-post") + 1] == "baker"
    assert command[command.index("--hessian-calc-mode") + 1] == "FiniteDifference"
    assert "--preopt" not in command and "--no-preopt" not in command
    assert "--tsopt" in command
    app["all_mode"].value = "mep"
    assert not app["w_ts"].value and not app["w_ts"].disabled
    app["w_th"].value = True
    assert app["w_ts"].value and app["w_ts"].disabled
    app["w_ts"].value = False
    assert app["w_ts"].value and app["w_ts"].disabled
    thermo_argv = app["_advanced_argv"]("all")
    assert "--hessian-calc-mode" in thermo_argv
    for flag in (
        "--irc-step-size", "--opt-mode-post", "--thresh-post",
        "--no-reject-uphill",
    ):
        assert flag in thermo_argv
    app["S"]["inputs"] = [str(primary), str(secondary)]
    app["w_charge_ok"].value = False
    app["w_charge_ok"].value = True
    thermo_cmd = app["build_cmd"]()
    assert thermo_cmd.count("--tsopt") == 1
    assert thermo_cmd.count("--thermo") == 1
    app["w_th"].value = False
    assert app["w_ts"].value and not app["w_ts"].disabled
    app["w_ts"].value = False
    assert "--tsopt" not in app["build_cmd"]()
    app["S"]["tsopt"] = False
    app["S"]["thermo"] = True
    with pytest.raises(ValueError, match="require TS optimization"):
        app["build_cmd"]()
    app["S"]["thermo"] = False

    result_json = tmp_path / "result.json"
    result_json.write_text('{"energy": -1.25}', encoding="utf-8")
    assert app["_artifact_kind"](str(result_json)) == "JSON"
    assert app["_artifact_kind"]("job.gjf") == "text"
    assert app["_artifact_kind"]("job.com") == "text"
    assert app["_artifact_kind"]("job.inp") == "text"
    preview = app["_text_preview_html"](str(result_json), "JSON")
    assert "&quot;energy&quot;" in preview and "-1.25" in preview
    summary = tmp_path / "summary.json"
    summary.write_text(json.dumps({
        "status": "success", "scientific_status": "partial",
        "scientific_status_reasons": ["IRC endpoint mismatch"],
        "segments": [{"index": 1, "barrier_kcal": 8.0, "delta_kcal": -1.0}],
        "post_segments": [{"index": 1, "mlip": {"barrier_kcal": 7.5, "delta_kcal": -1.2}}],
        "rate_limiting_step": {"barrier_kcal": 7.5, "method": "mlip"},
    }), encoding="utf-8")
    summary_html = app["_summary_html"](str(summary))
    assert "Highest local barrier: 7.5" in summary_html
    assert "IRC endpoint mismatch" in summary_html
    assert "MLIP" in summary_html and "partial result" in summary_html
    ts_only_summary = tmp_path / "ts_only_summary.json"
    ts_only_summary.write_text(json.dumps({
        "status": "success",
        "scientific_status": "success",
        "pipeline_mode": "tsopt-only",
        "n_images": 5,
        "mlip_backend": "mace",
        "mlip_model": "MACE-OMOL-0",
        "endpoint_assignment": {"chemical_direction_known": False},
        "segments": [{
            "index": 1, "kind": "tsopt",
            "barrier_kcal": 8.0, "delta_kcal": -1.0,
        }],
        "rate_limiting_step": {
            "segment": 1, "barrier_kcal": 8.0, "method": "MLIP",
        },
    }), encoding="utf-8")
    ts_only_html = app["_summary_html"](str(ts_only_summary))
    assert "raw MEP" not in ts_only_html
    assert "IRC frames: 5" in ts_only_html
    assert "energy order" in ts_only_html
    extra_artifacts = []
    for name in ("a.txt", "b.txt", "c.txt", "z.txt"):
        artifact = tmp_path / name
        artifact.write_text(name, encoding="utf-8")
        extra_artifacts.append(str(artifact))
    app["S"].update(
        _last_out_dir=str(tmp_path),
        _last_subcmd="all",
        _last_files=[str(summary), *extra_artifacts],
        _last_manifest={"status": "success", "exit_code": 0},
    )
    context_html = app["_result_context_html"](str(tmp_path))
    assert "IRC endpoint mismatch" not in context_html
    assert "<b>5 files</b>" in context_html
    assert context_html.count("<code>") == 3
    calls.clear()
    app["_structure_preview"](str(primary))
    assert any(
        isinstance(value, str) and 'class="rxmolstar-frame"' in value
        and "Mol*" in value
        for value in calls
    )

    # The advertised small-molecule example uses the same removable queue
    # transaction and opens a real XYZ viewer.
    calls.clear()
    example_dir = NOTEBOOK.parent / "aromatic_claisen"
    app["_example_file"] = lambda relpath: str(example_dir / Path(relpath).name)
    app["ex_choice"].value = "Aromatic Claisen rearrangement - small molecule"
    app["ex_btn"].click()
    assert len(app["S"]["inputs"]) == 2
    assert len(app["input_file_rows"].children) == 2
    assert app["all_mode"].value == "mep" and not app["b_clear_inputs"].disabled
    assert app["S"]["_view_format"] == "xyz"
    assert app["workspace"].layout.display == ""
    assert any(
        isinstance(value, str)
        and ('class="rxmolstar-frame"' in value or '"type":"rx-load-structure"' in value)
        for value in calls
    )
    # The other half of the routing: with no ligand-charge source the notebook
    # emits -q explicitly, because the CLI requires it for non-.gjf input, and
    # the override checkbox stays dead because there is nothing to override.
    xyz_command = app["build_cmd"]()
    assert xyz_command[xyz_command.index("-q") + 1] == "0"
    assert app["w_charge_ok"].disabled
    assert "-r" not in xyz_command
    app["b_clear_inputs"].click()
    assert app["S"]["inputs"] == []

    # A precomputed scan3d surface is a first-class plot-only input.
    surface = tmp_path / "surface.csv"
    surface.write_text("d1,d2,d3,energy\n1,1,1,0\n", encoding="utf-8")
    assert app["_ingest_saved_files"]([str(surface)], "test csv")
    csv_command = app["build_cmd"]()
    assert csv_command[:2] == ["pdb2reaction", "scan3d"]
    assert "--csv" in csv_command and "-i" not in csv_command

    app["load_pdb"]([str(primary)], center=["LIG"], mode="mmcif")
    original_inputs = list(app["S"]["inputs"])
    app["dd_subcmd"].value = "sp"
    app["prep_radius"].value = 4.2
    assert app["adv_radius"].value == 4.2
    app["charge_rows"]["LIG"]["use"].value = True
    extract_commands = []

    def fake_extract(command, **_kwargs):
        extract_commands.append(list(command))
        out_start = command.index("-o") + 1
        out_end = command.index("-c")
        for output in command[out_start:out_end]:
            Path(output).write_text(primary.read_text(encoding="utf-8"), encoding="utf-8")
            metadata[str(output)] = [dict(row) for row in metadata[str(primary)]]
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(app["subprocess"], "run", fake_extract)
    app["b_extract"].click()
    assert extract_commands
    assert extract_commands[-1][extract_commands[-1].index("-r") + 1] == "4.2"
    assert extract_commands[-1][extract_commands[-1].index("-l") + 1] == "LIG:0,MG:2"
    assert app["dd_subcmd"].value == "sp" and app["S"]["subcmd"] == "sp"
    assert app["S"]["_pre_extract"]["inputs"] == original_inputs
    assert app["b_revert"].layout.display == ""
    assert app["S"]["inputs"] != original_inputs
    app["b_revert"].click()
    assert app["dd_subcmd"].value == "sp" and app["S"]["subcmd"] == "sp"
    assert app["S"]["inputs"] == original_inputs
    assert app["S"]["mode"] == "mmcif"
    assert app["S"]["_pre_extract"] is None
    assert app["b_revert"].layout.display == "none"

    small = tmp_path / "small.xyz"
    small.write_text(
        "3\nsmall\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\nH 1.0 0.0 0.0\n",
        encoding="utf-8",
    )
    app["S"].update(mode="small", inputs=[str(small)])
    app["build_selection"]()
    app["_sync_capability_controls"]()
    assert app["workspace"].layout.display == ""
    assert app["selection_help"].layout.display == "none"
    assert app["selection_route"].layout.display == "none"
    assert app["center_panel"].layout.display == "none"
    assert app["charge_panel"].layout.display == "none"
    assert app["extract_panel"].layout.display == "none"
    assert app["S"]["_view_format"] == "xyz"

    app["S"]["_pre_extract"] = {"inputs": [str(primary)]}
    app["b_revert"].layout.display = ""
    app["_clear_structure_bound_state"]()
    assert app["S"]["_pre_extract"] is None
    assert app["b_revert"].layout.display == "none"


def test_colab_exercises_every_workflow_and_advanced_flag_widget(
    tmp_path: Path, monkeypatch,
) -> None:
    """Operate every workflow selector and every editable live Click option."""
    app, _ = _execute_app(monkeypatch, tmp_path)
    # The dropdown deliberately offers only the scientific workflows (author
    # decision, 2026-07-27). The utility subcommands stay in SUBS so the
    # command editor can still validate them, but they are not selectable, and
    # restoring a session that names one silently falls back to 'all'.
    option_values = {
        item[1] if isinstance(item, tuple) else item
        for item in app["dd_subcmd"].options
    }
    assert option_values == set(app["BASIC_SUBS"])
    utility_subs = set(app["SUBS"]) - set(app["BASIC_SUBS"])
    assert utility_subs and not utility_subs & option_values
    assert (
        "dd_subcmd.value = d['subcmd'] if d['subcmd'] in BASIC_SUBS else 'all'"
        in _notebook()["cells"][2]["source"]
    )

    for subcommand in sorted(option_values):
        app["dd_subcmd"].value = subcommand
        assert app["S"]["subcmd"] == subcommand
        assert app["subreq"].value and app["outputs_html"].value
        app["_render_advanced_rows"]()

    root = app["PRODUCT_CLI"]
    root_context = app["click"].Context(root)
    for subcommand in root.list_commands(root_context):
        rendered = {
            param.name for param in app["_advanced_options"](subcommand)
            if app["_advanced_status"](subcommand, param) == "rendered"
        }
        exercised: set[str] = set()
        scenarios = [None]
        if subcommand == "all":
            scenarios = [
                ("mep", False, False, False),
                ("scan", False, False, False),
                ("tsonly", True, False, False),
                ("mep", True, True, True),
            ]
        for scenario in scenarios:
            if scenario is not None:
                mode, tsopt, thermo, dft = scenario
                app["all_mode"].value = mode
                app["w_ts"].value = tsopt
                app["w_th"].value = thermo
                app["adv_dft"].value = dft
            for param in app["_advanced_options"](subcommand):
                if (
                    app["_advanced_status"](subcommand, param) != "rendered"
                    or not app["_advanced_semantic_applicable"](
                        subcommand, param.name
                    )
                ):
                    continue
                app["S"].setdefault("advanced_overrides", {})[subcommand] = {}
                row = app["_advanced_widget"](subcommand, param)
                widget = row.children[0]
                for value in _advanced_widget_values(app, param):
                    widget.value = value
                    assert (
                        app["S"]["advanced_overrides"][subcommand][param.name]
                        == value
                    )
                    argv = app["_advanced_argv"](subcommand)
                    _assert_advanced_argv_parses(app, subcommand, argv)
                widget.value = None if hasattr(widget, "options") else ""
                assert (
                    param.name
                    not in app["S"]["advanced_overrides"][subcommand]
                )
                exercised.add(param.name)
        assert exercised == rendered, (
            subcommand,
            sorted(rendered - exercised),
        )


def test_colab_operates_every_workflow_through_validate_run_and_results(
    tmp_path: Path, monkeypatch,
) -> None:
    """Build, parse, validate, run, and render every GUI workflow state."""
    app, _ = _execute_app(monkeypatch, tmp_path)
    pdb_a = tmp_path / "reactant.pdb"
    pdb_b = tmp_path / "product.pdb"
    trajectory = tmp_path / "path.xyz"
    pdb_text = (
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O1  LIG A   1       1.200   0.000   0.000  1.00  0.00           O\n"
        "HETATM    3  N1  LIG A   1       0.000   1.200   0.000  1.00  0.00           N\nEND\n"
    )
    pdb_a.write_text(pdb_text, encoding="utf-8")
    pdb_b.write_text(pdb_text.replace("1.200", "1.300"), encoding="utf-8")
    trajectory.write_text(
        "2\n0.0\nH 0 0 0\nH 0 0 1\n2\n1.0\nH 0 0 0\nH 0 0 1.1\n",
        encoding="utf-8",
    )

    app["DFT_READY"] = True
    if "dft" not in app["SUBS"]:
        app["SUBS"].append("dft")
    app["dd_subcmd"].options = app["_sub_options"](app["SUBS"])

    streamed: list[list[str]] = []
    results: list[tuple[str, str]] = []
    app["_stream"] = lambda argv: (
        streamed.append(list(argv)) or (0, "operation-matrix success")
    )
    app["_output_scope"] = lambda argv: {
        "target": str(tmp_path / ("out-" + argv[1])),
        "root": str(tmp_path / ("out-" + argv[1])),
        "shallow": False,
        "stdout_only": False,
        "direct_current": False,
    }
    app["_output_scope_collision"] = lambda _scope: False
    app["_snapshot_output_scope"] = lambda _scope: {}
    app["_begin_results_attempt"] = lambda root, sub: results.append((root, sub))
    app["_results"] = lambda root: results.append((root, "rendered"))

    atoms = [
        {"index": index, "chain": "A", "resn": "LIG", "resi": "1",
         "atom": name, "xyz": xyz}
        for index, (name, xyz) in enumerate((
            ("C1", (0.0, 0.0, 0.0)),
            ("O1", (1.2, 0.0, 0.0)),
            ("N1", (0.0, 1.2, 0.0)),
        ))
    ]
    primary_atom_meta = [
        {
            "index": atom["index"], "chain": atom["chain"],
            "resname": atom["resn"], "resseq": atom["resi"], "icode": "",
            "atom": atom["atom"], "xyz": atom["xyz"],
        }
        for atom in atoms
    ]
    cases = [
        ("all", "mep"), ("all", "scan"), ("all", "tsonly"),
        ("opt", None), ("sp", None), ("tsopt", None), ("freq", None),
        ("irc", None), ("dft", None), ("scan", None), ("scan2d", None),
        ("scan3d", None), ("path-opt", None), ("path-search", None),
        ("extract", None), ("fix-altloc", None), ("add-elem-info", None),
        ("energy-diagram", None), ("bond-summary", None), ("trj2fig", None),
    ]

    for index, (subcommand, all_kind) in enumerate(cases):
        streamed.clear()
        app["S"].update(
            mode="pdb", inputs=[str(pdb_a)], ref_pdbs=[], parm=None,
            model_pdb=None, center=[], center_ids=[], lcharge={},
            scan_atoms=[None, None], scan_stages=[], scan_axes=[],
            freeze_buf=[None, None], freeze_pairs=[], freeze_atoms=[],
            advanced_overrides={}, _session_file_identities={},
            _installed_backend="mace", backend="mace",
            _primary_atom_meta=primary_atom_meta,
        )
        app["w_ts"].value = False
        app["w_th"].value = False
        app["adv_dft"].value = False
        app["w_out"].value = "matrix-%02d-%s" % (index, subcommand)
        app["dd_subcmd"].value = subcommand
        app["all_mode"].value = all_kind or "mep"

        if subcommand in {"all", "path-opt", "path-search", "bond-summary"}:
            app["S"]["inputs"] = [str(pdb_a), str(pdb_b)]
        if subcommand == "all" and all_kind in {"scan", "tsonly"}:
            app["S"]["inputs"] = [str(pdb_a)]
        if subcommand == "all" and all_kind == "scan":
            bond = {"a": atoms[0], "b": atoms[1], "t": 1.5}
            app["S"]["scan_atoms"] = [atoms[0], atoms[1]]
            app["S"]["scan_stages"] = [[bond]]
        if subcommand == "scan":
            bond = {"a": atoms[0], "b": atoms[1], "t": 1.5}
            app["S"]["scan_atoms"] = [atoms[0], atoms[1]]
            app["S"]["scan_stages"] = [[bond]]
        if subcommand in {"scan2d", "scan3d"}:
            axes = [
                {"a": atoms[0], "b": atoms[1], "lo": 1.0, "hi": 2.0},
                {"a": atoms[0], "b": atoms[2], "lo": 1.0, "hi": 2.0},
                {"a": atoms[1], "b": atoms[2], "lo": 1.0, "hi": 2.0},
            ]
            app["S"]["scan_axes"] = axes[:2 if subcommand == "scan2d" else 3]
        if subcommand == "extract":
            app["S"]["center_ids"] = ["A:LIG:1"]
        elif subcommand == "energy-diagram":
            app["S"]["inputs"] = []
            app["S"]["mode"] = "utility"
            app["S"]["advanced_overrides"][subcommand] = {
                "input_values": "0 5", "output_path": str(tmp_path / "energy.png"),
            }
        elif subcommand == "trj2fig":
            app["S"]["inputs"] = [str(trajectory)]
            app["S"]["mode"] = "xyz"

        if subcommand in app["COMPUTE"]:
            app["w_q"].value = 0
            app["w_charge_ok"].value = False
            app["w_charge_ok"].value = True
        command = app["build_cmd"]()
        assert command[:2] == ["pdb2reaction", subcommand]
        click_command = app["PRODUCT_CLI"].get_command(
            app["click"].Context(app["PRODUCT_CLI"]), subcommand,
        )
        parser_argv = _root_normalized_subcommand_argv(
            app, subcommand, command[2:],
        )
        with click_command.make_context(
            subcommand, parser_argv, resilient_parsing=False,
        ) as context:
            expected_extra = (
                app["S"]["inputs"][1:]
                if subcommand in {"all", "path-search"} else []
            )
            assert list(context.args) == expected_extra, (
                subcommand, command, context.args,
            )

        app["_auto"]["on"] = False
        app["cmd_box"].value = shlex.join(command)
        if subcommand in app["COMPUTE"]:
            assert not app["b_validate"].disabled
            app["b_validate"].click()
            assert streamed and "--dry-run" in streamed[0]
        else:
            before_validate = len(streamed)
            app["b_validate"].click()
            assert len(streamed) == before_validate
        assert not app["b_run"].disabled, (subcommand, app["run_status"].value)
        app["b_run"].click()
        assert streamed[-1] == command
        assert app["S"]["_last_manifest"]["status"] == "success"
        assert app["S"]["_last_manifest"]["subcommand"] == subcommand
        assert results[-1][1] == "rendered"


def test_colab_results_keep_missing_energies_unknown(
    tmp_path: Path, monkeypatch,
) -> None:
    app, rendered = _execute_app(monkeypatch, tmp_path)
    missing_first = tmp_path / "missing_first_trj.xyz"
    missing_first.write_text(
        "2\nenergy unavailable\nH 0 0 0\nH 0.7 0 0\n"
        "2\nenergy=-0.990000 unit=hartree\nH 0.2 0 0\nH 0.9 0 0\n",
        encoding="utf-8",
    )
    app["S"]["_last_subcmd"] = "path-opt"
    rendered.clear()
    app["_load_trajectory"](str(missing_first), str(tmp_path))
    assert app["_TRAJ"]["energies"] == [None, -0.99]
    assert app["_rel_kcal"]() is None
    assert "profile not re-referenced" in app["frame_state"].value
    assert "Energy profile unavailable" in app["plot_out"].value

    missing_middle = tmp_path / "missing_middle_trj.xyz"
    missing_middle.write_text(
        "2\nenergy=-1.000000 unit=hartree\nH 0 0 0\nH 0.7 0 0\n"
        "2\nenergy unavailable\nH 0.2 0 0\nH 0.9 0 0\n",
        encoding="utf-8",
    )
    rendered.clear()
    app["_load_trajectory"](str(missing_middle), str(tmp_path))
    assert app["_rel_kcal"]() == [0.0, None]
    app["frame_slider"].value = 1
    assert "ΔE unavailable" in app["frame_state"].value
    assert app["plot_out"].value.count('class="rxenergy-frame"') == 1


def test_colab_operates_scientific_selectors_and_remaining_buttons(
    tmp_path: Path, monkeypatch,
) -> None:
    app, _rendered = _execute_app(monkeypatch, tmp_path)
    first = tmp_path / "first.pdb"
    second = tmp_path / "second.pdb"
    pdb_text = (
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O1  LIG A   1       1.200   0.000   0.000  1.00  0.00           O\n"
        "HETATM    3  N1  LIG A   1       0.000   1.200   0.000  1.00  0.00           N\n"
        "HETATM    4  O   HOH A   2       3.000   0.000   0.000  1.00  0.00           O\n"
        "END\n"
    )
    first.write_text(pdb_text, encoding="utf-8")
    second.write_text(pdb_text.replace("3.000", "3.100"), encoding="utf-8")
    app["load_pdb"]([str(first), str(second)])

    app["_set_running"](True)
    assert app["dd_subcmd"].disabled
    assert app["w_out"].disabled
    assert app["upl"].disabled
    app["_set_running"](False)
    assert not app["dd_subcmd"].disabled
    assert not app["w_out"].disabled
    assert not app["upl"].disabled

    app["dd_subcmd"].options = app["_sub_options"](app["SUBS"])
    app["dd_subcmd"].value = "energy-diagram"
    for control in (
        app["pick_action"], app["last_pick_info"], app["viewer_more"],
        app["last_pick_row"],
    ):
        assert control.layout.display == "none"
    app["dd_subcmd"].value = "scan"
    assert app["pick_action"].layout.display == ""
    assert app["last_pick_row"].layout.display == ""

    # Reaction-order arrows dispatch their row callbacks and preserve a usable
    # primary viewer after each reorder.
    _widget_with_description(app["input_file_rows"], "Move later").click()
    assert app["S"]["inputs"] == [str(second), str(first)]
    _widget_with_description(app["input_file_rows"], "Move earlier").click()
    assert app["S"]["inputs"] == [str(first), str(second)]
    app["load_pdb"]([str(first)])

    def pick(index: int) -> None:
        app["on_click"](str(index), live_marked=True)

    def scan_pair(first_index: int, second_index: int) -> None:
        _widget_with_description(app["scan_panel"], "1  Pick atom A").click()
        assert app["pick_action"].value == "scanA"
        pick(first_index)
        _widget_with_description(app["scan_panel"], "2  Pick atom B").click()
        assert app["pick_action"].value == "scanB"
        pick(second_index)

    # Operate 1D staged-scan controls rather than seeding their final state.
    app["dd_subcmd"].value = "scan"
    scan_pair(0, 1)
    target = _widget_with_description(app["scan_panel"], "3  Target Å")
    target.value = 1.8
    _widget_with_description(app["scan_panel"], "Add sequential stage").click()
    assert len(app["S"]["scan_stages"]) == 1
    assert app["S"]["scan_atoms"] == [None, None]   # adding a stage resets the pair
    scan_pair(0, 1)
    _widget_with_description(app["scan_panel"], "Add to current stage").click()
    assert len(app["S"]["scan_stages"][0]) == 1
    assert "already in the current stage" in app["S"]["_last_pick_message"]
    scan_pair(1, 0)
    _widget_with_description(app["scan_panel"], "Add to current stage").click()
    assert len(app["S"]["scan_stages"][0]) == 1
    scan_pair(0, 2)
    _widget_with_description(app["scan_panel"], "Add to current stage").click()
    assert len(app["S"]["scan_stages"][0]) == 2
    app["b_clear_scan"].click()
    assert app["S"]["scan_stages"] == []
    _widget_with_description(app["scan_panel"], "Clear current pair").click()
    assert app["S"]["scan_atoms"] == [None, None]

    # The shared multidimensional controls are exercised at both cardinalities.
    for subcommand, pairs in (
        ("scan2d", ((0, 1), (0, 2))),
        ("scan3d", ((0, 1), (0, 2), (1, 2))),
    ):
        app["dd_subcmd"].value = subcommand
        for pair_index, (left, right) in enumerate(pairs):
            scan_pair(left, right)
            low = _widget_with_description(app["scan_panel"], "3  Low Å")
            high = _widget_with_description(app["scan_panel"], "High Å")
            low.value = 1.0 + pair_index * 0.1
            high.value = 2.0 + pair_index * 0.1
            _widget_with_description(app["scan_panel"], "4  Add axis").click()
        assert len(app["S"]["scan_axes"]) == len(pairs)
        axis_removals = [
            widget for widget in _widget_descendants(app["scan_panel"])
            if getattr(widget, "tooltip", "") == "Remove axis"
        ]
        assert len(axis_removals) == len(pairs)
        axis_removals[-1].click()
        assert len(app["S"]["scan_axes"]) == len(pairs) - 1
        app["b_clear_scan"].click()
        assert app["S"]["scan_axes"] == []

    # Freeze-pair and frozen-atom actions go through the visible selector and
    # the actual add/clear buttons.
    app["dd_subcmd"].value = "opt"
    app["pick_action"].value = "freezeA"
    pick(0)
    app["pick_action"].value = "freezeB"
    pick(1)
    _widget_with_description(app["freeze_panel"], "add freeze pair").click()
    assert len(app["S"]["freeze_pairs"]) == 1
    _widget_with_description(app["freeze_panel"], "clear pairs").click()
    assert app["S"]["freeze_pairs"] == []
    app["pick_action"].value = "freezeatom"
    pick(2)
    assert app["S"]["freeze_atoms"] == [3]
    _widget_with_description(app["freeze_panel"], "clear atoms").click()
    assert app["S"]["freeze_atoms"] == []

    # Manual exact-atom entry and each removable selection chip dispatch.
    app["pick_action"].value = "center"
    app["exact_atom"].value = "1"
    app["exact_atom_btn"].click()
    assert app["S"]["_last_pick"]["index"] == 0
    app["center_widget"].value = ("LIG",)
    app["S"]["center_ids"] = ["A:LIG:1"]
    app["S"]["freeze_atoms"] = [3]
    app["_render_chips"]()
    for description in ("LIG ✕", "A:LIG:1 ✕", "⚓3 ✕"):
        _widget_with_description(app["chips_box"], description).click()
    assert app["center_widget"].value == ()
    assert app["S"]["center_ids"] == [] and app["S"]["freeze_atoms"] == []
    app["S"]["center_ids"] = ["A:LIG:1"]
    app["S"]["freeze_atoms"] = [2]
    app["b_clear"].click()
    assert app["S"]["center_ids"] == []
    _widget_with_description(app["freeze_panel"], "clear atoms").click()
    assert app["S"]["freeze_atoms"] == []

    # Search and every notebook-owned disclosure button make a full round trip.
    app["dd_subcmd"].options = app["_sub_options"](app["SUBS"])
    app["dd_subcmd"].value = "all"
    app["adv_search"].value = "cycles"
    assert "total" in app["adv_count"].value
    app["adv_search"].value = ""
    folds = [
        widget for widget in _widget_descendants(app["rootbox"])
        if "rxfold" in getattr(widget, "_dom_classes", ())
    ]
    assert folds
    for fold in folds:
        fold._rx_set_open(False)
        fold._rx_button.click()
        assert fold._rx_body.layout.display == ""
        fold._rx_button.click()
        assert fold._rx_body.layout.display == "none"

    # Clipboard and session controls use the same Colab imports as production.
    clipboard: list[str] = []
    downloads: list[str] = []
    google = types.ModuleType("google")
    colab = types.ModuleType("google.colab")
    output = types.ModuleType("google.colab.output")
    files = types.ModuleType("google.colab.files")
    output.eval_js = lambda script: clipboard.append(str(script))
    files.download = lambda path: downloads.append(str(path))
    colab.output = output
    colab.files = files
    google.colab = colab
    monkeypatch.setitem(sys.modules, "google", google)
    monkeypatch.setitem(sys.modules, "google.colab", colab)
    monkeypatch.setitem(sys.modules, "google.colab.output", output)
    monkeypatch.setitem(sys.modules, "google.colab.files", files)

    app["dd_subcmd"].value = "fix-altloc"
    app["b_rebuild"].click()
    app["b_copy"].click()
    assert clipboard and json.dumps(app["cmd_box"].value) in clipboard[-1]
    saved_out = app["w_out"].value
    app["b_save"].click()
    session_path = Path(downloads[-1])
    session_bytes = session_path.read_bytes()
    app["w_out"].value = "changed-after-save"
    app["up_sess"].value = ({
        "name": "session.json", "type": "application/json",
        "size": len(session_bytes), "content": memoryview(session_bytes),
        "last_modified": datetime.datetime.now(datetime.timezone.utc),
    },)
    assert app["w_out"].value == saved_out
    assert app["up_sess"].value == ()

    # Real Results discovery, both selectors, navigation, preview and ZIP.
    trajectory_a = tmp_path / "path_a_trj.xyz"
    trajectory_b = tmp_path / "path_b_trj.xyz"
    trajectory_text = (
        "2\nenergy=-1.000000 unit=hartree\nH 0 0 0\nH 0.7 0 0\n"
        "2\nenergy=-0.990000 unit=hartree\nH 0.2 0 0\nH 0.9 0 0\n"
    )
    trajectory_a.write_text(trajectory_text, encoding="utf-8")
    trajectory_b.write_text(
        trajectory_text.replace("-0.990000", "-0.980000"), encoding="utf-8",
    )
    artifact_a = tmp_path / "artifact_a.pdb"
    artifact_b = tmp_path / "artifact_b.pdb"
    artifact_a.write_text(pdb_text, encoding="utf-8")
    artifact_b.write_text(pdb_text.replace("1.200", "1.300"), encoding="utf-8")
    current_files = [trajectory_a, trajectory_b, artifact_a, artifact_b]
    app["S"].update(
        _last_out_dir=str(tmp_path), _last_subcmd="path-search",
        _last_files=[str(path) for path in current_files],
        _last_manifest={"status": "success", "exit_code": 0},
        _last_log="button coverage log",
    )
    app["res_btn"].disabled = False
    app["res_btn"].click()
    assert len(app["traj_choice"].options) == 2
    app["traj_choice"].value = str(trajectory_b)
    app["frame_next"].click()
    assert app["frame_slider"].value == 1
    app["frame_prev"].click()
    assert app["frame_slider"].value == 0
    app["frame_play"].value = 1
    app["artifact_choice"].value = str(artifact_b)
    app["artifact_fold"]._rx_set_open(False)
    app["artifact_fold"]._rx_button.click()
    assert app["artifact_fold"]._rx_body.layout.display == ""
    app["dl_btn"].click()
    archive_path = Path(downloads[-1])
    with zipfile.ZipFile(archive_path) as archive:
        members = set(archive.namelist())
    assert {"colab_run.json", "colab_run.log"} <= members
    assert {path.name for path in current_files} <= members

    # Operate every built-in example branch with local release-matched assets.
    app["_example_file"] = lambda relpath: str(NOTEBOOK.parent / relpath)
    for choice in app["ex_choice"].options:
        app["ex_choice"].value = choice
        app["ex_btn"].click()
        assert app["S"]["inputs"], choice
        assert "⚠️" not in app["example_msg"].value


def test_colab_adversarial_state_transactions_and_editor_ownership(
    tmp_path: Path, monkeypatch,
) -> None:
    app, calls = _execute_app(monkeypatch, tmp_path)
    primary = tmp_path / "primary.pdb"
    compatible = tmp_path / "product.pdb"
    incompatible = tmp_path / "incompatible.pdb"
    primary.write_text(
        "HETATM    1  C1  LIG A  10       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O1  LIG A  10       1.200   0.000   0.000  1.00  0.00           O\nEND\n",
        encoding="utf-8",
    )
    compatible.write_text(
        "HETATM    1  C1  LIG A  10      10.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O1  LIG A  10      12.400   0.000   0.000  1.00  0.00           O\nEND\n",
        encoding="utf-8",
    )
    incompatible.write_text(
        "HETATM    1  C1  ALT B  20       4.000   0.000   0.000  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    metadata = {
        str(primary): [
            {"serial": 1, "chain": "A", "resname": "LIG", "resseq": 10, "icode": "", "name": "C1"},
            {"serial": 2, "chain": "A", "resname": "LIG", "resseq": 10, "icode": "", "name": "O1"},
        ],
        str(compatible): [
            {"serial": 1, "chain": "A", "resname": "LIG", "resseq": 10, "icode": "", "name": "C1"},
            {"serial": 2, "chain": "A", "resname": "LIG", "resseq": 10, "icode": "", "name": "O1"},
        ],
        str(incompatible): [
            {"serial": 1, "chain": "B", "resname": "ALT", "resseq": 20, "icode": "", "name": "C1"},
        ],
    }
    real_loader = app["_load_view_structure"]

    def load_view(path):
        path = str(path)
        return Path(path).read_text(encoding="utf-8"), [dict(row) for row in metadata[path]], path

    app["_load_view_structure"] = load_view
    assert app["b_run"].disabled and app["b_validate"].disabled
    app["cmd_box"].value = "pdb2reaction --version"
    assert not app["b_run"].disabled
    assert app["b_validate"].disabled
    app["cmd_box"].value = "# incomplete"
    assert app["b_run"].disabled

    advanced_row = next(row for row in app["adv_rows_box"].children if hasattr(row, "_rx_search"))
    advanced_info = advanced_row.children[1]
    assert '<details class="rxinfo-details"' in advanced_info.value
    assert 'role="note"' in advanced_info.value
    app["_render_advanced_rows"]()
    assert " open" not in advanced_info.value

    app["load_pdb"]([str(primary), str(compatible)], center=["LIG"])
    assert "reaction order shown above" in app["input_order_note"].value
    app["S"]["scan_atoms"] = [
        {"index": 0, "chain": "A", "resn": "LIG", "resi": "10", "atom": "C1", "xyz": (0.0, 0.0, 0.0)},
        {"index": 1, "chain": "A", "resn": "LIG", "resi": "10", "atom": "O1", "xyz": (1.2, 0.0, 0.0)},
    ]
    app["view_input"].value = 1
    assert app["S"]["scan_atoms"][0]["xyz"] == pytest.approx((10.0, 0.0, 0.0))
    assert app["S"]["scan_atoms"][1]["xyz"] == pytest.approx((12.4, 0.0, 0.0))
    calls.clear()
    app["render_viewer"]()
    assert any(
        isinstance(value, str) and '"type":"rx-load-structure"' in value
        and '"generation":%d' % app["_VIEWER_GENERATION"]["value"] in value
        for value in calls
    )

    app["load_pdb"]([str(primary), str(incompatible)], center=["LIG"], lcharge={"LIG": 1})
    app["S"]["scan_atoms"] = [
        {"index": 0, "chain": "A", "resn": "LIG", "resi": "10", "atom": "C1", "xyz": (0.0, 0.0, 0.0)},
        {"index": 1, "chain": "A", "resn": "LIG", "resi": "10", "atom": "O1", "xyz": (1.2, 0.0, 0.0)},
    ]
    primary_owned = json.dumps(app["S"]["scan_atoms"], sort_keys=True)
    app["view_input"].value = 1
    assert json.dumps(app["S"]["scan_atoms"], sort_keys=True) == primary_owned
    before_selection = (
        list(app["S"]["center"]), dict(app["S"]["lcharge"]),
        list(app["S"]["scan_atoms"]), list(app["S"]["freeze_atoms"]),
    )
    assert app["b_clear"].disabled and app["center_widget"].disabled
    app["_clear_sel"](None)
    assert (
        app["S"]["center"], app["S"]["lcharge"],
        app["S"]["scan_atoms"], app["S"]["freeze_atoms"],
    ) == before_selection
    app["view_input"].value = 0
    assert not app["b_clear"].disabled
    app["_clear_sel"](None)
    assert app["S"]["center"] == [] and app["S"]["lcharge"] == {}

    app["S"]["center_ids"] = ["A:LIG:10"]
    app["S"]["_last_manifest"] = {"status": "success"}
    before = json.dumps(app["_session_dict"](), sort_keys=True)
    with pytest.raises(ValueError, match="one object"):
        app["_apply_session"]([])
    assert json.dumps(app["_session_dict"](), sort_keys=True) == before
    assert app["S"]["_last_manifest"] == {"status": "success"}
    wrong_tool = app["_session_dict"]()
    wrong_tool["tool"] = "mlmm-toolkit"
    with pytest.raises(ValueError, match="belong"):
        app["_apply_session"](wrong_tool)
    assert app["S"]["center_ids"] == ["A:LIG:10"]

    saved = app["_session_dict"]()
    saved.update(
        inputs=[str(primary), str(compatible)], mode="pdb", subcmd="all", all_mode="mep",
        tsopt=False, thermo=False, backend="uma", model="uma-m-1p1",
        rep="sticks", color="spectrum",
    )
    saved["file_identities"]["inputs"] = [
        app["_file_identity"](str(primary)),
        app["_file_identity"](str(compatible)),
    ]
    saved["primary_atom_signatures"] = [
        list(item) for item in app["_atom_signatures"](metadata[str(primary)])
    ]
    saved["advanced_overrides"] = {"all": {"verbose": 3}}
    app["all_mode"].value = "tsonly"
    app["_ALL_MODE_STATE"]["tsopt_before_tsonly"] = True
    stale_command = "pdb2reaction sp -i stale.pdb -q 99"
    app["cmd_box"].value = stale_command
    assert app["_auto"]["on"] is False
    with pytest.raises(ValueError, match="not installed"):
        app["_apply_session"](saved)
    saved["backend"] = app["BACKEND"]
    saved["model"] = app["DEFAULT_MODEL"][app["BACKEND"]]
    assert app["_apply_session"](saved) == []
    assert app["all_mode"].value == "mep" and app["w_ts"].value is False
    assert app["dd_backend"].value == app["BACKEND"]
    assert app["dd_model"].value == app["DEFAULT_MODEL"][app["BACKEND"]]
    assert app["dd_rep"].value == "sticks" and app["dd_col"].value == "spectrum"
    assert app["S"]["advanced_overrides"]["all"]["verbose"] == 3
    assert app["_auto"]["on"] is True and app["cmd_box"].value != stale_command
    bad_verbose = app["_session_dict"]()
    bad_verbose["advanced_overrides"] = {"all": {"verbose": "3"}}
    before_bad_verbose = json.dumps(app["_session_dict"](), sort_keys=True)
    with pytest.raises(ValueError, match="integer from 0 to 3"):
        app["_apply_session"](bad_verbose)
    assert json.dumps(app["_session_dict"](), sort_keys=True) == before_bad_verbose
    untrusted_confirmation = app["_session_dict"]()
    untrusted_confirmation["charge_explicit"] = True
    assert app["_apply_session"](untrusted_confirmation) == []
    assert app["w_charge_ok"].value is False
    assert app["S"]["charge_explicit"] is False
    assert app["S"]["_charge_scope"] is None

    dependent = app["_session_dict"]()
    dependent["tsopt"] = False
    dependent["thermo"] = True
    dependent["advanced"]["dft"] = False
    assert app["_apply_session"](dependent) == []
    assert app["w_ts"].value is True and app["w_ts"].disabled is True
    assert app["_session_dict"]()["tsopt"] is True
    app["w_th"].value = False
    assert app["w_ts"].value is True and app["w_ts"].disabled is False
    app["w_ts"].value = False

    app["all_mode"].value = "mep"
    app["w_ts"].value = False
    app["all_mode"].value = "tsonly"
    ts_only_saved = app["_session_dict"]()
    assert ts_only_saved["tsopt"] is False
    assert app["_apply_session"](ts_only_saved) == []
    app["all_mode"].value = "mep"
    assert app["w_ts"].value is False

    legacy = app["_session_dict"]()
    legacy["schema_version"] = 1
    before_legacy = json.dumps(app["_session_dict"](), sort_keys=True)
    with pytest.raises(ValueError, match="Save the session again"):
        app["_apply_session"](legacy)
    assert json.dumps(app["_session_dict"](), sort_keys=True) == before_legacy

    tampered_signature = app["_session_dict"]()
    tampered_signature["primary_atom_signatures"][0][4] = "WRONG"
    before_tamper = json.dumps(app["_session_dict"](), sort_keys=True)
    with pytest.raises(ValueError, match="atom identity/order"):
        app["_apply_session"](tampered_signature)
    assert json.dumps(app["_session_dict"](), sort_keys=True) == before_tamper

    occupied = tmp_path / "occupied-out"
    occupied.mkdir()
    (occupied / "stale.txt").write_text("stale", encoding="utf-8")
    collision_argv = [
        "pdb2reaction", "opt", "-i", str(primary), "-o", str(occupied),
    ]
    rerouted = app["_preflight_output_scope"](collision_argv)
    assert rerouted is not None
    assert os.path.abspath(rerouted["root"]) != os.path.abspath(str(occupied))
    assert (occupied / "stale.txt").read_text(encoding="utf-8") == "stale"
    assert collision_argv[collision_argv.index("-o") + 1] == rerouted["root"]

    missing = tmp_path / "missing-dir" / "missing.pdb"
    pending = app["_session_dict"]()
    pending.update(
        backend="mace", model="MACE-OMOL-0", inputs=[str(missing)], mode="pdb",
        subcmd="opt", all_mode="mep", center=[], center_ids=[], lcharge={},
    )
    pending["file_identities"]["inputs"] = [
        app["_file_identity"](str(primary)),
    ]
    pending["primary_atom_signatures"] = [
        list(item) for item in app["_atom_signatures"](metadata[str(primary)])
    ]
    assert app["_apply_session"](pending) == [str(missing)]
    uploaded = tmp_path / "missing.pdb"
    uploaded.write_text(incompatible.read_text(encoding="utf-8"), encoding="utf-8")
    metadata[str(uploaded)] = [dict(row) for row in metadata[str(incompatible)]]
    with pytest.raises(ValueError, match="SHA-256 identity"):
        app["_ingest_saved_files"]([str(uploaded)], "wrong session re-upload")
    assert app["S"]["inputs"] == [str(missing)]
    uploaded.write_text(primary.read_text(encoding="utf-8"), encoding="utf-8")
    metadata[str(uploaded)] = [dict(row) for row in metadata[str(primary)]]
    assert app["_ingest_saved_files"]([str(uploaded)], "test re-upload")
    assert app["S"]["inputs"] == [str(uploaded)]
    assert str(missing) not in app["S"]["inputs"]

    bad = tmp_path / "empty.pdb"
    bad.write_text("REMARK no atoms\nEND\n", encoding="utf-8")
    old_inputs = list(app["S"]["inputs"])
    app["_load_view_structure"] = real_loader
    assert not app["_ingest_saved_files"]([str(bad)], "test invalid")
    assert app["S"]["inputs"] == old_inputs
    assert "not attached" in app["input_msg"].value
    app["_load_view_structure"] = load_view

    duplicate = {"index": 0, "chain": "A", "resn": "LIG", "resi": "10", "atom": "C1", "xyz": (0, 0, 0)}
    app["dd_subcmd"].value = "scan"
    app["S"]["scan_atoms"] = [dict(duplicate), dict(duplicate)]
    with pytest.raises(ValueError, match="different atoms"):
        app["build_cmd"]()

    class DummyUpload:
        def __init__(self, value):
            self.value = value
            self._counter = 2

    old_upload = DummyUpload({"file": object()})
    app["_reset_file_upload"](old_upload)
    assert old_upload.value == {} and old_upload._counter == 0
    new_upload = DummyUpload((object(),))
    app["_reset_file_upload"](new_upload)
    assert new_upload.value == ()

    class BrokenStdout:
        def __iter__(self):
            return self

        def __next__(self):
            raise RuntimeError("reader failed")

        def close(self):
            pass

    class FakeProcess:
        pid = 4242

        def __init__(self):
            self.stdout = BrokenStdout()
            self.returncode = None
            self.waited = False

        def poll(self):
            return self.returncode

        def wait(self, timeout=None):
            self.waited = True
            self.returncode = -15
            return self.returncode

    process = FakeProcess()
    signals = []
    monkeypatch.setattr(app["subprocess"], "Popen", lambda *args, **kwargs: process)
    monkeypatch.setattr(app["os"], "getpgid", lambda pid: pid)
    monkeypatch.setattr(app["os"], "killpg", lambda pid, sig: signals.append((pid, sig)))
    rc, transcript = app["_stream"](["pdb2reaction", "--version"])
    assert rc == 125 and process.waited and signals
    assert "reader failed" in transcript

    vanished = FakeProcess()
    monkeypatch.setattr(
        app["os"], "killpg",
        lambda *_args: (_ for _ in ()).throw(ProcessLookupError()),
    )
    app["_stop_child"](vanished)
    assert vanished.waited


def test_colab_gui_guards_state_capabilities_and_current_run_results() -> None:
    setup = _notebook()["cells"][1]["source"]
    app = _notebook()["cells"][2]["source"]

    assert "INSTALL_DFT = install_dft" in setup
    assert "def _clear_structure_bound_state(preserve_model=None):" in app
    for key in (
        "scan_stages", "scan_axes", "freeze_buf", "freeze_pairs",
        "freeze_atoms", "measure_atoms", "_pdb_text", "_atom_meta",
    ):
        assert key in app
    assert "_INITIAL_MODEL" in app
    assert "_bk0 = BACKEND if BACKEND in MODELS else 'mace'" in app
    assert "DFT_READY" in app and "DMF_READY" in app
    assert "OUT_JSON_SUBS" in app and "if sub in OUT_JSON_SUBS: cmd += ['--out-json']" in app
    assert "class _DropUpload(anywidget.AnyWidget):" in app
    assert "upl = _DropUpload(accept=_acc, formats=_drop_formats)" in app
    assert "pdb2reaction_gui.on_drop" not in app
    assert "def _ingest_saved_files(" in app
    assert "_reset_file_upload(upl)" in app
    assert "if _UPLOAD_MODE == 'anywidget': upl.on_msg(_on_drop_upload)" in app
    assert "pts.append((k, 'TS'))" not in app
    assert "points[k] = 'TS candidate'" not in app
    assert "points[k] = 'minimum candidate'" not in app
    assert "'extrema': True" not in app
    assert "except KeyboardInterrupt:" in app
    assert "start_new_session=(os.name == 'posix')" in app
    assert "os.killpg(" in app and "signal.SIGTERM" in app and "signal.SIGKILL" in app
    assert "process.terminate()" in app and "process.kill()" in app and "process.wait(" in app
    assert "cancelled" in app
    assert "_MANUAL" in app
    assert "current_output_paths" in app and "key_output_files" in app
    assert "def _select_status_json(current, sub):" in app
    assert "if sub in ('all', 'path-search')" in app
    assert "Use machine-readable claims when present" in app
    assert "if real_run:" in app and "'exit_code': rc" in app
    assert "'status': 'success' if rc == 0" in app
    assert "Invalid command line:" in app
    assert "def _command_input_option_flags(argv):" in app
    assert "def _click_parsed_input_files(argv):" in app
    assert "isinstance(value_type, click.Path)" in app
    assert "flags.intersection({'-s', '--scan-lists'})" in app
    assert "_RUN_STATE['validation_log']" in app
    assert "tryCreateComponentFromExpression" not in app
    assert "sub == 'all' and all_kind == 'scan'" in app
    assert "bond table or JSON on stdout (--json)" in app
    assert "colab_run.json" in app and "zipfile.ZipFile" in app
    assert "shutil.make_archive(" not in app
    assert "W.Button(description='Show information'" not in app
    assert "<details class=\"rxinfo-details\" data-revision=" in app
    assert "submit(event.dataTransfer ? event.dataTransfer.files : []);" in app
    assert "_tab_body.children = [_TAB_PAGES[i][1]]" not in app
    assert "layout=W.Layout(width='560px')" not in app
    assert app.count("effective = _normalized_scope_argv(a)") == 2


def test_colab_output_scope_executes_cli_grammar_and_utility_defaults(
    tmp_path: Path, monkeypatch,
) -> None:
    contract = _output_contract()
    scope_for = contract["_output_scope"]
    monkeypatch.chdir(tmp_path)

    from pdb2reaction.cli import cli as product_cli
    import click

    root_context = click.Context(product_cli)
    compute_subcommands = {
        "all", "scan", "opt", "path-opt", "path-search", "tsopt", "freq",
        "irc", "dft", "sp", "scan2d", "scan3d",
    }
    for subcommand in compute_subcommands:
        command = product_cli.get_command(root_context, subcommand)
        output_param = next(
            param for param in command.params
            if {"-o", "--out", "--out-dir", "--output"}.intersection(param.opts)
        )
        assert scope_for(["pdb2reaction", subcommand])["root"] == str(output_param.default)
    assert scope_for(["pdb2reaction", "opt", "-oattached"])["root"] == "attached"
    assert scope_for(["pdb2reaction", "opt", "--OUT-DIR", "Upper"])["root"] == "Upper"
    assert scope_for(["pdb2reaction", "opt", "--OUT-DIR=UpperEq"])["root"] == "UpperEq"
    assert scope_for(["pdb2reaction", "-i", "input.pdb"])["root"] == "result_all"

    for info_argv in (
        ["pdb2reaction", "energy-diagram", "--help"],
        ["pdb2reaction", "all", "--help-advanced"],
        ["pdb2reaction", "--version"],
    ):
        info_scope = scope_for(info_argv)
        assert info_scope["stdout_only"] is True
        assert info_scope["targets"] == []
        assert not contract["_output_scope_collision"](info_scope)

    first = tmp_path / "plots-a" / "same.png"
    second = tmp_path / "plots-b" / "same.png"
    positional = tmp_path / "plots-c" / "report.csv"
    scope = scope_for([
        "pdb2reaction", "trj2fig", str(positional), "-i", "trajectory.xyz",
        "-o", str(first), f"--out={second}",
        "--out-json", "--no-out-json", "--out-json",
    ])
    assert scope["targets"] == [
        str(first.resolve()), str(second.resolve()), str(positional.resolve()),
    ]
    assert scope["direct_current"] is True
    assert str(first.parent / "result.json") in scope["exact_targets"]
    assert str(first.parent / "summary.json") in scope["exact_targets"]
    assert len([p for p in scope["targets"] if Path(p).name == "same.png"]) == 2

    second.parent.mkdir(parents=True)
    second.write_text("stale", encoding="utf-8")
    assert contract["_output_scope_collision"](scope)
    before = contract["_snapshot_output_scope"](scope)
    first.parent.mkdir(parents=True)
    positional.parent.mkdir(parents=True)
    first.write_text("new", encoding="utf-8")
    positional.write_text("new", encoding="utf-8")
    unrelated = first.parent / "unrelated.txt"
    unrelated.write_text("ignore", encoding="utf-8")
    after = contract["_snapshot_output_scope"](scope)
    changed = {path for path, stat in after.items() if before.get(path) != stat}
    assert changed == {str(first.resolve()), str(positional.resolve())}

    default_plot = scope_for(["pdb2reaction", "trj2fig", "-i", "trajectory.xyz"])
    assert default_plot["targets"] == [str((tmp_path / "energy.png").resolve())]
    no_json = scope_for([
        "pdb2reaction", "trj2fig", "-i", "trajectory.xyz",
        "--out-json", "--no-out-json",
    ])
    assert not any(Path(path).name == "result.json" for path in no_json["exact_targets"])

    legacy_false = scope_for([
        "pdb2reaction", "trj2fig", "-i", "trajectory.xyz",
        "--out-json", "False",
    ])
    assert legacy_false["targets"] == [str((tmp_path / "energy.png").resolve())]
    assert not any(Path(path).name == "False" for path in legacy_false["exact_targets"])
    assert not any(Path(path).name == "result.json" for path in legacy_false["exact_targets"])
    legacy_negative_false = scope_for([
        "pdb2reaction", "trj2fig", "-i", "trajectory.xyz",
        "--no-out-json", "False",
    ])
    assert legacy_negative_false["targets"] == [str((tmp_path / "energy.png").resolve())]
    assert str((tmp_path / "result.json").resolve()) in legacy_negative_false["exact_targets"]
    uppercase_inline_false = scope_for([
        "pdb2reaction", "energy-diagram", "-i", "0", "-i", "1",
        "--OUT-JSON=False",
    ])
    assert not any(Path(path).name == "result.json" for path in uppercase_inline_false["exact_targets"])

    normalize = contract["_normalized_scope_argv"]
    enabled = contract["_flag_enabled"]
    assert not enabled(normalize([
        "pdb2reaction", "all", "--dry-run", "--no-dry-run",
    ]), "--dry-run", "--no-dry-run")
    assert not enabled(normalize([
        "pdb2reaction", "all", "--dry-run", "False",
    ]), "--dry-run", "--no-dry-run")
    assert enabled(normalize([
        "pdb2reaction", "all", "--no-dry-run", "False",
    ]), "--dry-run", "--no-dry-run")
    forced = contract["_force_dry_run"]([
        "pdb2reaction", "all", "--no-dry-run", "--",
    ])
    assert forced == [
        "pdb2reaction", "all", "--no-dry-run", "--dry-run", "--",
    ]
    assert enabled(normalize(forced), "--dry-run", "--no-dry-run")

    energy = scope_for(["pdb2reaction", "energy-diagram", "-i", "0", "-i", "1"])
    assert energy["targets"] == [str((tmp_path / "energy_diagram.png").resolve())]
    energy_named = scope_for([
        "pdb2reaction", "energy-diagram", "-i", "0", "-i", "1",
        "-o", "plots/profile", "--out-json",
    ])
    assert energy_named["targets"] == [str((tmp_path / "plots/profile.png").resolve())]
    assert str((tmp_path / "plots/result.json").resolve()) in energy_named["exact_targets"]
    fake_home = tmp_path / "home"
    fake_home.mkdir()
    monkeypatch.setenv("HOME", str(fake_home))
    tilde_energy = scope_for([
        "pdb2reaction", "energy-diagram", "-i", "0", "-i", "1",
        "-o", "~/profile",
    ])
    assert tilde_energy["targets"] == [str((tmp_path / "~" / "profile.png").resolve())]
    tilde_trj = scope_for([
        "pdb2reaction", "trj2fig", "-i", "trajectory.xyz", "-o", "~/profile.png",
    ])
    assert tilde_trj["targets"] == [str((fake_home / "profile.png").resolve())]

    add_elem = scope_for([
        "pdb2reaction", "add-elem-info", "-i", "inputs/enzyme.pdb",
    ])
    assert add_elem["targets"] == [
        str((tmp_path / "inputs/enzyme_add_elem.pdb").resolve())
    ]
    add_elem_other = scope_for([
        "pdb2reaction", "add-elem-info", "-i", "inputs/enzyme.ent",
    ])
    assert add_elem_other["targets"] == [
        str((tmp_path / "inputs/enzyme.ent_add_elem.pdb").resolve())
    ]
    fixed = scope_for(["pdb2reaction", "fix-altloc", "-i", "inputs/enzyme.pdb"])
    assert fixed["targets"] == [
        str((tmp_path / "inputs/enzyme_clean.pdb").resolve())
    ]

    input_file = tmp_path / "inputs/enzyme.pdb"
    input_file.parent.mkdir(exist_ok=True)
    input_file.write_text("ATOM\n", encoding="utf-8")
    add_elem_inplace = scope_for([
        "pdb2reaction", "add-elem-info", "-i", str(input_file), "--overwrite",
    ])
    assert add_elem_inplace["targets"] == [str(input_file.resolve())]
    assert contract["_output_scope_collision"](add_elem_inplace)
    inplace = scope_for([
        "pdb2reaction", "fix-altloc", "-i", str(input_file), "--inplace",
    ])
    assert set(inplace["targets"]) == {
        str(input_file.resolve()), str(input_file.with_suffix(".pdb.bak").resolve()),
    }
    assert contract["_output_scope_collision"](inplace)

    input_dir = tmp_path / "batch"
    input_dir.mkdir()
    directory = scope_for(["pdb2reaction", "fix-altloc", "-i", str(input_dir)])
    assert directory["root"] == str(tmp_path / "batch_clean")
    assert directory["direct_current"] is False

    extract_default = scope_for([
        "pdb2reaction", "extract", "-i", "reactant.pdb", "-c", "LIG",
    ])
    assert extract_default["targets"] == [str((tmp_path / "model.pdb").resolve())]
    assert str((tmp_path / "model.cif").resolve()) not in extract_default["exact_targets"]
    extract_single_extra_outputs = scope_for([
        "pdb2reaction", "extract", "-i", "reactant.pdb", "-c", "LIG",
        "-o", "first.pdb", "ignored.pdb",
    ])
    assert extract_single_extra_outputs["targets"] == [
        str((tmp_path / "first.pdb").resolve())
    ]
    extract_multi = scope_for([
        "pdb2reaction", "extract", "-i", "reactant.cif", "product.cif",
        "-c", "LIG", "-o", "cluster.ent", "--out-json",
    ])
    assert extract_multi["targets"] == [str((tmp_path / "cluster.ent").resolve())]
    assert str((tmp_path / "cluster.cif").resolve()) in extract_multi["exact_targets"]
    duplicate_input = scope_for([
        "pdb2reaction", "extract", "-i", "same.pdb", "same.pdb", "-c", "LIG",
    ])
    assert duplicate_input["targets"] == [str((tmp_path / "model_same.pdb").resolve())]
    mixed = scope_for([
        "pdb2reaction", "extract", "-i", "normal.pdb", "wide.cif", "-c", "LIG",
        "-o", "normal.ent", "wide.ent",
    ])
    assert str((tmp_path / "normal.cif").resolve()) not in mixed["exact_targets"]
    assert str((tmp_path / "wide.cif").resolve()) in mixed["exact_targets"]
    combined_first_regular = scope_for([
        "pdb2reaction", "extract", "-i", "normal.pdb", "wide.cif", "-c", "LIG",
        "-o", "combined.ent",
    ])
    assert str((tmp_path / "combined.cif").resolve()) not in combined_first_regular["exact_targets"]

    (tmp_path / "result").mkdir()
    (tmp_path / "result" / "stale.txt").write_text("stale", encoding="utf-8")
    stdout_scope = scope_for([
        "pdb2reaction", "bond-summary", "-i", "reactant.pdb", "product.pdb",
    ])
    assert stdout_scope["stdout_only"] is True
    assert contract["_snapshot_output_scope"](stdout_scope) == {}
    assert contract["_output_scope_collision"](stdout_scope) is False


def test_colab_release_state_and_linked_results_regressions(
    tmp_path: Path, monkeypatch,
) -> None:
    app, rendered = _execute_app(monkeypatch, tmp_path)

    assert "never" in {value for _label, value in app["adv_thresh"].options}
    assert app["adv_dftfb"].placeholder == "wb97m-v/def2-tzvpd"

    reactant = tmp_path / "reactant.xyz"
    product = tmp_path / "product.xyz"
    reactant.write_text(
        "3\nenergy=-1.000000 unit=hartree\n"
        "C 0.0 0.0 0.0\nO 1.2 0.0 0.0\nH -0.5 0.8 0.0\n",
        encoding="utf-8",
    )
    product.write_text(
        "3\nenergy=-0.990000 unit=hartree\n"
        "C 0.1 0.0 0.0\nO 1.3 0.0 0.0\nH -0.4 0.8 0.0\n",
        encoding="utf-8",
    )
    app["load_pdb"]([str(reactant)], mode="small")
    app["dd_subcmd"].value = "scan"
    app["S"]["scan_atoms"] = [
        {"index": 0, "chain": "", "resn": "MOL", "resi": "1", "icode": "",
         "atom": "C1", "xyz": (0.0, 0.0, 0.0)},
        {"index": 1, "chain": "", "resn": "MOL", "resi": "1", "icode": "",
         "atom": "O2", "xyz": (1.2, 0.0, 0.0)},
    ]
    literal = app["scan_literals"]()[0]
    assert ast.literal_eval(literal) == [(1, 2, 1.6)]
    from pdb2reaction.workflows.scan_common import parse_staged_scan_request
    parsed_scan = parse_staged_scan_request(
        [literal], one_based=True, atom_meta=[],
    )
    assert parsed_scan.stages[0][0][:2] == (0, 1)
    assert app["_aspec"]({
        "index": 4, "chain": "", "resn": "LIG", "resi": "12A",
        "icode": "A", "atom": "C1",
    }) == 5

    app["S"]["scan_preset"] = "[('OLD 1 A','OLD 1 B',1.4)]"
    app["S"]["scan_atoms"] = [None, None]
    app["pick_action"].value = "scanA"
    app["on_click"]("0", "MOL", "1", "", "C1", "1", "", live_marked=True)
    assert app["S"]["scan_preset"] == ""
    assert app["scan_literals"]() == []

    center_file = tmp_path / "centers.pdb"
    center_file.write_text("MODEL        1\nENDMDL\n", encoding="utf-8")
    center_argv = [
        "pdb2reaction", "all", "-i", str(reactant), "-c", str(center_file),
    ]
    assert str(center_file.resolve()) in app["_command_input_files"](center_argv)
    first_center_fingerprint = app["_validation_fingerprint"](center_argv)
    center_file.write_text("MODEL        2\nENDMDL\n", encoding="utf-8")
    assert app["_validation_fingerprint"](center_argv) != first_center_fingerprint

    claimed_root = tmp_path / "claimed"
    claimed_root.mkdir()
    mode_path = claimed_root / "mode_0001_-100.00cm-1_trj.xyz"
    mode_path.write_text(reactant.read_text(encoding="utf-8"), encoding="utf-8")
    result_path = claimed_root / "result.json"
    result_path.write_text(json.dumps({
        "files": {"frequencies_txt": "frequencies_cm-1.txt",
                  "mode_files": [mode_path.name]},
    }), encoding="utf-8")
    assert mode_path.resolve() in {
        Path(path).resolve() for path in app["_structured_current_paths"](
            str(claimed_root), [str(result_path), str(mode_path)],
        )
    }
    old_mode = claimed_root / "old_mode_trj.xyz"
    old_mode.write_text(reactant.read_text(encoding="utf-8"), encoding="utf-8")
    new_mode = claimed_root / "new_mode_trj.xyz"
    new_mode.write_text(product.read_text(encoding="utf-8"), encoding="utf-8")
    result_path.write_text(json.dumps({
        "files": {"mode_files": [old_mode.name, new_mode.name]},
    }), encoding="utf-8")
    reused_paths = {
        Path(path).resolve() for path in app["_structured_current_paths"](
            str(claimed_root), [str(result_path), str(new_mode)],
        )
    }
    assert new_mode.resolve() in reused_paths
    assert old_mode.resolve() not in reused_paths
    assert '"mode_files": _mode_output_files' in (
        NOTEBOOK.parents[1] / "pdb2reaction" / "workflows" / "freq.py"
    ).read_text(encoding="utf-8")

    good = tmp_path / "good_trj.xyz"
    good.write_text(
        "2\nenergy=-1.000000 unit=hartree\nH 0 0 0\nH 0.7 0 0\n"
        "2\nenergy=-0.990000 unit=hartree\nH 0.2 0 0\nH 0.9 0 0\n",
        encoding="utf-8",
    )
    bad = tmp_path / "bad_trj.xyz"
    bad.write_text("not an xyz trajectory\n", encoding="utf-8")
    app["_result_pick_guard"]["active"] = True
    try:
        app["traj_choice"].options = [
            ("bad", str(bad)), ("good", str(good)),
        ]
        app["traj_choice"].value = str(bad)
    finally:
        app["_result_pick_guard"]["active"] = False
    app["_load_trajectory"](str(bad), str(tmp_path))
    assert app["trajectory_box"].layout.display == ""
    assert app["trajectory_content"].layout.display == "none"
    assert not app["traj_choice"].disabled

    rendered.clear()
    app["_load_trajectory"](str(good), str(tmp_path))
    assert app["trajectory_content"].layout.display == ""
    assert app["plot_out"].value.count('class="rxenergy-frame"') == 1
    assert "plotly_click" in html.unescape(app["plot_out"].value)
    generation = app["_TRAJ"]["generation"]
    app["_set_frame_from_browser"](generation, 1)
    assert app["frame_slider"].value == 1
    app["_set_frame_from_browser"](generation - 1, 0)
    assert app["frame_slider"].value == 1

    semantics = {
        "start": "R", "end": "P", "extrema": True,
        "bridge_ranges": [(2, 3)],
    }
    assert app["_stationary"]([0.0, 4.0, 10.0, 4.0, 0.0], semantics) == [
        (0, "R"), (4, "P"),
    ]

    leaf = tmp_path / "leaf-summary.json"
    leaf.write_text(json.dumps({
        "status": "completed", "scientific_status": "success",
        "n_modes": 6, "n_imaginary": 1,
    }), encoding="utf-8")
    leaf_html = app["_summary_html"](str(leaf))
    assert "modes: <b>6</b>" in leaf_html
    assert "images:" not in leaf_html and "None" not in leaf_html

    app["trajectory_box"].layout.display = ""
    app["artifact_fold"].layout.display = ""
    app["_begin_results_attempt"](str(tmp_path / "new-result"), "scan")
    assert app["trajectory_box"].layout.display == "none"
    assert app["artifact_fold"].layout.display == "none"
    assert "Run in progress" in app["results_empty"].value

    pdb_text = (
        "HETATM    1  C1  LIG A  10       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O1  LIG A  10       1.200   0.000   0.000  1.00  0.00           O\n"
        "END\n"
    )
    first_pdb = tmp_path / "first.pdb"
    second_pdb = tmp_path / "second.pdb"
    first_pdb.write_text(pdb_text, encoding="utf-8")
    second_pdb.write_text(pdb_text, encoding="utf-8")
    app["load_pdb"]([str(first_pdb), str(second_pdb)], center=["LIG"])
    app["charge_rows"]["LIG"]["use"].value = True
    app["charge_rows"]["LIG"]["val"].value = 1
    app["w_charge_ok"].value = True
    app["view_input"].value = 1
    app["view_input"].value = 0
    assert app["w_charge_ok"].value is True
    assert app["S"]["lcharge"] == {"LIG": 1.0}

    app["all_mode"].value = "tsonly"
    assert app["_ALL_MODE_STATE"]["user"] is True
    app["_queue_change"](clear=True)
    app["load_pdb"]([str(first_pdb)])
    assert app["all_mode"].value == "scan"
    app["load_pdb"]([str(second_pdb)], append=True)
    assert app["all_mode"].value == "mep"


def test_colab_document_iframes_survive_srcdoc_stripping(monkeypatch, tmp_path) -> None:
    """Colab strips `srcdoc`, so each document iframe also ships a blob delivery."""
    app, _rendered = _execute_app(monkeypatch, tmp_path)

    frame = app["_document_iframe"](
        "<b>x</b><script>void 0;</script>",
        'class="rxprobe" data-rx-channel="viewer"', "width:100%;",
    )
    assert 'srcdoc="' in frame                        # JupyterLab honours srcdoc
    assert 'class="rxprobe"' in frame and 'data-rx-channel="viewer"' in frame
    assert "URL.createObjectURL(new Blob([" in frame  # Colab needs a blob src instead
    assert 'frame.getAttribute("srcdoc")' in frame    # ... and only when srcdoc is gone
    script = frame[frame.index("<script>(function()"):]
    assert script.count("</script>") == 1             # the payload cannot close the script

    viewer = app["_molstar_iframe"]("HETATM\n", "pdb", interactive=True, generation=3)
    assert 'srcdoc="' in viewer and "URL.createObjectURL" in viewer
    assert 'data-rx-generation="3"' in viewer and 'class="rxmolstar-frame"' in viewer

    energy = app["_energy_plot_iframe"]([0.0, 1.0], {}, 4)
    assert 'srcdoc="' in energy and "URL.createObjectURL" in energy
    assert 'data-rx-channel="trajectory"' in energy and 'data-rx-generation="4"' in energy

    pdf = app["_binary_iframe"]("application/pdf", "AAAA", 'title="PDF result"', "width:100%;")
    assert 'src="data:application/pdf;base64,AAAA"' in pdf
    assert "atob(" in pdf and 'frame.getAttribute("src")' in pdf


def test_colab_uma_login_accepts_a_colab_secret(monkeypatch) -> None:
    """A stored Colab secret signs in without showing the token prompt."""
    import importlib.metadata

    setup = _notebook()["cells"][1]["source"].replace(
        'backend = "mace"', 'backend = "uma"', 1,
    ).replace("install_dft = True", "install_dft = False", 1)
    logins: list[tuple] = []
    fake_hf = types.ModuleType("huggingface_hub")
    fake_hf.login = lambda **kwargs: logins.append(("token", kwargs))
    fake_hf.notebook_login = lambda: logins.append(("notebook", {}))
    fake_hf.get_token = lambda: None
    userdata = types.ModuleType("google.colab.userdata")
    userdata.get = lambda name: "secret-token" if name == "HF_TOKEN" else None
    colab = types.ModuleType("google.colab")
    colab.userdata = userdata
    google = types.ModuleType("google")
    google.colab = colab
    monkeypatch.setitem(sys.modules, "huggingface_hub", fake_hf)
    monkeypatch.setitem(sys.modules, "google", google)
    monkeypatch.setitem(sys.modules, "google.colab", colab)
    monkeypatch.setitem(sys.modules, "google.colab.userdata", userdata)
    monkeypatch.setattr(
        subprocess, "run",
        lambda argv, **_kwargs: types.SimpleNamespace(stdout="GPU 0", stderr="", returncode=0),
    )
    monkeypatch.setattr(
        importlib.metadata, "version",
        lambda name: "0.4.12" if name == "pdb2reaction" else "test",
    )
    monkeypatch.delenv("HF_TOKEN", raising=False)

    exec(compile(setup, str(NOTEBOOK), "exec"), {})

    assert [entry[0] for entry in logins] == ["token"]
    assert logins[0][1]["token"] == "secret-token"


def test_colab_setup_cell_is_frozen() -> None:
    """The Setup cell is frozen for this release (user decision, 2026-07-25).

    Behaviour contracts live in the tests above; this digest additionally freezes
    everything else in the cell, including its printed output. Update the digest
    only together with a deliberate decision to change Setup.
    """
    setup = _notebook()["cells"][1]["source"]
    digest = hashlib.sha256(setup.encode("utf-8")).hexdigest()

    assert digest == "449a6fe5554cc0f5f3bf590816cfe66da6104aba3c9fe185f1b922256abf8e1a", (
        "the Colab Setup cell changed; it is frozen for this release. Re-read the "
        "Setup contracts above, then update this digest deliberately. Got: " + digest
    )


def test_colab_navigation_names_resolve_to_real_targets() -> None:
    """Every string that sends the user somewhere must name a place that exists.

    p2r renamed tab 2 "Viewer" -> "Setup" and cell 1 "Setup" -> "Installation",
    but eleven strings kept the old names: six told a blocked user to pick atoms
    "in Viewer" (no such tab), and five told them to "rerun Setup", which had
    become the name of a tab that cannot install anything. Nothing caught it
    because each string is individually valid English.
    """
    notebook = _notebook()
    setup = notebook["cells"][1]["source"]
    app = notebook["cells"][2]["source"]

    tab_labels = {
        label
        for label in ("1 Input", "2 Input", "2 Viewer", "2 Setup", "3 Options", "4 Results")
        if f"('{label}'" in app
    }
    assert "2 Setup" in tab_labels and "2 Viewer" not in tab_labels

    # Prose that points at a tab must use a label the tab strip actually has.
    for phrase in ("in Viewer", "in Options", "in Results"):
        target = f"2 {phrase.split()[1]}" if phrase != "in Options" else "3 Options"
        if phrase in app:
            assert target in tab_labels or phrase == "in Options", (
                f"{phrase!r} names a tab that does not exist; tabs are {sorted(tab_labels)}"
            )
    assert "in Viewer" not in app

    # The install cell is titled "Installation"; prose must send users there by
    # that name, never by "Setup", which is now tab 2.
    assert "Installation —" in setup or "Installation —" in setup
    for stale in ("rerun Setup", "restart Setup"):
        assert stale not in app, f"{stale!r} names tab 2, which cannot install anything"
