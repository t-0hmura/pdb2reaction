"""Static contracts for the release-matched Colab GUI notebook."""

from __future__ import annotations

import ast
import csv
import glob
import html
import json
import os
import shlex
import sys
import types
from pathlib import Path

import pytest


NOTEBOOK = Path(__file__).parents[1] / "examples" / "pdb2reaction_colab.ipynb"


class _RecordingViewer:
    def __init__(self, calls: list, width=None, height=None):
        self.calls = calls
        calls.append(("view", width, height))

    def _record(self, name, *args, **kwargs):
        self.calls.append((name, args, kwargs))
        return {"shape": name}

    def __getattr__(self, name):
        return lambda *args, **kwargs: self._record(name, *args, **kwargs)


def _execute_app(monkeypatch, tmp_path: Path) -> tuple[dict, list]:
    """Execute the complete app cell with real widgets and a recording viewer."""
    calls: list = []
    fake = types.ModuleType("py3Dmol")
    fake.VDW = "VDW"
    fake.view = lambda width=None, height=None: _RecordingViewer(calls, width, height)
    monkeypatch.setitem(sys.modules, "py3Dmol", fake)
    import IPython.display as ipd
    monkeypatch.setattr(ipd, "display", lambda *args, **kwargs: None)
    monkeypatch.setattr(ipd, "clear_output", lambda *args, **kwargs: None)
    monkeypatch.setattr(ipd, "HTML", lambda value: value)
    monkeypatch.setattr(ipd, "Image", lambda *args, **kwargs: (args, kwargs))
    monkeypatch.chdir(tmp_path)
    namespace = {"TOOL": "pdb2reaction", "BACKEND": "mace", "REPO_DIR": "unused"}
    source = _notebook()["cells"][2]["source"]
    exec(compile(source, str(NOTEBOOK), "exec"), namespace)
    return namespace, calls


def _notebook() -> dict:
    return json.loads(NOTEBOOK.read_text(encoding="utf-8"))


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
        "_pick_residue_indices",
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
        "_trajectory_semantics",
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
            ".pdb": "structure", ".ent": "structure", ".cif": "structure",
            ".mmcif": "structure",
        },
    }
    exec(compile(module, str(NOTEBOOK), "exec"), namespace)
    assert wanted <= namespace.keys()
    return namespace


def test_colab_notebook_has_valid_code_cells_and_gpu_metadata() -> None:
    notebook = _notebook()

    assert notebook["nbformat"] == 4
    assert notebook["metadata"]["accelerator"] == "GPU"
    assert len(notebook["cells"]) == 3
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            ast.parse(cell["source"])


def test_colab_setup_is_pinned_to_matching_release_and_one_backend() -> None:
    setup = _notebook()["cells"][1]["source"]

    assert 'pdb2reaction_ref = "v0.4.12"' in setup
    assert "first run ~5 min" in setup
    assert "first run ~5-10 min" not in setup
    # The release notebook installs the pinned wheel from PyPI, the same way a
    # normal user does, so the version guard compares what pip actually resolved
    # against the requested tag.
    assert "pip('pdb2reaction' + ('[dft]' if install_dft else '') + '==' + pdb2reaction_ref.lstrip('v'))" in setup
    assert "git clone" not in setup
    assert "v0.4.4" not in setup
    assert "Restart the Colab runtime first" in setup
    assert "mace-torch>=0.3.8" in setup
    assert "HF_TOKEN" in setup
    assert "installed_version != pdb2reaction_ref[1:]" in setup
    assert "install_dft = False" in setup
    assert "INSTALL_DFT = install_dft" in setup
    assert "[%%d/5] %%s".replace("%%", "%") in setup
    assert "time.monotonic()" in setup
    assert ".pysisyphusrc" in setup
    assert setup.index(".pysisyphusrc") < setup.index("pip('pdb2reaction'")
    assert "nglview" not in setup


def test_colab_gui_tracks_current_structure_and_execution_contracts() -> None:
    app = _notebook()["cells"][2]["source"]

    assert isinstance(app, str)
    assert ".pdb,.ent,.cif,.mmcif,.xyz,.gjf" in app
    assert "prepare_input_structure" in app
    assert "load_pdb_atom_metadata" in app
    assert "CHAIN:RESNAME:RESSEQ" in app
    assert "String(atom.index)" in app
    # One occurrence registers the Python callback and one invokes it from
    # py3Dmol JavaScript; keeping both names identical is the round-trip gate.
    assert app.count("pdb2reaction_gui.on_click") == 2
    assert app.count("pdb2reaction_gui.on_view_change") == 2
    assert "demo.on_click" not in app
    assert "return shlex.split(line)" in app
    assert "a = _force_dry_run(a)" in app
    assert "real_run = not _flag_enabled(effective, '--dry-run', '--no-dry-run')" in app
    assert "b_run = W.Button(description='Run'" in app
    assert "RUN exact command" not in app
    assert "reuse non-empty out dir" in app
    assert "uma-s-1p2" in app
    assert "MODELS = {'mace': ['MACE-OMOL-0']" in app
    assert "options=['auto', 'fp32', 'fp64']" in app
    assert "Appended to the command" not in app
    assert "Every remaining option" not in app
    assert "py3Dmol.view(width=int(S.get('viewer_width', 720)), height=int(S.get('viewer_height', 540)))" in app
    assert "rxworkspace" in app
    assert "rxviewer" in app
    assert "rxinspector" in app
    assert ":focus-visible" in app
    assert ".rxapp .widget-button .fa, .rxapp .widget-upload .fa { display:none !important; }" in app
    assert "overflow-wrap:anywhere" in app
    assert "max_width='100%'" in app
    assert "flex:0 0 auto; min-width:0" in app
    assert "layout=W.Layout(width='300px', max_width='100%')" in app
    assert "view_controls = W.HBox([dd_rep, dd_col, cb_water, dd_size]," in app
    assert "W.HBox([btn_reset, btn_zoomsel, btn_zoompick, btn_clear_pick]," in app
    assert "sel_lang" not in app
    assert "No structure loaded" in app
    # Colab renders ipywidgets' Tab as an empty block, so the tab strip is Buttons
    # + a swapping VBox. Pin the four pages, the navigation entry point, and guard
    # against a Tab regression.
    assert "_TAB_PAGES = [('1 Input', input_box), ('2 Viewer', viewer_box)," in app
    assert "('3 Options', options_box), ('4 Results', results_box)]" in app
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
    # -r/--radius is an extraction option: only all and extract accept it.
    assert "if sub in ('all', 'extract') and r and r > 0: cmd += ['-r', str(r)]" in app
    assert "adv_radius.disabled = sub not in FLAG_SUBS['adv_radius']" in app
    assert "_set_flag_visible(adv_radius, sub in FLAG_SUBS['adv_radius'])" in app
    assert "_set_flag_visible(adv_dftfb, _dftfb_applicable)" in app
    assert "req='1 file + scan-lists" in app
    assert "No results yet." in app
    assert "frame_slider.disabled = False" in app
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
    # Help uses a real click-to-toggle button rather than Colab's unreliable
    # native/hover tooltip path.
    assert "def _info_control(tip, target=None):" in app
    assert "def _set_info_text(control, tip):" in app
    assert "def _close_info_target(target):" in app
    assert "def _hdr(content, tip):" in app
    assert "def _flag_row(widget, tip, info_target=None):" in app
    assert "W.Button(description='Show information'" in app
    assert "body.layout.display = '' if state['open'] else 'none'" in app
    assert "_OPEN_INFO = {'control': None}" in app
    assert ".rxapp .rxinfo-popover { position:fixed !important" in app
    assert "width:min(360px,calc(100vw - 24px))" in app
    assert ".rxinfo-popover .widget-html-content small { color:#f8fafc !important; }" in app
    assert "📥 needs" not in app
    assert "run_log_fold = _collapsible('Run log', logbox)" in app
    assert "command_editor = _collapsible('Command line', W.VBox([cmd_box" in app
    assert "W.HBox([view_input, pick_action, last_pick_info])" in app
    assert "workflow_contract_row = W.HBox([subreq, outputs_html]" in app
    assert "workflow_controls.add_class('rxworkflow-controls')" in app
    assert "grid-template-columns:minmax(280px,320px) minmax(150px,1fr) 240px" in app
    assert ".rxviewer { flex:2 1 520px !important; min-width:300px; position:sticky" in app
    assert ".rxviewer { position:static; }" in app
    assert "role=\"tooltip\"" not in app
    # Standalone opt/tsopt/path-opt users need a cluster model first.
    assert "b_extract = W.Button(description='Extract cluster model'" in app
    assert "prep_radius = W.FloatText(" in app
    assert "keep_subcmd=True" in app
    # The wheel ships no examples, so Load example resolves them from the git
    # tag matching the installed release (a source checkout is used when present).
    assert "def _example_file(relpath):" in app
    assert "raw.githubusercontent.com/t-0hmura/pdb2reaction" in app
    assert "Run Setup first" not in app
    assert "cmd = [CLI, 'extract', '-i', *S['inputs'], '-o', *outs, '-c', ','.join(cen)]" in app
    assert "cen = _center_cli_selectors()" in app
    assert "utility_bar = W.HBox([b_clear])" in app
    assert "ps.get('mlip')" in app
    assert "ps.get('gibbs_mlip')" in app
    assert "ps.get('uma')" not in app
    assert "ps.get('gibbs_uma')" not in app
    assert "def _effective_out_dir(argv=None):" in app
    assert "flags = ('-o', '--out', '--out-dir', '--output')" in app
    assert "def _trj2fig_output_targets(argv):" in app
    assert "command.make_parser(sub_ctx).parse_args" in app
    assert "def _snapshot_output_scope(scope):" in app
    assert "def _output_scope_collision(scope):" in app
    assert "scope = _output_scope(a)" in app
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
        "'_last_pick': None", "def _pick_residue_indices",
        "def _draw_last_pick", "v.addSphere(",
        "'wireframe': True", "viewer.__rxPickHalo=viewer.addSphere(",
        "'opacity': 0.12 if S.get('_last_pick') else 0.35",
        "viewer.getView()", "v.setView(list(S['_viewer_view']))",
        "zoom to click", "clear last click", "aria-live=\"polite\"",
        "water_sel = {'resn': sorted(_WATER)}",
        "visible_atoms = {} if S['show_water'] else {'not': water_sel}",
        "exact atom", "set current pick", "view_input", "_view_mapping_ok",
        "last_pick_info", "Generated file preview", "Download current run (.zip)",
        "results_box.add_class('rxresults')", "overflow-x:auto",
        "colab_run.log",
        "energy unavailable", "Command was cancelled", "Command failed",
        "_frame_link = W.jslink", "linked structure + energy",
        "ax.axvline(xs[i]", "artifact_fold._rx_set_open",
        "if center: S['_zoomsel'] = True", "if indices: v.zoom(0.78)",
    ):
        assert marker in app
    assert "('✓ ' if _j == i else '')" not in app
    draw = app[app.index("def _draw_last_pick"):app.index("def render_viewer")]
    halo = draw.index("v.addSphere({'center': center")
    assert draw.index("v.addStyle({'index': indices}") < halo
    assert "'colorscheme': 'default'" in draw
    assert "'radius': 0.92" in draw
    assert halo < draw.index("v.addLabel(")
    assert "v.addBox(" not in app
    assert "_pick_box_edges" not in app
    assert "addSurface" not in draw
    assert "v.addSurface(py3Dmol.VDW, {'color': '#f59e0b'" not in app
    assert "pick_action.value = 'scanB'" in app
    assert "pick_action.value = 'freezeB'" in app
    assert "if(atom.icode)sel.icode=atom.icode" in app
    assert "color:'#111827',opacity:0.90,wireframe:true" in app
    assert "String(atom.serial),(atom.icode||'').trim()" in app
    assert "def _resolve_click_meta(" in app
    assert "'viewer_index': viewer_index" in app
    assert "v.addModel(text, 'pdb', {'keepH': True, 'altLoc': '*'})" in app
    assert "v.setViewChangeCallback(_VIEW_CHANGE_JS)" in app
    assert "py3Dmol.view(width=int(S.get('viewer_width', 720)), height=int(S.get('viewer_height', 540)))" in app
    assert "view = py3Dmol.view(width=720, height=540)" in app
    assert "view.setStyle({}, {'cartoon': {'opacity': 0.55, 'colorscheme': 'default'}})" in app
    assert "'greenCarbon'" not in app
    assert "'line': {'opacity': 0.35" not in app
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
    contract["S"].update(_atom_meta=metadata, _last_pick={"index": 0})
    assert contract["_pick_residue_indices"]() == [0, 1]
    serial_metadata = [dict(row, serial=100 + i) for i, row in enumerate(metadata)]
    contract["S"]["_atom_meta"] = serial_metadata
    # A wrong WebGL index must not beat the stable PDB serial/signature.
    assert contract["_resolve_click_meta"](
        0, "101", "LIG", "10", "A", "C1", "",
    )["index"] == 1
    contract["S"]["_atom_meta"] = metadata
    assert contract["_resolve_atom_query"]("A:LIG:10:C1") == 1
    assert contract["_resolve_atom_query"]("2") == 1
    contract["S"].update(center=["LIG"], center_ids=["B:LIG:10"],
                         _primary_atom_meta=metadata)
    assert contract["_center_cli_selectors"]() == ["A:LIG:10", "B:LIG:10", "A:LIG:10A"]
    contract["S"].update(center=[], center_ids=["A:LIG:10"])
    with pytest.raises(ValueError, match="insertion-code sibling"):
        contract["_center_cli_selectors"]()
    contract["SPEC"].update({"bond-summary": {"n_in": (2, None)},
                              "trj2fig": {"n_in": (1, 1)}})
    assert contract["_input_count_error"]("bond-summary", 1)
    assert not contract["_input_count_error"]("bond-summary", 2)
    assert contract["_input_count_error"]("trj2fig", 2)
    assert contract["_artifact_kind"]("profile.html") == "interactive HTML"
    assert contract["_artifact_kind"]("plot.jpeg") == "image"
    assert contract["_artifact_kind"]("final.pdb") == "structure"
    assert contract["_artifact_kind"]("final.cif") == "structure"
    assert contract["_artifact_kind"]("final.xyz") == "structure"
    assert contract["_artifact_kind"]("optimization_trj.xyz") is None
    opt = contract["_trajectory_semantics"]("opt", "optimization_trj.xyz")
    path = contract["_trajectory_semantics"]("path-opt", "mep_trj.xyz")
    vibration = contract["_trajectory_semantics"]("tsopt", "vib/imag_120i_trj.xyz")
    assert vibration["title"] == "Vibrational-mode animation"
    assert vibration["x"] == "phase frame" and not vibration["extrema"]
    assert contract["_stationary"]([0.0, 2.0, 0.0], opt) == [
        (0, "initial"), (2, "optimized"),
    ]
    assert (1, "peak candidate") in contract["_stationary"]([0.0, 2.0, 0.0], path)


def test_colab_app_executes_atomic_view_and_result_transitions(
    tmp_path: Path, monkeypatch,
) -> None:
    app, calls = _execute_app(monkeypatch, tmp_path)
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
    indexed_styles = [
        call for call in calls
        if call[0] == "addStyle" and call[1] and isinstance(call[1][0], dict)
        and "index" in call[1][0]
    ]
    assert indexed_styles == []
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
    app["S"]["_last_log"] = "old log"
    app["_stream"] = lambda argv: (0, "validation transcript")
    app["_do_validate"](None)
    assert app["S"]["_last_log"] == "old log"
    assert app["_RUN_STATE"]["validation_log"] == "validation transcript"
    app["_stream"] = lambda argv: (2, "invalid options")
    app["_do_validate"](None)
    assert app["run_log_fold"].children[1].layout.display == ""
    app["_invalidate_last_run"]("Input identity changed; run again.")
    assert app["S"]["_last_manifest"] == {} and app["S"]["_last_files"] == []
    assert app["artifact_choice"].disabled and app["dl_btn"].disabled
    assert "Input identity changed" in app["results_empty"].value


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

    # Information starts closed and opens only on an explicit click.
    info = app["last_pick_info"]
    assert info._rx_info_body.layout.display == "none"
    info._rx_info_button.click()
    assert info._rx_info_body.layout.display == ""
    second_info = app["_info_control"]("second")
    second_info._rx_info_button.click()
    assert info._rx_info_body.layout.display == "none"
    assert info._rx_info_button.description == "Show information"
    assert second_info._rx_info_body.layout.display == ""
    second_info._rx_info_button.click()
    info._rx_info_button.click()
    assert info._rx_info_body.layout.display == ""
    info._rx_info_button.click()
    assert info._rx_info_body.layout.display == "none"
    assert app["key_opts_box"].children[1].layout.display == "none"
    assert app["command_editor"].children[1].layout.display == "none"
    assert app["run_log_fold"].children[1].layout.display == "none"

    # The name picker contains only chemically useful hetero residues; waters
    # and canonical amino acids remain available through exact 3D clicks.
    assert set(app["_center_values"]()) == {"LIG", "MG"}
    mg_row = app["charge_rows"]["MG"]
    assert mg_row["auto"] and mg_row["use"].disabled and mg_row["val"].disabled
    assert mg_row["val"].value == 2
    assert "MG" not in app["S"]["lcharge"]

    app["cb_water"].value = True
    for representation in ("cartoon", "sticks", "line"):
        calls.clear()
        app["dd_rep"].value = representation
        app["render_viewer"]()
        assert any(
            call[0] == "addModel" and call[1][2] == {"keepH": True, "altLoc": "*"}
            for call in calls
        )
        assert any(
            call[0] == "addStyle"
            and call[1][0].get("elem") == "O"
            and call[1][1].get("sphere", {}).get("radius") == 0.50
            for call in calls
        )
        if representation == "cartoon":
            assert not any(
                call[0] == "addStyle" and "line" in call[1][1]
                for call in calls
            )

    # Serial/signature mapping wins over a drifted WebGL index. The persistent
    # render thickens the exact residue, then marks the clicked atom without
    # replacing its element colour.
    calls.clear()
    app["pick_action"].value = "center"
    app["on_click"]("1", "LIG", "3", "A", "C1", "5", "")
    assert app["S"]["_last_pick"]["index"] == 4
    assert app["S"]["_last_pick"]["viewer_index"] == 1
    assert app["_pick_residue_indices"]() == [4, 5]
    residue_style = next(
        i for i, call in enumerate(calls)
        if call[0] == "addStyle" and call[1][0] == {"index": [4, 5]}
    )
    halo = next(
        i for i, call in enumerate(calls)
        if i > residue_style and call[0] == "addSphere" and call[1][0].get("radius") == 0.92
    )
    label = next(i for i, call in enumerate(calls) if i > halo and call[0] == "addLabel")
    assert residue_style < halo < label
    assert calls[residue_style][1][1]["stick"]["colorscheme"] == "default"
    assert not any(
        call[0] == "addStyle" and call[1][0] == {"index": [4]}
        and "sphere" in call[1][1]
        for call in calls
    )
    assert calls[halo][1][0]["center"] == {"x": 4.0, "y": 0.0, "z": 0.0}

    # A row-local × removes only that file and keeps primary-bound selections.
    close_secondary = app["input_file_rows"].children[1].children[1].children[-1]
    close_secondary.click()
    assert app["S"]["inputs"] == [str(primary)]
    assert app["S"]["center_ids"] == ["A:LIG:3"]

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
    assert app["_advanced_coverage"]("dft")
    assert app["_advanced_coverage"]("sp")["hessian_calc_mode"] == "rendered"

    app["cb_advsub"].value = True
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

    app["dd_subcmd"].value = "sp"
    assert app["key_opts_box"].layout.display == ""
    assert [value for _label, value in app["pick_action"].options] == [
        "center", "ligand", "measure",
    ]
    assert app["center_panel"].layout.display == ""
    assert app["charge_panel"].layout.display == ""
    assert app["extract_panel"].layout.display == ""
    assert app["adv_radius"]._rx_flag_row.layout.display == "none"
    assert app["adv_dftfb"]._rx_flag_row.layout.display == "none"
    app["dd_subcmd"].value = "freq"
    assert app["key_opts_box"].layout.display == ""
    app["cb_advsub"].value = True
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
    app["w_ts"].value = False
    app["w_th"].value = False
    app["adv_dft"].value = False
    off_argv = app["_advanced_argv"]("all")
    for flag in (
        "--irc-step-size", "--opt-mode-post", "--thresh-post",
        "--hessian-calc-mode", "--no-reject-uphill",
    ):
        assert flag not in off_argv
    app["all_mode"].value = "tsonly"
    assert app["adv_refine"].layout.display == "none"
    assert app["w_ts"].value and app["w_ts"].disabled
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
    thermo_argv = app["_advanced_argv"]("all")
    assert "--hessian-calc-mode" in thermo_argv
    for flag in (
        "--irc-step-size", "--opt-mode-post", "--thresh-post",
        "--no-reject-uphill",
    ):
        assert flag not in thermo_argv
    app["w_th"].value = False
    app["S"]["inputs"] = [str(primary), str(secondary)]
    assert "--tsopt" not in app["build_cmd"]()

    result_json = tmp_path / "result.json"
    result_json.write_text('{"energy": -1.25}', encoding="utf-8")
    assert app["_artifact_kind"](str(result_json)) == "JSON"
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
    assert "Provisional barrier" in summary_html
    assert "IRC endpoint mismatch" in summary_html
    assert "refined MLIP" in summary_html and "⚡" not in summary_html
    calls.clear()
    app["_structure_preview"](str(primary))
    assert any(call[0] == "addModel" and call[1][1] == "pdb" for call in calls)
    assert any(call[0] == "show" for call in calls)

    # The advertised fast example uses the same removable queue transaction as uploads.
    app["ex_choice"].value = "HCN -> HNC - small molecule (fast)"
    app["ex_btn"].click()
    assert len(app["S"]["inputs"]) == 2
    assert len(app["input_file_rows"].children) == 2
    assert app["all_mode"].value == "mep" and not app["b_clear_inputs"].disabled
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
    assert extract_commands and extract_commands[-1][-2:] == ["-r", "4.2"]
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

    app["S"].update(mode="small", inputs=[str(primary)])
    app["_sync_capability_controls"]()
    assert app["workspace"].layout.display == "none"
    assert app["selection_help"].layout.display == "none"
    assert app["selection_route"].layout.display == ""
    assert "No residue selection is needed" in app["selection_route"].value
    for panel in ("center_panel", "charge_panel", "scan_panel", "extract_panel", "freeze_acc"):
        assert app[panel].layout.display == "none"

    app["S"]["_pre_extract"] = {"inputs": [str(primary)]}
    app["b_revert"].layout.display = ""
    app["_clear_structure_bound_state"]()
    assert app["S"]["_pre_extract"] is None
    assert app["b_revert"].layout.display == "none"


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
    advanced_info._rx_info_button.click()
    assert app["advanced_help"].layout.display == ""
    assert app["_OPEN_INFO"]["control"] is advanced_info
    app["_render_advanced_rows"]()
    assert app["advanced_help"].layout.display == "none"
    assert app["_OPEN_INFO"]["control"] is None
    assert advanced_info._rx_info_button.description == "Show information"

    app["load_pdb"]([str(primary), str(compatible)], center=["LIG"])
    assert "reaction order shown above" in app["input_order_note"].value
    app["S"]["scan_atoms"] = [
        {"index": 0, "chain": "A", "resn": "LIG", "resi": "10", "atom": "C1", "xyz": (0.0, 0.0, 0.0)},
        {"index": 1, "chain": "A", "resn": "LIG", "resi": "10", "atom": "O1", "xyz": (1.2, 0.0, 0.0)},
    ]
    camera = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 25.0]
    app["S"]["_viewer_view"] = list(camera)
    app["view_input"].value = 1
    assert app["S"]["scan_atoms"][0]["xyz"] == pytest.approx((10.0, 0.0, 0.0))
    assert app["S"]["scan_atoms"][1]["xyz"] == pytest.approx((12.4, 0.0, 0.0))
    assert app["S"]["_viewer_view"] == camera
    calls.clear()
    app["render_viewer"]()
    assert any(
        call[0] == "addSphere" and call[1][0].get("color") == "red"
        and call[1][0].get("wireframe") is True
        for call in calls
    )
    assert not any(
        call[0] == "addStyle"
        and call[1][1].get("sphere", {}).get("color") in {"red", "cyan", "blue"}
        for call in calls
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
        tsopt=False, backend="uma", model="uma-m-1p1", rep="sticks", color="spectrum",
    )
    app["all_mode"].value = "tsonly"
    app["_ALL_MODE_STATE"]["tsopt_before_tsonly"] = True
    stale_command = "pdb2reaction sp -i stale.pdb -q 99"
    app["cmd_box"].value = stale_command
    assert app["_auto"]["on"] is False
    assert app["_apply_session"](saved) == []
    assert app["all_mode"].value == "mep" and app["w_ts"].value is False
    assert app["dd_backend"].value == "uma" and app["dd_model"].value == "uma-m-1p1"
    assert app["dd_rep"].value == "sticks" and app["dd_col"].value == "spectrum"
    assert app["_auto"]["on"] is True and app["cmd_box"].value != stale_command

    app["all_mode"].value = "mep"
    app["w_ts"].value = False
    app["all_mode"].value = "tsonly"
    ts_only_saved = app["_session_dict"]()
    assert ts_only_saved["tsopt"] is False
    assert app["_apply_session"](ts_only_saved) == []
    app["all_mode"].value = "mep"
    assert app["w_ts"].value is False

    app["w_reuse"].value = True
    app["_invalidate_last_run"]("new system")
    assert app["w_reuse"].value is False

    missing = tmp_path / "missing-dir" / "missing.pdb"
    pending = app["_session_dict"]()
    pending.update(
        backend="mace", model="MACE-OMOL-0", inputs=[str(missing)], mode="pdb",
        subcmd="opt", all_mode="mep", center=[], center_ids=[], lcharge={},
    )
    assert app["_apply_session"](pending) == [str(missing)]
    uploaded = tmp_path / "missing.pdb"
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


def test_colab_gui_guards_state_capabilities_and_current_run_results() -> None:
    setup = _notebook()["cells"][1]["source"]
    app = _notebook()["cells"][2]["source"]

    assert "INSTALL_DFT = install_dft" in setup
    assert "def _clear_structure_bound_state():" in app
    for key in (
        "scan_stages", "scan_axes", "freeze_buf", "freeze_pairs",
        "freeze_atoms", "measure_atoms", "_pdb_text", "_atom_meta",
    ):
        assert key in app
    assert "_INITIAL_MODEL" in app
    assert "_bk0 = BACKEND if BACKEND in MODELS else 'mace'" in app
    assert "DFT_READY" in app and "DMF_READY" in app
    assert "OUT_JSON_SUBS" in app and "if sub in OUT_JSON_SUBS: cmd += ['--out-json']" in app
    assert "out = new Array(fl.length)" in app
    assert "out[index]" in app
    assert "out.push(" not in app
    assert "r.onerror = function()" in app
    assert "pts.append((k, 'TS'))" not in app
    assert "peak candidate" in app and "minimum candidate" in app
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
    assert "mapping_ok = bool(S.get('_view_mapping_ok', True))" in app
    assert "sub == 'all' and all_kind == 'scan'" in app
    assert "bond table or JSON on stdout (--json)" in app
    assert "colab_run.json" in app and "zipfile.ZipFile" in app
    assert "shutil.make_archive(" not in app
    assert "W.Button(description='Show information'" in app
    assert "body.layout.display = '' if state['open'] else 'none'" in app
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
