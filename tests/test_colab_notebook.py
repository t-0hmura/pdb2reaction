"""Static contracts for the release-matched Colab GUI notebook."""

from __future__ import annotations

import ast
import json
from pathlib import Path


NOTEBOOK = Path(__file__).parents[1] / "examples" / "pdb2reaction_colab.ipynb"


def _notebook() -> dict:
    return json.loads(NOTEBOOK.read_text(encoding="utf-8"))


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
    assert "v0.4.4" not in setup
    assert "Restart the Colab runtime first" in setup
    assert "mace-torch>=0.3.8" in setup
    assert "HF_TOKEN" in setup
    assert "installed_version != pdb2reaction_ref[1:]" in setup


def test_colab_gui_tracks_current_structure_and_execution_contracts() -> None:
    app = _notebook()["cells"][2]["source"]

    assert ".pdb,.ent,.cif,.mmcif,.xyz,.gjf" in app
    assert "prepare_input_structure" in app
    assert "load_pdb_atom_metadata" in app
    assert "CHAIN:RESNAME:RESSEQ" in app
    assert "String(atom.index)" in app
    # One occurrence registers the Python callback and one invokes it from
    # py3Dmol JavaScript; keeping both names identical is the round-trip gate.
    assert app.count("pdb2reaction_gui.on_click") == 2
    assert "demo.on_click" not in app
    assert "return shlex.split(line)" in app
    assert "if '--dry-run' not in a: a = a + ['--dry-run']" in app
    assert "RUN exact command" in app
    assert "reuse non-empty out dir" in app
    assert "uma-s-1p2" in app
    assert "MODELS = {'mace': ['MACE-OMOL-0']" in app
    assert "options=['auto', 'fp32', 'fp64']" in app
    assert "unknown options; use Validate" in app
    assert "py3Dmol.view(width='100%', height=380)" in app
    assert "rxworkspace" in app
    assert "rxviewer" in app
    assert "rxinspector" in app
    assert ":focus-visible" in app
    assert "overflow-wrap:anywhere" in app
    assert "max_width='100%'" in app
    assert "No structure loaded" in app
    # Colab renders ipywidgets' Tab as an empty block, so the tab strip is Buttons
    # + a swapping VBox. Pin the four pages, the navigation entry point, and guard
    # against a Tab regression.
    assert "_TAB_PAGES = [('1 Input', input_box), ('2 Select', select_box)," in app
    assert "('3 Options', options_box), ('4 Results', results_box)]" in app
    assert "def _tab_go(i):" in app
    assert "W.Tab(" not in app
    # Accordion has the same Colab rendering defect, so every collapsible (NGLView,
    # freezing/measurement, advanced flags) is the same Button + VBox pattern.
    assert "def _collapsible(title, child, on_open=None):" in app
    assert "W.Accordion(" not in app
    assert "adv_acc = _collapsible(" in app
    assert "freeze_acc = _collapsible(" in app
    assert "ngl_acc = _collapsible(" in app
    # "all workflow" mode and the depth switches only apply to `all`.
    assert "all_mode_box = W.VBox([W.HTML('<b>all workflow</b>" in app
    assert "depth_box = W.VBox([W.HTML('<b>Depth</b> (off = fast MEP only):')," in app
    assert "W.HBox([w_ts, w_th, adv_dft])" in app
    assert "for _name in ('all_mode_box', 'depth_box'):" in app
    # -r/--radius is an extraction option: only all and extract accept it.
    assert "if sub in ('all', 'extract') and r and r > 0: cmd += ['-r', str(r)]" in app
    assert "adv_radius.disabled = sub not in FLAG_SUBS['adv_radius']" in app
    assert "req='1 file + scan-lists" in app
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
    assert "if 'freeze' in SPEC.get(sub, {}).get('panels', ()) and S['freeze_atoms']:" in app
    assert "_panels = SPEC.get(sub, {}).get('panels', ())" in app
    # Surfaced key flags + their emission.
    assert "key_opts_box = W.VBox([" in app
    assert "cmd += ['--flatten']" in app
    assert "cmd += ['--refine-path']" in app
    assert "cmd += ['--max-cycles', str(int(mc))]" in app
    # Hover help is the HTML `title` global attribute (renders in Colab, where
    # widget-level tooltips on an HTML box do not), with a ⓘ affordance.
    assert "def _hdr(html, tip):" in app
    assert "&#9432;" in app
    assert 'title="%s"' in app
    # Standalone opt/tsopt/path-opt users need a cluster model first.
    assert "b_extract = W.Button(description='Extract cluster model'" in app
    assert "cmd = [CLI, 'extract', '-i', *S['inputs'], '-o', out, '-c', ','.join(cen)]" in app
    assert "utility_bar = W.HBox([b_rebuild, b_copy, b_clear])" in app
    assert "ps.get('mlip')" in app
    assert "ps.get('gibbs_mlip')" in app
    assert "ps.get('uma')" not in app
    assert "ps.get('gibbs_uma')" not in app
    assert "def _effective_out_dir(argv=None):" in app
    assert "out = _effective_out_dir(a)" in app
    assert "_results(out)" in app
    assert "'adv_mep':     {'all', 'path-opt', 'path-search'}," in app
    assert "sub in TOOL_CAPABILITIES['threshold']" in app
    assert "elif sub == 'dft':" in app
    assert "cmd += ['--func-basis', fb]" in app
    assert "adv_mep.disabled = sub not in TOOL_CAPABILITIES['mep_mode']" in app
    assert "d['all_mode'] = _wv('all_mode', 'mep')" in app
    assert "all_mode.value = saved_all_mode" in app
    assert "bytes(c).decode('utf-8')" in app
