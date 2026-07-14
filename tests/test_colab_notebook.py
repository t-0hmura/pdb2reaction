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
    assert "py3Dmol.view(width='100%', height=440)" in app
    assert "rxworkspace" in app
    assert "rxviewer" in app
    assert "rxinspector" in app
    assert ":focus-visible" in app
    assert "overflow-wrap:anywhere" in app
    assert "max_width='100%'" in app
    assert "No structure loaded" in app
    assert "['1 Input', '2 Select', '3 Options', '4 Results']" in app
    assert "utility_bar = W.HBox([b_rebuild, b_copy, b_clear])" in app
    assert "ps.get('mlip')" in app
    assert "ps.get('gibbs_mlip')" in app
    assert "ps.get('uma')" not in app
    assert "ps.get('gibbs_uma')" not in app
    assert "def _effective_out_dir(argv=None):" in app
    assert "out = _effective_out_dir(a)" in app
    assert "_results(out)" in app
