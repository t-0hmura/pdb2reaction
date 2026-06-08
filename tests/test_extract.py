"""Unit tests for pdb2reaction.extract pure helper functions."""

from __future__ import annotations



class TestFormatEchoMessage:
    def test_no_args(self):
        from pdb2reaction.workflows.extract import _format_echo_message

        assert _format_echo_message("hello") == "hello"

    def test_with_args(self):
        from pdb2reaction.workflows.extract import _format_echo_message

        assert _format_echo_message("value=%d", 42) == "value=42"

    def test_format_error_fallback(self):
        from pdb2reaction.workflows.extract import _format_echo_message

        # Wrong format specifier should fall back gracefully
        result = _format_echo_message("value=%d", "not_a_number")
        assert isinstance(result, str)


class TestParseResTokens:
    def test_single_resseq(self):
        from pdb2reaction.workflows.extract import _parse_res_tokens

        result = _parse_res_tokens("123")
        assert len(result) == 1
        chain, resseq, icode = result[0]
        assert resseq == 123

    def test_chain_resseq(self):
        from pdb2reaction.workflows.extract import _parse_res_tokens

        result = _parse_res_tokens("A:123")
        assert len(result) == 1
        chain, resseq, icode = result[0]
        assert chain == "A"
        assert resseq == 123

    def test_chain_resseq_icode(self):
        from pdb2reaction.workflows.extract import _parse_res_tokens

        result = _parse_res_tokens("A:123A")
        assert len(result) == 1
        chain, resseq, icode = result[0]
        assert chain == "A"
        assert resseq == 123
        assert icode == "A"

    def test_multiple(self):
        from pdb2reaction.workflows.extract import _parse_res_tokens

        result = _parse_res_tokens("A:10,B:20")
        assert len(result) == 2


class TestMaxSerialFromPdbText:
    def test_basic(self):
        from pdb2reaction.workflows.extract import _max_serial_from_pdb_text

        pdb_text = (
            "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N\n"
            "ATOM     42  CA  ALA A   1       4.000   5.000   6.000  1.00  0.00           C\n"
        )
        assert _max_serial_from_pdb_text(pdb_text) == 42

    def test_empty(self):
        from pdb2reaction.workflows.extract import _max_serial_from_pdb_text

        assert _max_serial_from_pdb_text("") == 0


class TestStripTrailingEND:
    def test_strips_end(self):
        from pdb2reaction.workflows.extract import _strip_trailing_END

        text = "ATOM      1  N   ALA A   1\nEND\n"
        result = _strip_trailing_END(text)
        assert "END" not in result.split("\n")[-1] or result.strip().endswith("ALA A   1")

    def test_no_end(self):
        from pdb2reaction.workflows.extract import _strip_trailing_END

        text = "ATOM      1  N   ALA A   1\n"
        result = _strip_trailing_END(text)
        assert result == text


class TestParseLigandChargeOption:
    def test_none(self):
        from pdb2reaction.workflows.extract import _parse_ligand_charge_option

        total, mapping = _parse_ligand_charge_option(None)
        assert total is None
        assert mapping is None

    def test_numeric_int(self):
        from pdb2reaction.workflows.extract import _parse_ligand_charge_option

        total, mapping = _parse_ligand_charge_option(-3)
        assert total == -3.0
        assert mapping is None

    def test_numeric_string(self):
        from pdb2reaction.workflows.extract import _parse_ligand_charge_option

        total, mapping = _parse_ligand_charge_option("-3")
        assert total == -3.0
        assert mapping is None

    def test_mapping_string(self):
        from pdb2reaction.workflows.extract import _parse_ligand_charge_option

        total, mapping = _parse_ligand_charge_option("GPP:-3,SAM:1")
        assert total is None
        assert mapping == {"GPP": -3.0, "SAM": 1.0}

    def test_dict_input(self):
        from pdb2reaction.workflows.extract import _parse_ligand_charge_option

        total, mapping = _parse_ligand_charge_option({"GPP": -3, "SAM": 1})
        assert total is None
        assert mapping == {"GPP": -3.0, "SAM": 1.0}


class TestFormatLinkHBlock:
    def test_basic(self):
        from pdb2reaction.workflows.extract import _format_linkH_block

        coords = [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0)]
        result = _format_linkH_block(coords, start_serial=100)
        lines = [l for l in result.strip().split("\n") if l.strip()]
        assert len(lines) == 2
        assert "HETATM" in lines[0]
        # start_serial=100 means first atom gets serial 101
        assert "101" in lines[0] or " 101" in lines[0]


class TestIonDictCaseFolding:
    """Regression test: lookup in compute_charge_summary is `rn = res.get_resname().upper()`,
    so any mixed-case ION key is unreachable and silently produces charge 0."""

    def test_all_keys_uppercase(self):
        from pdb2reaction.workflows.extract import ION

        non_upper = [k for k in ION if k != k.upper()]
        assert non_upper == [], (
            f"ION dict must use only uppercase keys (case-folded lookup); "
            f"found mixed-case keys: {non_upper}"
        )

    def test_no_duplicate_uppercase_keys(self):
        """If two source-level keys differ only in case, the dict literal collapses to one
        and the other charge is silently lost. Guard against future reintroduction."""
        import ast

        from pdb2reaction.workflows import extract as _extract_mod

        src = open(_extract_mod.__file__).read()
        tree = ast.parse(src)
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.AnnAssign)
                and getattr(node.target, "id", None) == "ION"
                and isinstance(node.value, ast.Dict)
            ):
                raw_keys = [
                    k.value for k in node.value.keys if isinstance(k, ast.Constant)
                ]
                upper_keys = [k.upper() for k in raw_keys]
                dupes = {k for k in upper_keys if upper_keys.count(k) > 1}
                assert not dupes, (
                    f"ION dict has source-level duplicate uppercase keys: {sorted(dupes)}"
                )
                return
        raise AssertionError("ION dict literal not found in extract.py via AST")


def test_compute_charge_summary_terminus_cap_charges():
    """A kept terminal cap carries the ionized-terminus formal charge the internal
    residue charge omits: C-terminus carboxylate (OXT) -> -1, N-terminus ammonium
    (H1/H2/H3) -> +1. No correction unless the cap is kept (keep_ncap/ccap)."""
    import io as _io
    from Bio import PDB
    from pdb2reaction.workflows.extract import compute_charge_summary

    def _atom(serial, atom, resname, chain, resseq, x, y, z, element):
        return (
            f"ATOM  {serial:>5} {atom:<4} {resname:>3} {chain:1}{resseq:>4}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{20.00:>6.2f}          {element:>2}\n"
        )

    lines = [
        _atom(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        _atom(2, "CA", "ALA", "A", 1, 1.5, 0.0, 0.0, "C"),
        _atom(3, "C", "ALA", "A", 1, 2.0, 1.4, 0.0, "C"),
        _atom(4, "O", "ALA", "A", 1, 1.3, 2.4, 0.0, "O"),
        _atom(5, "OXT", "ALA", "A", 1, 3.3, 1.4, 0.0, "O"),
        _atom(6, "CB", "ALA", "A", 1, 2.1, -1.2, 0.0, "C"),
        "TER\n",
        _atom(7, "N", "ALA", "B", 1, 10.0, 0.0, 0.0, "N"),
        _atom(8, "H1", "ALA", "B", 1, 9.5, 0.8, 0.0, "H"),
        _atom(9, "H2", "ALA", "B", 1, 9.5, -0.8, 0.0, "H"),
        _atom(10, "H3", "ALA", "B", 1, 10.8, 0.0, 0.0, "H"),
        _atom(11, "CA", "ALA", "B", 1, 11.5, 0.0, 0.0, "C"),
        _atom(12, "C", "ALA", "B", 1, 12.0, 1.4, 0.0, "C"),
        _atom(13, "O", "ALA", "B", 1, 11.3, 2.4, 0.0, "O"),
        _atom(14, "CB", "ALA", "B", 1, 12.1, -1.2, 0.0, "C"),
        "TER\n", "END\n",
    ]
    structure = PDB.PDBParser(QUIET=True).get_structure("term", _io.StringIO("".join(lines)))
    res = list(structure.get_residues())
    cterm = next(r.get_full_id() for r in res if r.get_parent().id == "A")
    nterm = next(r.get_full_id() for r in res if r.get_parent().id == "B")
    sel = {r.get_full_id() for r in res}

    assert compute_charge_summary(structure, sel, set())["protein_charge"] == 0.0
    assert compute_charge_summary(structure, sel, set(), keep_ccap_ids={cterm})["protein_charge"] == -1.0
    assert compute_charge_summary(structure, sel, set(), keep_ncap_ids={nterm})["protein_charge"] == 1.0
    assert compute_charge_summary(
        structure, sel, set(), keep_ncap_ids={nterm}, keep_ccap_ids={cterm}
    )["protein_charge"] == 0.0
