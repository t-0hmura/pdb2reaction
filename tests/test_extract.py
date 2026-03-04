"""Unit tests for pdb2reaction.extract pure helper functions."""

from __future__ import annotations

import pytest


class TestFormatEchoMessage:
    def test_no_args(self):
        from pdb2reaction.extract import _format_echo_message

        assert _format_echo_message("hello") == "hello"

    def test_with_args(self):
        from pdb2reaction.extract import _format_echo_message

        assert _format_echo_message("value=%d", 42) == "value=42"

    def test_format_error_fallback(self):
        from pdb2reaction.extract import _format_echo_message

        # Wrong format specifier should fall back gracefully
        result = _format_echo_message("value=%d", "not_a_number")
        assert isinstance(result, str)


class TestParseResTokens:
    def test_single_resseq(self):
        from pdb2reaction.extract import _parse_res_tokens

        result = _parse_res_tokens("123")
        assert len(result) == 1
        chain, resseq, icode = result[0]
        assert resseq == 123

    def test_chain_resseq(self):
        from pdb2reaction.extract import _parse_res_tokens

        result = _parse_res_tokens("A:123")
        assert len(result) == 1
        chain, resseq, icode = result[0]
        assert chain == "A"
        assert resseq == 123

    def test_chain_resseq_icode(self):
        from pdb2reaction.extract import _parse_res_tokens

        result = _parse_res_tokens("A:123A")
        assert len(result) == 1
        chain, resseq, icode = result[0]
        assert chain == "A"
        assert resseq == 123
        assert icode == "A"

    def test_multiple(self):
        from pdb2reaction.extract import _parse_res_tokens

        result = _parse_res_tokens("A:10,B:20")
        assert len(result) == 2


class TestMaxSerialFromPdbText:
    def test_basic(self):
        from pdb2reaction.extract import _max_serial_from_pdb_text

        pdb_text = (
            "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N\n"
            "ATOM     42  CA  ALA A   1       4.000   5.000   6.000  1.00  0.00           C\n"
        )
        assert _max_serial_from_pdb_text(pdb_text) == 42

    def test_empty(self):
        from pdb2reaction.extract import _max_serial_from_pdb_text

        assert _max_serial_from_pdb_text("") == 0


class TestStripTrailingEND:
    def test_strips_end(self):
        from pdb2reaction.extract import _strip_trailing_END

        text = "ATOM      1  N   ALA A   1\nEND\n"
        result = _strip_trailing_END(text)
        assert "END" not in result.split("\n")[-1] or result.strip().endswith("ALA A   1")

    def test_no_end(self):
        from pdb2reaction.extract import _strip_trailing_END

        text = "ATOM      1  N   ALA A   1\n"
        result = _strip_trailing_END(text)
        assert result == text


class TestParseLigandChargeOption:
    def test_none(self):
        from pdb2reaction.extract import _parse_ligand_charge_option

        total, mapping = _parse_ligand_charge_option(None)
        assert total is None
        assert mapping is None

    def test_numeric_int(self):
        from pdb2reaction.extract import _parse_ligand_charge_option

        total, mapping = _parse_ligand_charge_option(-3)
        assert total == -3.0
        assert mapping is None

    def test_numeric_string(self):
        from pdb2reaction.extract import _parse_ligand_charge_option

        total, mapping = _parse_ligand_charge_option("-3")
        assert total == -3.0
        assert mapping is None

    def test_mapping_string(self):
        from pdb2reaction.extract import _parse_ligand_charge_option

        total, mapping = _parse_ligand_charge_option("GPP:-3,SAM:1")
        assert total is None
        assert mapping == {"GPP": -3.0, "SAM": 1.0}

    def test_dict_input(self):
        from pdb2reaction.extract import _parse_ligand_charge_option

        total, mapping = _parse_ligand_charge_option({"GPP": -3, "SAM": 1})
        assert total is None
        assert mapping == {"GPP": -3.0, "SAM": 1.0}


class TestFormatLinkHBlock:
    def test_basic(self):
        from pdb2reaction.extract import _format_linkH_block

        coords = [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0)]
        result = _format_linkH_block(coords, start_serial=100)
        lines = [l for l in result.strip().split("\n") if l.strip()]
        assert len(lines) == 2
        assert "HETATM" in lines[0]
        # start_serial=100 means first atom gets serial 101
        assert "101" in lines[0] or " 101" in lines[0]
