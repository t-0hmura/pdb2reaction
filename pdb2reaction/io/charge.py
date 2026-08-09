"""Active-site model charge engine (structure -> charge summary).

Product-local service consumed by both the extraction workflow
(``workflows.extract``) and CLI-input preparation (``core.utils``). It lives
below ``core`` and depends only on leaf data and identity modules.

Pure of ``core``/``workflows``: it depends only on the canonical residue tables
(``domain.residue_data``) and structure-identity helpers (``io.structure_formats``).
Verbosity-gated logging of a summary (``log_charge_summary``) lives in
``core.utils`` because it reads the CLI console-gating state.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Set, Tuple

import click

from pdb2reaction.domain.residue_data import (
    AMINO_ACIDS,
    C_TERMINAL_RESNAMES,
    ION,
    N_TERMINAL_RESNAMES,
    ResidueKey,
    WATER_RES,
)
from pdb2reaction.io.structure_formats import residue_auth_identity


def _format_echo_message(msg: str, *args: Any) -> str:
    if not args:
        return str(msg)
    try:
        return str(msg) % args
    except Exception:
        tail = " ".join(str(x) for x in args)
        return f"{msg} {tail}".strip()


def _echo_warning(msg: str, *args: Any) -> None:
    click.echo(f"WARNING: {_format_echo_message(msg, *args)}", err=True)


def _fmt_res_id(res) -> str:
    """
    Return a compact residue tag like 'A:123A SER' or '123 SER'.
    """
    chain, resseq, icode_txt, resname = residue_auth_identity(res)
    chain_txt = f"{chain}:" if chain else ""
    return f"{chain_txt}{resseq}{icode_txt} {resname}"


def _fmt_fid(structure, fid: Tuple) -> str:
    """
    Format a full-id into a human-friendly residue tag.
    """
    res = structure[fid[1]][fid[2]].child_dict[fid[3]]
    return _fmt_res_id(res)


def _residue_key_from_res(res) -> ResidueKey:
    """
    Build a cross-structure residue key from a residue.
    """
    chain_id, resseq, icode, resname = residue_auth_identity(res)
    hetflag = str(res.id[0])
    return (chain_id, hetflag, str(resseq), icode, resname)


def _residue_key_from_fid(structure, fid: Tuple) -> ResidueKey:
    """
    Build a cross-structure residue key from a full-id.
    """
    res = structure[fid[1]][fid[2]].child_dict[fid[3]]
    return _residue_key_from_res(res)


def _residue_name(res) -> str:
    return residue_auth_identity(res)[3].upper()


# ---- helper for parsing --ligand-charge (number or 'RES:Q' mapping) ----
def _parse_ligand_charge_option(ligand_charge: float | str | Dict[str, float] | None
                                ) -> Tuple[Optional[float], Optional[Dict[str, float]]]:
    """
    Returns
    -------
    (total_charge, mapping)
      total_charge : float | None
      mapping      : dict[RESNAME -> float] | None
    """
    if ligand_charge is None:
        return None, None
    if isinstance(ligand_charge, (int, float)):
        return float(ligand_charge), None
    if isinstance(ligand_charge, dict):
        mapping = {str(k).upper(): float(v) for k, v in ligand_charge.items()}
        return None, mapping
    if isinstance(ligand_charge, str):
        s = ligand_charge.strip()
        if not s:
            return None, None
        # try numeric
        try:
            return float(s), None
        except ValueError:
            pass
        # mapping: tokens "RES:Q"
        tokens = [t for t in re.split(r"[,\s]+", s) if t]
        mapping: Dict[str, float] = {}
        for tok in tokens:
            if ":" not in tok:
                raise click.BadParameter(
                    f"Invalid --ligand-charge token '{tok}'. "
                    "Use 'RES:Q' (e.g., GPP:-3) or a bare number (e.g., -3).",
                    param_hint="-l / --ligand-charge",
                )
            res, qtxt = tok.split(":", 1)
            resname = res.strip().upper()
            if not resname:
                raise click.BadParameter(
                    f"Invalid --ligand-charge token '{tok}': empty residue name.",
                    param_hint="-l / --ligand-charge",
                )
            try:
                qval = float(qtxt.strip())
            except ValueError:
                raise click.BadParameter(
                    f"Invalid --ligand-charge token '{tok}': "
                    f"'{qtxt}' is not a number.",
                    param_hint="-l / --ligand-charge",
                )
            mapping[resname] = qval
        if not mapping:
            raise ValueError("Empty --ligand-charge mapping.")
        return None, mapping
    raise TypeError(f"Unsupported type for ligand_charge: {type(ligand_charge)!r}")


def _sorted_fids_by_file_order(structure, fids) -> List[Tuple]:
    """
    Sort full-ids by file order using a residue index map.
    """
    order: Dict[Tuple, int] = {}
    idx = 0
    for model in structure:
        for chain in model:
            for res in chain.get_residues():
                order[res.get_full_id()] = idx
                idx += 1
    return sorted(set(fids), key=lambda fid: order.get(fid, 10**12))


def infer_present_terminal_cap_ids(
    structure,
    selected_ids: Set[Tuple],
) -> Tuple[Set[Tuple], Set[Tuple]]:
    """Return selected residues whose explicit atoms encode ionized termini.

    This helper is for full/model structures that were not produced by the
    extraction truncation logic.  It deliberately infers only from explicit
    OXT or H1/H2/H3 atoms; it does not guess terminal protonation from chain
    position alone.
    """

    keep_ncap_ids: Set[Tuple] = set()
    keep_ccap_ids: Set[Tuple] = set()
    for fid in selected_ids:
        res = structure[fid[1]][fid[2]].child_dict[fid[3]]
        if _residue_name(res) not in AMINO_ACIDS:
            continue
        atom_names = {atom.get_name().strip().upper() for atom in res}
        if "OXT" in atom_names:
            keep_ccap_ids.add(fid)
        if {"H1", "H2", "H3"} <= atom_names:
            keep_ncap_ids.add(fid)
    return keep_ncap_ids, keep_ccap_ids


def compute_charge_summary(structure,
                           selected_ids: Set[Tuple],
                           substrate_ids: Set[Tuple],
                           ligand_charge: float | str | Dict[str, float] | None = None,
                           keep_ncap_ids: Set[Tuple] | None = None,
                           keep_ccap_ids: Set[Tuple] | None = None) -> Dict[str, Any]:
    """
    Compute model charge summary.

    Args
    ----
    structure : Bio.PDB.Structure.Structure
        The (first) structure to evaluate.
    selected_ids : set[tuple]
        Residues included in the model.
    substrate_ids : set[tuple]
        Residues designated as substrate.
    ligand_charge : float | str | dict[str,float] | None
        - float: total charge to assign across **unknown residues** (preferring unknown substrate).
        - str  : numeric string (total) or mapping like 'GPP:-3,SAM:1' (per‑resname).
        - dict : mapping {RESNAME: charge}. In mapping mode, other unknown residues remain 0.

    Returns
    -------
    dict with keys:
      - total_charge : float
      - protein_charge : float
      - ligand_total_charge : float
      - ion_total_charge : float
      - ion_charges : list[(str tag, float)]
      - unknown_residue_charges : dict[str -> float]  # for concise per‑resname log
    """
    keep_ncap_ids = keep_ncap_ids or set()
    keep_ccap_ids = keep_ccap_ids or set()
    terminal_corrections: List[Tuple[str, Any, int]] = []

    per_map: Dict[ResidueKey, float] = {}
    aa_charge = 0.0
    total = 0.0

    fids_in_order = _sorted_fids_by_file_order(structure, selected_ids)

    # First pass: dictionary/ion/water charges; collect unknowns and ions
    unknown_fids: List[Tuple] = []
    unknown_substrate_fids: List[Tuple] = []
    ion_entries: List[Tuple[str, float]] = []
    selected_ion_charges: Dict[str, float] = {}

    for fid in fids_in_order:
        res = structure[fid[1]][fid[2]].child_dict[fid[3]]
        rn = _residue_name(res)
        key = _residue_key_from_res(res)
        if rn in WATER_RES:
            q = 0.0
        elif rn in AMINO_ACIDS:
            q = float(AMINO_ACIDS[rn])
            # Terminal-cap charge: a residue at a true chain terminus whose
            # terminal cap is KEPT (keep_ncap_ids / keep_ccap_ids from the
            # truncation logic) carries an ionized-terminus formal charge that
            # the internal AMINO_ACIDS value omits. Mirror the cap atoms the
            # truncation preserves (N-cap N/H*; C-cap C/O/OXT):
            #   C-terminus carboxylate (COO-, has OXT)       -> -1
            #   N-terminus ammonium    (NH3+, has H1/H2/H3)  -> +1
            # TER-break caps have neither OXT nor NH3+ Hs, so they stay neutral.
            if keep_ccap_ids or keep_ncap_ids:
                atom_names = {a.get_name() for a in res}
                if fid in keep_ccap_ids and "OXT" in atom_names \
                        and rn not in C_TERMINAL_RESNAMES:
                    q -= 1.0
                    terminal_corrections.append(
                        (residue_auth_identity(res)[3], res.id[1], -1)
                    )
                if fid in keep_ncap_ids and {"H1", "H2", "H3"} <= atom_names \
                        and rn not in N_TERMINAL_RESNAMES:
                    q += 1.0
                    terminal_corrections.append(
                        (residue_auth_identity(res)[3], res.id[1], 1)
                    )
            aa_charge += q
        elif rn in ION:
            q = float(ION[rn])
            ion_entries.append((_fmt_fid(structure, fid), q))
            selected_ion_charges[rn] = q
        else:
            q = 0.0
            unknown_fids.append(fid)
            if fid in substrate_ids:
                unknown_substrate_fids.append(fid)
        per_map[key] = q
        total += q

    # Apply --ligand-charge if provided
    total_spec, mapping_spec = _parse_ligand_charge_option(ligand_charge)

    if total_spec is not None:
        # Distribute total across unknown substrate if present, else across all unknowns
        targets = unknown_substrate_fids if unknown_substrate_fids else unknown_fids
        if targets:
            per_res_val = float(total_spec) / float(len(targets))
            for fid in targets:
                key = _residue_key_from_fid(structure, fid)
                per_map[key] = per_res_val
            # recompute totals
            total = sum(per_map.values())
            aa_charge = sum(q for k, q in per_map.items() if k[4] in AMINO_ACIDS)
        else:
            _echo_warning(
                "[extract] --ligand-charge %s was provided but no unknown "
                "(non-dictionary) residues were found to apply it to; the value "
                "is ignored — check the substrate/ligand resname.", total_spec)
    elif mapping_spec is not None:
        # Per‑resname mapping. Unspecified unknown residues remain 0.
        matched_resnames: Set[str] = set()
        for fid in unknown_fids:
            res = structure[fid[1]][fid[2]].child_dict[fid[3]]
            rn = _residue_name(res)
            if rn in mapping_spec:
                key = _residue_key_from_fid(structure, fid)
                per_map[key] = float(mapping_spec[rn])
                matched_resnames.add(rn)
        # The Colab charge editor shows locked, authoritative ion charges in
        # the generated -l mapping so the complete charge choice is visible.
        # An exact table value is therefore an acknowledged no-op, not an
        # unmatched ligand warning.  Conflicting ion values remain unmatched
        # and are reported below; the internal ion table still wins.
        matched_resnames.update(
            rn
            for rn, q in selected_ion_charges.items()
            if rn in mapping_spec and abs(float(mapping_spec[rn]) - q) <= 1e-9
        )
        unmatched_resnames = sorted(set(mapping_spec) - matched_resnames)
        if unmatched_resnames:
            _echo_warning(
                "[extract] --ligand-charge mapping entr%s %s matched no "
                "unknown (non-dictionary) selected residue%s and %s ignored. "
                "Standard/modified amino acids, ions, and water use the "
                "internal charge tables; otherwise check the input and "
                "residue selector. Use -q to override the derived total charge.",
                "y" if len(unmatched_resnames) == 1 else "ies",
                ", ".join(unmatched_resnames),
                "" if len(unmatched_resnames) == 1 else "s",
                "was" if len(unmatched_resnames) == 1 else "were",
            )
        # recompute totals
        total = sum(per_map.values())
        aa_charge = sum(q for k, q in per_map.items() if k[4] in AMINO_ACIDS)

    # Net ligand and ion charges
    unknown_keys = {_residue_key_from_fid(structure, fid) for fid in unknown_fids}
    ligand_total = sum(per_map[k] for k in unknown_keys)
    ion_total = sum(q for _, q in ion_entries)

    # Build per‑resname mapping for unknown residues (after applying any overrides)
    unknown_residue_charges: Dict[str, float] = {}
    for fid in unknown_fids:
        res = structure[fid[1]][fid[2]].child_dict[fid[3]]
        rn = _residue_name(res)
        key = _residue_key_from_fid(structure, fid)
        unknown_residue_charges[rn] = float(per_map[key])

    return {
        "total_charge": float(total),
        "protein_charge": float(aa_charge),
        "ligand_total_charge": float(ligand_total),
        "ion_total_charge": float(ion_total),
        "ion_charges": [(tag, float(q)) for tag, q in ion_entries],
        "unknown_residue_charges": unknown_residue_charges,
        "terminal_corrections": terminal_corrections,
    }
