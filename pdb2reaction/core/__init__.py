"""L5 Foundation — default values, shared utilities, logging.

Modules:
- ``defaults`` — single source of truth for keyword arguments
  (``CALC_KW_DEFAULT``, ``UMA_CALC_KW``, ``OPT_BASE_KW``, ``LBFGS_KW``, ``RFO_KW``,
  ``BIAS_KW``, ``IRC_KW``, ``FREQ_KW``, ``SEARCH_KW``, etc.) and
  output-directory constants.
- ``utils`` — pure helpers (YAML parse, scan-list parsing, freeze-atom resolution,
  PDB metadata, format helpers, pretty-printers).
- ``logging`` — ``setup_logging(verbose)`` driven by the ``-v`` / ``-vv`` root flag;
  configures stdlib logging level (WARNING / INFO / DEBUG).

Dependency direction (design intent): ``L1 -> L2 -> {L3, L4} -> L5``; the
foundation is meant to be the leaf. In practice ``utils`` still reaches DOWN into
lower data/service modules — element inference (``domain.add_elem_info``), the
charge engine and structure I/O (``io.charge`` / ``io.structure_formats``) — but
it imports **no** workflow (``workflows/*``): the ``core -> workflows`` back-edge
was removed in the C12 refactor, so the product module graph is acyclic. This is
enforced by ``.github/scripts/check_import_graph.py`` (no product cycle; no
``core``/``domain`` -> ``workflows`` edge), not by a marker checker.
"""
