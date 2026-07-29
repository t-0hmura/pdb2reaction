"""L5 Foundation — default values, shared utilities, logging.

Modules:
- ``defaults`` — single source of truth for keyword arguments
  (``CALC_KW_DEFAULT``, ``UMA_CALC_KW``, ``OPT_BASE_KW``, ``LBFGS_KW``, ``RFO_KW``,
  ``BIAS_KW``, ``IRC_KW``, ``FREQ_KW``, ``SEARCH_KW``, etc.) and
  output-directory constants.
- ``utils`` — pure helpers (YAML parse, scan-list parsing, freeze-atom resolution,
  PDB metadata, format helpers, pretty-printers).
- ``logging`` — ``setup_logging(verbose)`` configures the stdlib logging level
  for the per-subcommand ``-v/--verbose LEVEL`` control.

Dependency direction (design intent): ``L1 -> L2 -> {L3, L4} -> L5``; the
foundation is meant to be the leaf. In practice ``utils`` still reaches DOWN into
lower data/service modules — element inference (``domain.add_elem_info``), the
charge engine and structure I/O (``io.charge`` / ``io.structure_formats``) — but
it imports **no** workflow (``workflows/*``), so the product module graph is
acyclic. ``.github/scripts/check_import_graph.py`` checks for product cycles and
``core``/``domain`` imports from ``workflows``.
"""
