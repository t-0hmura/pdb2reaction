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

Dependency direction: this layer **does not import any other in-package layer**
(L1 / L2 / L3 / L4). It is the leaf.
"""
