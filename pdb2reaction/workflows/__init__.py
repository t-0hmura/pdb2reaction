"""L2 Application — per-subcommand orchestration (one module per ``pdb2reaction`` subcommand).

Each module exposes a ``cli`` callable wired to the corresponding CLI subcommand
(``pdb2reaction <name>``). The ``_LAZY_SUBCOMMANDS`` registry in
``pdb2reaction.cli.app`` lazy-loads these modules on demand to keep
``pdb2reaction --help`` startup fast.

Modules:
- ``all`` — end-to-end pipeline (MEP → TS → IRC → freq → DFT).
- ``opt`` — optimization (LBFGS grad / RFO hess).
- ``tsopt`` — TS optimization (RFO / Dimer / Bofill flatten loop).
- ``freq`` — vibrational analysis + thermochemistry.
- ``irc`` — IRC integration (EulerPC).
- ``path_opt`` / ``path_search`` — reaction path search and optimization.
- ``scan`` / ``scan2d`` / ``scan3d`` — bond-length scans (1D/2D/3D).
- ``extract`` — pocket / model-region extractor.
- ``dft`` — single-point DFT on the ML region (PySCF / GPU4PySCF).
- ``align_freeze`` — coordinate alignment + freeze-atom selection.
- ``scan_common`` — shared scan CLI options + helpers (factory pattern).
"""
