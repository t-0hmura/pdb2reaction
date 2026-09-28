"""L4 Infra — output / I/O modules (no chemistry logic).

Modules:
- ``energy_diagram`` — render energy diagrams from numeric values (Plotly).
- ``hessian_cache`` — in-process Hessian cache shared by the stages of ``all``.
- ``pdb_fix`` — ``pdb2reaction fix-altloc`` subcommand backend (PDB altloc resolution).
- ``summary`` — per-run ``summary.json`` / ``summary.log`` writer.
- ``trj2fig`` — trajectory-to-figure plotting (``pdb2reaction trj2fig``).

Note: harmonic-restraint setup lives in ``pdb2reaction/workflows/restraints.py``
(L2 layer), not here.
"""
