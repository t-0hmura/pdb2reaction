"""MCP tool definitions — 1 tool per pdb2reaction CLI subcommand.

Each `@mcp.tool()`-decorated function:
- Declares typed parameters → FastMCP auto-generates JSON schema for the agent.
- Translates kwargs into the corresponding pdb2reaction CLI argv.
- Runs the CLI via subprocess and parses `<out_dir>/summary.json`.
- Returns a structured dict via `SubcmdResult.to_dict()`.

Each tool also accepts:
- `out_dir`: explicit output directory (= CLI's `--out-dir`). Tools default to
  a unique-per-call temp dir if the agent does not supply one, so concurrent
  agent calls do not collide.
- `extra_args`: list of additional CLI flags the agent wants to pass that are
  not surfaced as named parameters here. Provides forward-compat without a
  schema change.
- `timeout_seconds`: subprocess timeout (default None = no timeout). Useful
  for the agent to cap a runaway opt before its own request budget expires.
"""
from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any, Optional

from pdb2reaction.mcp._runner import run_subcmd


def _resolve_out_dir(out_dir: Optional[str], prefix: str) -> Path:
    """Pick a usable output directory: explicit > unique temp."""
    if out_dir:
        path = Path(out_dir).expanduser().resolve()
        path.mkdir(parents=True, exist_ok=True)
        return path
    return Path(tempfile.mkdtemp(prefix=f"p2r_mcp_{prefix}_"))


def _shared_calc_flags(
    *,
    backend: Optional[str],
    solvent: Optional[str],
    solvent_model: Optional[str],
    precision: Optional[str],
    workers: Optional[int],
    workers_per_node: Optional[int],
) -> list[str]:
    """Common backend / solvent / precision / workers flags shared by most subcmds."""
    args: list[str] = []
    if backend:
        args.extend(["-b", str(backend)])
    if solvent:
        args.extend(["--solvent", str(solvent)])
    if solvent_model:
        args.extend(["--solvent-model", str(solvent_model)])
    if precision:
        args.extend(["--precision", str(precision)])
    if workers is not None:
        args.extend(["--workers", str(workers)])
    if workers_per_node is not None:
        args.extend(["--workers-per-node", str(workers_per_node)])
    return args


def register_all(mcp) -> None:
    """Register every pdb2reaction subcommand as an MCP tool on `mcp`."""

    # Core stage runners (opt / tsopt / irc / freq)

    @mcp.tool()
    def optimize_geometry(
        input_pdb: str,
        charge: int,
        *,
        ligand_charge: Optional[str] = None,
        multiplicity: Optional[int] = None,
        opt_mode: str = "grad",
        max_cycles: int = 10000,
        thresh: Optional[str] = None,
        coord_type: Optional[str] = None,
        backend: Optional[str] = None,
        solvent: Optional[str] = None,
        solvent_model: Optional[str] = None,
        precision: Optional[str] = None,
        workers: Optional[int] = None,
        workers_per_node: Optional[int] = None,
        print_every: Optional[int] = None,
        dump_trajectory: bool = False,
        ref_pdb: Optional[str] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Optimize a single molecular geometry (CLI: `pdb2reaction opt`).

        Default `opt_mode="grad"` matches the CLI default; pass `"hess"` for
        RFO Hessian-based optimization. Returns: SubcmdResult dict with
        summary.json contents (final_energy, n_cycles, converged, files map).
        """
        od = _resolve_out_dir(out_dir, "opt")
        argv: list[str] = ["pdb2reaction", "opt", "-i", input_pdb, "-q", str(charge)]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if multiplicity is not None:
            argv.extend(["-m", str(multiplicity)])
        argv.extend(["--opt-mode", opt_mode, "--max-cycles", str(max_cycles)])
        if thresh:
            argv.extend(["--thresh", thresh])
        if coord_type:
            argv.extend(["--coord-type", coord_type])
        if print_every is not None:
            argv.extend(["--print-every", str(print_every)])
        if dump_trajectory:
            argv.append("--dump")
        if ref_pdb:
            argv.extend(["--ref-pdb", ref_pdb])
        argv.extend(_shared_calc_flags(
            backend=backend, solvent=solvent, solvent_model=solvent_model,
            precision=precision, workers=workers, workers_per_node=workers_per_node,
        ))
        argv.append("--out-json")
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def find_transition_state(
        ts_pdb: str,
        charge: int,
        *,
        ligand_charge: Optional[str] = None,
        multiplicity: Optional[int] = None,
        opt_mode: str = "hess",
        max_cycles: int = 10000,
        thresh: Optional[str] = None,
        coord_type: Optional[str] = None,
        backend: Optional[str] = None,
        solvent: Optional[str] = None,
        solvent_model: Optional[str] = None,
        precision: Optional[str] = None,
        workers: Optional[int] = None,
        workers_per_node: Optional[int] = None,
        flatten: bool = False,
        hessian_calc_mode: Optional[str] = None,
        print_every: Optional[int] = None,
        ref_pdb: Optional[str] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Optimize a transition-state candidate (CLI: `pdb2reaction tsopt`).

        opt_mode = one of {"grad", "hess", "dimer", "rsirfo", "trim", "rsprfo"}.
        - "hess" / "rsprfo" (default): RS-P-RFO (Banerjee 1985).
        - "grad" / "dimer": Hessian-guided Dimer.
        - "trim" (Helgaker 1991) / "rsirfo": alternative TS optimizers.

        Returns: SubcmdResult dict with summary.json (ts_energy, n_imaginary_modes,
        imaginary_frequencies_cm, n_opt_cycles, files map).
        """
        od = _resolve_out_dir(out_dir, "tsopt")
        argv: list[str] = ["pdb2reaction", "tsopt", "-i", ts_pdb, "-q", str(charge)]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if multiplicity is not None:
            argv.extend(["-m", str(multiplicity)])
        argv.extend(["--opt-mode", opt_mode, "--max-cycles", str(max_cycles)])
        if thresh:
            argv.extend(["--thresh", thresh])
        if coord_type:
            argv.extend(["--coord-type", coord_type])
        if flatten:
            argv.append("--flatten")
        if hessian_calc_mode:
            argv.extend(["--hessian-calc-mode", hessian_calc_mode])
        if print_every is not None:
            argv.extend(["--print-every", str(print_every)])
        if ref_pdb:
            argv.extend(["--ref-pdb", ref_pdb])
        argv.extend(_shared_calc_flags(
            backend=backend, solvent=solvent, solvent_model=solvent_model,
            precision=precision, workers=workers, workers_per_node=workers_per_node,
        ))
        argv.append("--out-json")
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def run_irc(
        ts_pdb: str,
        charge: int,
        *,
        ligand_charge: Optional[str] = None,
        multiplicity: Optional[int] = None,
        max_cycles: Optional[int] = None,
        step_size: Optional[float] = None,
        root: Optional[int] = None,
        forward: Optional[bool] = None,
        backward: Optional[bool] = None,
        backend: Optional[str] = None,
        solvent: Optional[str] = None,
        solvent_model: Optional[str] = None,
        precision: Optional[str] = None,
        workers: Optional[int] = None,
        workers_per_node: Optional[int] = None,
        hessian_calc_mode: Optional[str] = None,
        irc_pos_def: Optional[bool] = None,
        ref_pdb: Optional[str] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Run IRC integration from a TS geometry (CLI: `pdb2reaction irc`).

        Opt-in:
        - irc_pos_def=True: convergence additionally requires positive-definite
          mass-weighted Hessian, blocking the IRC 'shoulder' false-convergence.

        Returns: SubcmdResult dict with summary.json (forward / backward
        energies, n_steps, files map).
        """
        od = _resolve_out_dir(out_dir, "irc")
        argv: list[str] = ["pdb2reaction", "irc", "-i", ts_pdb, "-q", str(charge)]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if multiplicity is not None:
            argv.extend(["-m", str(multiplicity)])
        if max_cycles is not None:
            argv.extend(["--max-cycles", str(max_cycles)])
        if step_size is not None:
            argv.extend(["--step-size", str(step_size)])
        if root is not None:
            argv.extend(["--root", str(root)])
        if forward is not None:
            argv.append("--forward" if forward else "--no-forward")
        if backward is not None:
            argv.append("--backward" if backward else "--no-backward")
        if hessian_calc_mode:
            argv.extend(["--hessian-calc-mode", hessian_calc_mode])
        if irc_pos_def is not None:
            argv.append("--irc-pos-def" if irc_pos_def else "--no-irc-pos-def")
        if ref_pdb:
            argv.extend(["--ref-pdb", ref_pdb])
        argv.extend(_shared_calc_flags(
            backend=backend, solvent=solvent, solvent_model=solvent_model,
            precision=precision, workers=workers, workers_per_node=workers_per_node,
        ))
        argv.append("--out-json")
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def compute_frequencies(
        input_pdb: str,
        charge: int,
        *,
        ligand_charge: Optional[str] = None,
        multiplicity: Optional[int] = None,
        temperature: Optional[float] = None,
        pressure_atm: Optional[float] = None,
        amplitude_ang: Optional[float] = None,
        n_frames: Optional[int] = None,
        backend: Optional[str] = None,
        solvent: Optional[str] = None,
        solvent_model: Optional[str] = None,
        precision: Optional[str] = None,
        workers: Optional[int] = None,
        workers_per_node: Optional[int] = None,
        hessian_calc_mode: Optional[str] = None,
        ref_pdb: Optional[str] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Vibrational analysis + thermochemistry (CLI: `pdb2reaction freq`).

        Returns: SubcmdResult dict with summary.json (frequencies_cm,
        ZPE, thermal corrections, n_imaginary, mode trajectory files).
        """
        od = _resolve_out_dir(out_dir, "freq")
        argv: list[str] = ["pdb2reaction", "freq", "-i", input_pdb, "-q", str(charge)]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if multiplicity is not None:
            argv.extend(["-m", str(multiplicity)])
        if temperature is not None:
            argv.extend(["--temperature", str(temperature)])
        if pressure_atm is not None:
            argv.extend(["--pressure", str(pressure_atm)])
        if amplitude_ang is not None:
            argv.extend(["--amplitude-ang", str(amplitude_ang)])
        if n_frames is not None:
            argv.extend(["--n-frames", str(n_frames)])
        if hessian_calc_mode:
            argv.extend(["--hessian-calc-mode", hessian_calc_mode])
        if ref_pdb:
            argv.extend(["--ref-pdb", ref_pdb])
        argv.extend(_shared_calc_flags(
            backend=backend, solvent=solvent, solvent_model=solvent_model,
            precision=precision, workers=workers, workers_per_node=workers_per_node,
        ))
        argv.append("--out-json")
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    # Scans (1D / 2D / 3D)

    @mcp.tool()
    def scan_1d(
        input_pdb: str,
        charge: int,
        scan_lists: str,
        *,
        ligand_charge: Optional[str] = None,
        multiplicity: Optional[int] = None,
        max_step_size: Optional[float] = None,
        relax_max_cycles: Optional[int] = None,
        opt_mode: Optional[str] = None,
        thresh: Optional[str] = None,
        preopt: Optional[bool] = None,
        endopt: Optional[bool] = None,
        backend: Optional[str] = None,
        solvent: Optional[str] = None,
        solvent_model: Optional[str] = None,
        precision: Optional[str] = None,
        workers: Optional[int] = None,
        workers_per_node: Optional[int] = None,
        print_every: Optional[int] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """1D distance scan with harmonic restraints (CLI: `pdb2reaction scan`).

        scan_lists: list-of-tuples literal, e.g. "[(1,5,1.4)]" = pull atom 1↔5
        target distance 1.4 Å.
        """
        od = _resolve_out_dir(out_dir, "scan")
        argv: list[str] = ["pdb2reaction", "scan", "-i", input_pdb, "-q", str(charge),
                           "--scan-lists", scan_lists]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if multiplicity is not None:
            argv.extend(["-m", str(multiplicity)])
        if max_step_size is not None:
            argv.extend(["--max-step-size", str(max_step_size)])
        if relax_max_cycles is not None:
            argv.extend(["--relax-max-cycles", str(relax_max_cycles)])
        if opt_mode:
            argv.extend(["--opt-mode", opt_mode])
        if thresh:
            argv.extend(["--thresh", thresh])
        if preopt is not None:
            argv.append("--preopt" if preopt else "--no-preopt")
        if endopt is not None:
            argv.append("--endopt" if endopt else "--no-endopt")
        if print_every is not None:
            argv.extend(["--print-every", str(print_every)])
        argv.extend(_shared_calc_flags(
            backend=backend, solvent=solvent, solvent_model=solvent_model,
            precision=precision, workers=workers, workers_per_node=workers_per_node,
        ))
        argv.append("--out-json")
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def scan_2d(
        input_pdb: str,
        scan_lists: str,
        *,
        charge: Optional[int] = None,
        ligand_charge: Optional[str] = None,
        multiplicity: Optional[int] = None,
        max_step_size: Optional[float] = None,
        relax_max_cycles: Optional[int] = None,
        thresh: Optional[str] = None,
        backend: Optional[str] = None,
        precision: Optional[str] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """2D distance scan (CLI: `pdb2reaction scan2d`)."""
        od = _resolve_out_dir(out_dir, "scan2d")
        argv: list[str] = ["pdb2reaction", "scan2d", "-i", input_pdb, "--scan-lists", scan_lists]
        if charge is not None:
            argv.extend(["-q", str(charge)])
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if multiplicity is not None:
            argv.extend(["-m", str(multiplicity)])
        if max_step_size is not None:
            argv.extend(["--max-step-size", str(max_step_size)])
        if relax_max_cycles is not None:
            argv.extend(["--relax-max-cycles", str(relax_max_cycles)])
        if thresh:
            argv.extend(["--thresh", thresh])
        if backend:
            argv.extend(["-b", backend])
        if precision:
            argv.extend(["--precision", precision])
        argv.append("--out-json")
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def scan_3d(
        input_pdb: str,
        scan_lists: str,
        *,
        charge: Optional[int] = None,
        ligand_charge: Optional[str] = None,
        multiplicity: Optional[int] = None,
        max_step_size: Optional[float] = None,
        relax_max_cycles: Optional[int] = None,
        thresh: Optional[str] = None,
        backend: Optional[str] = None,
        precision: Optional[str] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """3D distance scan (CLI: `pdb2reaction scan3d`)."""
        od = _resolve_out_dir(out_dir, "scan3d")
        argv: list[str] = ["pdb2reaction", "scan3d", "-i", input_pdb, "--scan-lists", scan_lists]
        if charge is not None:
            argv.extend(["-q", str(charge)])
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if multiplicity is not None:
            argv.extend(["-m", str(multiplicity)])
        if max_step_size is not None:
            argv.extend(["--max-step-size", str(max_step_size)])
        if relax_max_cycles is not None:
            argv.extend(["--relax-max-cycles", str(relax_max_cycles)])
        if thresh:
            argv.extend(["--thresh", thresh])
        if backend:
            argv.extend(["-b", backend])
        if precision:
            argv.extend(["--precision", precision])
        argv.append("--out-json")
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()


    @mcp.tool()
    def optimize_path(
        reactant_pdb: str,
        product_pdb: str,
        charge: int,
        *,
        ligand_charge: Optional[str] = None,
        multiplicity: Optional[int] = None,
        max_nodes: Optional[int] = None,
        max_cycles: Optional[int] = None,
        mep_mode: Optional[str] = None,
        preopt: Optional[bool] = None,
        climb: Optional[bool] = None,
        backend: Optional[str] = None,
        precision: Optional[str] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Optimize a reaction-path segment (CLI: `pdb2reaction path-opt`).

        Two-endpoint MEP optimization (GSM / NEB family).
        """
        od = _resolve_out_dir(out_dir, "path_opt")
        argv: list[str] = ["pdb2reaction", "path-opt", "-i", reactant_pdb, product_pdb,
                           "-q", str(charge)]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if multiplicity is not None:
            argv.extend(["-m", str(multiplicity)])
        if max_nodes is not None:
            argv.extend(["--max-nodes", str(max_nodes)])
        if max_cycles is not None:
            argv.extend(["--max-cycles", str(max_cycles)])
        if mep_mode:
            argv.extend(["--mep-mode", mep_mode])
        if preopt is not None:
            argv.append("--preopt" if preopt else "--no-preopt")
        if climb is not None:
            argv.append("--climb" if climb else "--no-climb")
        if backend:
            argv.extend(["-b", backend])
        if precision:
            argv.extend(["--precision", precision])
        argv.append("--out-json")
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def search_paths(
        input_pdb: str,
        charge: int,
        *,
        ligand_charge: Optional[str] = None,
        multiplicity: Optional[int] = None,
        max_nodes: Optional[int] = None,
        mep_mode: Optional[str] = None,
        refine_mode: Optional[str] = None,
        backend: Optional[str] = None,
        precision: Optional[str] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Recursive reaction-pathway search (CLI: `pdb2reaction path-search`)."""
        od = _resolve_out_dir(out_dir, "path_search")
        argv: list[str] = ["pdb2reaction", "path-search", "-i", input_pdb, "-q", str(charge)]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if multiplicity is not None:
            argv.extend(["-m", str(multiplicity)])
        if max_nodes is not None:
            argv.extend(["--max-nodes", str(max_nodes)])
        if mep_mode:
            argv.extend(["--mep-mode", mep_mode])
        if refine_mode:
            argv.extend(["--refine-mode", refine_mode])
        if backend:
            argv.extend(["-b", backend])
        if precision:
            argv.extend(["--precision", precision])
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()


    @mcp.tool()
    def run_full_pipeline(
        reactant_pdb: str,
        product_pdb: Optional[str] = None,
        *,
        charge: Optional[int] = None,
        ligand_charge: Optional[str] = None,
        multiplicity: Optional[int] = None,
        extract_ligand: Optional[str] = None,
        extract_radius: Optional[float] = None,
        scan_lists: Optional[str] = None,
        max_cycles: Optional[int] = None,
        thresh: Optional[str] = None,
        thresh_post: Optional[str] = None,
        backend: Optional[str] = None,
        precision: Optional[str] = None,
        refine_path: Optional[bool] = None,
        do_tsopt: Optional[bool] = None,
        do_dft: Optional[bool] = None,
        do_thermo: Optional[bool] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """End-to-end workflow extract → MEP → TS → IRC → freq → DFT (CLI: `pdb2reaction`).

        This is the `all` subcommand (registered as the root cli group default).
        Returns: SubcmdResult dict with the full pipeline summary.json including
        per-stage results, activation energies, reaction energies.
        """
        od = _resolve_out_dir(out_dir, "all")
        argv: list[str] = ["pdb2reaction", "-i", reactant_pdb]
        if product_pdb:
            argv.append(product_pdb)
        if charge is not None:
            argv.extend(["-q", str(charge)])
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if multiplicity is not None:
            argv.extend(["-m", str(multiplicity)])
        if extract_ligand:
            argv.extend(["-c", extract_ligand])
        if extract_radius is not None:
            argv.extend(["-r", str(extract_radius)])
        if scan_lists:
            argv.extend(["--scan-lists", scan_lists])
        if max_cycles is not None:
            argv.extend(["--max-cycles", str(max_cycles)])
        if thresh:
            argv.extend(["--thresh", thresh])
        if thresh_post:
            argv.extend(["--thresh-post", thresh_post])
        if backend:
            argv.extend(["-b", backend])
        if precision:
            argv.extend(["--precision", precision])
        if refine_path is not None:
            argv.extend(["--refine-path", "true" if refine_path else "false"])
        if do_tsopt is not None:
            argv.extend(["--tsopt", "true" if do_tsopt else "false"])
        if do_dft is not None:
            argv.extend(["--dft", "true" if do_dft else "false"])
        if do_thermo is not None:
            argv.extend(["--thermo", "true" if do_thermo else "false"])
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def run_single_point_dft(
        input_pdb: str,
        charge: int,
        *,
        ligand_charge: Optional[str] = None,
        multiplicity: Optional[int] = None,
        func_basis: Optional[str] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Single-point DFT energy via gpu4pyscf subprocess (CLI: `pdb2reaction dft`).

        ``func_basis`` is the single combined CLI argument (``FUNC/BASIS``,
        e.g. ``"wb97m-v/def2-tzvpd"``).
        """
        od = _resolve_out_dir(out_dir, "dft")
        argv: list[str] = ["pdb2reaction", "dft", "-i", input_pdb, "-q", str(charge)]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if multiplicity is not None:
            argv.extend(["-m", str(multiplicity)])
        if func_basis:
            argv.extend(["--func-basis", func_basis])
        argv.append("--out-json")
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def run_single_point(
        input_pdb: str,
        charge: int,
        *,
        ligand_charge: Optional[str] = None,
        multiplicity: Optional[int] = None,
        do_hess: Optional[bool] = None,
        hessian_calc_mode: Optional[str] = None,
        backend: Optional[str] = None,
        solvent: Optional[str] = None,
        solvent_model: Optional[str] = None,
        precision: Optional[str] = None,
        workers: Optional[int] = None,
        workers_per_node: Optional[int] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Single-point MLIP energy / forces / optional Hessian (CLI: `pdb2reaction sp`).

        ``do_hess=True`` additionally writes the full Hessian to
        ``hessian.npy``. ``hessian_calc_mode`` selects the Hessian
        backend ("Analytical" only works for UMA; other backends fall
        back to FiniteDifference). NOTE: analytical Hessian raises a
        RuntimeError when ``workers > 1``; pass
        ``hessian_calc_mode="FiniteDifference"`` explicitly in that case.
        """
        od = _resolve_out_dir(out_dir, "sp")
        argv: list[str] = ["pdb2reaction", "sp", "-i", input_pdb, "-q", str(charge)]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if multiplicity is not None:
            argv.extend(["-m", str(multiplicity)])
        if do_hess is not None:
            argv.append("--hess" if do_hess else "--no-hess")
        if hessian_calc_mode:
            argv.extend(["--hessian-calc-mode", hessian_calc_mode])
        argv.extend(_shared_calc_flags(
            backend=backend, solvent=solvent, solvent_model=solvent_model,
            precision=precision, workers=workers, workers_per_node=workers_per_node,
        ))
        argv.append("--out-json")
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    # Structure / I/O helpers (no out_dir, no summary.json)

    @mcp.tool()
    def extract_active_site(
        complex_pdb: str,
        ligand_id: str,
        radius_angstrom: float,
        output_pdb: str,
        *,
        ligand_charge: Optional[str] = None,
        exclude_backbone: Optional[bool] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Extract an active-site model (CLI: `pdb2reaction extract`).

        Cuts a sphere of radius `radius_angstrom` around the ligand `ligand_id`
        (PDB residue 3-letter code). Returns the written PDB path.
        `exclude_backbone=None` (default) defers to the CLI default
        (`--no-exclude-backbone`); pass True / False to override explicitly.
        """
        argv: list[str] = ["pdb2reaction", "extract",
                           "-i", complex_pdb, "-c", ligand_id,
                           "-r", str(radius_angstrom), "-o", output_pdb]
        if ligand_charge:
            argv.extend(["--ligand-charge", ligand_charge])
        if exclude_backbone is True:
            argv.append("--exclude-backbone")
        elif exclude_backbone is False:
            argv.append("--no-exclude-backbone")
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=None, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def add_element_info(
        input_pdb: str,
        output_pdb: str,
        *,
        overwrite: bool = False,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Repair / add PDB element columns (CLI: `pdb2reaction add-elem-info`)."""
        argv: list[str] = ["pdb2reaction", "add-elem-info", "-i", input_pdb, "-o", output_pdb]
        if overwrite:
            argv.append("--overwrite")
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=None, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def fix_altloc(
        input_pdb: str,
        output_pdb: str,
        *,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Resolve PDB alternate locations (CLI: `pdb2reaction fix-altloc`)."""
        argv: list[str] = ["pdb2reaction", "fix-altloc", "-i", input_pdb, "-o", output_pdb]
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=None, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def plot_trajectory(
        input_trj_xyz: str,
        output_png: str,
        *,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Plot an energy profile from a trajectory (CLI: `pdb2reaction trj2fig`)."""
        argv: list[str] = ["pdb2reaction", "trj2fig", "-i", input_trj_xyz, "-o", output_png]
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=None, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def plot_energy_diagram(
        energies: str,
        output_png: str,
        *,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Draw an energy diagram (CLI: `pdb2reaction energy-diagram`).

        `energies`: Python-literal list of floats, e.g. "[0, 12.5, 4.3, 18.7, -1.2]".
        """
        argv: list[str] = ["pdb2reaction", "energy-diagram", "-i", energies, "-o", output_png]
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=None, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def detect_bond_changes(
        reactant_pdb: str,
        product_pdb: str,
        *,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Detect bond changes between two PDB structures (CLI: `pdb2reaction bond-summary`)."""
        argv: list[str] = ["pdb2reaction", "bond-summary", "-i", reactant_pdb, product_pdb]
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=None, timeout=timeout_seconds).to_dict()
