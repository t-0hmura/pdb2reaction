"""Shared Click option decorators across pdb2reaction subcommands.

Factories define the common
`--charge/-q + --ligand-charge/-l + --multiplicity/-m` option set used by
`freq / irc / opt / tsopt`.

Subcommands not covered by this factory:
  * `dft` — `--multiplicity` help text references .gjf inheritance, so the
    text differs per call site.
  * `all / extract / path_opt / path_search / scan / scan2d / scan3d` — the
    help text or `required=` value of one of the three options differs, so those
    commands define the options locally.
"""

from __future__ import annotations

from typing import Callable, Sequence

import click


def add_print_every_option() -> Callable[[Callable], Callable]:
    """Attach `--print-every N` (debug verbosity throttle).

    Routes to pysisyphus `Optimizer.__init__(print_every=N)`. N=1 (pysis
    default) prints every cycle; larger N prints every N-th macro cycle.
    Useful when running long opts and only the periodic summary is needed.
    Default ``None`` so omission falls through to defaults / YAML.
    """
    def decorator(func: Callable) -> Callable:
        return click.option(
            "--print-every",
            "print_every",
            type=click.IntRange(min=1),
            default=None,
            # Declared None so `cli_param_overridden` can tell an explicit value
            # from an omission; the string reports the effective default, so
            # --help, the generated reference and the Colab option list stop
            # reading as unset.
            show_default="100",
            help="Print optimizer status every N cycles (debug knob).",
        )(func)
    return decorator


def add_irc_pos_def_option() -> Callable[[Callable], Callable]:
    """Attach `--irc-pos-def/--no-irc-pos-def`.

    When enabled, IRC convergence additionally requires `eigvalsh(mw_hessian)[0] > 0`
    on top of rms(grad) <= threshold. Blocks the "shoulder" false-convergence
    where the IRC walker calls success on a downhill descent before reaching
    the local minimum. Default ``None`` falls through to the rms-only criterion.
    """
    def decorator(func: Callable) -> Callable:
        return click.option(
            "--irc-pos-def/--no-irc-pos-def",
            "irc_pos_def",
            default=None,
            show_default="no-irc-pos-def (rms-only criterion)",
            help="Require pos-def Hessian at IRC convergence (blocks shoulder false-convergence).",
        )(func)
    return decorator


def add_coord_type_option(
    choices: Sequence[str] = ("cart", "redund", "dlc", "tric"),
) -> Callable[[Callable], Callable]:
    """Attach `--coord-type` to a Click command.

    Selects the optimization coordinate system passed through to pysisyphus'
    Geometry constructor. ``cart`` is the default. ``redund`` and ``tric`` are
    accepted for single-structure optimizers
    (opt / tsopt / scan / freq) but NOT for Chain-of-States engines —
    ``path-opt`` and ``path-search`` pass ``choices=("cart", "dlc")`` here
    because pysisyphus' ChainOfStates only honours those two coordinate
    systems. Subcommands hard-coupled to Cartesian (irc, dft) skip the
    decorator entirely.

    Default is ``None`` so omission falls through to
    ``GEOM_KW_DEFAULT['coord_type']`` (cart) via the standard YAML override
    chain.
    """
    options_str = "|".join(choices)
    def decorator(func: Callable) -> Callable:
        # NOTE: dest is `cli_coord_type`, NOT `coord_type`, because every
        # downstream cli body already has a local `coord_type` variable
        # (`coord_type = geom_cfg.get("coord_type", ...)` right before the
        # geom_loader call). If we bound the CLI param to `coord_type`,
        # Python's local-scope inference would shadow the closure variable
        # and the assemble-block reference would UnboundLocalError. Keeping
        # the dest unique prevents the collision; the user-facing flag
        # remains `--coord-type`.
        return click.option(
            "--coord-type",
            "cli_coord_type",
            type=click.Choice(list(choices), case_sensitive=False),
            default=None,
            show_default="cart",
            help=(
                f"Optimization coordinate system ({options_str})."
            ),
        )(func)
    return decorator


def add_precision_option() -> Callable[[Callable], Callable]:
    """Attach `--precision fp32|fp64` to a Click command.

    Backend-agnostic precision flag. The CLI body routes the value into
    the backend-specific configuration key via
    ``pdb2reaction.backends.apply_precision_to_calc_cfg``:

    - ``uma``  -> ``precision`` ('fp32' | 'fp64')
    - ``orb``  -> ``precision`` ('float32-high' | 'float64')
    - ``mace`` -> ``default_dtype`` ('float32' | 'float64')
    - ``aimnet2`` -> fp32 is a no-op; fp64 is rejected (its model inputs
      are cast to float32 upstream, so fp64 cannot be honoured)

    Unset resolves per backend (``backends._BACKEND_DEFAULT_PRECISION``):
    UMA fp32 (its upstream fairchem baseline), ORB and MACE fp64. Precision
    affects model evaluation and finite-difference derivatives, so frequency
    and IRC results should be validated for the selected backend and precision.

    Wire targets: every subcommand that constructs a backend calculator
    (opt, tsopt, freq, irc, sp, scan / scan2d / scan3d, path-opt,
    path-search, all).
    """
    def decorator(func: Callable) -> Callable:
        return click.option(
            "--precision",
            "precision",
            type=click.Choice(["fp32", "fp64"], case_sensitive=False),
            default=None,
            show_default="per backend: uma fp32; orb, mace fp64",
            help=(
                "MLIP backend precision: fp32 or fp64. Unset defaults per "
                "backend (uma: fp32; orb, mace: fp64). Routed to "
                "backend-specific kwargs (UMA precision / ORB precision / "
                "MACE default_dtype). aimnet2: fp32 no-op; fp64 rejected."
            ),
        )(func)
    return decorator


def add_backend_model_option() -> Callable[[Callable], Callable]:
    """Attach ``--backend-model NAME`` to a Click command.

    Backend-agnostic model-variant override. The CLI body routes the value into
    the active backend's ``model`` kwarg via
    ``pdb2reaction.backends.apply_backend_model_to_calc_cfg``. Unset keeps the
    backend's built-in default model (uma-s-1p2 / orb_v3_conservative_omol /
    MACE-OMOL-0 / aimnet2). Same wire targets as ``add_precision_option``.
    """
    def decorator(func: Callable) -> Callable:
        return click.option(
            "--backend-model",
            "backend_model",
            type=str,
            default=None,
            show_default="the selected backend's own model",
            help=(
                "Model variant for the selected --backend (e.g. "
                "uma-s-1p2 / uma-m-1p1 for uma, orb_v3_conservative_omol for orb, "
                "MACE-OMOL-0 / off:small for mace). "
                "Default: the backend's built-in model."
            ),
        )(func)
    return decorator


def add_calc_file_option() -> Callable[[Callable], Callable]:
    """Attach ``--calc-file PATH`` (+ ``--calc-file-func-name NAME``) to a Click command.

    When set, pdb2reaction loads an arbitrary ASE Calculator from the user
    Python file and uses it as the energy/gradient backend (overriding
    ``--backend``). The CLI body routes the value via
    ``pdb2reaction.backends.apply_calc_file_to_calc_cfg``, which switches the
    backend to ``custom``. The file must expose a factory
    ``get_calculator(charge, spin, device, **kwargs) -> ase Calculator``
    (rename via ``--calc-file-func-name``). Lets users couple GFN-xTB (tblite),
    DFTB+, ORCA, or any ASE-compatible engine without modifying pdb2reaction.
    Same wire targets as ``add_precision_option``.
    """
    def decorator(func: Callable) -> Callable:
        func = click.option(
            "--calc-file-func-name",
            "calc_factory",
            type=str,
            default=None,
            show_default="get_calculator",
            help=(
                "Name of the callable in --calc-file that returns an ASE "
                "Calculator (or a module-level Calculator instance). "
                "CLI overrides config YAML; otherwise defaults to get_calculator."
            ),
        )(func)
        func = click.option(
            "--calc-file",
            "calc_file",
            type=click.Path(exists=True, dir_okay=False),
            default=None,
            help=(
                "Python file exposing get_calculator(...) -> an ASE Calculator, "
                "used as the energy/gradient backend (overrides --backend). "
                "Couples GFN-xTB / DFTB+ / any ASE engine. See --calc-file-func-name."
            ),
        )(func)
        return func
    return decorator


def _deterministic_callback(ctx, param, value):
    """Eager callback: activate strict-deterministic mode when --deterministic
    is set. Process-global, so it covers every backend used in the run and
    every in-process child stage of ``all``. ``expose_value=False`` keeps it
    out of the command function signature (no body changes needed)."""
    if ctx.resilient_parsing:
        return value
    if value:
        from pdb2reaction.backends._determinism import setup_deterministic
        setup_deterministic()
        return value
    # --no-deterministic cannot switch OFF what the environment turned on: the
    # env var is read independently downstream. Say so, or the user reads the
    # run as non-deterministic when it is not.
    from pdb2reaction.backends._determinism import is_deterministic_requested
    if is_deterministic_requested():
        source = ctx.get_parameter_source(param.name)
        suffix = (
            " despite --no-deterministic"
            if source is click.core.ParameterSource.COMMANDLINE
            else ""
        )
        click.echo(
            "[determinism] NOTE: PDB2REACTION_STRICT_DETERMINISTIC=1 is set in "
            f"the environment, so this run stays strictly deterministic{suffix}. "
            "Unset the variable to disable it.",
            err=True,
        )
    return value


def add_deterministic_option() -> Callable[[Callable], Callable]:
    """Attach ``--deterministic/--no-deterministic`` to a Click command.

    Strict PyTorch determinism: turns on ``torch.use_deterministic_algorithms`` plus
    an ``index_reduce_`` shim to support same-stack repeatability. It is a
    process-global side effect applied via an eager, value-less callback, so it
    propagates to all backends and to the in-process child stages of ``all``
    without per-stage forwarding. Slower than the default, and raises (rather
    than silently degrading) if the torch build cannot honour strict mode.
    Default off. This cannot guarantee identity across package versions or
    hardware, nor control an arbitrary custom ASE calculator. The env var
    ``PDB2REACTION_STRICT_DETERMINISTIC=1`` is the equivalent entry point for
    CI / the direct Python API.
    """
    def decorator(func: Callable) -> Callable:
        return click.option(
            "--deterministic/--no-deterministic",
            default=False,
            show_default=True,
            is_eager=True,
            expose_value=False,
            callback=_deterministic_callback,
            help=(
                "Request strict same-stack PyTorch determinism (deterministic "
                "algorithms + index_reduce_ shim). Slower; raises for detected "
                "unsupported ops; custom calculators are outside its scope."
            ),
        )(func)
    return decorator


def _allow_charge_mult_mismatch_callback(ctx, param, value):
    """Eager callback: when --allow-charge-mult-mismatch is set, disable the cluster
    electron-parity check process-globally (covers every backend + every in-process child
    stage of ``all`` without per-stage forwarding, like --deterministic). ``expose_value=False``
    keeps it out of the command signature."""
    if ctx.resilient_parsing:
        return value
    if value:
        from pdb2reaction.core.utils import set_allow_charge_mult_mismatch
        set_allow_charge_mult_mismatch(True)
    return value


def add_allow_charge_mult_mismatch_option() -> Callable[[Callable], Callable]:
    """Attach ``--allow-charge-mult-mismatch`` to a Click command.

    Skips the cluster charge/multiplicity electron-parity check (``validate_charge_spin``)
    and logs that it was skipped. Every integer-electron spin state obeys the parity
    relation, so an open-shell cluster needs a matching multiplicity rather than this
    option. Reserved for users who intentionally submit a nonstandard electron count.
    Process-global via an eager, value-less callback, so it propagates to every backend
    and child stage without per-stage forwarding.
    """
    def decorator(func: Callable) -> Callable:
        return click.option(
            "--allow-charge-mult-mismatch",
            is_flag=True,
            default=False, show_default=True,
            is_eager=True,
            expose_value=False,
            callback=_allow_charge_mult_mismatch_callback,
            help=(
                "Skip the cluster charge/multiplicity electron-parity check (logs that it was "
                "skipped). Open-shell clusters need a matching multiplicity instead; use this "
                "only for an intentionally nonstandard electron count."
            ),
        )(func)
    return decorator


def add_ml_charge_spin_options(
    *, allow_ref_pdb: bool = True
) -> Callable[[Callable], Callable]:
    """Attach the standard ML region charge/spin triple to a Click command.

    Options: -q/--charge, -l/--ligand-charge, -m/--multiplicity (spin).
    The ligand-charge help mentions ``--ref-pdb`` only on commands that
    actually expose that option.
    """
    options = [
        click.option(
            "-q", "--charge",
            type=int,
            required=False,
            help=(
                "Total charge. Required for non-.gjf inputs unless --ligand-charge is provided "
                "(.gjf templates inherit the charge automatically)."
            ),
        ),
        click.option(
            "-l", "--ligand-charge",
            type=str,
            default=None,
            show_default=False,
            help=(
                "Total charge or per-resname mapping (e.g., GPP:-3,SAM:1) used to derive "
                "charge when -q is omitted "
                + (
                    "(requires PDB/mmCIF input or --ref-pdb)."
                    if allow_ref_pdb
                    else "(requires PDB/mmCIF input)."
                )
            ),
        ),
        click.option(
            "-m", "--multiplicity",
            "spin",
            type=click.IntRange(min=1),
            default=None,
            show_default="1",
            help="Spin multiplicity (2S+1).",
        ),
    ]

    def decorator(func: Callable) -> Callable:
        for opt in reversed(options):
            func = opt(func)
        return func

    return decorator
