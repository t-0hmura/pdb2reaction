"""Resolved configuration for the optional PySCF/GPU4PySCF backend."""

from __future__ import annotations

import os
import re
from copy import deepcopy
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Tuple

import click


DFT_DEFAULT_FUNC_BASIS = "wb97m-v/def2-svp"
DFT_CLI_META_KEY = "pdb2reaction.dft_cli"
PCM_DIELECTRIC = {
    "water": 78.3553,
    "methanol": 32.613,
    "ethanol": 24.852,
    "acetonitrile": 35.688,
    "dmso": 46.826,
    "thf": 7.4257,
    "toluene": 2.3741,
}

_PYSCF_SECTIONS = {
    "mol",
    "mf",
    "grids",
    "density_fit",
    "with_df",
    "with_solvent",
}
_MANAGED_PYSCF_FIELDS = {
    "atom",
    "charge",
    "spin",
    "unit",
    "output",
    "stdout",
    "chkfile",
}
_TOP_LEVEL_FIELDS = {
    "func_basis",
    "func",
    "functional",
    "basis",
    "charge",
    "multiplicity",
    "engine",
    "conv_tol",
    "max_cycle",
    "grid_level",
    "verbose",
    "lowmem",
    "density_fit",
    "auxbasis",
    "solvent",
    "solvent_model",
    "save_scf_checkpoint",
    "checkpoint_path",
    "embedcharge",
    "embedcharge_cutoff",
    "nprocs",
    "nprocs_source",
    "memory",
    "memory_mb",
    "memory_source",
    "pyscf",
}

_MEMORY_RE = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*([kmgt]?i?b)?\s*$", re.IGNORECASE)


def parse_memory_mb(value: Any) -> int:
    """Normalize a PySCF host-memory value to decimal MB."""

    if isinstance(value, bool):
        raise click.BadParameter("DFT memory must be a size such as 64000MB or 180GB.")
    match = _MEMORY_RE.match(str(value))
    if match is None:
        raise click.BadParameter("DFT memory must be a size such as 64000MB or 180GB.")
    amount = float(match.group(1))
    unit = (match.group(2) or "mb").lower()
    factors = {
        "kb": 1.0e-3, "kib": 1024.0 / 1.0e6,
        "mb": 1.0, "mib": 1024.0**2 / 1.0e6,
        "gb": 1.0e3, "gib": 1024.0**3 / 1.0e6,
        "tb": 1.0e6, "tib": 1024.0**4 / 1.0e6,
    }
    memory_mb = int(amount * factors[unit])
    if memory_mb < 256:
        raise click.BadParameter("DFT memory must be at least 256 MB.")
    return memory_mb


def _positive_env_int(name: str) -> Optional[int]:
    value = os.environ.get(name)
    if value in (None, ""):
        return None
    try:
        parsed = int(str(value).split("(", 1)[0])
    except ValueError:
        return None
    return parsed if parsed > 0 else None


def _affinity_count() -> Optional[int]:
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return None


def _pbs_slot_count() -> Optional[int]:
    """Return a portable PBS/Torque CPU allocation when one is exported."""

    total = _positive_env_int("PBS_NP")
    if total is not None:
        return total
    nodefile = os.environ.get("PBS_NODEFILE")
    if nodefile:
        try:
            slots = sum(
                1
                for line in Path(nodefile).read_text(encoding="utf-8").splitlines()
                if line.strip()
            )
        except OSError:
            slots = 0
        if slots:
            return slots
    return _positive_env_int("PBS_NUM_PPN")


def _available_memory_mb() -> Optional[int]:
    candidates = []
    slurm_node = os.environ.get("SLURM_MEM_PER_NODE")
    if slurm_node:
        try:
            candidates.append(parse_memory_mb(slurm_node))
        except click.BadParameter:
            pass
    slurm_cpu = os.environ.get("SLURM_MEM_PER_CPU")
    slurm_cpus = _positive_env_int("SLURM_CPUS_PER_TASK")
    if slurm_cpu and slurm_cpus:
        try:
            candidates.append(parse_memory_mb(slurm_cpu) * slurm_cpus)
        except click.BadParameter:
            pass
    try:
        meminfo = Path("/proc/meminfo").read_text(encoding="utf-8")
        match = re.search(r"^MemAvailable:\s+(\d+)\s+kB", meminfo, re.MULTILINE)
        if match:
            candidates.append(int(match.group(1)) // 1000)
    except OSError:
        pass
    for limit_path, current_path in (
        (Path("/sys/fs/cgroup/memory.max"), Path("/sys/fs/cgroup/memory.current")),
        (Path("/sys/fs/cgroup/memory/memory.limit_in_bytes"), Path("/sys/fs/cgroup/memory/memory.usage_in_bytes")),
    ):
        try:
            limit_text = limit_path.read_text(encoding="utf-8").strip()
            if limit_text != "max":
                remaining = int(limit_text) - int(current_path.read_text(encoding="utf-8").strip())
                if remaining > 0:
                    candidates.append(remaining // 1_000_000)
        except (OSError, ValueError):
            pass
    return min((value for value in candidates if value > 0), default=None)


def resolve_dft_resources(raw: Mapping[str, Any]) -> Tuple[int, str, Optional[int], str]:
    """Resolve PySCF threads and host RAM without treating them as PES identity."""

    affinity = _affinity_count()
    explicit_nprocs = raw.get("nprocs")
    if explicit_nprocs not in (None, "", "auto"):
        try:
            nprocs = int(explicit_nprocs)
        except (TypeError, ValueError) as exc:
            raise click.BadParameter("DFT nprocs must be a positive integer.") from exc
        if nprocs < 1:
            raise click.BadParameter("DFT nprocs must be >= 1.")
        source = str(raw.get("nprocs_source") or "explicit")
    else:
        nprocs = _positive_env_int("OMP_NUM_THREADS") or _positive_env_int(
            "SLURM_CPUS_PER_TASK"
        ) or _pbs_slot_count() or _positive_env_int("NSLOTS") or affinity or os.cpu_count() or 1
        source = "environment"
    if affinity is not None:
        nprocs = min(nprocs, affinity)

    explicit_memory = raw.get("memory", raw.get("memory_mb"))
    if explicit_memory not in (None, "", "auto"):
        memory_mb = parse_memory_mb(explicit_memory)
        memory_source = str(raw.get("memory_source") or "explicit")
    elif os.environ.get("PYSCF_MAX_MEMORY"):
        memory_mb = parse_memory_mb(os.environ["PYSCF_MAX_MEMORY"])
        memory_source = "PYSCF_MAX_MEMORY"
    elif os.environ.get("PYSCF_CONFIG_FILE") or (Path.home() / ".pyscf_conf.py").is_file():
        memory_mb = None
        memory_source = "pyscf-config"
    else:
        available_mb = _available_memory_mb()
        memory_mb = None if available_mb is None else max(256, int(available_mb * 0.8))
        memory_source = "environment" if memory_mb is not None else "pyscf-default"
    return int(nprocs), source, memory_mb, memory_source


def parse_func_basis(value: str) -> Tuple[str, str]:
    """Parse ``FUNCTIONAL/BASIS`` without importing PySCF."""

    text = str(value or "").strip()
    if "/" not in text:
        raise click.BadParameter(
            "Expected FUNCTIONAL/BASIS, for example 'pbe0/def2-svp'.",
            param_hint="--func-basis",
        )
    functional, basis = (part.strip() for part in text.split("/", 1))
    if not functional or not basis:
        raise click.BadParameter(
            "Functional and basis must both be non-empty.",
            param_hint="--func-basis",
        )
    return functional, basis


def _mapping(value: Any, path: str) -> Dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise click.BadParameter(f"{path} must be a mapping.")
    return deepcopy(dict(value))


def _same_value(left: Any, right: Any) -> bool:
    return left == right or str(left) == str(right)


def _strict_bool(value: Any, path: str) -> bool:
    if not isinstance(value, bool):
        raise click.BadParameter(f"{path} must be true or false.")
    return value


@dataclass(frozen=True)
class DFTSettings:
    functional: str
    basis: str
    engine: str = "gpu"
    charge: int = 0
    multiplicity: int = 1
    conv_tol: float = 1.0e-9
    max_cycle: int = 100
    grid_level: int = 3
    verbose: int = 0
    lowmem: bool = True
    density_fit: bool = False
    auxbasis: Optional[str] = None
    solvent: str = "none"
    solvent_model: str = "none"
    save_scf_checkpoint: bool = False
    checkpoint_path: Optional[str] = None
    embedcharge: bool = False
    embedcharge_cutoff: Optional[float] = 12.0
    nprocs: int = 1
    nprocs_source: str = "environment"
    memory_mb: Optional[int] = None
    memory_source: str = "pyscf-default"
    pyscf: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    @property
    def func_basis(self) -> str:
        return f"{self.functional}/{self.basis}"

    @property
    def is_hf(self) -> bool:
        return self.functional.strip().casefold() in {"hf", "rhf", "uhf"}

    # CHEMISTRY-RULE:4 The low-memory GPU driver is valid only for closed-shell RKS.
    @property
    def use_rks_lowmem(self) -> bool:
        """Whether the effective method can use GPU4PySCF's closed-shell driver."""
        return self.lowmem and self.engine == "gpu" and self.multiplicity == 1

    @property
    def memory_mode(self) -> str:
        if self.use_rks_lowmem:
            return "gpu4pyscf_rks_lowmem"
        if self.lowmem:
            return "direct_jk"
        if self.density_fit:
            return "density_fit"
        return "standard_direct"

    def to_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        data["func_basis"] = self.func_basis
        return data

    def scientific_identity(self) -> Dict[str, Any]:
        """Return settings that can change the electronic Hamiltonian/PES."""

        scientific_pyscf = deepcopy(self.pyscf)
        mol_cfg = scientific_pyscf.get("mol")
        if isinstance(mol_cfg, dict):
            mol_cfg.pop("max_memory", None)
            if not mol_cfg:
                scientific_pyscf.pop("mol", None)
        return {
            "functional": self.functional,
            "basis": self.basis,
            "charge": self.charge,
            "multiplicity": self.multiplicity,
            "conv_tol": self.conv_tol,
            "max_cycle": self.max_cycle,
            "grid_level": self.grid_level,
            "lowmem": self.lowmem,
            "density_fit": self.density_fit,
            "auxbasis": self.auxbasis,
            "solvent": self.solvent,
            "solvent_model": self.solvent_model,
            "embedcharge": self.embedcharge,
            "pyscf": scientific_pyscf,
        }


def resolve_dft_settings(
    calc_cfg: Mapping[str, Any],
    *,
    cli_values: Optional[Mapping[str, Any]] = None,
    allow_embedding: bool = False,
) -> DFTSettings:
    """Resolve ``calc.dft`` and explicit flat CLI values once."""

    nested = _mapping(calc_cfg.get("dft"), "calc.dft")
    serialized = calc_cfg.get("dft_settings")
    if serialized is not None:
        serialized_map = _mapping(serialized, "calc.dft_settings")
        nested = {**serialized_map, **nested}

    unknown = sorted(set(nested) - _TOP_LEVEL_FIELDS)
    if unknown:
        raise click.BadParameter(
            "Unknown calc.dft key(s): " + ", ".join(f"calc.dft.{key}" for key in unknown)
        )

    explicit = {key: value for key, value in dict(cli_values or {}).items() if value is not None}
    raw = {**nested, **explicit}

    func_basis = raw.get("func_basis")
    functional = raw.get("functional", raw.get("func"))
    basis = raw.get("basis")
    if func_basis is not None:
        parsed_functional, parsed_basis = parse_func_basis(str(func_basis))
        if functional is not None and not _same_value(functional, parsed_functional):
            raise click.BadParameter("calc.dft.func_basis conflicts with calc.dft.functional.")
        if basis is not None and not _same_value(basis, parsed_basis):
            raise click.BadParameter("calc.dft.func_basis conflicts with calc.dft.basis.")
        functional, basis = parsed_functional, parsed_basis
    elif functional is None and basis is None:
        functional, basis = parse_func_basis(DFT_DEFAULT_FUNC_BASIS)
    elif functional is None or basis is None:
        raise click.BadParameter(
            "calc.dft.functional and calc.dft.basis must be supplied together."
        )

    engine = str(raw.get("engine", "gpu")).strip().lower()
    if engine not in {"cpu", "gpu"}:
        raise click.BadParameter("DFT engine must be 'cpu' or 'gpu'.", param_hint="--engine")

    pyscf_cfg = _mapping(raw.get("pyscf"), "calc.dft.pyscf")
    bad_sections = sorted(set(pyscf_cfg) - _PYSCF_SECTIONS)
    if bad_sections:
        raise click.BadParameter(
            "Unknown calc.dft.pyscf section(s): " + ", ".join(bad_sections)
        )
    for section, values in tuple(pyscf_cfg.items()):
        if section == "density_fit" and isinstance(values, bool):
            pyscf_cfg[section] = {"enabled": bool(values)}
        else:
            pyscf_cfg[section] = _mapping(values, f"calc.dft.pyscf.{section}")
        managed = sorted(_MANAGED_PYSCF_FIELDS & set(pyscf_cfg[section]))
        if managed:
            raise click.BadParameter(
                f"calc.dft.pyscf.{section} cannot set workflow-owned field(s): "
                + ", ".join(managed)
            )

    convenience = {
        "conv_tol": ("mf", "conv_tol", 1.0e-9, float),
        "max_cycle": ("mf", "max_cycle", 100, int),
        "grid_level": ("grids", "level", 3, int),
        "verbose": ("mol", "verbose", 0, int),
    }
    resolved_simple: Dict[str, Any] = {}
    for name, (section, attribute, default, cast) in convenience.items():
        direct = raw.get(name)
        detailed = pyscf_cfg.get(section, {}).get(attribute)
        if direct is not None and detailed is not None and not _same_value(direct, detailed):
            raise click.BadParameter(
                f"calc.dft.{name} conflicts with calc.dft.pyscf.{section}.{attribute}."
            )
        resolved_simple[name] = cast(direct if direct is not None else detailed if detailed is not None else default)

    lowmem = _strict_bool(raw.get("lowmem", True), "calc.dft.lowmem")
    density_cfg = pyscf_cfg.get("density_fit", {})
    direct_density = raw.get("density_fit")
    if isinstance(direct_density, Mapping):
        density_cfg = {**density_cfg, **dict(direct_density)}
        direct_density = None
    detailed_density = density_cfg.get("enabled")
    if direct_density is not None:
        direct_density = _strict_bool(
            direct_density, "calc.dft.density_fit"
        )
    if detailed_density is not None:
        detailed_density = _strict_bool(
            detailed_density, "calc.dft.pyscf.density_fit.enabled"
        )
    if (
        direct_density is not None
        and detailed_density is not None
        and direct_density != detailed_density
    ):
        raise click.BadParameter(
            "calc.dft.density_fit conflicts with calc.dft.pyscf.density_fit.enabled."
        )
    density_fit = bool(
        direct_density
        if direct_density is not None
        else density_cfg.get("enabled", not lowmem)
    )
    if lowmem and density_fit:
        raise click.BadParameter(
            "calc.dft.lowmem=true conflicts with density fitting. Use "
            "--no-lowmem or set calc.dft.density_fit=false."
        )
    direct_auxbasis = raw.get("auxbasis")
    detailed_auxbasis = density_cfg.get("auxbasis")
    if (
        direct_auxbasis is not None
        and detailed_auxbasis is not None
        and str(direct_auxbasis) != str(detailed_auxbasis)
    ):
        raise click.BadParameter(
            "calc.dft.auxbasis conflicts with calc.dft.pyscf.density_fit.auxbasis."
        )
    auxbasis = direct_auxbasis if direct_auxbasis is not None else detailed_auxbasis
    if auxbasis is not None:
        auxbasis = str(auxbasis)

    direct_solvent = raw.get("solvent")
    detailed_solvent = pyscf_cfg.get("with_solvent", {}).get("solvent")
    if (
        direct_solvent is not None
        and detailed_solvent is not None
        and str(direct_solvent).strip().casefold()
        != str(detailed_solvent).strip().casefold()
    ):
        raise click.BadParameter(
            "calc.dft.solvent conflicts with calc.dft.pyscf.with_solvent.solvent."
        )
    solvent = str(
        direct_solvent if direct_solvent is not None else detailed_solvent or "none"
    ).strip()
    if solvent.casefold() in {"", "none", "gas", "vacuum"}:
        solvent = "none"
        solvent_model = "none"
    else:
        solvent_model = str(raw.get("solvent_model", "smd") or "smd").strip().lower()
        if solvent_model == "auto":
            solvent_model = "smd"
        if solvent_model not in {"pcm", "smd"}:
            raise click.BadParameter(
                "DFT solvent_model must be 'pcm' or 'smd'.",
                param_hint="--solvent-model",
            )

    embedcharge = _strict_bool(
        raw.get("embedcharge", False), "calc.dft.embedcharge"
    )
    if embedcharge and not allow_embedding:
        raise click.BadParameter("Electrostatic embedding is available only in mlmm-toolkit.")

    max_cycle = int(resolved_simple["max_cycle"])
    if max_cycle < 1:
        raise click.BadParameter("calc.dft.max_cycle must be >= 1.")
    conv_tol = float(resolved_simple["conv_tol"])
    if conv_tol <= 0:
        raise click.BadParameter("calc.dft.conv_tol must be positive.")

    save_checkpoint = _strict_bool(
        raw.get("save_scf_checkpoint", False),
        "calc.dft.save_scf_checkpoint",
    )
    checkpoint_path = raw.get("checkpoint_path")
    detailed_memory = pyscf_cfg.get("mol", {}).get("max_memory")
    direct_memory = raw.get("memory", raw.get("memory_mb"))
    if direct_memory not in (None, "", "auto") and detailed_memory is not None:
        if parse_memory_mb(direct_memory) != int(detailed_memory):
            raise click.BadParameter(
                "calc.dft.memory conflicts with calc.dft.pyscf.mol.max_memory."
            )
    resource_raw = dict(raw)
    if direct_memory in (None, "", "auto") and detailed_memory is not None:
        resource_raw["memory_mb"] = int(detailed_memory)
    nprocs, nprocs_source, memory_mb, memory_source = resolve_dft_resources(resource_raw)

    return DFTSettings(
        functional=str(functional),
        basis=str(basis),
        engine=engine,
        charge=int(calc_cfg.get("charge", raw.get("charge", 0))),
        multiplicity=int(calc_cfg.get("spin", raw.get("multiplicity", 1))),
        conv_tol=conv_tol,
        max_cycle=max_cycle,
        grid_level=int(resolved_simple["grid_level"]),
        verbose=int(resolved_simple["verbose"]),
        lowmem=lowmem,
        density_fit=density_fit,
        auxbasis=auxbasis,
        solvent=solvent,
        solvent_model=solvent_model,
        save_scf_checkpoint=save_checkpoint,
        checkpoint_path=None if checkpoint_path in (None, "") else str(checkpoint_path),
        embedcharge=embedcharge,
        embedcharge_cutoff=(
            None
            if raw.get("embedcharge_cutoff", 12.0) is None
            else float(raw.get("embedcharge_cutoff", 12.0))
        ),
        nprocs=nprocs,
        nprocs_source=nprocs_source,
        memory_mb=memory_mb,
        memory_source=memory_source,
        pyscf=pyscf_cfg,
    )


def finalize_dft_calculator_config(
    ctx: click.Context,
    calc_cfg: Dict[str, Any],
    *,
    output_dir: Optional[Any] = None,
) -> None:
    """Attach one serialized :class:`DFTSettings` to a resolved calc config."""

    backend = str(calc_cfg.get("backend", "uma")).strip().lower()
    cli_values = dict(ctx.meta.get(DFT_CLI_META_KEY, {}))
    nested_present = bool(calc_cfg.get("dft"))
    if backend != "dft":
        if cli_values or nested_present:
            raise click.BadParameter("DFT-only options require --backend dft.")
        solvent = str(calc_cfg.get("solvent", "none") or "none").strip().casefold()
        solvent_model = str(calc_cfg.get("solvent_model", "alpb") or "alpb").strip().casefold()
        if solvent not in {"", "none", "gas", "vacuum"} and solvent_model in {"pcm", "smd"}:
            raise click.BadParameter(
                "PCM/SMD solvent models require --backend dft; MLIP backends use ALPB/CPCMx."
            )
        return

    for name in ("solvent", "solvent_model"):
        try:
            source = ctx.get_parameter_source(name)
        except Exception:
            source = None
        if source not in (None, click.core.ParameterSource.DEFAULT):
            cli_values[name] = ctx.params.get(name)

    settings = resolve_dft_settings(calc_cfg, cli_values=cli_values)
    if (
        settings.save_scf_checkpoint
        and settings.checkpoint_path is None
        and ctx.info_name != "all"
        and (output_dir is not None or ctx.params.get("out_dir") is not None)
    ):
        settings = replace(
            settings,
            checkpoint_path=str(
                Path(
                    output_dir
                    if output_dir is not None
                    else ctx.params["out_dir"]
                )
                / "_work"
                / "dft_scf"
                / "state.chk"
            ),
        )
    calc_cfg.pop("dft", None)
    calc_cfg["dft_settings"] = settings.to_dict()
    for key in (
        "model",
        "task_name",
        "precision",
        "default_dtype",
        "compile_model",
        "workers",
        "workers_per_node",
    ):
        calc_cfg.pop(key, None)


__all__ = [
    "DFTSettings",
    "DFT_CLI_META_KEY",
    "DFT_DEFAULT_FUNC_BASIS",
    "finalize_dft_calculator_config",
    "parse_func_basis",
    "parse_memory_mb",
    "resolve_dft_resources",
    "PCM_DIELECTRIC",
    "resolve_dft_settings",
]
