# Two-tier PEP 562 lazy import infrastructure.
#
# `_LAZY_SYMBOLS` rescues `from pdb2reaction import <symbol>` (a name re-exported
# from a sub-module) and `_LAZY_MODULES` rescues `from pdb2reaction import
# <module>` (a sub-module accessed as a root attribute) without forcing eager
# imports at package import time. Dotted paths
# (`from pdb2reaction.utils import X`) cannot be intercepted by PEP 562 and are
# served by separate path-shim files instead.

from __future__ import annotations

# Release metadata and runtime imports share the committed ``_version.py``
# module; ``pyproject.toml`` reads the same attribute when building artifacts.
try:
    from pdb2reaction._version import __version__, __version_tuple__
except ImportError:
    __version__ = "0.0.0.dev0"
    __version_tuple__ = (0, 0, 0, "dev0")

# Root symbol → owning module path, e.g. {"UMACalculator": "pdb2reaction.backends.uma"}.
_LAZY_SYMBOLS: dict[str, str] = {}

# Root module attribute → real module path. Symbol re-export cannot return a
# module object; this map serves `from pdb2reaction import <module>` directly.
_LAZY_MODULES: dict[str, str] = {}

__all__ = [
    "__version__",
    "__version_tuple__",
    *_LAZY_SYMBOLS.keys(),
]


def __getattr__(name: str):
    import importlib

    if name in _LAZY_SYMBOLS:
        mod = importlib.import_module(_LAZY_SYMBOLS[name])
        value = getattr(mod, name)
        globals()[name] = value  # cache on hot path; subsequent lookups skip __getattr__
        return value
    if name in _LAZY_MODULES:
        module = importlib.import_module(_LAZY_MODULES[name])
        globals()[name] = module
        return module
    raise AttributeError(f"module 'pdb2reaction' has no attribute {name!r}")


def __dir__():
    return sorted(set(__all__) | set(_LAZY_MODULES) | set(globals()))
