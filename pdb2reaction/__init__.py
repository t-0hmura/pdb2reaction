# pdb2reaction/__init__.py

try:
    from pdb2reaction._version import __version__, __version_tuple__
except ImportError:
    __version__ = "0.0.0.dev0"
    __version_tuple__ = (0, 0, 0, "dev0")

__all__ = [
    "__version__",
    "__version_tuple__",
]
