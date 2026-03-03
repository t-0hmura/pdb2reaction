# pdb2reaction/uma_pysis.py

"""
Backward-compatible shim — delegates to ``pdb2reaction.backends.uma``.

Existing code using ``from pdb2reaction.uma_pysis import uma_pysis`` will
continue to work unchanged.
"""

from .backends.uma import UMACalculator as uma_pysis, UMAASECalculator as uma_ase  # noqa: F401
from .defaults import GEOM_KW_DEFAULT, UMA_CALC_KW  # noqa: F401

# Alias for code that does ``from .uma_pysis import CALC_KW``
CALC_KW = UMA_CALC_KW
