"""Shared alternate-location scoring policy for structure readers.

The typed structure bridge and the text-preserving PDB fixer intentionally
render records differently, but they must choose the same deposited
conformer.  This module keeps that choice independent of either renderer.
"""

from __future__ import annotations

from collections.abc import Iterable
import math
from typing import Optional, TypeVar


Label = TypeVar("Label")


def parsed_occupancy(value: object) -> Optional[float]:
    """Return a finite parsed occupancy, or ``None`` for missing/invalid data."""

    if value is None:
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def occupancy_rank(value: Optional[float], first_index: int) -> tuple[int, float, int]:
    """Return the canonical score for one occupancy observation.

    Any parsed value, including zero, ranks above a missing value.  Equal
    parsed values and all-missing cases are resolved by first appearance.
    """

    parsed = parsed_occupancy(value)
    if parsed is None:
        return (0, 0.0, -int(first_index))
    return (1, parsed, -int(first_index))


def choose_altloc_label(
    observations: Iterable[tuple[Label, Optional[float], int]],
) -> Label:
    """Choose one label by parsed mean occupancy, then first appearance.

    ``observations`` contains ``(label, occupancy-or-None, order)`` triples.
    A label with at least one parsed occupancy outranks a missing-only label;
    missing values are excluded from its mean rather than replaced by 1.0.
    """

    stats: dict[Label, list[float | int]] = {}
    for label, occupancy, order in observations:
        if label not in stats:
            stats[label] = [0.0, 0, int(order)]
        parsed = parsed_occupancy(occupancy)
        if parsed is not None:
            stats[label][0] = float(stats[label][0]) + parsed
            stats[label][1] = int(stats[label][1]) + 1

    if not stats:
        raise ValueError("Cannot choose an alternate-location label from no observations.")

    def score(item: tuple[Label, list[float | int]]) -> tuple[int, float, int]:
        _label, (total, count, first_index) = item
        if int(count) == 0:
            return occupancy_rank(None, int(first_index))
        return occupancy_rank(float(total) / int(count), int(first_index))

    return max(stats.items(), key=score)[0]


__all__ = ["choose_altloc_label", "occupancy_rank", "parsed_occupancy"]
