# Tests Directory Compatibility

`pdb2reaction` historically stored tests under `test/`.

To align with common tooling and with `mlmm_toolkit` (`tests/`), pytest is now
configured to discover both:

- `tests/`
- `test/`

Keep existing fixtures under `test/` unless there is a concrete need to move
them. New tests may be added in either location during the transition period.
