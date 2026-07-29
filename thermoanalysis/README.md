# `thermoanalysis/` (bundled fork)

> **This is a repo-internal fork of [thermoanalysis](https://github.com/eljost/thermoanalysis), NOT the upstream PyPI package. Do not `pip install thermoanalysis` alongside this package — it will silently overwrite the bundled copy.**

The fork is shipped inside this repository; treat it as part of `pdb2reaction` rather than a swappable upstream.

## Why a fork?

The bundled fork carries a small but important divergence from upstream that breaks side-by-side installation:

- **`QCData.py` branding / IO diff** — the bundled file emits output identifiers and quantum-chemistry data format expected by `pdb2reaction/freq.py`. Upstream's signatures do not match.

The rest of the package (`thermo.py`, `config.py`, `constants.py`) is close to upstream but rebuilt against the bundled `pysisyphus/` fork's unit conventions.

## Divergent files (do NOT replace with upstream)

| file | divergence | rule |
|------|------------|------|
| `QCData.py` | branding + IO signature differences | freq stage consumer contract |

## Change policy

Logic edits to `QCData.py` require a demonstrated I/O or numerical need and
verification against the thermochemistry golden tests. Use a
`[CHEMISTRY-RULE:N]` prefix only if a marked chemistry rule is actually
changed; the bundled-fork location alone does not imply that prefix.
Preserve the public symbols and I/O signatures consumed by
`pdb2reaction/workflows/freq.py`.

## Upstream compatibility

If you `pip install thermoanalysis` into the same environment as `pdb2reaction`, Python's import machinery may resolve to the bundled copy or the upstream copy depending on `sys.path` order. The flat-top placement at `<repo-top>/thermoanalysis/` is required for the bundled copy to take precedence in the editable-install path; do not move it under `pdb2reaction/`.

## See also

- [`../docs/architecture.md`](../docs/architecture.md) §5.3, §6 — repo-internal fork policy
- [`../CONTRIBUTING.md`](../CONTRIBUTING.md) §4.3 — validation policy for logic edits in bundled forks
- `THIRD_PARTY_NOTICES.txt` — thermoanalysis upstream attribution
