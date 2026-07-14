# xTB / ALPB solvent layer (xtb.md)

`pdb2reaction` can apply an xTB implicit-solvent **delta correction** to
its MLIP/custom calculator. It evaluates xTB with and without implicit solvent
on the same coordinates, then adds their difference to the base calculator.

## When to use

| Use it when | Skip it when |
|---|---|
| You need a quick, approximate solvent correction for a bare cluster | Full QM/MM embedding or explicit waters (out of scope) |
| Comparing barrier heights between gas-phase and solvated approximation | Reporting absolute solvent free energies (xTB-ALPB is empirical) |
| Substrate is in a common solvent (water, DMSO, methanol, …) | Solvent isn't in the ALPB parameter set (unusual organics) |

## Install

`pdb2reaction` invokes the `xtb` binary via `subprocess.run`; it does
not import any Python xtb package. The supported install path is the
conda-forge binary:

```bash
conda install -c conda-forge xtb
```

(Verify with `xtb --version`.) `pip install xtb` provides Python
bindings that `pdb2reaction` does not use.

**Site-installed binary:**

If your site already has `xtb` as a module or in `$PATH`:

```bash
xtb --version
which xtb
```

`pdb2reaction` invokes `xtb` via subprocess; make sure it is on `$PATH`
in any PBS / SLURM job that needs it.

## CLI usage

The solvent correction wraps any other calculator transparently:

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    --solvent water        # ALPB water on top of MLIP
```

Available solvents (pass to `--solvent`): the keyword is forwarded
verbatim to xTB's ALPB parameter set, so use the spelling xTB
recognises. Common entries include `water`, `methanol`, `acetone`,
`acetonitrile`, `dmso`, `dmf`, `chcl3`, `ch2cl2`, `hexane`,
`benzene`, `toluene`, `thf`, `nhexan`, `phenol`, `octanol`,
`woctanol`, `aniline`, `furane`, `ether`, `noctane`, `co2`. The
authoritative list comes from xTB's ALPB parameter set — run
`xtb --help` and consult the xTB docs for the version you have
installed (the set has expanded over xTB releases).

To turn off: omit `--solvent` or pass `--solvent none`.

`--solvent-model` selects `alpb` (default, conda-forge binary) or
`cpcmx` (requires a source build with CPCM-X enabled).

## How it composes with the base calculator

The corrected energy at each step is:

```
E_total = E_base + [E_xTB(solvent) - E_xTB(vacuum)]
```

The force and Hessian corrections use the same subtraction. Thus optimizers,
IRC, and frequency calculations see a position-consistent corrected surface;
this is not the full xTB energy added on top of MLIP. Each correction requires
both vacuum and solvent xTB evaluations, so it materially increases cost.

This wrapper is used by MLIP/custom-calculator workflows. The standalone
`pdb2reaction dft` command has no `--solvent` option and is not automatically
combined with this correction.

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `xtb` binary not on PATH | `pdb2reaction` calls `xtb` via `subprocess.run`; install via `conda install -c conda-forge xtb` (binary; no Python bindings needed). |
| `xtb` binary too old for current ALPB parameter set | ALPB parameter set has been updated across xtb releases; upgrade if results disagree with documentation. |
| Different barrier vs literature | Could be `--solvent` mismatch or the literature used a different solvation model (CPCM / SMD). State the model in any comparison. |

## See also

- `core.md` — `pdb2reaction` install (xTB layer is plumbed via the
  bundled `backends/xtb_alpb_correction.py` and `backends/solvent.py`).
- `dft.md` — the standalone DFT command is a separate refinement and does not
  accept this xTB correction.
- `pdb2reaction-cli/SKILL.md` — `--solvent` is accepted by `all`,
  `tsopt`, `freq`, `irc`, `opt`, `sp`, `path-search`, `path-opt`, `scan`,
  `scan2d`, and `scan3d`. `trj2fig` also accepts it, but only uses it when
  `-q` and/or `-m` triggers MLIP recomputation of every frame; comment-energy
  plotting does not run a calculator. It is not accepted by `dft`, `extract`,
  or the structure-only utilities.
