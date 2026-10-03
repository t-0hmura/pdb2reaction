# `pdb2reaction all`

`all` runs extraction (with `-c`), the MEP search or a staged scan, and, when
asked, TS optimization, IRC, frequencies, and DFT in one command. It succeeded
when the console prints `[tsopt] Converged (n_imag=1).` for each TS and
`Scientific status: success` under the last `====== Pipeline summary ======`.

## When to use

Use `all` when one job should produce R, TS, P (and IM) structures and barrier
candidates for one or more path segments. The MEP stage runs single-pass
`path-opt` by default; `--refine-path` runs the recursive `path-search`
instead. Frequency and IRC checks decide which candidates are validated
elementary steps. To inspect each stage before the next, run the stages one by
one ([overview](../pdb2reaction-overview/SKILL.md#run-stage-by-stage-and-judge-each-stage)).

## Pick the mode

| Input | Mode (`Pipeline mode` in `summary.log`) | Page |
|---|---|---|
| Two or more structures in reaction order | Multi-structure MEP search (`MEP`) | [all-endpoint-mep.md](all-endpoint-mep.md) |
| One structure with `-s` | Single structure + scan (`Scan`) | [all-scan-list.md](all-scan-list.md) |
| One structure with `--tsopt` and no `-s` | TS-only mode (`TS-only`) | [all-ts-only.md](all-ts-only.md) |

One structure without `-s` or `--tsopt` stops with `BadParameter` ("Provide at
least two structures with -i/--input in reaction order, or use a single
structure with --scan-lists, or a single structure with --tsopt."). `-s` with
two or more structures also stops with an error. One structure with both `-s`
and `--tsopt` runs the scan mode.

## Minimal run

```bash
pdb2reaction all -i <inputs> [-c <centers>] [-l 'RES:Q,...'] [-s '...'] \
    [--tsopt] [--thermo] [--dft] [-b uma|orb|mace|aimnet2|dft] [-o result_all/]
```

Add `--dry-run` first to check the inputs, the charge, and the planned stages
without a calculation. Each mode page has a complete command.

## Judge success

- **Console**: `[tsopt] Converged (n_imag=1).` per TS, and `Scientific status: success` under the last `====== Pipeline summary ======`. `[time] Elapsed` only marks the end of the run.
- **`summary.json`**: `scientific_status` is `success`, `partial`, or `failed`, with `scientific_status_reasons`. `success` means every requested stage converged and, with `--tsopt`, every TS has n_imag = 1. A TS with n_imag ≥ 2 gives `partial`; n_imag = 0 stops before IRC and is never `success`.
- **TS**: a successful TS optimization gives one imaginary mode along the reaction coordinate. Read `post_segments[].tsopt.n_imaginary_modes` and `.imaginary_frequencies_cm`, and play `segments/seg_NN/ts/vib/imag_*_trj.xyz` to see that the mode moves the bonds that form or break.
- **Endpoints**: whether they are the intended R and P is for you to check. Compare `segments/seg_NN/reactant.*` and `product.*`, and the bond changes in section [2] of `summary.log`, with the intended states. Even if the IRC does not converge, the result is usable when the endpoint optimizations reach the intended R and P.
- **Barriers**: `segments[].barrier_kcal` is the barrier on the MEP before TS optimization (TS − R in TS-only mode). After `--tsopt`, read `post_segments[].mlip.barrier_kcal`; `gibbs_mlip` (`--thermo`), `dft` (`--dft`), and `gibbs_dft_mlip` (both) carry the same keys. `rate_limiting_step` is the highest local barrier at the highest method available for every segment, not a microkinetic assignment.

```python
import json
d = json.load(open("result_all/summary.json"))
print(d["execution_status"], d["scientific_status"], d.get("scientific_status_reasons"))
print(d["charge"], d["spin"], d["rate_limiting_step"])
for seg in d["segments"]:
    print(seg["index"], seg["kind"], seg["barrier_kcal"], seg["delta_kcal"])
for post in d.get("post_segments", []):   # match to segments by "index"
    print(post["index"], post.get("tsopt", {}).get("n_imaginary_modes"))
```

`segments` holds the MEP records; requested post-processing is in
`post_segments`. Match the two by `index`: list positions differ when a segment
was skipped or `segments` has a `kind` other than `"seg"`. Key details:
[outputs](../pdb2reaction-overview/outputs.md#per-segment-keys).

## Resume a failed segment

Repeat the original command with the same inputs, extraction, path, and
calculator options and the same `--out-dir`, add `--resume-segment N`, and
change only post-processing options (for example `--tsopt-max-cycles`). The
saved MEP is checked, earlier segments are kept, and segment N onward, the
summary, and the diagrams are written again. The run stops with an error when
the saved inputs and MEP do not match the command.

## Pitfalls and recovery

- **`-s` given more than once.** Use one `-s` occurrence; each following Python literal is one sequential stage. Repeating the flag is rejected. Quote each literal with single quotes outside and double quotes inside.
- **Charge not checked before a long job.** `--dry-run` validates inputs and prints the plan, then exits before scan, MEP, TSOPT, IRC, freq, and DFT. With `-c/--center`, it runs extraction in a temporary directory to validate the derived charge and electron parity, then deletes it.
- **`--refine-path` (off by default).** Refinement can improve a poor TS seed but may split a bad path into unnecessary segments and greatly increase cost. Extra segments are candidates, not proof of hidden intermediates.
- **Cycle limits.** `--max-cycles-gsm` and `--dmf-max-iterations` bound only the MEP stage; scan, TS, IRC, freq, and DFT have their own limits (for example `--tsopt-max-cycles`, `--irc-max-cycles`).
- **Status is not `success`.** Read `scientific_status_reasons`, then the matching block of `summary.log`. Stage `result.json` files are in `segments/seg_NN/ts/` and `irc/`, and under `_work/` for the scan and `path-opt`; freq and DFT write none.
- **A segment directory exists but the stage failed.** `segments/seg_NN/` is created when post-processing starts and can be partial; check `summary.json` and the stage `result.json`, not the directory.
- **TS stops before IRC.** IRC starts only when the TS optimization converged, its final Hessian was computed, and n_imag ≥ 1. With n_imag ≥ 2, IRC runs with a warning along the mode closest to the MEP direction; it is a diagnostic, not a first-order TS. With n_imag = 0, a cycle limit (no final Hessian, so no n_imag), a plateau stop (`--stop-plateau`; the Hessian still gives n_imag), or a skipped or failed final Hessian, the run stops before IRC, keeps the TS files in `segments/seg_NN/ts/`, and does not post-process later segments. Next: [Wrong n_imag after tsopt](../pdb2reaction-overview/ts-strategy.md#wrong-n_imag-after-tsopt) and [When the TS does not come out](../pdb2reaction-overview/ts-strategy.md#when-the-ts-does-not-come-out).
- **`--dft` with `-b dft`.** The run stops at startup. Run `pdb2reaction dft` or `pdb2reaction sp -b dft` as a separate job.
- **UMA with `--uma-workers` above 1 and an explicit `Analytical` Hessian.** This raises `BackendError`; use one worker or `FiniteDifference`.
- **`_work/`.** It holds intermediate files, including the TS candidates (HEI); keep it while you use them.

## Outputs

Cite `segments/seg_NN/reactant.*`, `ts.*`, and `product.*` (written with
`--tsopt`). `summary.json` and `summary.log` sit at the top of `--out-dir` once
the run reaches its summary; an early input error leaves neither. The top level
also has `mep_trj.xyz` (and `.pdb`), `energy_diagram_MEP.png`, and
`energy_diagram_*_all.png`. Each `seg_NN/` has `ts/`, `irc/`,
`freq/{R,TS,P}/` (`--thermo`), and `dft/{R,TS,P}/` (`--dft`). `_work/` has
`models/` (with `-c`), `scan/` (with `-s`), and `path_opt/` (`path_search/` with
`--refine-path`) with the TS candidates `hei_seg_NN.*` and MEP scratch
`seg_NNN_<tag>/`. Full tree: [Output tree](../pdb2reaction-overview/outputs.md#output-tree).

## Next step

- Mode pages: [all-endpoint-mep.md](all-endpoint-mep.md), [all-scan-list.md](all-scan-list.md), [all-ts-only.md](all-ts-only.md).
- The stages it runs: [extract.md](extract.md), [path.md](path.md), [tsopt.md](tsopt.md), [irc.md](irc.md), [freq.md](freq.md), [dft.md](dft.md).
- [Reading outputs](../pdb2reaction-overview/outputs.md): `summary.json` keys and R/TS/P paths.
- Defaults (`OUT_DIR_ALL` and the per-stage `*_KW`): [Where flags and defaults live](SKILL.md#where-flags-and-defaults-live).
