# `pdb2reaction trj2fig`

```text
Usage: pdb2reaction trj2fig [OPTIONS] [EXTRA_OUTS]...

  Plot ΔE or E from an XYZ trajectory and export figure/CSV.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+detailed step logging
                                  and deliverable paths; 3=everything (full
                                  config blocks, per-file paths, DEBUG logging).
                                  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                XYZ trajectory file.  [required]
  -o, --output, --out FILE        Output file(s). You can repeat -o and/or list
                                  extra filenames after options
                                  (.png/.jpg/.jpeg/.html/.svg/.pdf/.csv). If
                                  nothing is given, defaults to energy.png.
  --unit [kcal|hartree]           Target unit for plotted and exported energy
                                  values.  [default: kcal]
  -r, --reference TEXT            Reference: 'init' (initial frame; last frame
                                  if --reverse-x), 'None' (absolute E), or a
                                  zero-based integer frame index.  [default:
                                  init]
  -q, --charge INTEGER            Total charge. Triggers energy recomputation
                                  when supplied.
  -m, --multiplicity INTEGER RANGE
                                  Spin multiplicity (2S+1). Triggers energy
                                  recomputation when supplied.  [default: (1);
                                  x>=1]
  --reverse-x / --no-reverse-x    Reverse the x-axis (last frame on the left).
                                  [default: no-reverse-x]
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend.  [default: uma]
  --backend-model TEXT            Model variant for the selected --backend (e.g.
                                  uma-s-1p2 / uma-m-1p1 for uma,
                                  orb_v3_conservative_omol for orb, MACE-OMOL-0
                                  / off:small for mace).  [default: (the
                                  selected backend's own model)]
  --precision [fp32|fp64]         MLIP backend precision: fp32 or fp64. Unset
                                  defaults per backend (uma: fp32; orb, mace:
                                  fp64). Routed to backend-specific kwargs (UMA
                                  precision / ORB precision / MACE
                                  default_dtype). aimnet2: fp32 no-op; fp64
                                  rejected.  [default: (per backend: uma fp32;
                                  orb, mace fp64)]
  --solvent TEXT                  Computationally expensive xTB solvent delta
                                  correction. Examples: water, methanol,
                                  acetonitrile, dmso, thf, toluene. 'none'
                                  disables it.  [default: none]
  --solvent-model [alpb|cpcmx]    xTB solvent model.  [default: alpb]
  --solvent-xtb-cmd TEXT          Command for the xTB solvent correction,
                                  including optional xTB arguments. If xTB SCC
                                  convergence is poor, increasing --etemp may
                                  help (for example: 'xtb --etemp 1000').
  --out-json / --no-out-json      Write machine-readable result.json to the
                                  output directory.  [default: no-out-json]
  -h, --help                      Show this message and exit.
```
