# `pdb2reaction trj2fig`

```text
Usage: pdb2reaction trj2fig [OPTIONS] [EXTRA_OUTS]...

  Plot ΔE or E from an XYZ trajectory and export figure/CSV.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+detailed step logging
                                  and output paths; 3=everything (full config
                                  blocks, per-file paths, DEBUG logging).
                                  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                XYZ trajectory file.  [required]
  -o, --output, --out FILE        Output file(s)
                                  (.png/.jpg/.jpeg/.html/.svg/.pdf/.csv). You
                                  can repeat -o and/or list extra filenames
                                  after options. Defaults to energy.png.
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
                                  uma-s-1p2p1 / uma-s-1p2p1-ft-d4rev / uma-m-1p1
                                  for uma, orbmol-v2 / orb_v3_conservative_omol
                                  for orb, MACE-OMOL-0 / MACE-POLAR-1-M for
                                  mace, aimnet2 / aimnet2-nse for aimnet2).
                                  [default: (the selected backend's own model)]
  --precision [tf32|fp32|fp64]    MLIP backend precision: tf32, fp32 (IEEE) or
                                  fp64. tf32 uses TF32 matmul for speed (uma,
                                  orb only); avoid it for Hessians. aimnet2:
                                  fp32 only.  [default: (per backend: uma fp32;
                                  orb, mace fp64)]
  --compile / --no-compile        torch.compile the UMA or ORB model: slower
                                  first call, slightly faster afterwards; worth
                                  it only for long runs.  [default: (off)]
  --solvent TEXT                  Optional xTB or AnySolv solvent correction.
                                  'none' disables it. AnySolv MoE is recommended
                                  for accuracy-sensitive reaction calculations;
                                  see docs for model access and setup.
                                  [default: none]
  --solvent-model [alpb|cpcmx|anysolv]
                                  Solvent correction: ALPB/CPCMx or AnySolv. For
                                  AnySolv, MoE is recommended; select its
                                  weights with calc.anysolv_checkpoint (see docs
                                  for model access).  [default: alpb]
  --solvent-xtb-cmd TEXT          Command for the xTB solvent correction,
                                  including optional xTB arguments. If xTB SCC
                                  convergence is poor, increasing --etemp may
                                  help (for example: 'xtb --etemp 1000').
  --out-json / --no-out-json      Write result.json to the output directory.
                                  [default: no-out-json]
  -h, --help                      Show this message and exit.
```
