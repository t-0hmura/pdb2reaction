# YAML 設定リファレンス

YAML 設定ファイル（`--config`）に書けるキーと既定値を、セクションごとに引くページです。セクションの一覧、優先順位、CLI フラグと YAML キーの対応も最初にまとめています。

| セクション | 説明 | 使用されるコマンド |
|---------|-------------|---------|
| [`geom`](#geom) | ジオメトリと座標設定 | opt, scan, scan2d, scan3d, tsopt, freq, irc, path-opt, path-search, dft, sp |
| [`calc`](#calc) | 機械学習原子間ポテンシャル（MLIP）のバックエンドの設定 | opt, scan, scan2d, scan3d, tsopt, freq, irc, path-opt, path-search, sp, dft（`charge` と `spin` のみ） |
| [`opt`](#opt) | 最適化の共通設定 | opt, scan, scan2d, scan3d, tsopt, path-opt, path-search |
| [`lbfgs`](#lbfgs) | L-BFGS の設定 | opt, scan, scan2d, scan3d, path-search, path-opt |
| [`rfo`](#rfo) | RFO の設定 | opt, scan, scan2d, scan3d, path-search, path-opt |
| [`gs`](#gs) | GSM（Growing String Method）のストリング: ノード数、climbing image、再パラメータ化 | path-opt, path-search |
| [`dmf`](#dmf) | DMF（Direct Max Flux）設定 | path-opt, path-search |
| [`stopt`](#stopt) | GSM のストリングを動かす StringOptimizer。`thresh` と `max_cycles` が GSM の止まり方を決める | path-opt, path-search |
| [`irc`](#ja-irc-section) | IRC 積分設定 | irc |
| [`freq`](#ja-freq-section) | 振動解析設定 | freq（`zero_cutoff_cm` は opt、tsopt も） |
| [`thermo`](#thermo) | 熱化学設定 | freq |
| [`dft`](#ja-dft-section) | `dft` コマンドの DFT 計算設定 | dft |
| [`bias`](#bias) | 調和バイアス設定 | scan, scan2d, scan3d |
| [`bond`](#bond) | 結合変化検出設定 | scan, path-search |
| [`search`](#search) | 再帰的経路探索設定 | path-search |
| [`hessian_dimer`](#hessian_dimer) | Hessian Guided Dimer TS 最適化 | tsopt |
| [`rsirfo`](#rsirfo) | RS-P-RFO / RS-I-RFO TS 最適化 | tsopt |
| `sp` | 一点計算の設定（`hess`（既定 `false`）・`hessian_calc_mode`。`--hess`・`--hessian-calc-mode` と同じ）。[sp](sp.md) を参照 | sp |

(ja-yaml-configuration-precedence)=
## 設定の優先順位

設定は以下の順序で適用されます（後のものが前のものを上書き）:

```
組み込みデフォルト  <  --config (YAML)  <  CLI フラグ
```

1. **組み込みデフォルト** — `pdb2reaction <subcmd> --help-advanced` と [コマンドリファレンス（英語のみ）](../reference/commands/index.md) の `[default: …]` に出る値。
2. **`--config`** — デフォルトを上書きする YAML ファイル（例: `--config my_settings.yaml`）。
3. **CLI フラグ** — コマンドラインで明示的に指定したオプション（例: `-q -1`, `--thresh gau_loose`）。CLI デフォルトのままのオプションは YAML の値を上書きしません。

例: YAML で `charge: 0` を設定し、CLI で `-q -1` を渡した場合、電荷は `-1` になります。

この優先順位は `all`, `opt`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, `dft`, `sp` に共通です。

実行で使われる値は {ref}`-v 3 <ja-verbosity-levels>` で確かめられます。各セクションを、名前、`-` の下線、実際に使う値の順に表示します:

```text
opt
---
thresh: gau
max_cycles: 100000
…
```

`--show-config`（`scan`・`scan2d`・`scan3d` には無い）は、読み込んだ YAML ファイルと最上位のキーを表示します。

セクション名を書き間違えると `[config] WARNING: YAML section(s) … are not recognized and were ignored.` が出ます。セクションの中のキーの書き間違いの扱いはセクションで決まります。

- `calc` は `[backend] WARNING: … ignored calc setting(s) …` を出して続けます。
- `freq`、`thermo`、`bias`、`bond`、`search`、`sp`、`hessian_dimer` の直下のキー、scan と path のコマンドでの `geom` は、何も出さずに無視します。`-v 3` ではそのセクションの表示に余分な行として出ます。
- ほかのセクションは、そのキーの名前を示すエラーで止まります。

(ja-common-cli-to-yaml-mapping)=
## 主要な CLI→YAML マッピング

| CLI フラグ | YAML キー | セクション |
|----------|----------|---------|
| `-q` / `--charge` | `charge` | `calc` |
| `-m` / `--multiplicity` | `spin` | `calc` |
| `-b` / `--backend` | `backend` | `calc` |
| `--backend-model` | `model` | `calc` |
| `--solvent` | `solvent` | `calc` |
| _(YAML のみ)_ | `device` | `calc` |
| `--thresh` | `thresh` | `opt` |
| `--max-cycles` | `max_cycles` | コマンド別: `opt`/`tsopt` は `opt`、`irc` は `irc` |
| `--max-cycles-gsm` | `max_cycles` | `stopt`（`stopt.stop_in_when_full` も設定） |
| `--dmf-max-iterations` | `max_cycles` | `dmf` |
| `--gsm-param` | `param` | `gs` |
| `--dump` | `dump` | `opt`（opt、tsopt、scan）、`stopt`（path-opt、path-search）、`thermo`（freq） |
| `--step-size`（irc） | `step_length` | `irc` |
| `--opt-mode` | _(CLI のみ)_ | — |
| `--freeze-atoms` | `freeze_atoms` | `geom` |
| `--coord-type` | `coord_type` | `geom` |
| `--temperature`（freq、`all --freq-temperature`） | `temperature` | `thermo` |
| `--pressure`（freq、`all --freq-pressure`） | `pressure_atm` | `thermo` |
| `--dft-engine` | `engine` | `dft` |

```{note}
**名前不一致 — `--pressure` vs `pressure_atm`.** どちらも atm で、内部で Pa に変換されます。単位を名前に含むのは YAML キーだけです。
```

### サブコマンド別の `--thresh` デフォルト

| サブコマンド | デフォルト `--thresh` |
|------------|---------------------|
| `opt` | `gau` |
| `tsopt`（Hessian Dimer） | `baker` |
| `tsopt`（RS-P-RFO / RS-I-RFO） | `baker` |
| `scan` | `gau` |
| `scan2d`, `scan3d` | `baker` |
| `path-opt`、`path-search`（単一構造の最適化） | `gau` |
| `path-opt`、`path-search`（GSM のストリング: `--thresh-gsm`、`stopt.thresh`） | `gau_loose` |
| `all`（事前最適化、後処理の極小化） | `gau` |
| `all`（後処理の TS の段） | `baker` |

受け付ける値: `gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`。実行ごとに `--thresh <preset>` または YAML の `opt.thresh` で上書きできます。

```{note}
**`--thresh` を持たないサブコマンド。** `irc`、`freq`、`dft` には `--thresh` が**ありません**:

- `irc` — 収束は `irc.rms_grad_thresh`、`irc.energy_thresh`、`irc.max_cycles` で制御されます。IRC は予測子–修正子積分器に従うため、力ベース極小化用のプリセット群は適用されません。
- `freq` — 最適化ステップが無いため `--thresh` は存在しません。数値精度は `--hessian-calc-mode` と MLIP 自体の精度で決まります。
- `dft` — SCF 収束は `dft.conv_tol` と `dft.max_cycle` で制御されます。`gau`/`baker` プリセットファミリは使用しません。
```

## 共通セクション

### `geom`

ジオメトリ読み込みと座標系の設定。

```yaml
geom:
 coord_type: cart # "cart"（デカルト座標）、"redund"（冗長内部座標）、"dlc"（非局在化内部座標）、"tric"（並進・回転を含む内部座標）。opt・tsopt・scan・scan2d・scan3d は 4 つとも、all・path-opt・path-search は cart と dlc のみ
 freeze_atoms: [] # 1 始まりの凍結原子番号。--freeze-links が有効なとき（PDB/mmCIF 入力、または --ref-pdb つきの XYZ/GJF）、自動検出した cap 水素の親原子の番号を合わせる
```

**注記:**
- 凍結原子の力はゼロ化されます。デフォルトの `return_partial_hessian: true` では、動かせる原子の部分の Hessian だけを返します。false にすると、凍結原子の行・列をゼロにした全体の行列を返します。
- デカルト座標の PHVA（部分 Hessian 振動解析）では、凍結原子を動かさない全系の剛体運動だけを除きます。詳細は [freq](freq.md#凍結境界での剛体モード) を参照してください。
- `irc` では、YAML や CLI の指定によらず `geom.coord_type` は `cart` です。

---

### `calc`

エネルギーと力を計算する calculator の設定。

```yaml
calc:
 backend: uma           # uma, orb, mace, aimnet2, dft, auto
 precision: auto # auto (uma/aimnet2 fp32、orb/mace fp64) | fp32 | fp64; aimnet2 は auto/fp32 のみ受理し fp64 を拒否
 charge: 0 # 全電荷。ここに書いたときだけ使う（組み込みの既定値は無い）。-q と -l が優先
 spin: 1 # Spin multiplicity 2S+1 (overridden by CLI -m)
 model: uma-s-1p2 # UMA: uma-s-1p2 | uma-m-1p1。model を書かずに backend を orb / mace / aimnet2 にすると orb_v3_conservative_omol / MACE-OMOL-0 / aimnet2
 task_name: omol # Task tag recorded in UMA batches
 device: auto # Device: "cuda", "cpu", or "auto"
 max_neigh: null # Maximum neighbors for graph construction
 radius: null # Cutoff radius for neighbor search
 r_edges: false # Store radial edges
 workers: 1 # UMA inference workers
 workers_per_node: 1 # Workers per node for parallel predictor
 out_hess_torch: true # Return Hessian as torch.Tensor
 hessian_double: true # Assemble/return Hessian in float64
 # freeze_atoms: null # geom.freeze_atoms から継承されるため直接指定しない
 hessian_calc_mode: FiniteDifference # Hessian mode: "Analytical" or "FiniteDifference"
 return_partial_hessian: true  # 動かせる原子の部分の Hessian だけを返す
 print_timing: true # Hessian計算のタイミング内訳を表示
 print_vram: true # Hessian計算中の CUDA VRAM 使用量を表示 (UMA バックエンドのみ)
 # xTB 溶媒補正（計算コスト大）
 solvent: none           # none, water, methanol, acetonitrile, dmso, thf, toluene
 solvent_model: alpb     # xTB solvent model: "alpb" or "cpcmx"
 xtb_cmd: xtb            # xTBコマンドと追加オプションを指定可能（例: "xtb --etemp 1000"）
 xtb_acc: 0.2            # xTB accuracy parameter
 # backend: dft の場合だけ使用
 dft:
  func_basis: wb97m-v/def2-svp
  engine: gpu             # gpu (GPU4PySCF) | cpu (PySCF)
  lowmem: true             # DF tensorを保持しないdirect JK
  density_fit: false       # 既定は lowmem の逆
  nprocs: auto             # scheduler/affinityからPySCF thread数を決定
  memory: auto             # host RAM上限（例64GB、GPU VRAMではない）
  solvent: none
  solvent_model: smd      # pcm | smd
  save_scf_checkpoint: false
  checkpoint_path: null   # 有効時の既定（all 以外のコマンド）: <out-dir>/_work/dft_scf/state.chk
  pyscf:
   mol: {}
   mf: {}
   grids: {}
   density_fit: {}
   with_df: {}
   with_solvent: {}
```

`backend: dft` は、上の `calc.dft` ブロックの設定で、エネルギー・力・Hessian をすべて DFT（PySCF または GPU4PySCF）で計算します。詳しくは [MLIP の TS を DFT で確かめる](dft-backend.md) を参照してください。後述の最上位の [`dft` セクション](#ja-dft-section) は、別の `dft` コマンド（`all --dft` も実行する）の設定です。どちらも同じ SCF のキー（`conv_tol`・`max_cycle`・`grid_level` など）を受け付け、`save_scf_checkpoint` と `checkpoint_path` は `calc.dft` だけにあります。`hessian_calc_mode: Analytical` は有限変位の誤差を避けられますが、計算時間とメモリはバックエンドと系によって変わるので、先に自分の系で試してください。`charge` と `spin` が CLI や `.gjf` テンプレートとどう組み合わさるかは {ref}`電荷の指定 <ja-charge-specification>` を参照してください。

---

### `opt`

L-BFGS/RFO で共通の単一構造最適化設定。

```yaml
opt:
 thresh: gau # Convergence preset: gau_loose, gau, gau_tight, gau_vtight, baker, never
 max_cycles: 100000 # Maximum optimizer iterations
 print_every: 100 # Logging stride
 min_step_norm: 1.0e-08 # Minimum step norm for acceptance
 assert_min_step: true # Stop if steps fall below threshold
 rms_force: null # Explicit RMS force target
 rms_force_only: false # Rely only on RMS force convergence
 max_force_only: false # Rely only on max force convergence
 force_only: false # Skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # RMS threshold when converging to reference geometry
 overachieve_factor: 0.0 # 0.0 = off; >0: converge when forces < thresh/factor, ignoring step (not used by baker)
 check_eigval_structure: false # Validate Hessian eigenstructure
 line_search: true # Enable line search
 energy_plateau: false # opt-in: エネルギー地形が平坦になった場合に stalled として停止（下記注記を参照）
 energy_plateau_thresh: 1.0e-04 # au (~0.06 kcal/mol); 平坦判定のレンジ閾値
 energy_plateau_window: 50 # 平坦判定に用いる直近ステップ数
 dump: false # Dump trajectory/restart data
 dump_restart: false # Dump restart checkpoints
 prefix: "" # Filename prefix
 out_dir: ./result_opt/ # Output directory
```

**平坦なエネルギー地形による停止（デフォルト無効）:** `opt` / `tsopt` / `all` の `--stop-plateau` で `energy_plateau` を有効にし、`--stop-plateau-thresh` / `--stop-plateau-window` で下記の 2 つの値を設定します。有効にすると、直近 `energy_plateau_window` ステップのエネルギーレンジ（max − min）が `energy_plateau_thresh` を下回ったとき、オプティマイザは収束扱いにせず `stalled` として停止します。エネルギーが平坦になっても MLIP の力のノイズで力が `baker` の閾値を下回らないときに使ってください。`tsopt` がプラトーで止まったときも Hessian を計算して n_imag を出します。TS 最適化が収束しないまま `max_cycles` に達したときは Hessian を計算しません。GSM と DMF の経路最適化はこの停止を使いません。

**収束プリセット**（デカルト座標では力は Hartree/Bohr、ステップは Bohr。角度は Hartree/rad と rad）:

| Preset | Max Force | RMS Force | Max Step | RMS Step |
|--------|-----------|-----------|----------|----------|
| `gau_loose` | 2.5e-3 | 1.7e-3 | 1.0e-2 | 6.7e-3 |
| `gau` | 4.5e-4 | 3.0e-4 | 1.8e-3 | 1.2e-3 |
| `gau_tight` | 1.5e-5 | 1.0e-5 | 6.0e-5 | 4.0e-5 |
| `gau_vtight` | 2.0e-6 | 1.0e-6 | 6.0e-6 | 4.0e-6 |
| `baker` | 3.0e-4 | 2.0e-4 | 3.0e-4 | 2.0e-4 |

`baker` は 4 列すべてに加えて前サイクルとの `|delta E| < 1e-6` hartree を要求し、Bakken and Helgaker（*J. Chem. Phys.* **117**, 9160 (2002)）が示した Baker 基準（`max(|force|) <= 3e-4` **かつ**（`|delta E| < 1e-6` **または** `max(|step|) <= 3e-4`））より厳しい設定です。文献の形は力の RMS が残った構造も収束と判定しうるため、機械学習ポテンシャル上では高次の鞍点で止まることがあり、厳しい形を使っています。

---

### `lbfgs`

L-BFGS の設定（`opt` を拡張）。

```yaml
lbfgs:
  # Inherits all opt settings, plus:
 keep_last: 7 # History size for L-BFGS buffers
 beta: 1.0 # Initial damping beta
 gamma_mult: false # Multiplicative gamma update toggle
 max_step: 0.3 # Maximum step length
 control_step: true # Control step length adaptively
 double_damp: true # Double damping safeguard
 mu_reg: null # Regularization strength
 max_mu_reg_adaptions: 10 # Cap on mu adaptations
 reject_uphill: false # 許容値を超えるenergy上昇の拒否を明示的に有効化
 uphill_tolerance: 0.0001 # energy上昇許容値（Hartree）
 rejection_step_floor: 1.0e-07 # retry stepの下限
 max_rejections_at_floor: 3 # 下限での連続拒否後に停止
```

---

### `rfo`

RFO（Rational Function Optimizer）の設定（`opt` を拡張）。

```yaml
rfo:
  # Inherits all opt settings, plus:
 trust_radius: 0.10 # Trust-region radius
 trust_update: true # Enable trust-region updates
 trust_min: 0.0001 # Minimum trust radius
 trust_max: 0.10 # Maximum trust radius (bohr)
 max_energy_incr: null # Allowed energy increase per step
 reject_uphill: false # 許容値を超えるenergy上昇の拒否を明示的に有効化
 uphill_tolerance: 0.0001 # energy上昇許容値（Hartree）
 rejection_trust_floor: 1.0e-07 # retry trust radiusの下限
 max_rejections_at_floor: 3 # 下限での連続拒否後に停止
 hessian_update: ts_bfgs # Hessian update scheme: ts_bfgs, bfgs, bofill, etc.
 hessian_init: calc # Hessian initialization: calc, unit, etc.
 hessian_recalc: 500 # Rebuild Hessian every N steps
 hessian_recalc_adapt: null # Adaptive Hessian rebuild factor
 small_eigval_thresh: 1.0e-08 # Eigenvalue threshold for stability
 alpha0: 1.0 # Initial micro step
 max_micro_cycles: 50 # RS iteration limit per step
 rfo_overlaps: false # Enable RFO overlaps
 gediis: false # Enable GEDIIS
 gdiis: true # Enable GDIIS
 gdiis_thresh: 0.0025 # GDIIS acceptance threshold
 gediis_thresh: 0.01 # GEDIIS acceptance threshold
 gdiis_test_direction: true # Test descent direction before DIIS
 adapt_step_func: true # Adaptive step scaling
```

## 経路最適化セクション

### `gs`

Growing String Method（GSM）の設定。

```yaml
gs:
 fix_first: true # Keep first endpoint fixed
 fix_last: true # Keep last endpoint fixed
 max_nodes: 20 # Maximum string nodes (internal images); GSM では端点 2 つを加えた数が全イメージ数
 perp_thresh: 0.005 # Perpendicular displacement threshold
 reparam_check: rms # Reparameterization check metric
 reparam_every: 1 # Reparameterization stride
 reparam_every_full: 1 # Full reparameterization stride
 param: equi # Parameterization scheme
 max_micro_cycles: 10 # RS iteration limit per step
 reset_dlc: true # Rebuild delocalized coordinates each step
 climb: true # Enable climbing image
 climb_rms: 0.0005 # Climbing RMS threshold
 climb_lanczos: true # Lanczos refinement for climbing
 climb_lanczos_rms: 0.0005 # Lanczos RMS threshold
 climb_fixed: false # Keep climbing image fixed
 scheduler: null # Optional scheduler backend
```

```{note}
`gs.max_nodes` / `--max-nodes` は **GSM** と **DMF** のどちらでも、端点の間で動かすイメージの数です。両エンジンとも端点 2 つを加えるため、イメージは全部で `max_nodes + 2` 個です。詳細は [`path-opt`](path-opt.md) を参照してください。

`gs.param` は `equi` または `energy` を受け付けます。energy weighting はGSMストリングの完全成長後にのみ適用され、高エネルギー領域へノード密度を寄せます。
```

---

### `dmf`

Direct Max Flux（DMF）による MEP 最適化。DMF は初期経路を FB-ENM（flat-bottom elastic network model）で、`correlated: true` のときは CFB-ENM（correlated FB-ENM）で作ります。

```yaml
dmf:
 backend: gpu # gpu (dmf.torch / CUDA、default) | cpu (dmf / NumPy)
 max_cycles: 3000 # DMF/IPOPT の最大反復数（--dmf-max-iterations で上書き）
 tol: tight # IPOPT dual_inf_tol: tight(0.04) | middle(0.10) | loose(0.20) または正の float（--dmf-tol で上書き）
 correlated: true # 初期経路を FB-ENM ではなく CFB-ENM で作る
 sequential: true # Sequential DMF execution
 fbenm_only_endpoints: false # Run FB-ENM beyond endpoints
 fbenm_options:
   delta_scale: 0.2 # FB-ENM displacement scaling
   bond_scale: 1.25 # Bond cutoff scaling
   fix_planes: true # Enforce planar constraints
 cfbenm_options:
   bond_scale: 1.25 # CFB-ENM bond cutoff scaling
   corr0_scale: 1.1 # Correlation scale for corr0
   corr1_scale: 1.5 # Correlation scale for corr1
   corr2_scale: 1.6 # Correlation scale for corr2
   eps: 0.05 # Correlation epsilon
   pivotal: true # Pivotal residue handling
   single: true # Single-atom pivots
   remove_fourmembered: true # Prune four-membered rings
 dmf_options:
   remove_rotation_and_translation: false # Keep rigid-body motions
   mass_weighted: false # Toggle mass weighting
   parallel: false # Enable parallel DMF
   eps_vel: 0.01 # Velocity tolerance
   eps_rot: 0.01 # Rotational tolerance
   beta: 10.0 # Beta parameter for DMF
   update_teval: false # Update transition evaluation
 ipopt_options: {} # 生の IPOPT option（例: {dual_inf_tol: 0.04}）
 k_fix: 300.0 # Harmonic constant for restraints (dmf 直下、dmf_options 配下ではない)
```

`dmf.tol` は DMF ソルブが最後に適用する許容値なので、同じファイル内の `ipopt_options.dual_inf_tol` より優先されます。生の IPOPT オプションを固定したい場合は `dmf.tol` を書かず `ipopt_options.dual_inf_tol` のみを指定してください。`gau_tight` などの Gaussian プリセットはここでは拒否され、`--thresh` / `--thresh-gsm` の担当です。

---

### `search`

再帰的経路探索（path-search のみ）。

```yaml
search:
 max_depth: 10 # 許可する再帰分割の階層数（0 = 分割しない）
 stitch_rmsd_thresh: 0.0001 # RMSD threshold (Bohr) for stitching segments
 bridge_rmsd_thresh: 0.0001 # RMSD threshold (Bohr) for bridging nodes
 max_nodes_segment: 20 # Max nodes per segment
 max_nodes_bridge: 5 # Max nodes per bridge
 kink_max_nodes: 3 # Max nodes for kink optimizations
 max_seq_kink: 2 # Max sequential kinks
 refine_mode: null # Refinement strategy: peak, minima, or null (auto)
```

---

### `stopt`

chain-of-states 経路最適化の StringOptimizer 設定。`stopt.lbfgs` / `stopt.rfo` は、`opt.lbfgs` / `opt.rfo` と同じように単一構造の最適化を設定します。

```yaml
stopt:
 type: string # Optimizer type label
 thresh: gau_loose # StringOptimizer convergence preset
 stop_in_when_full: 300 # ストリングが伸びきった後に許すサイクル数。使い切ると未収束で止まる
 align: false # path-opt/path-search では常に false。イメージの重ね合わせは別に Kabsch 法で行う
 scale_step: global # Step scaling mode
 max_cycles: 300 # Maximum StringOptimizer iterations
 dump: false # Dump trajectory/restart data
 dump_restart: false # Dump restart checkpoints
 reparam_thresh: 0.0 # Reparameterization threshold
 coord_diff_thresh: 0.0 # Coordinate-difference threshold
 out_dir: ./result_path_opt/ # Output directory
 print_every: 10 # Logging stride
```

## TS 最適化セクション

TS 最適化は `--opt-mode` で**2 つのセクション**のどちらかを使います:
- `--opt-mode dimer`（または `grad`）→ `hessian_dimer` セクション
- `--opt-mode rsprfo`（または `hess`、デフォルト）、`rsirfo`、`trim` → `rsirfo` セクション

`opt` と使用中のセクションの片方だけに書いたキーは両方で使い、両方に異なる値を書くとエラーで止まります。どちらにも `thresh` が無ければ、TS 最適化は `opt` に示した `gau` ではなく `baker` を使います。

### `hessian_dimer`

Hessian Guided Dimer TS 最適化。

```yaml
hessian_dimer:
 thresh_loose: gau_loose # Loose convergence preset
 thresh: baker # Main convergence preset
 update_interval_hessian: 500 # Hessian rebuild cadence
 neg_freq_thresh_cm: 5.0 # n_imag を数える閾値（cm⁻¹）。freq セクションを参照
 flatten_amp_ang: 0.1 # Flattening amplitude (Å)
 flatten_max_iter: 0 # flatten の回数（0 = 無効）。0 のとき --flatten は 50 を使う
 flatten_sep_cutoff: 0.0 # Minimum distance between representative atoms
 flatten_k: 10 # Representative atoms sampled per mode
 flatten_loop_bofill: false # Bofill update for flatten displacements
 mem: 100000 # Memory limit for solver
 device: auto # Device selection for eigensolver
 root: 0 # Targeted TS root index
 dimer:
   length: 0.0189 # Dimer separation (Bohr)
   rotation_max_cycles: 15 # Max rotation iterations
   rotation_method: fourier # Rotation optimizer method
   rotation_thresh: 0.0001 # Rotation convergence threshold
   rotation_tol: 1 # Rotation tolerance factor
   rotation_max_element: 0.001 # Max rotation matrix element
   rotation_interpolate: true # Interpolate rotation steps
   rotation_disable: false # Disable rotations entirely
   rotation_disable_pos_curv: true # Disable when positive curvature detected
   rotation_remove_trans: true # 選択した剛体null成分を除去
   trans_force_f_perp: true # Project forces perpendicular to translation
   bonds: null # Bond list for constraints
   N_hessian: null # Hessian size override
   bias_rotation: false # Bias rotational search
   bias_translation: false # Bias translational search
   bias_gaussian_dot: 0.1 # Gaussian bias dot product
   seed: null # RNG seed for rotations
   write_orientations: false # 方向を出力（明示的な true も可）
   forward_hessian: true # Propagate Hessian forward
 lbfgs:                    # `hessian_dimer` 内で `dimer` と同階層
   # Same keys as lbfgs section
   thresh: baker
   line_search: false # 必須: Dimer の有効力は物理エネルギーと共役でない
```

内側の L-BFGS 固有設定は、最上位の `lbfgs` ではなく `hessian_dimer.lbfgs` に置きます。共通の `print_every` と `energy_plateau*` は上記の競合規則に従います。`line_search` は `false` 固定です。Dimer の射影・反転した有効力は表示する物理エネルギーの勾配ではないため、`true` は拒否されます。`max_cycles` はここでは設定しません。Hessian の更新の間に走る L-BFGS の 1 回ごとに、`opt.max_cycles` の残りのサイクル数までを使います。

```{note}
**`flatten_max_iter`。** `--flatten` は余分な虚振動を最大 `flatten_max_iter` 回（値が 0 のときは 50 回）の flatten で取り除き、`--no-flatten` は flatten を無効にします。どちらのフラグも無いときは、ここに正の値を書くと flatten が有効になります。
{ref}`--flatten を使うとき <ja-flatten-precedence-caveat>` を参照してください。
```

---

### `rsirfo`

RS-I-RFO / RS-P-RFO TS 最適化。

```yaml
rsirfo:
 thresh: baker # RS-IRFO convergence preset
 max_cycles: 100000 # opt.max_cycles と共有。異なる明示値はエラー
 print_every: 100 # Logging stride
 min_step_norm: 1.0e-08 # Minimum accepted step norm
 assert_min_step: true # Assert when steps stagnate
 roots: [0] # 対象root indexを1個だけ指定（一次鞍点のみ）
 hessian_ref: null # Reference Hessian
 rx_modes: null # Reaction-mode definitions
 prim_coord: null # Primary coordinates to monitor
 rx_coords: null # Reaction coordinates to monitor
 hessian_update: bofill # Hessian update scheme
 hessian_recalc: 500 # Rebuild exact Hessian every N macro steps (rfo から継承)
 hessian_recalc_reset: true # Reset recalc counter after exact Hessian
 max_micro_cycles: 50 # RS iteration limit per step
 augment_bonds: false # Augment reaction path based on bond analysis
 min_line_search: false # 常に false: RS-P-RFO は line search を使わない
 max_line_search: false # 常に false: RS-P-RFO は line search を使わない
 assert_neg_eigval: false # Require negative eigenvalue at convergence
 track_mode_by_overlap: false # 前回の Hessian との重なりで追跡対象 TS モードを選ぶ
 reject_mode_loss: false # TS モードが見つかった後、それを失うステップを棄却し、trust radius を小さくしてやり直す
 mode_loss_trust_floor: 1.0e-05 # そのやり直しでの trust radius の下限
 max_mode_loss_rejections: 5 # 下限到達後に許す棄却回数
 verify_saddle: true # 収束時に厳密な Hessian で n_imag を数える。n_imag = 0 は収束として受けない
 saddle_imaginary_threshold_cm: 5.0 # n_imag を数える閾値（cm⁻¹）。freq セクションを参照
 saddle_recovery_step: 0.01 # 極小（n_imag = 0）から抜けるための上り方向のステップ
 saddle_recovery_check_interval: 50 # その回復中に厳密な Hessian を確かめる間隔（ステップ数）
 saddle_recovery_max_cycles: 0 # 回復のステップ数の上限。0 で回復しない
 out_dir: ./result_tsopt/ # 出力ディレクトリ
 # Also inherits rfo-like settings: trust_radius, trust_update, etc.
```

RS-P-RFO では、`min_line_search` または `max_line_search` に `true` を書くと警告を出して `false` に戻します。RS-I-RFO と TRIM はどちらのキーも無視し、Dimer は `hessian_dimer.lbfgs.line_search` を使います。

```{note}
**`--flatten` の優先順位。** RS-P-RFO・RS-I-RFO・TRIM の flatten ループも `hessian_dimer.flatten_max_iter` を読み、規則は `hessian_dimer` の注記と同じです。
{ref}`--flatten を使うとき <ja-flatten-precedence-caveat>` を参照してください。
```

## IRC セクション

(ja-irc-section)=
### `irc` セクション

IRC 積分設定。

```yaml
irc:
 step_length: 0.1 # 積分のステップ長（Bohr、質量加重しないデカルト座標。--step-size）
 never_stop: false # 物理的な端点判定を無視してmax_cyclesまで追跡
 max_cycles: 125 # Maximum steps along IRC
 forward: true # Propagate in forward direction
 backward: true # Propagate in backward direction
 root: 0 # Normal-mode root index
 hessian_init: calc # Hessian initialization source
 hessian_update: bofill # Hessian update scheme
 hessian_recalc: null # Hessian rebuild cadence
 energy_increase_thresh: 0.0   # 通常modeでは1 stepでもenergyが上昇すれば停止
 dump_every: null # デフォルト無効。正の間隔では座標・energy・gradientのみをcheckpoint保存（Hessianなし）
 dump_fn: irc_data.h5 # dump_every指定時のcheckpointファイル名
 displ: energy # Displacement construction method
 displ_energy: 0.001 # Energy-based displacement scaling
 displ_length: 0.1 # Length-based displacement fallback
 rms_grad_thresh: 0.001 # RMS gradient convergence threshold
 hard_rms_grad_thresh: null # Hard RMS gradient stop
 energy_thresh: 0.000001 # Energy change threshold
 imag_below: 0.0 # 出発の mode（root）の ν がこの値以下のときだけ IRC を始める（cm⁻¹）
 force_inflection: true # Enforce inflection detection
 check_bonds: false # Check bonds during propagation
 out_dir: ./result_irc/ # Output directory
 prefix: "" # Filename prefix
 max_pred_steps: 500 # Predictor-corrector max steps
 loose_cycles: 3 # Loose cycles before tightening
 corr_func: mbs # EulerPC コレクタ関数
```

`corr_func` は、予測子–修正子法の IRC 積分器（EulerPC）が使う修正子ステップを選びます。登録されているのは `"mbs"`（Modified Bulirsch–Stoer）だけで、それ以外の値は構築時にエラーになります。

## 振動解析セクション

(ja-freq-section)=
### `freq` セクション

振動解析設定。

```yaml
freq:
 zero_cutoff_cm: 5.0 # ν < -zero_cutoff_cm を虚振動として数える
 amplitude_ang: 0.8 # Displacement amplitude for modes (Å)
 n_frames: 20 # モードtrajectoryのフレーム数
 max_write: 10 # Maximum number of modes to write
 sort: value # Sort order: "value" or "abs"
 out_dir: ./result_freq/ # Output directory
```

虚振動の既定の分類基準は ν < −5.00 cm⁻¹ です。`freq`、`opt`（flatten）、`tsopt` はこの `freq.zero_cutoff_cm` で n_imag を数え、`irc` はこの値を読みません。`tsopt` では `hessian_dimer.neg_freq_thresh_cm` か `rsirfo.saddle_imaginary_threshold_cm` でもこの閾値を設定でき、3 つのうち 2 つを異なる値で明示するとエラーで止まります。`n_negative_modes` は cutoff の内側の負の振動数も数えます。n_imag と `n_negative_modes` のどちらも、最適化の収束の判定には使いません。cutoff の値によらず、出力は符号付きの全振動数を残し、熱化学は正のモードをすべて使います。

---

### `thermo`

熱化学設定。

```yaml
thermo:
 temperature: 298.15 # Thermochemistry temperature (K)
 pressure_atm: 1.0 # Thermochemistry pressure (atm)
 symmetry_number: null # 自動判定。正整数は高度な上書き指定
 dump: false # Write thermoanalysis.yaml
```

## DFT セクション

(ja-dft-section)=
### `dft` セクション

DFT 計算設定。

```yaml
dft:
 func: wb97m-v # Exchange-correlation functional
 basis: def2-svp # Basis set name
 func_basis: null # Combined "FUNC/BASIS" string (overrides func/basis)
 conv_tol: 1.0e-09 # SCF convergence tolerance (hartree)
 max_cycle: 100 # Maximum SCF iterations
 grid_level: 3 # PySCF grid level
 engine: gpu # SCF backend: "gpu" (GPU4PySCF) or "cpu" (PySCF)
 solvent: none # PySCF native solvent名
 solvent_model: smd # pcm | smd
 pyscf: {} # PySCF object名attribute転送
 lowmem: true # 低memory direct JK。falseでdensity fitting
 nprocs: auto # scheduler/affinityからPySCF thread数を決定
 memory: auto # host RAM上限（例64GB、GPU VRAMではない）
 verbose: 0 # PySCF verbosity (0-9)。-v 0/1 で効く。既定の -v 2 と -v 3 では 4 以上に上がる
 out_dir: ./result_dft/ # Output directory root
```

## スキャン関連セクション

スキャン座標は `--config` の YAML ではなく `-s/--scan-lists` で指定します。構文は {ref}`スキャンリスト仕様 <ja-scan-list-spec>` を参照してください。

(ja-bias-section)=
### `bias`

`scan`・`scan2d`・`scan3d` で使用する調和バイアス設定。

```yaml
bias:
 k: 300 # Harmonic bias strength (eV·Å⁻²)
```

**サブコマンド間で共有されるばね定数。** 同じ物理的な調和ペナルティ（`k`、単位 eV·Å⁻²）が以下の箇所にデフォルト値 `300.0` で現れます:

| YAML キー | 使用元 | CLI フラグ |
|----------|-------|-----------|
| `bias.k` | `scan`, `scan2d`, `scan3d` | `--restraint-k` |
| `dmf.k_fix` | `path-opt` / `path-search` で `--mep-mode dmf` を使用する場合 | —（YAML 専用） |

`opt` も `--distance-restraint` の組に同じ既定値の `--restraint-k` を使いますが、CLI フラグからだけ読み、`bias:` セクションは読みません。

調和拘束の強さを調整したい場合はこれらのいずれかを上書きしてください。値を小さく（例: `20.0`）すると、柔らかい誘導項としてジオメトリが緩和しやすくなります。デフォルト値はほぼ剛体的に固定する値です。

---

### `bond`

元素の共有結合半径による結合変化検出。

```yaml
bond:
 device: auto # 距離計算に使うデバイス: "cuda"、"cpu"、"auto"
 bond_factor: 1.2 # Covalent-radius scaling for cutoff
 margin_fraction: 0.05 # Fractional tolerance for comparisons
 delta_fraction: 0.05 # Minimum relative change to flag bond formation/breaking
```

## 例: 設定ファイルの全体例

```yaml
# pdb2reaction configuration example

geom:
 coord_type: cart
 freeze_atoms: []

calc:
 backend: uma
 model: uma-s-1p2 # 選んだバックエンドのモデル名（UMA: uma-s-1p2 | uma-m-1p1）
 device: auto
 hessian_calc_mode: FiniteDifference # 移植性のあるデフォルト値。Analytical は事前検証して選択

gs:
 max_nodes: 12
 climb: true
 climb_lanczos: true

stopt:
 thresh: gau_loose
 max_cycles: 300
 dump: false

lbfgs:
 max_cycles: 100000

rfo:
 max_cycles: 100000

bond:
 bond_factor: 1.2
 delta_fraction: 0.05

search:
 max_depth: 10
 max_nodes_segment: 20

freq:
 max_write: 10
 amplitude_ang: 0.8

thermo:
 temperature: 298.15
 pressure_atm: 1.0
 symmetry_number: null

dft:
 func: wb97m-v
 basis: def2-svp
 grid_level: 3
```

## 使用上の注意点

- `workers` と `workers_per_node` は UMA バックエンドでだけ効きます。
- `workers > 1` の UMA は解析 Hessian を計算できず、`hessian_calc_mode: Analytical` を明示するとエラーで止まります。`workers: 1` にするか `FiniteDifference` を使ってください。詳細は {ref}`workers と解析 Hessian <ja-workers-analytical-error>` を参照してください。
- `freq` と `irc` は、`calc.return_partial_hessian` の指定によらず部分 Hessian を使います。
- `all` は実行する各段に同じファイルを渡し、各段は概要の表の「使用されるコマンド」の列でそのコマンドが挙がっているセクションを読みます。例えば TS の段は `opt`・`hessian_dimer`・`rsirfo` を含む `tsopt` のセクションを読み、`dft` セクションは `all --dft` のときに効きます。
- `opt.lbfgs` / `opt.rfo` は `lbfgs` / `rfo`、`freq.thermo` は `thermo` の別の書き方です。同じ設定に 2 か所で異なる値を書くとエラーで止まります。例えば `lbfgs.max_cycles` と `opt.lbfgs.max_cycles`、また L-BFGS を選んでいるときの `opt.max_cycles` と `lbfgs.max_cycles` です。`-o/--out-dir` は `out_dir` キーより優先し、`all` は各段の出力先を自分で決めます。

## 関連ドキュメント

- [all](all.md) - 一気通貫ワークフロー
- [opt](opt.md) - 単一構造最適化
- [tsopt](tsopt.md) - 遷移状態最適化
- [path-search](path-search.md) - 再帰的 MEP 探索
- [freq](freq.md) - 振動解析
- [dft](dft.md) - DFT 計算
- [バックエンド](backends.md) - MLIP バックエンドの詳細
- [トラブルシューティング](troubleshooting.md) - よくあるエラーと対処法
