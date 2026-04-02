# `tsopt`

## 概要

> **要約:** 遷移状態（TS）*候補*を、Dimer（`--opt-mode grad`）または RS-I-RFO（Restricted-Step Image Rational Function Optimization）（`--opt-mode hess`、デフォルト）で最適化します。`tsopt` は最後に自動で Hessian 計算と虚振動数チェックを行います。妥当な TS（一次鞍点）では虚振動数が **ちょうど 1 つ** です。端点の接続性は `irc` で確認してください。

### 要点
- **入力:** `path-opt` / `path-search` が出力する HEI、または自前の TS 初期構造（`geom_loader` が扱える形式）。
- **モード:** `hess`（`rsirfo`）= RS‑I‑RFO（デフォルト、一般的により堅牢）。`grad`（`dimer`）= Hessian Guided Dimer（1ステップあたりのコストが低いことが多い）。
- **品質確認:** `tsopt` は最終 Hessian 計算と虚振動数チェックを内部で実行します（出力の n=1 を確認）。結果はなお *候補* であり、[irc](irc.md) で端点の接続性を確認してください。完全な振動解析や熱化学補正が必要な場合のみ、別途 [freq](freq.md) を実行します。
- **任意の後処理:** `--flatten`（デフォルト無効）で余分な虚振動数モードの除去を制御します。
- **出力変換:** `--convert-files`（デフォルト）で、PDB 入力は（`--dump` のとき）`.pdb` を併記し、Gaussian テンプレートは最終構造の `.gjf` を書き出します。

### `--opt-mode` の選び方
- **`--opt-mode hess`（RS‑I‑RFO）**: デフォルト。ヘシアン計算のコストを許容でき、堅牢性を重視する場合に推奨。
- **`--opt-mode grad`（Dimer）**: 軽量な探索や、複数の TS 初期構造から素早く反復したい場合に有効。

> **命名規則の注意:** CLI は `grad|dimer`（= Dimer）および `hess|rsirfo`（= RS-I-RFO、デフォルト）を受け付けます。YAML では `dimer` または `rsirfo` を直接指定してください。

XYZ/GJF 入力では `--ref-pdb` で参照 PDB トポロジーを指定し、XYZ 座標を保持したまま PDB/GJF への変換が可能です。TS 初期構造が必要な場合は、2 端点なら [path-opt](path-opt.md)、2 構造以上なら [path-search](path-search.md) で HEI を取得してから `tsopt`（内部で虚振動数チェック済み）→ `irc` の順で検証してください。

## 最小例

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --out-dir ./result_tsopt
```

## 出力の見方

- `result_tsopt/final_geometry.pdb`（または `final_geometry.xyz`）
- `result_tsopt/vib/imag_*_trj.xyz`
- `result_tsopt/vib/imag_*.pdb`（PDB 入力の場合）

## よくある例

1. VRAM に余裕がある場合に dimer モード + 解析的ヘシアンで実行する。

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
 --opt-mode grad --hessian-calc-mode Analytical --out-dir ./result_tsopt_grad
```

2. 最適化軌跡を保存して確認する。

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --dump --out-dir ./result_tsopt_dump
```

3. rsirfo モードを YAML 上書きと併用する。

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
 --opt-mode hess --config tsopt.yaml --out-dir ./result_tsopt_hess
```

4. rsirfo モードで flatten を有効化して実行する。

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
 --opt-mode hess --flatten --out-dir ./result_tsopt_flatten
```

## 使用法
```bash
pdb2reaction tsopt -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m 2S+1] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [--opt-mode grad|hess|dimer|rsirfo] [--flatten/--no-flatten] \
 [--freeze-links/--no-freeze-links] [--max-cycles N] [--thresh PRESET] \
 [--hessian-calc-mode Analytical|FiniteDifference] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

### 例
```bash
# 推奨ベースライン: 電荷/多重度を指定し rsirfo を使う
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode hess --out-dir ./result_tsopt/

# YAML 上書き、有限差分ヘシアン、freeze-links処理を伴う dimer モード
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --freeze-links \
 --opt-mode grad --max-cycles 10000 --no-dump \
 --out-dir ./result_tsopt/ --config ./args.yaml \
 --hessian-calc-mode FiniteDifference

# YAMLのみで駆動される rsirfo モード（RS-I-RFO）
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode hess \
 --config ./args.yaml --out-dir ./result_tsopt/
```

## ワークフロー
- **電荷/スピン解決**: 電荷の解決順序の詳細は [CLI 規約: 電荷の指定](cli-conventions.md#電荷の指定) を参照してください。
- **構造ロードと freeze-links**: `--freeze-links` が有効な場合、リンク水素の親原子は自動的に凍結されます（[リンク水素と凍結原子](extract.md#リンク水素と凍結原子) を参照）。
- **MLIP ヘシアン（デフォルト: UMA）**: `--hessian-calc-mode` で解析的ヘシアンと有限差分ヘシアンを切り替えます。いずれも活性（PHVA）部分空間を考慮します。凍結原子が存在する場合、MLIP バックエンドは活性ブロックのみを返すことがあります。ヘシアン評価モードの詳細は [MLIP 計算機](uma-pysis.md#ヘシアンモード) を参照してください。
- **Dimerモード詳細**:
 - Hessian Guided Dimer段階は、正確ヘシアン（活性サブスペース、TR射影）を周期的に評価してダイマー方向を更新します。`root == 0` のときは最小固有対に `torch.lobpcg` を優先し、失敗時は `torch.linalg.eigh` にフォールバックします。
 - `--flatten` が有効な場合、フラット化ループはΔxとΔgを用い、Bofill（SR1/MS ↔ PSBブレンド; `hessian_dimer.flatten_loop_bofill` で切替）で活性ヘシアンを更新します。各ループは虚振動数モード推定 → 1回フラット化 → ダイマー方向再更新 → dimer+LBFGSマイクロ区間 → （任意で）Bofill更新を実行します。虚振動数モードが1つになったら最終的な正確ヘシアンで振動解析を行います。
 - `root != 0` の場合は初期ダイマー方向のみそのrootを使用し、以降の更新は最も負のモード（`root = 0`）に従います。
- **RS-I-RFOモード**: RS-I-RFOを実行し、任意のヘシアン参照やR+S分割セーフガード、マイクロサイクル制御は `rsirfo` セクションで設定します。`--flatten` が有効で収束後も虚振動数モードが複数残る場合、追加モードをフラット化してRS-I-RFOを再実行し、虚振動数モードが1つになるか上限に達するまで繰り返します。
- **モード出力と変換**: 検出された虚振動数モードはすべて `vib/imag_*_trj.xyz` に書き出されます。PDB 入力で変換が有効な場合は `.pdb` としても出力されます。最適化軌跡と最終構造は `--dump` 時に入力テンプレート経由で PDB に変換されます。Gaussian テンプレートでは最終構造のみ `.gjf` が生成されます。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | 総電荷。`.gjf` テンプレートまたは `--ligand-charge`（PDB 入力または `--ref-pdb` 付きXYZ/GJF）が提供しない限り必須。両方指定時は `-q` が優先 | テンプレート/導出が適用されない限り必須 |
| `-l, --ligand-charge TEXT` | 残基別電荷マッピング（例: `GPP:-3,SAM:1`）。PDB の残基電荷から全系の電荷を自動導出します（手動計算不要）。`-q` 省略時に使用（PDB 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `--workers INT` | MLIP予測器の並列度（workers > 1 で解析ヘシアン無効） | `1` |
| `--workers-per-node INT` | ノードあたりのワーカー数。並列予測器に渡されます | `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| `--freeze-links/--no-freeze-links` | PDBのみ。リンク水素の親を凍結（`geom.freeze_atoms` にマージ） | `True` |
| `--max-cycles INT` | `opt.max_cycles` に渡されるマクロサイクル上限 | `10000` |
| `--opt-mode TEXT` | 最適化モード: `grad`（`dimer`）または `hess`（`rsirfo`）。`dimer`/`rsirfo` も指定可。 | `hess` |
| `--dump/--no-dump` | 軌跡をダンプ | `False` |
| `-o, --out-dir TEXT` | 出力ディレクトリ | `./result_tsopt/` |
| `--thresh TEXT` | 収束プリセットの上書き（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`、`never`） | `baker` |
| `--flatten/--no-flatten` | 余分な虚振動数モードのフラット化ループを有効化（`False` は `flatten_max_iter=0` を強制）。dimer（dimerループ）/rsirfo（RS-IRFO後）に適用 | `False` |
| `--hessian-calc-mode CHOICE` | MLIPヘシアンモード（`Analytical` または `FiniteDifference`） | `FiniteDifference` |
| `--convert-files/--no-convert-files` | PDB または Gaussian 入力用の XYZ/TRJ → PDB/GJF コンパニオン出力を切り替え | `True` |
| `--ref-pdb FILE` | 入力がXYZ/GJFの場合に使用する参照 PDB トポロジー | _None_ |
| `--config FILE` | 明示 CLI オプションより前に適用するベース YAML 設定ファイル | _None_ |
| `--show-config/--no-show-config` | 解決後の設定レイヤーを表示して実行を継続 | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--solvent TEXT` | xTB 暗黙溶媒（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |
| `--dry-run/--no-dry-run` | 実行せずに入力/設定を検証し、実行計画を表示。 | `False` |

```{note}
**`--flatten` はデフォルトで無効です。** `defaults.py` では `flatten_max_iter=50` と定義されていますが、CLI の初期化で `--flatten` が明示的に渡されない限り `flatten_max_iter=0` が設定されます。TS 候補に複数の虚振動数がある場合は、`--flatten` を追加して余分なモードの除去ループを有効にしてください。
```

## 出力
```
out_dir/ (デフォルト:./result_tsopt/)
├─ final_geometry.xyz # 常に書き込み
├─ final_geometry.pdb # 入力がPDBの場合（変換有効時）
├─ final_geometry.gjf # 入力がGaussianの場合（変換有効時）
├─ optimization_all_trj.xyz # --dumpがTrueのときの dimer モードダンプ
├─ optimization_all.pdb # PDB 入力の dimer モードPDB コンパニオン（変換有効時、--dump）
├─ optimization_trj.xyz # --dumpがTrueのときの rsirfo モード軌跡
├─ optimization.pdb # rsirfo モードPDB コンパニオン（変換有効時、--dump）
├─ vib/
│ ├─ imag_±XXXX.Xcm-1_trj.xyz
│ └─ imag_±XXXX.Xcm-1.pdb
└─.dimer_mode.dat # dimer モード方向シード
```

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。
- `--opt-mode` はワークフロー選択用です（デフォルト: `hess` = `rsirfo`）。YAML キーを手動で変更するのではなく、目的のアルゴリズムに合ったモードを選択してください。
- 虚振動数モード検出の閾値はデフォルトで約 5 cm⁻¹（`hessian_dimer.neg_freq_thresh_cm` で変更可能）。複数残る場合は `root` がどの虚振動数モードを追跡するかに影響します。
- 設定の優先順位は [CLI 規約: 設定の優先順位](cli-conventions.md#設定の優先順位) を参照してください。
- PHVAの並進/回転射影は `freq` と同じ実装を使用し、GPU メモリ消費を抑えつつ、活性空間の正しい固有ベクトルを保持します。

共通セクションについては [YAML リファレンス](yaml-reference.md) を参照してください。必要な値だけ変更することを推奨します。

### 共通設定（両モード共通）

```yaml
geom:
 coord_type: cart # coordinate type: cartesian vs dlc internals
 freeze_atoms: [] # 1-based frozen atoms merged with CLI/link detection
calc:
 charge: 0 # total charge (CLI/template override)
 spin: 1 # spin multiplicity 2S+1
 model: uma-s-1p1 # uma-s-1p1 | uma-m-1p1
 task_name: omol # UMA task name
 device: auto # UMA device selection
 max_neigh: null # maximum neighbors for graph construction
 radius: null # cutoff radius for neighbor search
 r_edges: false # store radial edges
 out_hess_torch: true # request torch-form Hessian
 freeze_atoms: null # calculator-level frozen atoms
 hessian_calc_mode: FiniteDifference # Hessian mode selection
 return_partial_hessian: true # partial Hessian (active-DOF only)
opt:
 thresh: baker # convergence preset (Gaussian/Baker-style)
 max_cycles: 10000 # optimizer cycle cap
 print_every: 100 # logging stride
 min_step_norm: 1.0e-08 # minimum norm for step acceptance
 assert_min_step: true # stop if steps fall below threshold
 rms_force: null # explicit RMS force target
 rms_force_only: false # rely only on RMS force convergence
 max_force_only: false # rely only on max force convergence
 force_only: false # skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # geom RMS threshold when converging to ref
 overachieve_factor: 0.0 # factor to tighten thresholds
 check_eigval_structure: false # validate Hessian eigenstructure
 line_search: true # enable line search
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 prefix: "" # filename prefix
 out_dir: ./result_tsopt/ # output directory
```

### Dimer モード（`--opt-mode grad`）

`--opt-mode grad`（Hessian Guided Dimer + LBFGS）で使用します。

```yaml
hessian_dimer:
 thresh_loose: gau_loose # loose convergence preset
 thresh: baker # main convergence preset
 update_interval_hessian: 500 # Hessian rebuild cadence
 neg_freq_thresh_cm: 5.0 # negative frequency threshold (cm⁻¹)
 flatten_amp_ang: 0.1 # flattening amplitude (Å)
 flatten_max_iter: 50 # flattening iteration cap (disabled when --no-flatten)
 flatten_sep_cutoff: 0.0 # minimum distance between representative atoms (Å)
 flatten_k: 10 # representative atoms sampled per mode
 flatten_loop_bofill: false # Bofill update for flatten displacements
 mem: 100000 # memory limit for solver
 device: auto # device selection for eigensolver
 root: 0 # targeted TS root index
 dimer:
 length: 0.0189 # dimer separation (Bohr)
 rotation_max_cycles: 15 # max rotation iterations
 rotation_method: fourier # rotation optimizer method
 rotation_thresh: 0.0001 # rotation convergence threshold
 rotation_tol: 1 # rotation tolerance factor
 rotation_max_element: 0.001 # max rotation matrix element
 rotation_interpolate: true # interpolate rotation steps
 rotation_disable: false # disable rotations entirely
 rotation_disable_pos_curv: true # disable when positive curvature detected
 rotation_remove_trans: true # remove translational components
 trans_force_f_perp: true # project forces perpendicular to translation
 bonds: null # bond list for constraints
 N_hessian: null # Hessian size override
 bias_rotation: false # bias rotational search
 bias_translation: false # bias translational search
 bias_gaussian_dot: 0.1 # Gaussian bias dot product
 seed: null # RNG seed for rotations
 write_orientations: true # write rotation orientations
 forward_hessian: true # propagate Hessian forward
 lbfgs:
 thresh: baker # LBFGS convergence preset
 max_cycles: 10000 # iteration limit
 print_every: 100 # logging stride
 min_step_norm: 1.0e-08 # minimum accepted step norm
 assert_min_step: true # assert when steps stagnate
 rms_force: null # explicit RMS force target
 rms_force_only: false # rely only on RMS force convergence
 max_force_only: false # rely only on max force convergence
 force_only: false # skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # RMS threshold when targeting geometry
 overachieve_factor: 0.0 # tighten thresholds
 check_eigval_structure: false # validate Hessian eigenstructure
 line_search: true # enable line search
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 prefix: "" # filename prefix
 out_dir: ./result_tsopt/ # output directory (defaults.py default is ./result_opt/, overridden to ./result_tsopt/ by tsopt at runtime)
 keep_last: 7 # history size for LBFGS buffers
 beta: 1.0 # initial damping beta
 gamma_mult: false # multiplicative gamma update toggle
 max_step: 0.3 # maximum step length
 control_step: true # control step length adaptively
 double_damp: true # double damping safeguard
 mu_reg: null # regularization strength
 max_mu_reg_adaptions: 10 # cap on mu adaptations
```

### RS-I-RFO モード（`--opt-mode hess`、デフォルト）

`--opt-mode hess`（RS-I-RFO、デフォルト）で使用します。

```yaml
rsirfo:
 thresh: baker # RS-IRFO convergence preset
 max_cycles: 10000 # iteration cap
 trust_radius: 0.10 # trust-region radius
 trust_update: true # enable trust-region updates
 trust_min: 0.0001 # minimum trust radius
 trust_max: 0.20 # maximum trust radius
 print_every: 100 # logging stride
 min_step_norm: 1.0e-08 # minimum accepted step norm
 assert_min_step: true # assert when steps stagnate
 rms_force: null # explicit RMS force target
 rms_force_only: false # rely only on RMS force convergence
 max_force_only: false # rely only on max force convergence
 force_only: false # skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # RMS threshold when targeting geometry
 overachieve_factor: 0.0 # tighten thresholds
 check_eigval_structure: false # validate Hessian eigenstructure
 line_search: true # enable line search
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 prefix: "" # filename prefix
 out_dir: ./result_tsopt/ # output directory (defaults.py default is ./result_opt/, overridden to ./result_tsopt/ by tsopt at runtime)
 roots: [0] # target root indices
 hessian_ref: null # reference Hessian
 rx_modes: null # reaction-mode definitions for projection
 prim_coord: null # primary coordinates to monitor
 rx_coords: null # reaction coordinates to monitor
 hessian_update: bofill # Hessian update scheme override
 hessian_recalc_reset: true # reset recalc counter after exact Hessian
 max_micro_cycles: 50 # micro-iterations per macro cycle
 augment_bonds: false # augment reaction path based on bond analysis
 min_line_search: true # enforce minimum line-search step
 max_line_search: true # enforce maximum line-search step
 assert_neg_eigval: false # require a negative eigenvalue at convergence
```

---

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け

- [path-search](path-search.md) — TS 候補（HEI）を特定するMEP 探索
- [irc](irc.md) — 最適化されたTSからの反応経路追跡
- [freq](freq.md) — 完全な振動解析と熱化学補正（虚振動数チェックは `tsopt` が内部で実行済み）
- [all](all.md) — 抽出 → MEP → tsopt → IRC → freq を連鎖する一気通貫ワークフロー
- [YAML リファレンス](yaml-reference.md) — `hessian_dimer`（Hessian Guided Dimer）と `rsirfo` の完全な設定オプション
- [用語集](glossary.md) — TS、Dimer、RS-I-RFO、ヘシアンの定義
