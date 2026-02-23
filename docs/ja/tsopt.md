# `tsopt`

## 概要

> **要約:** 遷移状態（TS）*候補*を、Dimer（`--opt-mode light`）または RS‑I‑RFO（`--opt-mode heavy`、デフォルト）で最適化します。VRAM に余裕がある場合は `--hessian-calc-mode Analytical` により高速化できることが多いです。妥当な TS では虚振動数が **ちょうど 1 つ**であるべきです。必ず freq/IRC でモードと接続性を確認してください。

### 要点
- **入力:** `path-opt` / `path-search` が出力する HEI、または自前の TS 初期構造（`geom_loader` が扱える形式）。
- **モード:** `heavy` = RS‑I‑RFO（既定、一般的により堅牢）。`light` = Hessian Dimer（1ステップあたりのコストが低いことが多い）。
- **品質確認:** 最適化後も TS は *候補* です。[freq](freq.md) と [irc](irc.md) でモードと接続性を確認してください。
- **任意の後処理:** `--flatten-imag-mode` は、収束後に余分な虚数モードが残っている場合にその除去を試みます。
- **出力変換:** `--convert-files`（デフォルト）で、PDB 入力は（`--dump` のとき）`.pdb` を併記し、Gaussian テンプレートは最終構造の `.gjf` を書き出します。

### `--opt-mode` の選び方
- **`--opt-mode heavy`（RS‑I‑RFO）**: デフォルト。ヘシアン計算のコストを許容でき、堅牢性を重視する場合に推奨。
- **`--opt-mode light`（Dimer）**: 軽量な探索を行いたい場合や、複数の TS 初期構造から素早く反復したい場合に有効。

XYZ/GJF 入力では `--ref-pdb` により参照 PDB トポロジーを与え、XYZ 座標を保持したまま PDB/GJF へのフォーマット対応変換ができます。TS 初期構造が必要な場合は、2 端点なら [path-opt](path_opt.md)、2 構造以上なら [path-search](path_search.md) で HEI を得てから `tsopt` → `freq` → `irc` の順で検証してください。

## 最小例

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --out-dir ./result_tsopt
```

## 出力の見方

- `result_tsopt/final_geometry.pdb`（または `final_geometry.xyz`）
- `result_tsopt/vib/final_imag_mode_*_trj.xyz`
- `result_tsopt/vib/final_imag_mode_*.pdb`（PDB 入力の場合）

## よくある例

1. VRAM に余裕がある場合に light モード + 解析的ヘシアンで実行する。

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
 --opt-mode light --hessian-calc-mode Analytical --out-dir ./result_tsopt_light
```

2. 最適化軌跡を保存して確認する。

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --dump --out-dir ./result_tsopt_dump
```

3. heavy モードを YAML 上書きと併用する。

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
 --opt-mode heavy --config tsopt.yaml --out-dir ./result_tsopt_heavy
```

## 使用法
```bash
pdb2reaction tsopt -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m 2S+1] \
 [--opt-mode light|heavy] [--flatten-imag-mode/--no-flatten-imag-mode] \
 [--freeze-links/--no-freeze-links] [--max-cycles N] [--thresh PRESET] \
 [--hessian-calc-mode Analytical|FiniteDifference] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

### 例
```bash
# 推奨ベースライン: 電荷/多重度を指定しlightワークフローを選択
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode light --out-dir ./result_tsopt/

# YAML 上書き、有限差分ヘシアン、freeze-links処理を伴うlightモード
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --freeze-links \
 --opt-mode light --max-cycles 10000 --no-dump \
 --out-dir ./result_tsopt/ --config ./args.yaml \
 --hessian-calc-mode FiniteDifference

# YAMLのみで駆動されるheavyモード（RS-I-RFO）
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode heavy \
 --config ./args.yaml --out-dir ./result_tsopt/
```

## ワークフロー
- **電荷/スピン解決**: 入力が `.gjf` の場合、電荷と多重度はテンプレート値を継承。`-q` が省略され `--ligand-charge` が与えられている場合、構造は酵素–基質複合体として扱われ、PDB 入力（または `--ref-pdb` 付きXYZ/GJF）で `extract.py` の電荷サマリーから総電荷を導出。明示的な `-q` は常に優先されます。テンプレートや導出が適用されない場合は `-q/--charge` が必須で、多重度は省略時 `1` です。
- **構造ロード & freeze-links**: 構造は `pysisyphus.helpers.geom_loader` を介して読み込まれます。PDB 入力では `--freeze-links` がリンク水素を検出して親原子を凍結し、`geom.freeze_atoms` にマージしてログに表示します。凍結原子はUMAの `calc.freeze_atoms` にも伝播します。
- **UMAヘシアン**: `--hessian-calc-mode` は解析的評価と有限差分評価を切り替えます。凍結原子がある場合、UMA は活性ブロックのみを返すことがあります。VRAMが十分な場合は `--hessian-calc-mode` を `Analytical` に設定することを強く推奨します。
- **Lightモード詳細**:
 - Hessian Dimer段階は、正確ヘシアン（活性サブスペース、TR射影）を周期的に評価してダイマー方向を更新します。`root == 0` のときは最小固有対に `torch.lobpcg` を優先し、失敗時は `torch.linalg.eigh` にフォールバックします。
 - `--flatten-imag-mode` が有効な場合、フラットンループはΔxとΔgを用い、Bofill（SR1/MS ↔ PSBブレンド; `hessian_dimer.flatten_loop_bofill` で切替）で活性ヘシアンを更新します。各ループは虚数モード推定 → 1回フラットン → ダイマー方向再更新 → dimer+LBFGSマイクロ区間 → （任意で）Bofill更新を実行します。虚数モードが1つになったら最終的な正確ヘシアンで周波数解析を行います。
 - `root != 0` の場合は初期ダイマー方向のみそのrootを使用し、以降の更新は最も負のモード（`root = 0`）に従います。
- **Heavyモード（RS-I-RFO）**: RS-I-RFOを実行し、任意のヘシアン参照やR+S分割セーフガード、マイクロサイクル制御は `rsirfo` セクションで設定します。`--flatten-imag-mode` が有効で収束後も虚数モードが複数残る場合、追加モードをフラットンしてRS-I-RFOを再実行し、虚数モードが1つになるか上限に達するまで繰り返します。
- **モード出力 & 変換**: 収束した虚数モードは常に `vib/final_imag_mode_*_trj.xyz` に書き出され、PDB 入力で変換が有効な場合は `.pdb` にもミラーされます。最適化軌跡と最終構造は、`--dump` のとき入力テンプレート経由で PDB に変換されます。PDB 入力で変換が有効な場合は `.pdb` にもミラーされ、Gaussian テンプレートでは最終構造のみ `.gjf` が生成されます。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | 総電荷。`.gjf` テンプレートまたは `--ligand-charge`（PDB 入力または `--ref-pdb` 付きXYZ/GJF）が提供しない限り必須。両方指定時は `-q` が優先 | テンプレート/導出が適用されない限り必須 |
| `--ligand-charge TEXT` | `-q` 省略時に使用する総電荷または残基名ごとのマッピング。PDB 入力（または `--ref-pdb` 付きXYZ/GJF）でextract方式の電荷導出を有効化 | _None_ |
| `--workers`, `--workers-per-node` | UMA予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| `--freeze-links/--no-freeze-links` | PDBのみ。リンク水素の親を凍結（`geom.freeze_atoms` にマージ） | `True` |
| `--max-cycles INT` | `opt.max_cycles` に転送されるマクロサイクル上限 | `10000` |
| `--opt-mode TEXT` | 上記のLight/Heavyエイリアス | `heavy` |
| `--dump/--no-dump` | 軌跡をダンプ | `False` |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_tsopt/` |
| `--thresh TEXT` | 収束プリセットの上書き（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`、`never`） | `baker` |
| `--flatten-imag-mode/--no-flatten-imag-mode` | 余分な虚数モードフラットンループを有効化（`False` は `flatten_max_iter=0` を強制）。lightモード（dimerループ）とheavyモード（RS-IRFO後）の両方に適用 | `False` |
| `--hessian-calc-mode CHOICE` | UMAヘシアンモード（`Analytical` または `FiniteDifference`） | `FiniteDifference` |
| `--convert-files/--no-convert-files` | PDB または Gaussian 入力用の XYZ/TRJ → PDB/GJF コンパニオン出力を切り替え | `True` |
| `--ref-pdb FILE` | 入力がXYZ/GJFの場合に使用する参照 PDB トポロジー | _None_ |
| `--config FILE` | 明示 CLI オプションより前に適用するベース YAML 設定ファイル | _None_ |
| `--show-config/--no-show-config` | 解決後の設定レイヤーを表示して実行を継続 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに入力/設定を検証し、実行計画を表示 | `False` |

## 出力
```
out_dir/ (デフォルト:./result_tsopt/)
├─ final_geometry.xyz # 常に書き込み
├─ final_geometry.pdb # 入力がPDBの場合（変換有効時）
├─ final_geometry.gjf # 入力がGaussianの場合（変換有効時）
├─ optimization_all_trj.xyz # --dumpがTrueのときのLightモードダンプ
├─ optimization_all.pdb # PDB 入力のLightモードPDB コンパニオン（変換有効時、--dump）
├─ optimization_trj.xyz # --dumpがTrueのときのHeavyモード軌跡
├─ optimization.pdb # HeavyモードPDB コンパニオン（変換有効時、--dump）
├─ vib/
│ ├─ final_imag_mode_±XXXX.Xcm-1_trj.xyz
│ └─ final_imag_mode_±XXXX.Xcm-1.pdb
└─.dimer_mode.dat # Lightモード方向シード
```

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- `--opt-mode` のエイリアスは上記のワークフローに対応します。YAML キーを手動で変更するのではなく、目的のアルゴリズムに合ったモードを選択してください（デフォルト: `heavy`）。
- 虚数モード検出の閾値はデフォルトで約 5 cm⁻¹（`hessian_dimer.neg_freq_thresh_cm` で設定可能）。複数残る場合は `root` がどの虚数モードを出力するかに影響します。
- 設定マージ優先順位は `defaults < config < 明示 CLI < override` です。
- - PHVAの並進/回転射影は `freq` と同じ実装を使用し、GPU メモリ消費を抑えつつ、活性空間の正しい固有ベクトルを保持します。



`defaults < config < 明示 CLI < override`。

共通セクションについては [YAML リファレンス](yaml_reference.md) を参照してください。下記ブロックが既にワークフローに合っている場合は、必要な値だけ変更することを推奨します。

```yaml
geom:
 coord_type: cart # coordinate type: cartesian vs dlc internals
 freeze_atoms: [] # 0-based frozen atoms merged with CLI/link detection
calc:
 charge: 0 # total charge (CLI/template override)
 spin: 1 # spin multiplicity 2S+1
 model: uma-s-1p1 # UMA model tag
 task_name: omol # UMA task name
 device: auto # UMA device selection
 max_neigh: null # maximum neighbors for graph construction
 radius: null # cutoff radius for neighbor search
 r_edges: false # store radial edges
 out_hess_torch: true # request torch-form Hessian
 freeze_atoms: null # calculator-level frozen atoms
 hessian_calc_mode: FiniteDifference # Hessian mode selection
 return_partial_hessian: false # full Hessian (avoids shape mismatches)
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
 out_dir:./result_tsopt/ # output directory
hessian_dimer:
 thresh_loose: gau_loose # loose convergence preset
 thresh: baker # main convergence preset
 update_interval_hessian: 500 # Hessian rebuild cadence
 neg_freq_thresh_cm: 5.0 # negative frequency threshold (cm^-1)
 flatten_amp_ang: 0.1 # flattening amplitude (Å)
 flatten_max_iter: 50 # flattening iteration cap (disabled when --no-flatten-imag-mode)
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
 out_dir:./result_tsopt/ # output directory
 keep_last: 7 # history size for LBFGS buffers
 beta: 1.0 # initial damping beta
 gamma_mult: false # multiplicative gamma update toggle
 max_step: 0.3 # maximum step length
 control_step: true # control step length adaptively
 double_damp: true # double damping safeguard
 mu_reg: null # regularization strength
 max_mu_reg_adaptions: 10 # cap on mu adaptations
rsirfo:
 thresh: baker # RS-IRFO convergence preset
 max_cycles: 10000 # iteration cap
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
 out_dir:./result_tsopt/ # output directory
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

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け

- [path-search](path_search.md) — TS 候補（HEI）を特定するMEP 探索
- [irc](irc.md) — 最適化されたTSからの反応経路追跡
- [freq](freq.md) — 虚振動数が 1 つのみであることを確認（妥当な TS の条件）
- [all](all.md) — 抽出 → MEP → tsopt → IRC → freq を連鎖するエンドツーエンドワークフロー
- [YAML リファレンス](yaml_reference.md) — `hessian_dimer` と `rsirfo` の完全な設定オプション
- [用語集](glossary.md) — TS、Dimer、RS-I-RFO、ヘシアンの定義
