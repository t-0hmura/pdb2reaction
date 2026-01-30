# `tsopt` サブコマンド

## 概要
`pdb2reaction tsopt` は2つの補完的なワークフローを使用して遷移状態を最適化します:

- **light** モード: 定期的な正確ヘシアン更新を伴うHessian Dimer探索、余分な虚数モードを除去するためのオプションのメモリ効率的なフラットンループ（デフォルトで無効）、活性自由度のPHVA対応ヘシアン更新
- **heavy** モード: 設定可能な信頼領域セーフガードを持つRS-I-RFOオプティマイザー、収束後に余分な虚数モードが残る場合のオプションの後最適化フラットンループ

両モードはエネルギー/勾配/ヘシアンにUMA計算機を使用し、YAMLから `geom`/`calc`/`opt` 設定を継承し、最終的な虚数モードを常に `.trj` に書き込みます。`--convert-files`（デフォルト有効）を有効にすると、PDB入力は軌跡を `.pdb` コンパニオンにミラーし、Gaussianテンプレートは最終構造の `.gjf` を出力します（軌跡は `.gjf` に変換されません）。XYZ/GJF入力では `--ref-pdb` が参照PDBトポロジーを提供しXYZ座標を保持するため、PDB/GJFへのフォーマット対応変換が可能です。デフォルトの `--opt-mode` は **heavy**（RS-I-RFO）です; Hessian Dimerワークフローを実行するには `--opt-mode light` に切り替えてください。

## 使用法
```bash
pdb2reaction tsopt -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m 2S+1] \
                    [--opt-mode light|heavy] [--flatten-imag-mode {True\|False}] \
                    [--freeze-links {True\|False}] [--max-cycles N] [--thresh PRESET] \
                    [--dump {True\|False}] [--out-dir DIR] [--args-yaml FILE] \
                    [--hessian-calc-mode Analytical|FiniteDifference] \
                    [--convert-files {True\|False}] [--ref-pdb FILE]
```

### 例
```bash
# 推奨ベースライン: 電荷/多重度を指定しlightワークフローを選択
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode light --out-dir ./result_tsopt/

# YAMLオーバーライド、有限差分ヘシアン、freeze-links処理を伴うlightモード
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --freeze-links True \
    --opt-mode light --max-cycles 10000 --dump False \
    --out-dir ./result_tsopt/ --args-yaml ./args.yaml \
    --hessian-calc-mode FiniteDifference

# YAMLのみで駆動されるheavyモード（RS-I-RFO）
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode heavy \
    --args-yaml ./args.yaml --out-dir ./result_tsopt/
```

## ワークフロー
- **電荷/スピン解決**: 入力が `.gjf` の場合、電荷と多重度はテンプレート値を継承。`-q` が省略され `--ligand-charge` が与えられている場合、構造は酵素–基質複合体として扱われ、PDB入力（または `--ref-pdb` 付きXYZ/GJF）で `extract.py` の電荷サマリーから総電荷を導出。明示的な `-q` は常に優先されます。テンプレート/導出が使えない場合は `-q/--charge` が必須で、多重度は省略時 `1` です。
- **構造ロード & freeze-links**: 構造は `pysisyphus.helpers.geom_loader` を介して読み込まれます。PDB入力では `--freeze-links True` がリンク水素を検出して親原子を凍結し、`geom.freeze_atoms` にマージしてログに表示します。凍結原子はUMAの `calc.freeze_atoms` にも伝播します。
- **UMAヘシアン**: `--hessian-calc-mode` は解析的評価と有限差分評価を切り替えます。凍結原子がある場合、UMAは活性サブスペースの部分ヘシアンのみを返すことがあります。VRAMが十分な場合は `--hessian-calc-mode` を `Analytical` に設定することを強く推奨します。
- **Lightモード詳細**:
  - Hessian Dimer段階は、正確ヘシアン（活性サブスペース、TR射影）を周期的に評価してダイマー方向を更新します。`root == 0` のときは最小固有対に `torch.lobpcg` を優先し、失敗時は `torch.linalg.eigh` にフォールバックします。
  - `--flatten-imag-mode` が有効な場合、フラットンループはΔxとΔgを用い、Bofill（SR1/MS ↔ PSBブレンド; `hessian_dimer.flatten_loop_bofill` で切替）で活性ヘシアンを更新します。各ループは虚数モード推定 → 1回フラットン → ダイマー方向再更新 → dimer+LBFGSマイクロ区間 → （任意で）Bofill更新を実行します。虚数モードが1つになったら最終的な正確ヘシアンで周波数解析を行います。
  - `root != 0` の場合は初期ダイマー方向のみそのrootを使用し、以降の更新は最も負のモード（`root = 0`）に従います。
- **Heavyモード（RS-I-RFO）**: RS-I-RFOを実行し、任意のヘシアン参照やR+S分割セーフガード、マイクロサイクル制御は `rsirfo` セクションで設定します。`--flatten-imag-mode` が有効で収束後も虚数モードが複数残る場合、追加モードをフラットンしてRS-I-RFOを再実行し、虚数モードが1つになるか上限に達するまで繰り返します。
- **モード出力 & 変換**: 収束した虚数モードは常に `vib/final_imag_mode_*.trj` に書き出され、PDB入力で変換が有効な場合は `.pdb` にもミラーされます。最適化軌跡と最終構造は `--dump True` のときPDBに変換され、Gaussianテンプレートでは最終構造のみ `.gjf` が生成されます。

## CLIオプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | 総電荷。`.gjf` テンプレートまたは `--ligand-charge`（PDB入力または `--ref-pdb` 付きXYZ/GJF）が提供しない限り必須。両方指定時は `-q` が優先 | テンプレート/導出が適用されない限り必須 |
| `--ligand-charge TEXT` | `-q` 省略時に使用する総電荷または残基名ごとのマッピング。PDB入力（または `--ref-pdb` 付きXYZ/GJF）でextract方式の電荷導出を有効化 | _None_ |
| `--workers`, `--workers-per-node` | UMA予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| `--freeze-links {True\|False}` | PDBのみ。リンク水素の親を凍結（`geom.freeze_atoms` にマージ） | `True` |
| `--max-cycles INT` | `opt.max_cycles` に転送されるマクロサイクル上限 | `10000` |
| `--opt-mode TEXT` | 上記のLight/Heavyエイリアス | `heavy` |
| `--dump {True\|False}` | 軌跡をダンプ | `False` |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_tsopt/` |
| `--thresh TEXT` | 収束プリセットのオーバーライド | `baker` |
| `--flatten-imag-mode {True\|False}` | 余分な虚数モードフラットンループを有効化 | `False` |
| `--hessian-calc-mode CHOICE` | UMAヘシアンモード（`Analytical` または `FiniteDifference`） | `FiniteDifference` |
| `--convert-files {True\|False}` | PDBまたはGaussian入力用のXYZ/TRJ → PDB/GJFコンパニオンをトグル | `True` |
| `--ref-pdb FILE` | 入力がXYZ/GJFの場合に使用する参照PDBトポロジー | _None_ |
| `--args-yaml FILE` | YAMLオーバーライド（`geom`、`calc`、`opt`、`hessian_dimer`、`rsirfo`） | _None_ |

## 出力（& ディレクトリレイアウト）
```
out_dir/ (デフォルト: ./result_tsopt/)
├─ final_geometry.xyz            # 常に書き込み
├─ final_geometry.pdb            # 入力がPDBの場合（変換有効時）
├─ final_geometry.gjf            # 入力がGaussianの場合（変換有効時）
├─ optimization_all.trj          # --dumpがTrueのときのLightモードダンプ
├─ optimization_all.pdb          # PDB入力のLightモードPDBコンパニオン（変換有効時、--dump True）
├─ optimization.trj              # --dumpがTrueのときのHeavyモード軌跡
├─ optimization.pdb              # HeavyモードPDBコンパニオン（変換有効時、--dump True）
├─ vib/
│  ├─ final_imag_mode_±XXXX.Xcm-1.trj
│  └─ final_imag_mode_±XXXX.Xcm-1.pdb
└─ .dimer_mode.dat               # Lightモード方向シード
```

## 注意事項
- `--opt-mode` エイリアスは上記のワークフローに正確にマップされる; YAMLキーを手動で調整するよりも意図したアルゴリズム用に1つを選択（デフォルト: `heavy`）
- 虚数モード検出は〜5 cm⁻¹がデフォルト（`hessian_dimer.neg_freq_thresh_cm` で設定可能）。複数残る場合は `root` がどの虚数モードを出力するかに影響します。
- `--hessian-calc-mode` はYAMLマージ後に `calc.hessian_calc_mode` をオーバーライド
- PHVAの並進/回転射影は `freq` と同じ実装を使用し、GPUメモリ消費を抑えつつ活性空間の固有ベクトルを保持します。


## YAML設定（`--args-yaml`）
YAMLはマッピングで指定します。YAMLはCLIを上書きします。共通セクションは [YAMLリファレンス](yaml-reference.md) を再利用してください。下記ブロックが既にワークフローに合っている場合は、必要な値だけ変更することを推奨します。

```yaml
geom:
  coord_type: cart           # coordinate type: cartesian vs dlc internals
  freeze_atoms: []           # 0-based frozen atoms merged with CLI/link detection
calc:
  charge: 0                  # total charge (CLI/template override)
  spin: 1                    # spin multiplicity 2S+1
  model: uma-s-1p1           # UMA model tag
  task_name: omol            # UMA task name
  device: auto               # UMA device selection
  max_neigh: null            # maximum neighbors for graph construction
  radius: null               # cutoff radius for neighbor search
  r_edges: false             # store radial edges
  out_hess_torch: true       # request torch-form Hessian
  freeze_atoms: null         # calculator-level frozen atoms
  hessian_calc_mode: FiniteDifference   # Hessian mode selection
  return_partial_hessian: false         # full Hessian (avoids shape mismatches)
opt:
  thresh: baker              # convergence preset (Gaussian/Baker-style)
  max_cycles: 10000          # optimizer cycle cap
  print_every: 100           # logging stride
  min_step_norm: 1.0e-08     # minimum norm for step acceptance
  assert_min_step: true      # stop if steps fall below threshold
  rms_force: null            # explicit RMS force target
  rms_force_only: false      # rely only on RMS force convergence
  max_force_only: false      # rely only on max force convergence
  force_only: false          # skip displacement checks
  converge_to_geom_rms_thresh: 0.05   # geom RMS threshold when converging to ref
  overachieve_factor: 0.0    # factor to tighten thresholds
  check_eigval_structure: false   # validate Hessian eigenstructure
  line_search: true          # enable line search
  dump: false                # dump trajectory/restart data
  dump_restart: false        # dump restart checkpoints
  prefix: ""                 # filename prefix
  out_dir: ./result_tsopt/   # output directory
hessian_dimer:
  thresh_loose: gau_loose    # loose convergence preset
  thresh: baker              # main convergence preset
  update_interval_hessian: 500   # Hessian rebuild cadence
  neg_freq_thresh_cm: 5.0    # negative frequency threshold (cm^-1)
  flatten_amp_ang: 0.1       # flattening amplitude (Å)
  flatten_max_iter: 50       # flattening iteration cap (disabled when --flatten-imag-mode False)
  flatten_sep_cutoff: 0.0    # minimum distance between representative atoms (Å)
  flatten_k: 10              # representative atoms sampled per mode
  flatten_loop_bofill: false # Bofill update for flatten displacements
  mem: 100000                # memory limit for solver
  device: auto               # device selection for eigensolver
  root: 0                    # targeted TS root index
  dimer:
    length: 0.0189           # dimer separation (Bohr)
    rotation_max_cycles: 15  # max rotation iterations
    rotation_method: fourier # rotation optimizer method
    rotation_thresh: 0.0001  # rotation convergence threshold
    rotation_tol: 1          # rotation tolerance factor
    rotation_max_element: 0.001   # max rotation matrix element
    rotation_interpolate: true    # interpolate rotation steps
    rotation_disable: false   # disable rotations entirely
    rotation_disable_pos_curv: true   # disable when positive curvature detected
    rotation_remove_trans: true   # remove translational components
    trans_force_f_perp: true  # project forces perpendicular to translation
    bonds: null               # bond list for constraints
    N_hessian: null           # Hessian size override
    bias_rotation: false      # bias rotational search
    bias_translation: false   # bias translational search
    bias_gaussian_dot: 0.1    # Gaussian bias dot product
    seed: null                # RNG seed for rotations
    write_orientations: true  # write rotation orientations
    forward_hessian: true     # propagate Hessian forward
  lbfgs:
    thresh: baker              # LBFGS convergence preset
    max_cycles: 10000          # iteration limit
    print_every: 100           # logging stride
    min_step_norm: 1.0e-08     # minimum accepted step norm
    assert_min_step: true      # assert when steps stagnate
    rms_force: null            # explicit RMS force target
    rms_force_only: false      # rely only on RMS force convergence
    max_force_only: false      # rely only on max force convergence
    force_only: false          # skip displacement checks
    converge_to_geom_rms_thresh: 0.05   # RMS threshold when targeting geometry
    overachieve_factor: 0.0    # tighten thresholds
    check_eigval_structure: false   # validate Hessian eigenstructure
    line_search: true          # enable line search
    dump: false                # dump trajectory/restart data
    dump_restart: false        # dump restart checkpoints
    prefix: ""                 # filename prefix
    out_dir: ./result_opt/     # output directory
    keep_last: 7               # history size for LBFGS buffers
    beta: 1.0                  # initial damping beta
    gamma_mult: false          # multiplicative gamma update toggle
    max_step: 0.3              # maximum step length
    control_step: true         # control step length adaptively
    double_damp: true          # double damping safeguard
    mu_reg: null               # regularization strength
    max_mu_reg_adaptions: 10   # cap on mu adaptations
rsirfo:
  thresh: baker              # RS-IRFO convergence preset
  max_cycles: 10000          # iteration cap
  print_every: 100           # logging stride
  min_step_norm: 1.0e-08     # minimum accepted step norm
  assert_min_step: true      # assert when steps stagnate
  rms_force: null            # explicit RMS force target
  rms_force_only: false      # rely only on RMS force convergence
  max_force_only: false      # rely only on max force convergence
  force_only: false          # skip displacement checks
  converge_to_geom_rms_thresh: 0.05   # RMS threshold when targeting geometry
  overachieve_factor: 0.0    # tighten thresholds
  check_eigval_structure: false   # validate Hessian eigenstructure
  line_search: true          # enable line search
  dump: false                # dump trajectory/restart data
  dump_restart: false        # dump restart checkpoints
  prefix: ""                 # filename prefix
  out_dir: ./result_opt/     # output directory
  roots: [0]                 # target root indices
  hessian_ref: null          # reference Hessian
  rx_modes: null             # reaction-mode definitions for projection
  prim_coord: null           # primary coordinates to monitor
  rx_coords: null            # reaction coordinates to monitor
  hessian_update: bofill     # Hessian update scheme override
  hessian_recalc_reset: true # reset recalc counter after exact Hessian
  max_micro_cycles: 50       # micro-iterations per macro cycle
  augment_bonds: false       # augment reaction path based on bond analysis
  min_line_search: true      # enforce minimum line-search step
  max_line_search: true      # enforce maximum line-search step
  assert_neg_eigval: false   # require a negative eigenvalue at convergence
```
