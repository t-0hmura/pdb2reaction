# `path-opt` サブコマンド

## 概要
`pdb2reaction path-opt` は、`--mep-mode` で選択されるpysisyphusのGrowing String法（GSM）またはDirect Max Flux（DMF）を使用して、2つのエンドポイント構造間の最小エネルギー経路（MEP）を探索します。UMAはすべてのイメージにエネルギー/勾配/ヘシアンを提供し、外部の剛体アライメントルーチンがオプティマイザー開始前にストリングを整えます。設定は **デフォルト → CLI → `--args-yaml`** の優先順位で `geom`/`calc`/`gs`/`opt` に適用されます。`--convert-files`（デフォルト有効）を有効にすると、PDB参照がある場合は軌跡を `.pdb` コンパニオンへ、Gaussianテンプレートがある場合はXYZスナップショット（例: HEI）を `.gjf` コンパニオンへミラーします。XYZ/GJF入力では `--ref-pdb` が参照PDBトポロジーを提供しXYZ座標を保持するため、PDB変換が可能です。GSMがデフォルトの経路生成器であり、単一構造最適化は `light`（LBFGS）プリセットがデフォルトです。

## 使用法
```bash
pdb2reaction path-opt -i REACTANT.{pdb|xyz} PRODUCT.{pdb|xyz} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m MULT] \
                      [--workers N] [--workers-per-node N] \
                      [--mep-mode {gsm|dmf}] [--freeze-links {True\|False}] [--max-nodes N] [--max-cycles N] \
                      [--climb {True\|False}] [--dump {True\|False}] [--thresh PRESET] \
                      [--preopt {True\|False}] [--preopt-max-cycles N] [--opt-mode light|heavy] [--fix-ends {True\|False}] \
                      [--out-dir DIR] [--args-yaml FILE] \
                      [--convert-files {True\|False}] [--ref-pdb FILE]
```

## ワークフロー
1. **事前アライメント & 凍結解決**
   - 最初以降のすべてのエンドポイントは最初の構造にKabschアライメントされます。いずれかのエンドポイントで `freeze_atoms` が定義されている場合、その原子のみでRMSDフィットし、変換は全原子に適用されます。
   - `--freeze-links=True`（デフォルト）のPDB入力では、リンク水素の親原子が検出され `freeze_atoms` にマージされます。
   - 明示的なアライメント/リファインがあるため、`StringOptimizer.align` は無効のまま維持されます。

2. **ストリング成長とHEIエクスポート**
   - 経路が成長・精密化された後、最高エネルギー内部局所極大を優先的に探索します。内部局所極大がない場合は内部ノードの最大値へ、内部ノードが無い場合は全体最大へフォールバックします。
   - 最高エネルギーイメージ（HEI）は `.xyz` と、PDB参照がある場合は `.pdb` として書き込み、Gaussianテンプレートがある場合は `.gjf` も出力します（いずれも `--convert-files` を尊重）。

### 主要な挙動
- **エンドポイント**: 入力は2構造のみ。形式は `geom_loader` に準拠。PDB入力（または `--ref-pdb` 付きXYZ/GJF）で軌跡/HEIのPDB出力が有効。
- **電荷/スピン**: CLIが`.gjf`テンプレートを上書き。`-q` 省略時に `--ligand-charge` がある場合、エンドポイントは酵素–基質複合体として扱われ、PDB入力では `extract.py` の電荷サマリーで総電荷を導出。明示的な `-q` が常に優先。テンプレート/導出がない場合は電荷 `0`、スピン `1` にフォールバックします。正しい状態のために明示指定を推奨します。
- **MEPセグメント**: `--max-nodes` はGSM/DMFの内部ノード数を制御（GSMの総画像数は `max_nodes + 2`）。`--thresh` またはYAMLで収束プリセット（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`）を指定。
- **クライミングイメージ**: `--climb` は標準のクライミング手順とLanczos接線リファインの両方を切り替え。
- **ダンプ**: `--dump True` で StringOptimizer の `opt.dump=True` を有効化し、`out_dir` に軌跡/再開情報を出力。
- **終了コード**: `0` 成功、`3` 最適化失敗、`4` 軌跡書き込みエラー、`5` HEI出力エラー、`130` 割り込み、`1` 予期せぬエラー。

## CLIオプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH PATH` | 反応物と生成物構造 | 必須 |
| `-q, --charge INT` | 総電荷（`calc.charge`）。省略時は `.gjf` テンプレートまたは `--ligand-charge`（PDB入力）が供給し、なければ `0` にフォールバック。両方指定時は `-q` が優先 | テンプレート/`0` |
| `--ligand-charge TEXT` | `-q` 省略時に使用する総電荷または残基名ごとのマッピング。PDB入力でextract方式の電荷導出を有効化し、それ以外では電荷は `0` にフォールバック | _None_ |
| `--workers`, `--workers-per-node` | UMA予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度 | テンプレート/`1` |
| `--freeze-links {True\|False}` | PDBのみ: リンクH親を凍結（YAMLとマージ） | `True` |
| `--max-nodes INT` | 内部ノード数（ストリングイメージ = `max_nodes + 2`） | `10` |
| `--mep-mode {gsm\|dmf}` | GSM（ストリングベース）またはDMF（ダイレクトフラックス）経路生成器を選択 | `gsm` |
| `--max-cycles INT` | オプティマイザーマクロイテレーション上限 | `300` |
| `--climb {True\|False}` | クライミングイメージ精密化を有効化 | `True` |
| `--dump {True\|False}` | MEP軌跡/リスタートをダンプ | `False` |
| `--opt-mode TEXT` | エンドポイント事前最適化用の単一構造オプティマイザー | `light` |
| `--convert-files {True\|False}` | PDB/Gaussian入力用のXYZ/TRJ → PDB/GJFコンパニオンをトグル | `True` |
| `--ref-pdb FILE` | XYZ/GJF入力用の参照PDBトポロジー | _None_ |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_path_opt/` |
| `--thresh TEXT` | GSM/ストリングオプティマイザーの収束プリセットをオーバーライド | `gau` |
| `--args-yaml FILE` | YAMLオーバーライド（セクション `geom`、`calc`、`gs`、`opt`） | _None_ |
| `--preopt {True\|False}` | アライメント/MEP探索前に各エンドポイントを事前最適化 | `False` |
| `--fix-ends {True\|False}` | GSM成長/精密化中にエンドポイント構造を固定 | `False` |

## 出力
```
out_dir/
├─ final_geometries.trj        # XYZ経路（コメント行にエネルギーを保持）
├─ final_geometries.pdb        # PDB参照が利用可能で変換が有効な場合
├─ hei.xyz                     # 最高エネルギーイメージ
├─ hei.pdb                     # PDB参照が利用可能な場合のHEI（変換有効時）
├─ hei.gjf                     # Gaussianテンプレートを使用して書き込まれたHEI（変換有効時）
├─ align_refine/               # 剛体アライメント/リファイン段階の中間ファイル（アライメント実行時）
└─ <オプティマイザーダンプ/リスタート>
```
コンソールには解決済みYAMLブロックが出力され、GSM/DMFのMEP進行状況とタイミングが報告されます。


## YAML設定（`--args-yaml`）
YAML値はCLIを上書きし、CLIはデフォルトを上書きします。

### `geom`
- [`opt`](opt.md) と同じキー（`coord_type`, `freeze_atoms` など）。`--freeze-links` がPDB入力で `freeze_atoms` にマージされます。

### `calc`
- UMA計算機の設定（`model`, `device`, 近傍半径, ヘシアンなど）。

### `dmf`
- Direct Max Flux + (C)FB-ENM 補間の制御。CLIで露出している `dmf` ブロックと同じキー。

### `gs`
- Growing String表現の制御: `max_nodes`, `perp_thresh`, 再パラメータ化（`reparam_check`, `reparam_every`, `reparam_every_full`, `param`）、`max_micro_cycles`, DLCリセット、climb関連、scheduler。

### `opt`
- StringOptimizer設定: type, `stop_in_when_full`, `align=False`（固定）、`scale_step`, `max_cycles`, dump系、`reparam_thresh`, `coord_diff_thresh`, `out_dir`, `print_every`。

### YAML例（デフォルト値）
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
gs:
  fix_first: true            # keep the first endpoint fixed during optimization
  fix_last: true             # keep the last endpoint fixed during optimization
  max_nodes: 10              # maximum string nodes
  perp_thresh: 0.005         # perpendicular displacement threshold
  reparam_check: rms         # reparametrization check metric
  reparam_every: 1           # reparametrization stride
  reparam_every_full: 1      # full reparametrization stride
  param: equi                # parametrization scheme
  max_micro_cycles: 10       # micro-iteration limit
  reset_dlc: true            # rebuild delocalized coordinates each step
  climb: true                # enable climbing image
  climb_rms: 0.0005          # climbing RMS threshold
  climb_lanczos: true        # Lanczos refinement for climbing
  climb_lanczos_rms: 0.0005  # Lanczos RMS threshold
  climb_fixed: false         # keep climbing image fixed
  scheduler: null            # optional scheduler backend
opt:
  type: string               # optimizer type label
  stop_in_when_full: 300     # early stop threshold when string is full
  align: false               # alignment toggle (kept off)
  scale_step: global         # step scaling mode
  max_cycles: 300            # maximum optimization cycles
  dump: false                # dump trajectory/restart data
  dump_restart: false        # dump restart checkpoints
  reparam_thresh: 0.0        # reparametrization threshold
  coord_diff_thresh: 0.0     # coordinate difference threshold
  out_dir: ./result_path_opt/   # output directory
  print_every: 10            # logging stride
dmf:
  correlated: true           # correlated DMF propagation
  sequential: true           # sequential DMF execution
  fbenm_only_endpoints: false   # run FB-ENM beyond endpoints
  fbenm_options:
    delta_scale: 0.2         # FB-ENM displacement scaling
    bond_scale: 1.25         # bond cutoff scaling
    fix_planes: true         # enforce planar constraints
    two_hop_mode: sparse     # neighbor traversal strategy
  cfbenm_options:
    bond_scale: 1.25         # CFB-ENM bond cutoff scaling
    corr0_scale: 1.1         # Correlation scale for corr0
    corr1_scale: 1.5         # Correlation scale for corr1
    corr2_scale: 1.6         # Correlation scale for corr2
    eps: 0.05                # Correlation epsilon
    pivotal: true            # Pivotal residue handling
    single: true             # Single-atom pivots
    remove_fourmembered: true  # Prune four-membered rings
    two_hop_mode: sparse     # Neighbor traversal strategy
  dmf_options:
    remove_rotation_and_translation: false  # Keep rigid-body motions
    mass_weighted: false               # Toggle mass weighting
    parallel: false                    # Enable parallel DMF
    eps_vel: 0.01                      # Velocity tolerance
    eps_rot: 0.01                      # Rotational tolerance
    beta: 10.0                         # Beta parameter for DMF
    update_teval: false                # Update transition evaluation
  k_fix: 300.0                         # Harmonic constant for restraints
```
