# `path-opt`

## 概要

> **要約:** **2 つの構造**（R → P）から、GSM（デフォルト）または DMF（`--mep-mode dmf`）で MEP を最適化します。経路軌跡を書き出し、最高エネルギー画像（HEI）を TS 候補として出力します。

### 要点
- **想定場面:** 反応物と生成物の **2 端点**が揃っていて、まず MEP の初期推定を得たい場合に使います。
- **手法:** 既定は GSM。`--mep-mode dmf` で DMF に切り替え可能。
- **主な出力:** `final_geometries.trj`（経路）と `hei.xyz`（HEI）。変換が有効なら `.pdb`/`.gjf` コンパニオンも生成。
- **既定値:** `--opt-mode light`（LBFGS）、`--climb`、`--max-nodes 10`、`--thresh gau`。
- **次にやること:** HEI は **TS 候補**です。`tsopt` → `freq`（虚数振動数は **1 つ**） → `irc` で検証します。

`pdb2reaction path-opt` は 2 端点間の最小エネルギー経路（MEP）を探索し、最高エネルギー画像（HEI）を報告します。HEI は *候補* に過ぎないため、[freq](freq.md) と [irc](irc.md) によるモード/接続性の確認が必須です。**2 構造以上**を入力して反応領域だけを自動で精密化したい場合は、[path-search](path_search.md) を使用してください。

UMA 計算機で各イメージのエネルギー/勾配/ヘシアンを評価します。最適化の前に剛体アライメントを行い、ストリングが不安定になりにくいようにします。`freeze_atoms` を指定した場合、RMSD フィットにはその原子群のみを使用します（変換自体は全原子へ適用されます）。

設定の優先順位は **defaults < config < 明示指定 CLI < override** です（`geom`, `calc`, `gs`, `opt`, `dmf`, `sopt.*`）。`--config` を先に適用し、`--override-yaml` を最後に適用します。`--args-yaml` は `--override-yaml` の legacy alias です。`--convert-files` が有効（デフォルト）な場合、参照 PDB があるときは `.pdb`、Gaussian テンプレートがあるときは `.gjf` のコンパニオンが生成されます。

## 最小例

```bash
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
  --out-dir ./result_path_opt
```

## 出力の見方

- `result_path_opt/summary.md`
- `result_path_opt/key_mep.trj` / `result_path_opt/key_ts.xyz`
- `result_path_opt/final_geometries.trj`
- `result_path_opt/hei.xyz`
- `result_path_opt/hei.pdb`（PDB 変換が有効な場合）

## よくある例

1. MEP 探索前に端点を事前最適化する。

```bash
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
  --preopt --preopt-max-cycles 20000 --out-dir ./result_path_opt_preopt
```

2. GSM ではなく DMF モードで実行する。

```bash
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
  --mep-mode dmf --max-nodes 12 --out-dir ./result_path_opt_dmf
```

3. リンク親原子を凍結し、climb を切って短時間で確認する。

```bash
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
  --freeze-links --no-climb --out-dir ./result_path_opt_fast
```

## 使用法
```bash
pdb2reaction path-opt -i REACTANT.{pdb|xyz} PRODUCT.{pdb|xyz} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m MULT] \
                      [--workers N] [--workers-per-node N] \
                      [--mep-mode {gsm|dmf}] [--freeze-links/--no-freeze-links] [--max-nodes N] [--max-cycles N] \
                      [--climb/--no-climb] [--dump/--no-dump] [--thresh PRESET] \
                      [--preopt/--no-preopt] [--preopt-max-cycles N] [--opt-mode light|heavy] [--fix-ends/--no-fix-ends] \
                      [--out-dir DIR] [--config FILE] [--override-yaml FILE] [--args-yaml FILE] \
                      [--show-config/--no-show-config] [--dry-run/--no-dry-run] \
                      [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

## ワークフロー
1. **事前アライメント & 凍結解決**
   - 2 番目以降のすべてのエンドポイントは、最初の構造に対して Kabsch アライメントされます。いずれかのエンドポイントで `freeze_atoms` が定義されている場合、RMSD フィットにはその原子のみを使用しますが、算出された変換は全原子に適用されます。
   - `--freeze-links=True`（デフォルト）のPDB 入力では、リンク水素の親原子が検出され `freeze_atoms` にマージされます。

2. **ストリング成長とHEIエクスポート**
   - 経路が成長・精密化された後、内部ノード間の局所極大のうちエネルギーが最も高いものを優先的に探索します。内部の局所極大がない場合は内部ノードの最大値に、内部ノードもない場合は全体の最大値にフォールバックします。
   - 最高エネルギーイメージ（HEI）は `.xyz` と、PDB 参照がある場合は `.pdb` として書き込み、Gaussian テンプレートがある場合は `.gjf` も出力します（いずれも `--convert-files` を尊重）。

### 主要な挙動
- **エンドポイント**: 入力は2構造のみ。形式は `geom_loader` に準拠。PDB 入力（または `--ref-pdb` 付きXYZ/GJF）で軌跡/HEIのPDB 出力が有効。
- **電荷/スピン**: CLI は `.gjf` テンプレートのメタデータより優先されます。`-q` 省略時に `--ligand-charge` がある場合、エンドポイントは酵素–基質複合体として扱われ、PDB 入力（または `--ref-pdb` 付き XYZ/GJF）では `extract.py` の電荷サマリーで総電荷を導出します。明示的な `-q` は常に優先されます。`.gjf` 以外の入力で `-q` を省略すると、導出が成功しない限り中断します。`.gjf` 入力で電荷メタデータが無く `-q` も無い場合は中断します。多重度は省略時 `1` がデフォルトです。正しい状態を得るため、常に明示的に指定してください。
- **MEPセグメント**: `--max-nodes` はGSM/DMFの内部ノード数を制御（GSMの総画像数は `max_nodes + 2`）。`--thresh` またはYAMLで収束プリセット（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`）を指定。
- **クライミングイメージ**: `--climb` は標準のクライミング手順とLanczos接線リファインの両方を切り替え。
- **ダンプ**: `--dump` で StringOptimizer の `opt.dump=True` に対応し、`out_dir` 内に軌跡ダンプを出力します。リスタート YAML は YAML で有効化した場合のみ書き出されます。
- **終了コード**: `0` 成功、`3` 最適化失敗、`4` 軌跡書き込みエラー、`5` HEI 出力エラー、`130` キーボード割り込み、`1` 予期せぬエラー。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH PATH` | 反応物と生成物構造 | 必須 |
| `-q, --charge INT` | 総電荷（`calc.charge`）。`.gjf` 以外では `--ligand-charge` 導出が成功しない限り必須（PDB 入力または `--ref-pdb` 付きXYZ/GJF）。`.gjf` テンプレートがあればそれを使用し、電荷メタデータが無い `.gjf` 入力は `-q` が無いと中断。両方指定時は `-q` が優先 | テンプレート/導出がない限り必須 |
| `--ligand-charge TEXT` | `-q` 省略時に使用する総電荷または残基名ごとのマッピング。PDB 入力（または `--ref-pdb` 付きXYZ/GJF）でextract方式の電荷導出を有効化 | _None_ |
| `--workers`, `--workers-per-node` | UMA予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度 | テンプレート/`1` |
| `--freeze-links/--no-freeze-links` | PDBのみ: リンクH親を凍結（YAMLとマージ） | `True` |
| `--max-nodes INT` | 内部ノード数（ストリングイメージ = `max_nodes + 2`） | `10` |
| `--mep-mode {gsm\|dmf}` | GSM（ストリングベース）またはDMF（ダイレクトフラックス）経路生成器を選択 | `gsm` |
| `--max-cycles INT` | オプティマイザーマクロイテレーション上限 | `300` |
| `--climb/--no-climb` | クライミングイメージ精密化を有効化 | `True` |
| `--dump/--no-dump` | MEP軌跡/リスタートをダンプ | `False` |
| `--opt-mode TEXT` | エンドポイント事前最適化用の単一構造オプティマイザー（`light` = LBFGS、`heavy` = RFO） | `light` |
| `--convert-files/--no-convert-files` | PDB/Gaussian入力用のXYZ/TRJ → PDB/GJFコンパニオンをトグル | `True` |
| `--ref-pdb FILE` | XYZ/GJF 入力用の参照 PDB トポロジー | _None_ |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_path_opt/` |
| `--thresh TEXT` | GSM/ストリングオプティマイザーの収束プリセットを上書き | `gau` |
| `--config FILE` | 明示 CLI 指定より前に適用されるベース YAML | _None_ |
| `--override-yaml FILE` | 最後に適用される上書き YAML（YAML 最優先レイヤ） | _None_ |
| `--args-yaml FILE` | `--override-yaml` の legacy alias | _None_ |
| `--show-config/--no-show-config` | 解決済み設定（YAML レイヤ情報を含む）を表示して実行継続 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画表示のみを行う | `False` |
| `--preopt/--no-preopt` | アライメント/MEP 探索前に各エンドポイントを事前最適化（GSM/DMF） | `False` |
| `--preopt-max-cycles INT` | エンドポイント事前最適化サイクルの上限 | `10000` |
| `--fix-ends/--no-fix-ends` | GSM成長/精密化中にエンドポイント構造を固定 | `False` |

## 出力
```
out_dir/
├─ summary.md                  # 主要成果物へ移動しやすいナビゲーションページ
├─ key_mep.trj                 # 主要 MEP 軌跡へのショートカット（symlink/copy）
├─ key_mep.pdb                 # 主要 MEP PDB へのショートカット（symlink/copy）
├─ key_mep.gjf                 # 主要 MEP GJF へのショートカット（利用可能時）
├─ key_ts.xyz / key_ts.pdb     # TS 候補スナップショットへのショートカット（symlink/copy）
├─ key_ts.gjf                  # TS 候補 GJF へのショートカット（利用可能時）
├─ final_geometries.trj        # XYZ経路（コメント行にエネルギーを保持）
├─ final_geometries.pdb        # PDB 参照が利用可能で変換が有効な場合
├─ hei.xyz                     # 最高エネルギーイメージ
├─ hei.pdb                     # PDB 参照が利用可能な場合のHEI（変換有効時）
├─ hei.gjf                     # Gaussian テンプレートを使用して書き込まれたHEI（変換有効時）
├─ align_refine/               # 剛体アライメント/リファイン段階の中間ファイル（アライメント実行時）
└─ <オプティマイザーダンプ/リスタート>
```
コンソールには解決済みYAMLブロックが出力され、GSM/DMFのMEP進行状況とタイミングが報告されます。


## YAML 設定（`--config`, `--override-yaml`, `--args-yaml`）
マージ順は **defaults < config < 明示指定 CLI < override** です。
`--args-yaml` は `--override-yaml` の legacy alias として維持されています。

### `geom`
- [`opt`](opt.md) と同じキー（`coord_type`, `freeze_atoms` など）。`--freeze-links` がPDB 入力で `freeze_atoms` にマージされます。

### `calc`
- UMA 計算機の設定（`model`, `device`, 近傍半径, ヘシアンなど）。

### `dmf`
- Direct Max Flux + (C)FB-ENM 補間の制御。CLIで露出している `dmf` ブロックと同じキー。

### `gs`
- Growing String表現の制御: `max_nodes`, `perp_thresh`, 再パラメータ化（`reparam_check`, `reparam_every`, `reparam_every_full`, `param`）、`max_micro_cycles`, DLCリセット、climb関連、scheduler。

### `opt`
- StringOptimizer設定: type, `stop_in_when_full`, `scale_step`, `max_cycles`, dump系、`reparam_thresh`, `coord_diff_thresh`, `out_dir`, `print_every`。

### `sopt.lbfgs` / `sopt.rfo`
- エンドポイント事前最適化の単一構造オプティマイザー設定。キーは [YAML リファレンス](yaml_reference.md) の `lbfgs` / `rfo` と同等で、YAML が CLI の `--preopt-max-cycles` を上書きします。

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
  scale_step: global         # step scaling mode
  max_cycles: 300            # maximum optimization cycles
  dump: false                # dump trajectory/restart data
  dump_restart: false        # dump restart checkpoints
  reparam_thresh: 0.0        # reparametrization threshold
  coord_diff_thresh: 0.0     # coordinate difference threshold
  out_dir: ./result_path_opt/   # output directory
  print_every: 10            # logging stride
dmf:
  max_cycles: 300            # DMF/IPOPT の最大反復数
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

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [path-search](path_search.md) — 自動精密化を伴う再帰的MEP 探索（2+構造用）
- [tsopt](tsopt.md) — HEIをTS 候補として最適化（freq/IRCで検証）
- [extract](extract.md) — path-opt入力用のポケットPDBを生成
- [all](all.md) — エンドツーエンドワークフロー（デフォルトでpath-searchを使用）
- [YAML リファレンス](yaml_reference.md) — `gs`、`dmf`、`opt` の完全な設定オプション
- [用語集](glossary.md) — MEP、GSM、DMF、HEIの定義
