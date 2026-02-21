# `irc`

## 概要

> **要約:** 遷移状態（TS）から反応物・生成物方向へ固有反応座標（IRC）を追跡します。デフォルトで前方・後方の両方向を実行します。VRAM に余裕がある場合は `--hessian-calc-mode Analytical` が推奨されます。

### 要点
- **入力:** TS 構造（最適化・検証済みが望ましい）。
- **分岐:** 既定で両方向（`--forward`, `--backward`）。
- **主要パラメータ:** `--step-size`（質量重み付き座標でのステップ長）、`--max-cycles`（ステップ数）。
- **強制上書き:** IRC はマージ後に `geom.coord_type = cart` と `calc.return_partial_hessian = false` を強制します（YAML 設定より優先）。
- **主な出力:** `finished_irc.trj` と `forward_irc.trj`/`backward_irc.trj`（参照 PDB があり変換が有効なら `.pdb` も生成）。

`pdb2reaction irc` は UMA を用いて EulerPC ベースの固有反応座標（IRC）積分を実行します。CLI は意図的にシンプルに設計されています。CLI に出ていないパラメータは YAML で明示的に指定することで、再現性のある実行が可能です。

XYZ/GJF 入力では `--ref-pdb` が参照 PDB トポロジーを提供し、XYZ 座標を保持したまま PDB 出力変換が可能になります。一般的な手順は `tsopt` → `freq`（虚数振動数が **1 つ**であることを確認）→ `irc` です。

## 最小例

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc
```

## 出力の見方

- `result_irc/summary.md`
- `result_irc/key_irc.trj`
- `result_irc/key_irc_forward.trj`
- `result_irc/finished_irc.trj`
- `result_irc/forward_irc.trj`
- `result_irc/backward_irc.trj`

## よくある例

1. 正方向のみを実行する。

```bash
pdb2reaction irc -i ts.xyz -q -1 -m 2 --forward --no-backward \
  --out-dir ./result_irc_forward
```

2. ステップを大きくし、解析ヘシアンを使う。

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --step-size 0.20 \
  --hessian-calc-mode Analytical --out-dir ./result_irc_analytical
```

3. 両方向を維持したままステップ数上限を増やす。

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --max-cycles 150 \
  --out-dir ./result_irc_long
```

## 使用法
```bash
pdb2reaction irc -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] \
                 [--workers N] [--workers-per-node N] [-m 2S+1]
                 [--max-cycles N] [--step-size Δs] [--root k]
                 [--forward/--no-forward] [--backward/--no-backward]
                 [--freeze-links/--no-freeze-links]
                 [--out-dir DIR]
                 [--convert-files/--no-convert-files] [--ref-pdb FILE]
                 [--hessian-calc-mode Analytical|FiniteDifference]
                 [--config FILE] [--override-yaml FILE|--args-yaml FILE]
                 [--show-config] [--dry-run]
```

### 例
```bash
# 順方向のみ、有限差分ヘシアン、大きいステップサイズ
pdb2reaction irc -i ts.xyz -q -1 -m 2 --forward --no-backward \
                --step-size 0.2 --hessian-calc-mode FiniteDifference --out-dir ./irc_fd/

# PDB 入力: 完成軌跡と方向別軌跡もPDBとしてエクスポート
pdb2reaction irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc/
```

## ワークフロー
1. **入力準備** – `geom_loader` がサポートする任意のフォーマットを受け入れます。参照 PDB が利用可能な場合（`.pdb` 入力または `--ref-pdb` 指定時）、EulerPC 軌跡はそのトポロジーで PDB に変換されます。`--freeze-links` がリンク水素の親原子を凍結して `geom.freeze_atoms` にマージします。
2. **設定マージ** – デフォルト → `--config` → 明示CLI → `--override-yaml`（`geom`、`calc`、`irc`）。`--args-yaml` は `--override-yaml` の legacy alias です。電荷/多重度は `.gjf` テンプレートがあれば継承し、`.gjf` 以外では `-q/--charge` が必須（PDB 入力または `--ref-pdb` 付きXYZ/GJFに対する `--ligand-charge` がある場合を除く）です。明示的な `-q` は常に優先され、多重度は省略時 `1` がデフォルトです。**IRC は常に `geom.coord_type = cart` と `calc.return_partial_hessian = false` を強制します（YAML の設定に関わらず）。**
3. **IRC積分** – EulerPCが `irc.forward/backward`、`irc.step_length`、`irc.root` に従って順方向/逆方向分岐を積分します。ヘシアンはUMA設定（`calc.*`、`--hessian-calc-mode`）に従い、更新スキーム（既定 `bofill`）や再計算間隔を反映します。
4. **出力** – 軌跡（`finished`、`forward`、`backward`）は `.trj` として書き込まれ、参照 PDB が利用可能な場合は `.pdb` にもミラーリングされます。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる遷移状態構造 | 必須 |
| `-q, --charge INT` | 総電荷; YAML が `calc.charge` を指定していない場合に使用。`.gjf` テンプレートまたは `--ligand-charge`（PDB 入力または `--ref-pdb` 付きXYZ/GJF）が提供しない限り必須。明示的な `-q` は `--ligand-charge` より優先 | テンプレート/導出が適用されない限り必須 |
| `--ligand-charge TEXT` | `-q` が省略された場合に使用される総電荷または残基名ごとのマッピング。PDB 入力（または `--ref-pdb` 付きXYZ/GJF）でextract方式の電荷導出を有効化 | _None_ |
| `--workers`, `--workers-per-node` | UMA予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1）。YAML が `calc.spin` を指定していない場合に使用 | `.gjf` テンプレート値または `1` |
| `--max-cycles INT` | 最大IRCステップ（YAML が `irc.max_cycles` を指定していない場合に使用） | `125` |
| `--step-size FLOAT` | 質量重み付き座標でのステップ長（YAML が `irc.step_length` を指定していない場合に使用） | `0.10` |
| `--root INT` | 初期変位の虚数モードインデックス（YAML が `irc.root` を指定していない場合に使用） | `0` |
| `--forward/--no-forward` | 順方向分岐を実行（YAML が `irc.forward` を指定していない場合に使用） | `True` |
| `--backward/--no-backward` | 逆方向分岐を実行（YAML が `irc.backward` を指定していない場合に使用） | `True` |
| `--freeze-links/--no-freeze-links` | PDB 入力用、リンクH親を凍結（`geom.freeze_atoms` にマージ） | `True` |
| `--out-dir TEXT` | 出力ディレクトリ（YAML が `irc.out_dir` を指定していない場合に使用） | `./result_irc/` |
| `--convert-files/--no-convert-files` | 参照 PDB が利用可能な場合に XYZ/TRJ → PDB コンパニオンを出力するかどうか | `True` |
| `--ref-pdb FILE` | 入力がXYZ/GJFの場合に使用する参照 PDB トポロジー | _None_ |
| `--hessian-calc-mode CHOICE` | UMAヘシアンモード（YAML が `calc.hessian_calc_mode` を指定していない場合に使用） | `FiniteDifference` |
| `--config FILE` | 明示CLI適用前に読み込むベース YAML。 | _None_ |
| `--override-yaml FILE` | 最終 YAML 上書き（YAML レイヤー最優先）。 | _None_ |
| `--args-yaml FILE` | `--override-yaml` の legacy alias。 | _None_ |
| `--show-config/--no-show-config` | 解決済み YAML レイヤー/設定を表示して続行。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画のみ表示。 | `False` |

## 出力
```
out_dir/ (デフォルト: ./result_irc/)
├─ summary.md                # 主要成果物のインデックス
├─ key_irc.trj               # finished_irc.trj へのショートカット
├─ key_irc_forward.trj       # forward_irc.trj へのショートカット
├─ key_irc_backward.trj      # backward_irc.trj へのショートカット
├─ key_irc.pdb               # finished_irc.pdb へのショートカット（存在時）
├─ key_irc_data.h5           # irc_data.h5 へのショートカット（存在時）
├─ <prefix>finished_irc.trj   # 完全な IRC 軌跡
├─ <prefix>forward_irc.trj    # 順方向分岐が実行された場合
├─ <prefix>backward_irc.trj   # 逆方向分岐が実行された場合
└─ *.pdb                      # PDB 入力用の軌跡コンパニオン（変換有効時）
```
コンソールには確定済みの `geom`/`calc`/`irc` 設定と実行時間の要約が表示されます。

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- CLI のブール値（`--forward`, `--backward`）は toggle 形式（`--forward/--no-forward`, `--backward/--no-backward`）で指定します。
- UMAはIRC全体で再利用されます。`step_length` を大きくし過ぎると EulerPC が不安定になることがあります。
- VRAMが十分な場合は `--hessian-calc-mode` を `Analytical` に設定することを強く推奨します。
- `--freeze-links` はPDB 入力にのみ適用され、リンク水素の親原子を凍結したままヘシアンを構築します。


## YAML 設定（`--config` / `--override-yaml` / `--args-yaml`）
マッピング形式で指定し、マージ順は **デフォルト < config < 明示CLI < override** です。`--args-yaml` は `--override-yaml` の legacy alias です。共通セクションについては [YAML リファレンス](yaml-reference.md) を参照してください: PDB 入力では `--freeze-links` が `geom.freeze_atoms` にマージされ、`--hessian-calc-mode` とCLIの電荷/スピンが `calc` に反映されます。`irc` では `geom.coord_type` が `cart` に、`calc.return_partial_hessian` が `false` に強制されます（YAML/CLIより優先）。

`irc` キー（括弧内はデフォルト）:
- `step_length` (`0.10`), `max_cycles` (`125`): 主な積分制御（`--step-size`/`--max-cycles`）。
- `downhill` (`False`), `forward` (`True`), `backward` (`True`), `root` (`0`): 分岐選択と初期変位のroot（`--forward`, `--backward`, `--root`）。
- `hessian_init` (`"calc"`), `hessian_update` (`"bofill"`), `hessian_recalc` (`null`): ヘシアン初期化/更新サイクル。
- `displ`, `displ_energy`, `displ_length`: 変位生成の制御。通常はデフォルト推奨。
- 収束閾値: `rms_grad_thresh` (`1.0e-3`), `hard_rms_grad_thresh` (`null`), `energy_thresh` (`1.0e-6`), `imag_below` (`0.0`).
- 出力/診断: `force_inflection` (`True`), `check_bonds` (`False`), `out_dir` (`"./result_irc/"`), `prefix` (`""`), `max_pred_steps` (`500`), `loose_cycles` (`3`), `corr_func` (`"mbs"`).

```yaml
geom:
  coord_type: cart           # irc では cart に強制（YAML値は無視）
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
  return_partial_hessian: false         # irc では false に強制（完全ヘシアン）
irc:
  step_length: 0.1           # integration step length
  max_cycles: 125            # maximum steps along IRC
  downhill: false            # follow downhill direction only
  forward: true              # propagate in forward direction
  backward: true             # propagate in backward direction
  root: 0                    # normal-mode root index
  hessian_init: calc         # Hessian initialization source
  displ: energy              # displacement construction method
  displ_energy: 0.001        # energy-based displacement scaling
  displ_length: 0.1          # length-based displacement fallback
  rms_grad_thresh: 0.001     # RMS gradient convergence threshold
  hard_rms_grad_thresh: null # hard RMS gradient stop
  energy_thresh: 0.000001    # energy change threshold
  imag_below: 0.0            # imaginary frequency cutoff
  force_inflection: true     # enforce inflection detection
  check_bonds: false         # check bonds during propagation
  out_dir: ./result_irc/     # output directory
  prefix: ""                 # filename prefix
  dump_fn: irc_data.h5       # IRC data filename
  dump_every: 5              # dump stride
  hessian_update: bofill     # Hessian update scheme
  hessian_recalc: null       # Hessian rebuild cadence
  max_pred_steps: 500        # predictor-corrector max steps
  loose_cycles: 3            # loose cycles before tightening
  corr_func: mbs             # correlation function choice
```

---

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け

- [tsopt](tsopt.md) — IRC実行前にTSを最適化
- [freq](freq.md) — TS 候補の虚数振動数を確認し、IRC端点を解析
- [opt](opt.md) — IRC端点を真の極小に最適化
- [all](all.md) — tsopt後にIRCを実行するエンドツーエンドワークフロー
- [YAML リファレンス](yaml-reference.md) — `irc` の完全な設定オプション
- [用語集](glossary.md) — IRC（固有反応座標）の定義
