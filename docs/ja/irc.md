# `irc` サブコマンド

## 概要
UMAを使用してEulerPCベースの固有反応座標（IRC）積分を実行します。CLIは意図的に狭く設計されています: 以下に記載されていないものはすべてYAMLで提供する必要があります。XYZ/GJF入力では `--ref-pdb` が参照PDBトポロジーを提供しXYZ座標を保持するため、フォーマット対応のPDB/GJF変換が可能です。

## 使用法
```bash
pdb2reaction irc -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] \
                 [--workers N] [--workers-per-node N] [-m 2S+1]
                 [--max-cycles N] [--step-size Δs] [--root k]
                 [--forward True|False] [--backward True|False]
                 [--freeze-links True|False]
                 [--out-dir DIR]
                 [--convert-files {True\|False}] [--ref-pdb FILE]
                 [--hessian-calc-mode Analytical|FiniteDifference]
                 [--args-yaml FILE]
```

### 例
```bash
# 順方向のみ、有限差分ヘシアン、大きいステップサイズ
pdb2reaction irc -i ts.xyz -q -1 -m 2 --forward True --backward False \
                --step-size 0.2 --hessian-calc-mode FiniteDifference --out-dir ./irc_fd/

# PDB入力: 完成軌跡と方向別軌跡もPDBとしてエクスポート
pdb2reaction irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc/
```

## ワークフロー
1. **入力準備** – `geom_loader` がサポートする任意のフォーマットを受け入れます。`.pdb` 入力ではEulerPC軌跡を元のトポロジーでPDBに自動変換し、`--freeze-links` がリンク水素の親原子を凍結して `geom.freeze_atoms` にマージします。
2. **設定マージ** – デフォルト → CLI → YAML（`geom`、`calc`、`irc`）。電荷/多重度は `.gjf` テンプレートがあれば継承し、`.gjf` 以外では `-q/--charge` が必須（PDB入力または `--ref-pdb` 付きXYZ/GJFに対する `--ligand-charge` がある場合を除く）です。明示的な `-q` は常に優先され、複数度は省略時 `1` がデフォルトです。
3. **IRC積分** – EulerPCが `irc.forward/backward`、`irc.step_length`、`irc.root` に従って順方向/逆方向分岐を積分します。ヘシアンはUMA設定（`calc.*`、`--hessian-calc-mode`）に従い、更新スキーム（既定 `bofill`）や再計算間隔を反映します。
4. **出力** – 軌跡（`finished`、`forward`、`backward`）は `.trj` として書き込まれ、PDB入力の場合は `.pdb` にミラーリングされます。`dump_every` > 0 の場合はHDF5ダンプも作成されます。

## CLIオプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる遷移状態構造 | 必須 |
| `-q, --charge INT` | 総電荷; `calc.charge` をオーバーライド。`.gjf` テンプレートまたは `--ligand-charge`（PDB入力または `--ref-pdb` 付きXYZ/GJF）が提供しない限り必須 | テンプレート/導出が適用されない限り必須 |
| `--ligand-charge TEXT` | `-q` が省略された場合に使用される総電荷または残基名ごとのマッピング。PDB入力（または `--ref-pdb` 付きXYZ/GJF）でextract方式の電荷導出を有効化 | _None_ |
| `--workers`, `--workers-per-node` | UMA予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1）。`calc.spin` を上書き | `.gjf` テンプレート値または `1` |
| `--max-cycles INT` | 最大IRCステップ | `125` |
| `--step-size FLOAT` | 質量重み付き座標でのステップ長 | `0.10` |
| `--root INT` | 初期変位の虚数モードインデックス | `0` |
| `--forward {True\|False}` | 順方向分岐を実行 | `True` |
| `--backward {True\|False}` | 逆方向分岐を実行 | `True` |
| `--freeze-links {True\|False}` | PDB入力用、リンクH親を凍結（`geom.freeze_atoms` にマージ） | `True` |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_irc/` |
| `--convert-files {True\|False}` | PDB入力用のXYZ/TRJ → PDBコンパニオンをトグル | `True` |
| `--ref-pdb FILE` | 入力がXYZ/GJFの場合に使用する参照PDBトポロジー | _None_ |
| `--hessian-calc-mode CHOICE` | UMAヘシアンモード | `FiniteDifference` |
| `--args-yaml FILE` | YAMLオーバーライド | _None_ |

## 出力
```
out_dir/ (デフォルト: ./result_irc/)
├─ irc_data.h5                # dump_every > 0の場合に書き込み
├─ <prefix>finished_irc.trj   # 完全なIRC軌跡
├─ <prefix>forward_irc.trj    # 順方向分岐が実行された場合
├─ <prefix>backward_irc.trj   # 逆方向分岐が実行された場合
└─ *.pdb                      # PDB入力用の軌跡コンパニオン（変換有効時）
```
コンソールには解決済みの `geom`/`calc`/`irc` 設定と経過時間が表示されます。

## 注意事項
- CLIのブール値（`--forward`, `--backward`）は、YAMLに反映させたい場合に `True`/`False` を明示して指定する必要があります。
- UMAはIRC全体で再利用されます。`step_length` を大きくし過ぎるとEulerPCが不安定になる可能性があります。
- VRAMが十分な場合は `--hessian-calc-mode` を `Analytical` に設定することを強く推奨します。
- 電荷/スピンは `.gjf` テンプレートがあればそれを継承します。`.gjf` 以外では、`-q/--charge` が必須（PDB入力または `--ref-pdb` 付きXYZ/GJFに対する `--ligand-charge` がある場合を除く）で、明示的な `-q` が常に優先されます。多重度は省略時に `1` がデフォルトです。
- `--freeze-links` はPDB入力にのみ適用され、リンク水素の親原子を凍結したままヘシアンを構築します。


## YAML設定（`--args-yaml`）
YAMLはマッピングで指定します。YAMLはCLIを上書きします。共通セクションは [YAMLリファレンス](yaml-reference.md) のジオメトリ/計算機キーを再利用します: PDB入力では `--freeze-links` が `geom.freeze_atoms` にマージされ、`--hessian-calc-mode` とCLIの電荷/スピンが `calc` に反映されます。

`irc` キー（括弧内はデフォルト）:
- `step_length` (`0.10`), `max_cycles` (`125`): 主な積分制御（`--step-size`/`--max-cycles`）。
- `downhill` (`False`), `forward` (`True`), `backward` (`True`), `root` (`0`): 分岐選択と初期変位のroot（`--forward`, `--backward`, `--root`）。
- `hessian_init` (`"calc"`), `hessian_update` (`"bofill"`), `hessian_recalc` (`null`): ヘシアン初期化/更新サイクル。
- `displ`, `displ_energy`, `displ_length`: 変位生成の制御。通常はデフォルト推奨。
- 収束閾値: `rms_grad_thresh` (`1.0e-3`), `hard_rms_grad_thresh` (`null`), `energy_thresh` (`1.0e-6`), `imag_below` (`0.0`).
- 出力/診断: `force_inflection` (`True`), `check_bonds` (`False`), `out_dir` (`"./result_irc/"`), `prefix` (`""`), `dump_fn` (`"irc_data.h5"`), `dump_every` (`5`), `max_pred_steps` (`500`), `loose_cycles` (`3`), `corr_func` (`"mbs"`).

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
  max_pred_steps: 500        # predictor-corrector max steps
  loose_cycles: 3            # loose cycles before tightening
  corr_func: mbs             # correlation function choice
```
