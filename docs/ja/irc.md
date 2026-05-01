# `irc`

## 概要

> **要約:** 遷移状態（TS）から反応物・生成物方向へ EulerPC（Euler Predictor-Corrector）ベースの固有反応座標（IRC）積分を実行します。デフォルトで前方・後方の両方向を実行します。VRAM に余裕がある場合は `--hessian-calc-mode Analytical` が推奨されます。

### 要点
- **想定場面:** `tsopt` で最適化・検証済みの TS 構造を出発点に、固有反応座標を追跡して端点接続性（R ↔ TS ↔ P）を確認するケース。
- **手法:** EulerPC（Euler Predictor-Corrector）積分 + MLIP バックエンドのヘシアン（デフォルト UMA、ORB/MACE/AIMNet2 も選択可）。デフォルトでは前方・後方両方の分岐を実行します。
- **主な出力:** `finished_irc_trj.xyz`、`forward_irc_trj.xyz`、`backward_irc_trj.xyz`（参照 PDB が利用可能なら `.pdb` コンパニオンも）。
- **デフォルト:** `--max-cycles 125`、`--step-size 0.10` Bohr、`--root 0`、`--forward`/`--backward` 両方有効、`--hessian-calc-mode FiniteDifference`、バックエンド `uma`。**強制上書き:** IRC は YAML/CLI マージ後に `geom.coord_type = cart` および `calc.return_partial_hessian = true` を強制します。
- **次ステップ:** IRC 端点を [opt](opt.md) で真の極小に最適化、または [freq](freq.md) と組み合わせて熱化学量を取得。ヘシアン評価モードの詳細は {ref}`ja-hessian-evaluation` を参照。

XYZ/GJF 入力では `--ref-pdb` で参照 PDB トポロジーを指定し、XYZ 座標を保持したまま PDB 出力変換が可能になります。一般的な手順は `tsopt` → `irc` です。

## 最小例

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc
```

## 出力の見方

- `result_irc/finished_irc_trj.xyz`
- `result_irc/forward_irc_trj.xyz`
- `result_irc/backward_irc_trj.xyz`

## よくある例

1. 正方向のみを実行する。

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --no-backward \
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
pdb2reaction irc -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [--workers N] [--workers-per-node N] [-m 2S+1] \
 [--max-cycles N] [--step-size Δs] [--root k] \
 [--forward/--no-forward] [--backward/--no-backward] \
 [--freeze-links/--no-freeze-links] \
 [--out-dir DIR] [--config FILE] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE] \
 [--hessian-calc-mode Analytical|FiniteDifference] \
 [--show-config] [--dry-run]
```

### 例
```bash
# 順方向のみ、有限差分ヘシアン、大きいステップサイズ
pdb2reaction irc -i ts.pdb -q 0 -m 1 --no-backward \
 --step-size 0.2 --hessian-calc-mode FiniteDifference --out-dir ./irc_fd/

# PDB 入力: 完成軌跡と方向別軌跡もPDBとしてエクスポート
pdb2reaction irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc/
```

## ワークフロー
1. **入力準備** – `geom_loader` がサポートする任意のフォーマットを受け入れます。参照 PDB が利用可能な場合（PDB 入力時、または `--ref-pdb` で指定した場合）、EulerPC 軌跡はそのトポロジーで PDB に変換されます。PDB 入力に対して `--freeze-links` はリンク水素の親原子を凍結し、`geom.freeze_atoms` を拡張します。
2. **EulerPC 積分** – EulerPC 予測子-修正子積分器が遷移状態から IRC 経路をたどります。`--forward`/`--backward` フラグに従って順方向および/または逆方向の分岐が実行されます。各ステップではエネルギーベースの予測子と修正子ステップを使用します。
3. **軌跡出力** – 完了済み、順方向、逆方向の IRC 軌跡が XYZ ファイルとして書き込まれます。参照 PDB が利用可能な場合、PDB コンパニオンも生成されます（`--convert-files`）。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる遷移状態構造 | 必須 |
| `-q, --charge INT` | 総電荷; YAML が `calc.charge` を指定していない場合に使用。`.gjf` テンプレートまたは `--ligand-charge/-l`（PDB 入力または `--ref-pdb` 付き XYZ/GJF）が提供しない限り必須。両方が設定された場合でも、明示的な `-q` は `--ligand-charge/-l` より優先されます | テンプレート/導出が適用されない限り必須 |
| `-l, --ligand-charge TEXT` | 残基別電荷マッピング（例: `GPP:-3,SAM:1`）。PDB の残基電荷から全系の電荷を自動導出します（手動計算不要）。`-q` 省略時に使用（PDB 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `--workers INT` | MLIP 予測器の並列度（workers > 1 で解析ヘシアン無効）。診断上の注意は {ref}`ja-workers-fd-downgrade` を参照 | `1` |
| `--workers-per-node INT` | ノードあたりのワーカー数。並列予測器に渡されます | `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1）。YAML が `calc.spin` を指定していない場合に使用 | `.gjf` テンプレート値または `1` |
| `--max-cycles INT` | 最大 IRC ステップ（YAML が `irc.max_cycles` を指定していない場合に使用） | `125` |
| `--step-size FLOAT` | ステップ長（Bohr、非質量加重デカルト座標）（YAML が `irc.step_length` を指定していない場合に使用） | `0.10` |
| `--root INT` | 射影ヘシアンの固有値を**昇順**（最も負の値を先頭）に並べたときの**0 始まり**のインデックス。初期 IRC 変位に使用するモードを指定します。虚振動が 1 個だけの妥当な TS では `--root 0`（唯一の負の固有値）のままにしてください。`--root 1`、`--root 2` などは、活性な虚モードがより負のスプリアス（疑似）モードよりも上位にランクされていることが分かっている場合にのみ使用します。YAML が `irc.root` を指定していない場合に使用 | `0` |
| `--forward/--no-forward` | 順方向分岐を実行（YAML が `irc.forward` を指定していない場合に使用） | `True` |
| `--backward/--no-backward` | 逆方向分岐を実行（YAML が `irc.backward` を指定していない場合に使用） | `True` |
| `--freeze-links/--no-freeze-links` | PDB 入力用、リンク H 親を凍結（`geom.freeze_atoms` にマージ）。詳細は [extract](extract.md) を参照 | `True` |
| `--freeze-atoms TEXT` | 凍結する原子の 1 始まりインデックスをカンマ区切りで明示的に指定（例: `'1,3,5'`）。`--freeze-links` と併用可、任意の入力形式に適用 | _None_ |
| `-o, --out-dir TEXT` | 出力ディレクトリ（YAML が `irc.out_dir` を指定していない場合に使用） | `./result_irc/` |
| `--convert-files/--no-convert-files` | 参照 PDB が利用可能な場合に XYZ/TRJ → PDB コンパニオンを出力するかどうか | `True` |
| `--ref-pdb FILE` | 入力が XYZ/GJF の場合に使用する参照 PDB トポロジー | _None_ |
| `--hessian-calc-mode CHOICE` | MLIP ヘシアンモード（YAML が `calc.hessian_calc_mode` を指定していない場合に使用） | `FiniteDifference` |
| `--config FILE` | 明示 CLI 適用前に読み込むベース YAML | _None_ |
| `--show-config/--no-show-config` | 解決済み YAML レイヤー/設定を表示して続行 | `False` |
| `--out-json/--no-out-json` | `out_dir` に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照 | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--solvent TEXT` | xTB 暗黙溶媒（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画のみ表示 | `False` |

## 出力
```
out_dir/ (デフォルト:./result_irc/)
├─ <prefix>finished_irc_trj.xyz # 完全な IRC 軌跡
├─ <prefix>forward_irc_trj.xyz # 順方向分岐が実行された場合
├─ <prefix>backward_irc_trj.xyz # 逆方向分岐が実行された場合
└─ *.pdb # 参照 PDB が利用可能な場合の軌跡コンパニオン（変換有効時）
```
コンソールには確定済みの `geom`/`calc`/`irc` 設定と実行時間の要約が表示されます。

## 終了コード

終了コードは CLI 規約の {ref}`ja-exit-codes` を参照。

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- MLIP バックエンド（デフォルト: UMA）は IRC 全体で再利用されます。`step_length` を大きくし過ぎると EulerPC が不安定になることがあります。
- `--freeze-links` が有効な場合、リンク水素の親原子が自動的に凍結されます（{ref}`リンク水素と凍結原子 <ja-link-hydrogen-and-frozen-atoms>` を参照）。

設定の優先順位は {ref}`CLI 規約: 設定の優先順位 <ja-configuration-precedence>` を参照してください。

`geom`、`calc`、`irc` の各セクションは [YAML リファレンス](yaml-reference.md) の正規定義から変更ありません: [`geom`](yaml-reference.md#geom)、[`calc`](yaml-reference.md#calc)、[`irc`](yaml-reference.md#irc-section) を参照してください。`--freeze-links` は PDB 入力で `geom.freeze_atoms` を拡張し、`--hessian-calc-mode` と CLI の charge/spin 値はマージ済み `calc` ブロックを補完します。

**`irc` 固有の強制上書き**（YAML/CLI マージ後に YAML 値を無視して適用）:

```yaml
geom:
 coord_type: cart # irc では cart に強制（YAML 値は無視）
calc:
 return_partial_hessian: true # irc では true に強制（partial Hessian、active-DOF 処理）
```

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け

- [tsopt](tsopt.md) — IRC 実行前に TS を最適化
- [freq](freq.md) — 完全な振動解析と熱化学補正
- [opt](opt.md) — IRC 端点を真の極小に最適化
- [all](all.md) — tsopt 後に IRC を実行する一気通貫ワークフロー
- [YAML リファレンス](yaml-reference.md) — `irc` の完全な設定オプション
- [用語集](glossary.md) — IRC（固有反応座標）の定義
