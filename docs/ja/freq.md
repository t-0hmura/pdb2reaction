# `freq`

## 概要

> **要約:** MLIP（デフォルト: UMA、`-b/--backend` で ORB ・ MACE ・ AIMNet2 も選択可能）を用いて振動数と熱化学量（ZPE、ギブズ自由エネルギーなど）を計算します。VRAM に余裕がある場合、`--hessian-calc-mode Analytical` によりヘシアン計算を高速化できます。虚振動数は負の値で表示されます。

### 要点
- **想定場面:** 完全な振動解析（極小点に虚振動数がないこと、TS にちょうど 1 つあること等の確認）または熱化学補正（ZPE、ギブズ自由エネルギーなど）が必要な場合。注: `tsopt` には虚振動数チェックが内蔵されているため、別途 `freq` を実行するのは主に熱化学量の取得や振動モードの詳細検討のためです。
- **手法:** MLIP バックエンドのヘシアン（デフォルト UMA、`-b/--backend` で ORB/MACE/AIMNet2 も選択可）+ 凍結原子に対する PHVA（Partial Hessian Vibrational Analysis）。`thermoanalysis` を介したオプションの QRRHO 風熱化学補正。
- **主な出力:** `frequencies_cm-1.txt`、モードごとの `mode_*_trj.xyz` アニメーション（PDB 入力で変換有効なら `.pdb` も）、`--dump` 指定時で `thermoanalysis` 利用可能なら `thermoanalysis.yaml`。
- **デフォルト:** バックエンド `uma`、`--hessian-calc-mode FiniteDifference`、`--max-write 10`、`--amplitude-ang 0.8`、`--n-frames 20`、`--sort value`、`--temperature 298.15`、`--pressure 1.0`、`--dump False`。
- **次ステップ:** 収束した TS では虚振動数が **ちょうど 1 つ**（負の cm⁻¹）。5 cm⁻¹ 検出閾値と 100 cm⁻¹ 品質ゲートの違いは {ref}`ja-imaginary-mode-thresholds` を参照。ヘシアン評価モードの詳細は {ref}`ja-hessian-evaluation` を参照。

## 最小例

```bash
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 --out-dir ./result_freq
```

## 出力の見方

- `result_freq/frequencies_cm-1.txt`
- `result_freq/mode_*_trj.xyz`
- `result_freq/mode_*.pdb`（PDB 入力かつ変換有効時）

## よくある例

1. まずは出力モード数を絞って確認する。

```bash
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 --max-write 6 --out-dir ./result_freq_quick
```

2. freeze-links + 熱化学ダンプを有効化して実行する。

```bash
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 --freeze-links --dump --out-dir ./result_freq_phva
```

3. VRAM に余裕があるノードで解析的ヘシアンを使う。

```bash
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 \
 --hessian-calc-mode Analytical --out-dir ./result_freq_analytical
```

## 使用法
```bash
pdb2reaction freq -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m 2S+1] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [--workers N] [--workers-per-node N] \
 [--freeze-links/--no-freeze-links] \
 [--max-write N] [--amplitude-ang Å] [--n-frames N] [--sort value|abs] \
 [--out-dir DIR] [--config FILE] [--show-config] [--dry-run] \
 [--temperature K] [--pressure FLOAT] [--dump/--no-dump] \
 [--hessian-calc-mode Analytical|FiniteDifference] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

### 例
```bash
# 明示的な電荷とスピンでの最小実行
pdb2reaction freq -i a.pdb -q 0 -m 1

# YAML 設定ファイルとカスタム出力ディレクトリを使用した PHVA
pdb2reaction freq -i a.xyz -q -1 --config ./freq.yaml --out-dir ./result_freq/
```

## ワークフロー
- **構造の読み込みと凍結処理**: 構造は `pysisyphus.helpers.geom_loader` で読み込まれます。PDB 入力では `--freeze-links` によりリンク水素を検出して親原子を凍結し、その結果を `geom.freeze_atoms` にマージします。マージされたインデックスはログに表示され、MLIP バックエンドと PHVA に伝播されます。
- **MLIP バックエンド**: `--hessian-calc-mode` で解析的または有限差分ヘシアンを選択します。MLIP バックエンドは原子が凍結されている場合、部分（活性）ヘシアンブロックを返すことがあります。ヘシアン評価モードの詳細は {ref}`ja-hessian-evaluation` を参照してください。
- **PHVA と並進・回転射影**: 凍結原子がある場合、固有値解析は活性部分空間内で行われ、並進・回転モードはその空間内で射影されます。3N×3N ヘシアンと活性ブロックヘシアンの両方に対応し、振動数は cm⁻¹ で報告されます（負の値は虚振動数）。
- **モードのエクスポート**: `--max-write` でアニメーション化するモード数を制限できます。モードは値順 (`value`) でソートされ、`--sort abs` を指定すると絶対値順になります。正弦波アニメーションの振幅（`--amplitude-ang`）とフレーム数（`--n-frames`）は YAML のデフォルトに従います。すべての入力に対して `_trj.xyz` が出力され、PDB テンプレートが存在し `--convert-files` が有効な場合のみ `.pdb` も出力されます（ASE 変換がフォールバックとして使用されます）。
- **熱化学**: `thermoanalysis` がインストールされている場合、QRRHO に準じたサマリー（EE、ZPE、E/H/G 補正、熱容量、エントロピー）が PHVA 振動数に基づいて出力されます。CLI の圧力（atm）は内部で Pa に変換されます。`--dump` を指定すると `thermoanalysis.yaml` も書き込まれます。
- **性能と終了挙動**: GPU メモリ使用量を抑えるため、ヘシアンは 1 つだけ保持します。キーボード割り込みは終了コード 130、その他のエラーはトレースバックを出力して終了コード 1 で終了します。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | 総電荷。省略時は `--ligand-charge` から導出可能。明示的な `-q` は導出値より優先される | `.gjf` テンプレートまたは `--ligand-charge` が提供しない限り必須 |
| `-l, --ligand-charge TEXT` | 残基別電荷マッピング（例: `GPP:-3,SAM:1`）。PDB の残基電荷から全系の電荷を自動導出します（手動計算不要）。`-q` 省略時に使用（PDB 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `--workers INT` | MLIP 予測器の並列度（workers > 1 で解析ヘシアン無効）。診断上の注意は {ref}`ja-workers-fd-downgrade` を参照 | `1` |
| `--workers-per-node INT` | ノードあたりのワーカー数。並列予測器に渡されます | `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| `--freeze-links/--no-freeze-links` | PDB のみ。リンク水素の親を凍結し `geom.freeze_atoms` にマージ。リンク水素の詳細は [extract](extract.md) を参照 | `True` |
| `--freeze-atoms TEXT` | 凍結する原子の 1 始まりインデックスをカンマ区切りで明示的に指定（例: `'1,3,5'`）。`--freeze-links` と併用可、任意の入力形式に適用 | _None_ |
| `--max-write INT` | エクスポートするモード数 | `10` |
| `--amplitude-ang FLOAT` | モードアニメーション振幅（Å） | `0.8` |
| `--n-frames INT` | モードアニメーションのフレーム数 | `20` |
| `--sort CHOICE` | モード順序: `value`（cm⁻¹）または `abs` | `value` |
| `-o, --out-dir TEXT` | 出力ディレクトリ | `./result_freq/` |
| `--temperature FLOAT` | 熱化学計算の温度（K） | `298.15` |
| `--pressure FLOAT` | 熱化学計算の圧力（atm）。CLI では `--pressure` ですが、対応する YAML キー（`thermo:` 配下）は `pressure_atm`（単位接尾辞付き）です。いずれも atm で指定し、内部で Pa に変換されます | `1.0` |
| `--dump/--no-dump` | `thermoanalysis.yaml` を書き込み。単体の `freq` ではデフォルト `False` ですが、`pdb2reaction all --thermo` から呼び出された場合、ラッパー側で明示的に `--no-dump` を指定しない限り `True` に切り替わります | `False` |
| `--hessian-calc-mode CHOICE` | MLIP ヘシアンモード（`Analytical` または `FiniteDifference`） | `FiniteDifference` |
| `--convert-files/--no-convert-files` | PDB テンプレートが利用可能な場合に XYZ/TRJ → PDB コンパニオンを出力するかどうか（GJF は出力しない） | `True` |
| `--ref-pdb FILE` | 入力が XYZ/GJF の場合に使用する参照 PDB トポロジー（XYZ 座標は保持） | _None_ |
| `--config FILE` | 明示 CLI 適用前に読み込むベース YAML | _None_ |
| `--show-config/--no-show-config` | 解決済み YAML レイヤー/設定を表示して続行 | `False` |
| `--out-json/--no-out-json` | `out_dir` に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照 | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--solvent TEXT` | xTB 暗黙溶媒（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画のみ表示 | `False` |

## 出力
```
out_dir/ (デフォルト:./result_freq/)
├─ mode_XXXX_±freqcm-1_trj.xyz # モードごとのアニメーション
├─ mode_XXXX_±freqcm-1.pdb # PDB テンプレートが存在し変換が有効な場合のみ
├─ frequencies_cm-1.txt # 選択したソート順での全振動数リスト
└─ thermoanalysis.yaml # thermoanalysisがインポート可能で--dumpがTrueの場合
```
- コンソールには `geom`/`calc`/`freq`/熱化学の設定要約ブロックが出力されます。

## 終了コード

終了コードは CLI 規約の {ref}`ja-exit-codes` を参照。

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- 虚振動数モードは負の振動数として報告されます。`freq` は検出された虚振動数の個数を表示し、`--dump` で詳細を出力します。
- `--hessian-calc-mode` は **デフォルト < config < 明示 CLI** の優先順位で解決されます。CLI で明示的に指定した値は config YAML の `calc.hessian_calc_mode` より優先されます。

`geom`、`calc`、`freq`、`thermo` の各セクションは [YAML リファレンス](yaml-reference.md) の正規定義から変更ありません: [`geom`](yaml-reference.md#geom)、[`calc`](yaml-reference.md#calc)、[`freq`](yaml-reference.md#freq-section)、[`thermo`](yaml-reference.md#thermo) を参照してください。`freq` ではデフォルトで `calc.return_partial_hessian = true`（PHVA）が自動設定されます（YAML で上書き可）。

`freq` 固有のデフォルトとして異なるのは出力ディレクトリのみです:

```yaml
freq:
 out_dir: ./result_freq/ # freq のデフォルト
```

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け

- [tsopt](tsopt.md) — 遷移状態の最適化（内部で虚振動数チェック済み）。続けて IRC で端点の接続性を確認
- [irc](irc.md) — TS からの IRC（端点での freq と組み合わせることが多い）
- [dft](dft.md) — より高精度なエネルギー評価のための DFT 一点計算
- [all](all.md) — `--thermo` を含む一気通貫ワークフロー
- [YAML リファレンス](yaml-reference.md) — `freq` と `thermo` の設定オプション一覧
- [用語集](glossary.md) — ZPE、ギブズ自由エネルギー、エンタルピー、エントロピーの定義
