# `freq`

## 概要

> **要約:** UMA を用いて振動数と熱化学量（ZPE、ギブズエネルギーなど）を計算します。VRAM に余裕がある場合、`--hessian-calc-mode Analytical` によりヘシアン計算を高速化できます。虚振動数は負の値で表示されます。

### 要点
- **想定場面:** 構造が極小点か TS かを検証する場合や、UMA による熱化学補正を求める場合に使用します。
- **凍結原子:** PHVA（部分ヘシアン振動解析）として扱われます。
- **主な出力:** `frequencies_cm-1.txt`、モードアニメーション（`_trj.xyz`、条件により `.pdb`）、`thermoanalysis.yaml`（有効化/利用可能な場合）。
- **TS のチェック:** 適切に収束した TS では虚振動数が **1 つだけ**（負の cm⁻¹）であることが期待されます。
- **性能:** VRAM が十分なら `--hessian-calc-mode Analytical` を推奨します。

`pdb2reaction freq` は UMA 計算機で振動解析を行い、凍結原子がある場合は PHVA として活性部分空間で固有解析を行います。基準振動のアニメーションを `_trj.xyz` として出力し、PDB テンプレートがあり `--convert-files` が有効な場合は `.pdb` も生成します。`thermoanalysis` パッケージがインストールされていれば、Gaussian 風の熱化学サマリーも出力されます。


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
pdb2reaction freq -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m 2S+1] \
 [--freeze-links/--no-freeze-links] \
 [--max-write N] [--amplitude-ang Å] [--n-frames N] \
 [--show-config] [--dry-run] \
 [--temperature K] [--pressure atm] [--dump/--no-dump] \
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
- **構造の読み込みと凍結処理**: 構造は `pysisyphus.helpers.geom_loader` で読み込まれます。PDB 入力では `--freeze-links` によりリンク水素を検出して親原子を凍結し、その結果を `geom.freeze_atoms` にマージします。マージされたインデックスはログに表示され、UMA と PHVA に伝播されます。
- **UMA 計算機**: `--hessian-calc-mode` で解析的または有限差分ヘシアンを選択します。凍結原子がある場合、UMAは活性ブロックのみのヘシアンを返すことがあります。VRAMが十分な場合は `Analytical` を強く推奨します。
- **PHVA と並進・回転射影**: 凍結原子がある場合、固有値解析は活性部分空間内で行われ、並進・回転モードはその空間内で射影されます。3N×3N ヘシアンと活性ブロックヘシアンの両方に対応し、振動数は cm⁻¹ で報告されます（負の値は虚振動数）。
- **モードのエクスポート**: `--max-write` でアニメーション化するモード数を制限できます。`--sort abs` を指定すると絶対値順にソートされます。正弦波アニメーションの振幅（`--amplitude-ang`）とフレーム数（`--n-frames`）は YAML のデフォルトに従います。すべての入力に対して `_trj.xyz` が出力され、PDB テンプレートが存在し `--convert-files` が有効な場合のみ `.pdb` も出力されます（ASE 変換がフォールバックとして使用されます）。
- **熱化学**: `thermoanalysis` がインストールされている場合、QRRHO に準じたサマリー（EE、ZPE、E/H/G 補正、熱容量、エントロピー）が PHVA 振動数に基づいて出力されます。CLI の圧力（atm）は内部で Pa に変換されます。`--dump` を指定すると `thermoanalysis.yaml` も書き込まれます。
- **性能と終了挙動**: GPU メモリ使用量を最小化するため、ヘシアンは 1 つだけ保持し、上三角固有値分解（`UPLO="U"`）を優先します。キーボード割り込みは終了コード 130、その他のエラーはトレースバックを出力して終了コード 1 で終了します。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | 総電荷。省略時は `--ligand-charge` から導出可能。明示的な `-q` は導出値より優先される | `.gjf` テンプレートまたは `--ligand-charge` が提供しない限り必須 |
| `--ligand-charge TEXT` | 残基別電荷マッピング（例: `GPP:-3,SAM:1`）。PDB の残基電荷から全系の電荷を自動導出します（手動計算不要）。`-q` 省略時に使用（PDB 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `--workers`, `--workers-per-node` | UMA予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| `--freeze-links/--no-freeze-links` | PDBのみ。リンク水素の親を凍結し `geom.freeze_atoms` にマージ | `True` |
| `--max-write INT` | エクスポートするモード数 | `10` |
| `--amplitude-ang FLOAT` | モードアニメーション振幅（Å） | `0.8` |
| `--n-frames INT` | モードアニメーションのフレーム数 | `20` |
| `--sort CHOICE` | モード順序: `value`（cm⁻¹）または `abs` | `value` |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_freq/` |
| `--temperature FLOAT` | 熱化学計算の温度（K） | `298.15` |
| `--pressure FLOAT` | 熱化学計算の圧力（atm） | `1.0` |
| `--dump/--no-dump` | `thermoanalysis.yaml` を書き込み | `False` |
| `--hessian-calc-mode CHOICE` | UMAヘシアンモード | `FiniteDifference` |
| `--convert-files/--no-convert-files` | PDB テンプレートが利用可能な場合に XYZ/TRJ → PDB コンパニオンを出力するかどうか（GJF は出力しない） | `True` |
| `--ref-pdb FILE` | 入力がXYZ/GJFの場合に使用する参照 PDB トポロジー | _None_ |
| `--config FILE` | 明示CLI適用前に読み込むベース YAML。 | _None_ |
| `--show-config/--no-show-config` | 解決済み YAML レイヤー/設定を表示して続行。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画のみ表示。 | `False` |

## 出力
```
out_dir/ (デフォルト:./result_freq/)
├─ mode_XXXX_±freqcm-1_trj.xyz # モードごとのアニメーション
├─ mode_XXXX_±freqcm-1.pdb # PDB テンプレートが存在し変換が有効な場合のみ
├─ frequencies_cm-1.txt # 選択したソート順での全振動数リスト
└─ thermoanalysis.yaml # thermoanalysisがインポート可能で--dumpがTrueの場合
```
コンソールには確定済みの `geom`/`calc`/`freq` 設定と熱化学設定の要約が出力されます。

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- 虚振動モードは負の振動数として報告されます。`freq` は検出された虚振動数の個数を表示し、`--dump` で詳細を出力します。
- `--hessian-calc-mode` は **デフォルト < config < 明示CLI < override** の優先順位で解決されます。YAML で `calc.hessian_calc_mode` が指定されている場合、最終 override レイヤーが優先されます。


マッピング形式で指定し、マージ順は **デフォルト < config < 明示CLI < override** です。共通セクションについては [YAML リファレンス](yaml_reference.md) を参照してください。熱化学制御用に `thermo` セクションも利用できます。

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
 return_partial_hessian: true # allow partial Hessians
freq:
 amplitude_ang: 0.8 # displacement amplitude for modes (Å)
 n_frames: 20 # number of frames per mode
 max_write: 10 # maximum number of modes to write
 sort: value # sort order: value vs abs
thermo:
 temperature: 298.15 # thermochemistry temperature (K)
 pressure_atm: 1.0 # thermochemistry pressure (atm)
 dump: false # write thermoanalysis.yaml when true
```

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け

- [tsopt](tsopt.md) — 遷移状態の最適化（妥当な TS であれば虚振動数は 1 つのみ）
- [irc](irc.md) — TS からの IRC（端点での freq と組み合わせることが多い）
- [dft](dft.md) — より高精度なエネルギー評価のための DFT 一点計算
- [all](all.md) — `--thermo` を含むend-to-endワークフロー
- [YAML リファレンス](yaml_reference.md) — `freq` と `thermo` の設定オプション一覧
- [用語集](glossary.md) — ZPE、ギブズエネルギー、エンタルピー、エントロピーの定義
