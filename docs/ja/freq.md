# `freq`

## 概要

> **要約:** UMA を用いて振動数と熱化学量（ZPE、ギブズエネルギーなど）を計算します。VRAM に余裕がある場合、`--hessian-calc-mode Analytical` によりヘシアン計算を高速化できます。虚振動数は負の値で表示されます。

### 要点
- **想定場面:** 構造が極小点か TS かを検証する場合や、UMA による熱化学補正を求める場合に使用します。
- **凍結原子:** PHVA（部分ヘシアン振動解析）として扱われます。
- **主な出力:** `frequencies_cm-1.txt`、モードアニメーション（`.trj`、条件により `.pdb`）、`thermoanalysis.yaml`（有効化/利用可能な場合）。
- **TS のチェック:** 適切に収束した TS では虚振動数が **1 つだけ**（負の cm⁻¹）であることが期待されます。
- **性能:** VRAM が十分なら `--hessian-calc-mode Analytical` を推奨します。

`pdb2reaction freq` は UMA 計算機で振動解析を行い、凍結原子がある場合は PHVA として活性部分空間で固有解析を行います。基準振動のアニメーションを `.trj` として出力し、PDB テンプレートがあり `--convert-files` が有効な場合は `.pdb` も生成します。`thermoanalysis` パッケージがインストールされていれば、Gaussian 風の熱化学サマリーも出力されます。

設定は **デフォルト → CLI → `--args-yaml`**（`geom`, `calc`, `freq`, `thermo`）の優先順位で解決されます。XYZ/GJF 入力では `--ref-pdb` で参照 PDB トポロジーを指定でき、XYZ 座標を保持したまま PDB 出力変換が可能になります。

## 使用法
```bash
pdb2reaction freq -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m 2S+1] \
                  [--freeze-links {True\|False}] \
                  [--max-write N] [--amplitude-ang Å] [--n-frames N] \
                  [--sort value|abs] [--out-dir DIR] [--args-yaml FILE] \
                  [--temperature K] [--pressure atm] [--dump {True\|False}] \
                  [--hessian-calc-mode Analytical|FiniteDifference] \
                  [--convert-files {True\|False}] [--ref-pdb FILE]
```

### 例
```bash
# 明示的な電荷とスピンでの最小実行
pdb2reaction freq -i a.pdb -q 0 -m 1

# YAML 上書きとカスタム出力ディレクトリを使用したPHVA
pdb2reaction freq -i a.xyz -q -1 --args-yaml ./args.yaml --out-dir ./result_freq/
```

## ワークフロー
- **構造の読み込みと凍結処理**: 構造は `pysisyphus.helpers.geom_loader` で読み込まれます。PDB 入力では `--freeze-links True` によりリンク水素を検出して親原子を凍結し、その結果を `geom.freeze_atoms` にマージします。マージされたインデックスはログに表示され、UMA と PHVA に伝播されます。
- **UMA 計算機**: `--hessian-calc-mode` で解析的または有限差分ヘシアンを選択します。凍結原子がある場合、UMAは活性ブロックのみのヘシアンを返すことがあります。VRAMが十分な場合は `Analytical` を強く推奨します。
- **PHVA と並進・回転射影**: 凍結原子がある場合、固有値解析は活性部分空間内で行われ、並進・回転モードはその空間内で射影されます。3N×3N ヘシアンと活性ブロックヘシアンの両方に対応し、振動数は cm⁻¹ で報告されます（負の値は虚振動数）。
- **モードのエクスポート**: `--max-write` でアニメーション化するモード数を制限できます。`--sort abs` を指定すると絶対値順にソートされます。正弦波アニメーションの振幅（`--amplitude-ang`）とフレーム数（`--n-frames`）は YAML のデフォルトに従います。すべての入力に対して `.trj` が出力され、PDB テンプレートが存在し `--convert-files` が有効な場合のみ `.pdb` も出力されます（ASE 変換がフォールバックとして使用されます）。
- **熱化学**: `thermoanalysis` がインストールされている場合、QRRHO に準じたサマリー（EE、ZPE、E/H/G 補正、熱容量、エントロピー）が PHVA 振動数に基づいて出力されます。CLI の圧力（atm）は内部で Pa に変換されます。`--dump True` を指定すると `thermoanalysis.yaml` も書き込まれます。
- **性能と終了挙動**: GPU メモリ使用量を最小化するため、ヘシアンは 1 つだけ保持し、上三角固有値分解（`UPLO="U"`）を優先します。キーボード割り込みは終了コード 130、その他のエラーはトレースバックを出力して終了コード 1 で終了します。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | 総電荷。省略時は `--ligand-charge` から導出可能。明示的な `-q` は導出値より優先される | `.gjf` テンプレートまたは `--ligand-charge` が提供しない限り必須 |
| `--ligand-charge TEXT` | `-q` が省略された場合に使用される総電荷または残基名ごとのマッピング。PDB 入力（または `--ref-pdb` 付きXYZ/GJF）でextract方式の電荷導出を有効化 | _None_ |
| `--workers`, `--workers-per-node` | UMA予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| `--freeze-links {True\|False}` | PDBのみ。リンク水素の親を凍結し `geom.freeze_atoms` にマージ | `True` |
| `--max-write INT` | エクスポートするモード数 | `10` |
| `--amplitude-ang FLOAT` | モードアニメーション振幅（Å） | `0.8` |
| `--n-frames INT` | モードアニメーションのフレーム数 | `20` |
| `--sort CHOICE` | モード順序: `value`（cm⁻¹）または `abs` | `value` |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_freq/` |
| `--temperature FLOAT` | 熱化学計算の温度（K） | `298.15` |
| `--pressure FLOAT` | 熱化学計算の圧力（atm） | `1.0` |
| `--dump {True\|False}` | `thermoanalysis.yaml` を書き込み | `False` |
| `--hessian-calc-mode CHOICE` | UMAヘシアンモード | `FiniteDifference` |
| `--convert-files {True\|False}` | PDB テンプレートが利用可能な場合に XYZ/TRJ → PDB コンパニオンを出力するかどうか（GJF は出力しない） | `True` |
| `--ref-pdb FILE` | 入力がXYZ/GJFの場合に使用する参照 PDB トポロジー | _None_ |
| `--args-yaml FILE` | YAML 上書き（セクション: `geom`、`calc`、`freq`、`thermo`） | _None_ |

## 出力
```
out_dir/ (デフォルト: ./result_freq/)
├─ mode_XXXX_±freqcm-1.trj  # モードごとのアニメーション
├─ mode_XXXX_±freqcm-1.pdb  # PDB テンプレートが存在し変換が有効な場合のみ
├─ frequencies_cm-1.txt     # 選択したソート順での全振動数リスト
└─ thermoanalysis.yaml      # thermoanalysisがインポート可能で--dumpがTrueの場合
```
コンソールには確定済みの `geom`/`calc`/`freq` 設定と熱化学設定の要約が出力されます。

## 注意事項
- 虚振動モードは負の振動数として報告されます。`freq` は検出された虚振動数の個数を表示し、`--dump True` で詳細を出力します。
- `--hessian-calc-mode` は **デフォルト → CLI → YAML** の優先順位で解決されます。YAML で `calc.hessian_calc_mode` が指定されている場合、CLI の値より優先されます。
- 電荷/スピンは `.gjf` テンプレートがあればそれを継承します。`.gjf` 以外では、`-q/--charge` が必須（PDB 入力または `--ref-pdb` 付きXYZ/GJFに対する `--ligand-charge` がある場合を除く）で、明示的な `-q` が常に優先されます。多重度は省略時に `1` がデフォルトです。意図した状態を確実にするため、明示的に指定してください。


## YAML 設定（`--args-yaml`）
マッピング形式で指定します。YAML の値はデフォルトと CLI の両方を上書きします（最優先）。共通セクションについては [YAML リファレンス](yaml-reference.md) を参照してください。熱化学制御用に `thermo` セクションも利用できます。

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
  return_partial_hessian: true          # allow partial Hessians
freq:
  amplitude_ang: 0.8         # displacement amplitude for modes (Å)
  n_frames: 20               # number of frames per mode
  max_write: 10              # maximum number of modes to write
  sort: value                # sort order: value vs abs
thermo:
  temperature: 298.15        # thermochemistry temperature (K)
  pressure_atm: 1.0          # thermochemistry pressure (atm)
  dump: false                # write thermoanalysis.yaml when true
```

---

## 関連項目

- [tsopt](tsopt.md) — 遷移状態の最適化（妥当な TS であれば虚振動数は 1 つのみ）
- [irc](irc.md) — TS からの IRC（端点での freq と組み合わせることが多い）
- [dft](dft.md) — より高精度なエネルギー評価のための DFT 一点計算
- [all](all.md) — `--thermo True` を含むエンドツーエンドワークフロー
- [YAML リファレンス](yaml-reference.md) — `freq` と `thermo` の設定オプション一覧
- [用語集](glossary.md) — ZPE、ギブズエネルギー、エンタルピー、エントロピーの定義
