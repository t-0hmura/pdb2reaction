# `freq`

## 概要

> **要約:** UMA を使用して振動数と熱化学（ZPE、ギブズエネルギーなど）を計算します。VRAM に余裕がある場合は `--hessian-calc-mode Analytical` でヘシアン評価が高速化します。虚数振動数は負の値で表示されます。

`pdb2reaction freq` は、凍結原子を部分ヘシアン振動解析（PHVA）で考慮しながら、UMA 計算機を使用して振動解析を実行します。質量重み付き基準モードを `.trj`/`.pdb` アニメーションとしてエクスポートし、オプションの `thermoanalysis` パッケージがインストールされている場合はGaussianスタイルの熱化学サマリーを出力し、`--dump True` の場合はYAML サマリーを出力できます。設定はデフォルト → CLI → YAML（`geom`/`calc`/`freq`）の順に適用され、YAMLが最も優先されます。XYZ/GJF入力では `--ref-pdb` が参照 PDB トポロジーを提供しXYZ座標を保持するため、フォーマット対応のPDB出力変換が可能です。

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
- **構造ロード & 凍結処理**: 構造は `pysisyphus.helpers.geom_loader` を介して読み込みます。PDB 入力では `--freeze-links True` がリンク水素を検出して親原子を凍結し、`geom.freeze_atoms` にマージしてログに表示し、UMAとPHVAに伝播します。
- **UMA 計算機**: `--hessian-calc-mode` で解析的または有限差分ヘシアンを選択します。凍結原子がある場合、UMAは活性ブロックのみのヘシアンを返すことがあります。VRAMが十分な場合は `Analytical` を強く推奨します。
- **PHVA & TR射影**: 凍結原子がある場合、固有解析は活性部分空間内で並進/回転モードを射影して行われます。3N×3Nヘシアンと活性ブロックの両方に対応し、周波数はcm⁻¹で報告されます（負の値は虚数）。
- **モードエクスポート**: `--max-write` でアニメーション化するモード数を制限し、`--sort abs` で絶対値順に並べ替えます。`--amplitude-ang` と `--n-frames` はYAMLの既定値に一致します。すべての入力で `.trj` を出力し、PDB テンプレートが存在し `--convert-files` が有効な場合のみ `.pdb` を出力します（ASE変換がフォールバックとして使用されます）。
- **熱化学**: `thermoanalysis` がインストールされている場合、QRRHOライクなサマリー（EE, ZPE, E/H/G補正、熱容量、エントロピー）を出力します。CLIの圧力（atm）は内部でPaに変換され、`--dump True` で `thermoanalysis.yaml` を書き込みます。
- **性能 & 終了挙動**: 1つのヘシアンのみを保持してGPUメモリを抑え、上三角固有分解（`UPLO="U"`）を優先します。割り込みは終了コード130、その他の失敗は終了コード1で終了します。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | 総電荷。`--ligand-charge` による導出より優先される | `.gjf` テンプレートまたは `--ligand-charge` が提供しない限り必須 |
| `--ligand-charge TEXT` | `-q` が省略された場合に使用される総電荷または残基名ごとのマッピング。PDB 入力（または `--ref-pdb` 付きXYZ/GJF）でextract方式の電荷導出を有効化 | _None_ |
| `--workers`, `--workers-per-node` | UMA予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| `--freeze-links {True\|False}` | PDBのみ。リンク水素の親を凍結し `geom.freeze_atoms` にマージ | `True` |
| `--max-write INT` | エクスポートするモード数 | `10` |
| `--amplitude-ang FLOAT` | モードアニメーション振幅（Å） | `0.8` |
| `--n-frames INT` | モードアニメーションあたりのフレーム数 | `20` |
| `--sort CHOICE` | モード順序: `value`（cm⁻¹）または `abs` | `value` |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_freq/` |
| `--temperature FLOAT` | 熱化学温度（K） | `298.15` |
| `--pressure FLOAT` | 熱化学圧力（atm） | `1.0` |
| `--dump {True\|False}` | `thermoanalysis.yaml` を書き込み | `False` |
| `--hessian-calc-mode CHOICE` | UMAヘシアンモード | `FiniteDifference` |
| `--convert-files {True\|False}` | PDB テンプレートが利用可能な場合のXYZ/TRJ → PDBコンパニオンをトグル（GJFは出力しません） | `True` |
| `--ref-pdb FILE` | 入力がXYZ/GJFの場合に使用する参照 PDB トポロジー | _None_ |
| `--args-yaml FILE` | YAML 上書き（セクション: `geom`、`calc`、`freq`、`thermo`） | _None_ |

## 出力
```
out_dir/ (デフォルト: ./result_freq/)
├─ mode_XXXX_±freqcm-1.trj  # モードごとのアニメーション
├─ mode_XXXX_±freqcm-1.pdb  # PDB テンプレートが存在し変換が有効な場合のみ
├─ frequencies_cm-1.txt     # 選択したソート順での完全な周波数リスト
└─ thermoanalysis.yaml      # thermoanalysisがインポート可能で--dumpがTrueの場合
```
コンソールには解決済みの `geom`/`calc`/`freq` と熱化学設定の要約が出力されます。

## 注意事項
- 虚数モードは負の周波数として報告されます。`freq` は検出数を表示し、`--dump True` で詳細を出力します。
- `--hessian-calc-mode` はYAMLマージ後に `calc.hessian_calc_mode` を上書きします。
- 電荷/スピンは `.gjf` テンプレートがあればそれを継承します。`.gjf` 以外では、`-q/--charge` が必須（PDB 入力または `--ref-pdb` 付きXYZ/GJFに対する `--ligand-charge` がある場合を除く）で、明示的な `-q` が常に優先されます。多重度は省略時に `1` がデフォルトです。


## YAML 設定（`--args-yaml`）
YAMLはマッピングで指定します。YAMLはデフォルトとCLIの両方を上書きします（最優先）。共通セクションは [YAML リファレンス](yaml-reference.md) を再利用します。熱化学用に `thermo` セクションも利用できます。

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

- [tsopt](tsopt.md) — 遷移状態最適化（虚数振動数が1つ必要）
- [irc](irc.md) — TSからのIRC（端点でのfreqと組み合わせることが多い）
- [dft](dft.md) — より高精度なエネルギー精密化のためのDFT一点計算
- [all](all.md) — `--thermo True` を使用したエンドツーエンドワークフロー
- [YAML リファレンス](yaml-reference.md) — `freq` と `thermo` の完全な設定オプション
- [用語集](glossary.md) — ZPE、ギブズエネルギー、エンタルピー、エントロピーの定義