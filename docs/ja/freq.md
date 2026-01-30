# `freq` サブコマンド

## 概要
`pdb2reaction freq` は、凍結原子を部分ヘシアン振動解析（PHVA）で考慮しながら、UMA計算機を使用して振動解析を実行します。質量重み付き基準モードを `.trj`/`.pdb` アニメーションとしてエクスポートし、オプションの `thermoanalysis` パッケージがインストールされている場合はGaussianスタイルの熱化学サマリーを出力し、`--dump True` の場合はYAMLサマリーを出力できます。

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

# YAMLオーバーライドとカスタム出力ディレクトリを使用したPHVA
pdb2reaction freq -i a.xyz -q -1 --args-yaml ./args.yaml --out-dir ./result_freq/
```

## ワークフロー
- **構造ロード & 凍結処理**: 構造は `pysisyphus.helpers.geom_loader` を介して読み込み。PDB入力では、`--freeze-links True` がリンク水素を検出してその親原子を凍結
- **UMA計算機**: `--hessian-calc-mode` で解析的または有限差分ヘシアンを選択。VRAMが十分な場合は `Analytical` を強く推奨
- **PHVA & TR射影**: 凍結原子がある場合、固有解析は活性部分空間内で並進/回転モードが射影された状態で行われる
- **モードエクスポート**: `--max-write` でアニメーション化するモード数を制限
- **熱化学**: `thermoanalysis` がインストールされている場合、QRRHOライクなサマリーを出力

## CLIオプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | 総電荷 | テンプレート/`--ligand-charge` が提供しない限り必須 |
| `--ligand-charge TEXT` | `-q` が省略された場合に使用される総電荷または残基名ごとのマッピング | _None_ |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| `--freeze-links {True\|False}` | PDBのみ。リンク水素の親を凍結 | `True` |
| `--max-write INT` | エクスポートするモード数 | `10` |
| `--amplitude-ang FLOAT` | モードアニメーション振幅（Å） | `0.8` |
| `--n-frames INT` | モードアニメーションあたりのフレーム数 | `20` |
| `--sort CHOICE` | モード順序: `value`（cm⁻¹）または `abs` | `value` |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_freq/` |
| `--temperature FLOAT` | 熱化学温度（K） | `298.15` |
| `--pressure FLOAT` | 熱化学圧力（atm） | `1.0` |
| `--dump {True\|False}` | `thermoanalysis.yaml` を書き込み | `False` |
| `--hessian-calc-mode CHOICE` | UMAヘシアンモード | `FiniteDifference` |
| `--convert-files {True\|False}` | PDBテンプレートが利用可能な場合のXYZ/TRJ → PDBコンパニオンをトグル | `True` |
| `--ref-pdb FILE` | 入力がXYZ/GJFの場合に使用する参照PDBトポロジー | _None_ |
| `--args-yaml FILE` | YAMLオーバーライド（セクション: `geom`、`calc`、`freq`、`thermo`） | _None_ |

## 出力
```
out_dir/ (デフォルト: ./result_freq/)
├─ mode_XXXX_±freqcm-1.trj  # モードごとのアニメーション
├─ mode_XXXX_±freqcm-1.pdb  # PDBテンプレートが存在し変換が有効な場合のみ
├─ frequencies_cm-1.txt     # 選択したソート順での完全な周波数リスト
└─ thermoanalysis.yaml      # thermoanalysisがインポート可能で--dumpがTrueの場合
```

## 注意事項
- 虚数モードは負の周波数として報告される
- `--hessian-calc-mode` はYAMLマージ後に `calc.hessian_calc_mode` をオーバーライド
