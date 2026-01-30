# `dft` サブコマンド

## 概要
GPU（利用可能な場合はGPU4PySCF、そうでない場合はCPU PySCF）を使用して一点DFT計算を実行します。総エネルギーに加えて、Mulliken、meta-Löwdin、およびIAO原子電荷/スピン密度を報告します。

## 使用法
```bash
pdb2reaction dft -i INPUT.{pdb|xyz|gjf|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m MULTIPLICITY] \
                 [--func-basis 'FUNC/BASIS'] \
                 [--max-cycle N] [--conv-tol Eh] [--grid-level L] \
                 [--out-dir DIR] [--engine gpu|cpu|auto] [--convert-files {True\|False}] \
                 [--ref-pdb FILE] [--args-yaml FILE]
```

### 例
```bash
# 明示的な汎関数/基底を使用したデフォルトGPU優先ポリシー
pdb2reaction dft -i input.pdb -q 0 -m 1 --func-basis 'wb97m-v/6-31g**'

# より厳密な制御、大きい基底、CPUのみバックエンド
pdb2reaction dft -i input.pdb -q 1 -m 2 --func-basis 'wb97m-v/def2-tzvpd' --max-cycle 150 --conv-tol 1e-9 --engine cpu
```

## ワークフロー
1. **入力処理** – `geom_loader` でロード可能な任意のファイル（.pdb/.xyz/.trj/…）を受け入れ。座標は `input_geometry.xyz` として再エクスポート
2. **設定マージ** – デフォルト → CLI → YAML（`dft` ブロック）
3. **SCFビルド** – `--func-basis` を汎関数と基底に解析。密度フィッティングはPySCFデフォルトで自動的に有効化
4. **ポピュレーション解析 & 出力** – 収束後（または失敗後）、エネルギー（Hartree/kcal·mol⁻¹）、収束メタデータ、タイミング、バックエンド情報、および原子ごとのMulliken/meta-Löwdin/IAO電荷とスピン密度を要約する `result.yaml` を書き込み

## CLIオプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | PySCFに提供される総電荷 | テンプレート/導出が適用されない限り必須 |
| `--ligand-charge TEXT` | `-q` が省略された場合に使用される総電荷または残基名ごとのマッピング | _None_ |
| `-m, --multiplicity INT` | スピン多重度（2S+1）。PySCF用に `2S` に変換 | `.gjf` テンプレート値または `1` |
| `--func-basis TEXT` | `FUNC/BASIS` 形式の汎関数/基底ペア | `wb97m-v/def2-tzvpd` |
| `--max-cycle INT` | 最大SCF反復 | `100` |
| `--conv-tol FLOAT` | SCF収束許容値（Hartree） | `1e-9` |
| `--grid-level INT` | PySCF数値積分グリッドレベル | `3` |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_dft/` |
| `--engine [gpu\|cpu\|auto]` | バックエンドポリシー: GPU4PySCF優先、CPUのみ、または自動 | `gpu` |
| `--convert-files {True\|False}` | インターフェース一貫性のために受け入れ; `dft` はPDB/GJF出力を生成しない | `True` |
| `--ref-pdb FILE` | 原子数検証とXYZ/GJF入力のリガンド電荷導出を有効にする参照PDBトポロジー | _None_ |
| `--args-yaml FILE` | YAMLオーバーライド | _None_ |

## 出力
```
out_dir/ (デフォルト: ./result_dft/)
├─ input_geometry.xyz   # PySCFに送信された構造スナップショット
└─ result.yaml          # 収束/エンジンメタデータを含むエネルギー/電荷/スピンサマリー
```

## 注意事項
- GPU4PySCFは利用可能な場合は常に使用; そうでない場合はCPU PySCFがビルドされる（`--engine cpu` がCPUを強制しない限り）
- **Blackwellアーキテクチャ**GPUが検出された場合、現在のGPU4PySCFがサポートされていない可能性があるため警告が出力される
- 密度フィッティングは常にPySCFデフォルトで試行される
