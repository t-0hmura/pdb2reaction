# `dft` サブコマンド

## 概要
GPU（利用可能な場合はGPU4PySCF、そうでない場合はCPU PySCF）を使用して一点DFT計算を実行します。総エネルギーに加えて、Mulliken、meta-Löwdin、およびIAO原子電荷/スピン密度を報告します。XYZ/GJF入力では `--ref-pdb` が参照PDBトポロジーを提供し、原子数の検証や（`--ligand-charge` 使用時の）電荷導出に用いられますが、`dft` 自体はPDB/GJF出力を生成しません。

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
1. **入力処理** – `geom_loader` でロード可能な任意のファイル（.pdb/.xyz/.trj/…）を受け入れ、座標は `input_geometry.xyz` として再エクスポートします。XYZ/GJF入力では `--ref-pdb` が参照PDBトポロジーを提供し、原子数検証や（`--ligand-charge` 使用時の）電荷導出に使われます。DFT段階はPDB/GJF出力を行いません。
2. **設定マージ** – デフォルト → CLI → YAML（`dft` ブロック）。YAMLがCLIより優先されます。電荷/多重度は `.gjf` があればそのメタデータを継承し、`-q` が省略され `--ligand-charge` が与えられている場合は酵素–基質複合体として扱って `extract.py` の電荷サマリーから総電荷を導出します。明示的な `-q` は常に優先され、`.gjf` 以外で `-q` が無く `--ligand-charge` も使えない場合は中断します。多重度は省略時 `1` がデフォルトです。
3. **SCFビルド** – `--func-basis` を汎関数と基底に解析し、密度フィッティングはPySCFデフォルトで自動的に有効化します。`--engine` でGPU/CPUの優先度を制御します（`gpu` はGPU4PySCF優先、`cpu` はCPU固定、`auto` はGPU→CPUの順）。非局所補正（例: VV10）はバックエンドの既定以上には明示設定しません。
4. **電子密度解析 & 出力** – 収束後（または失敗後）、エネルギー（Hartree/kcal·mol⁻¹）、収束メタデータ、タイミング、バックエンド情報、および原子ごとのMulliken/meta-Löwdin/IAO電荷とスピン密度を要約する `result.yaml` を書き込みます。解析に失敗した列は `null` に設定され、警告が出力されます。

## CLIオプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | PySCFに提供される総電荷。`.gjf` テンプレートまたは `--ligand-charge`（PDB入力または `--ref-pdb` 付きXYZ/GJF）が提供しない限り必須。両方指定時は `-q` が優先 | テンプレート/導出が適用されない限り必須 |
| `--ligand-charge TEXT` | `-q` が省略された場合に使用される総電荷または残基名ごとのマッピング。PDB入力（または `--ref-pdb` 付きXYZ/GJF）でextract方式の電荷導出を有効化 | _None_ |
| `-m, --multiplicity INT` | スピン多重度（2S+1）。PySCF用に `2S` に変換 | `.gjf` テンプレート値または `1` |
| `--func-basis TEXT` | `FUNC/BASIS` 形式の汎関数/基底ペア | `wb97m-v/def2-tzvpd` |
| `--max-cycle INT` | 最大SCF反復 | `100` |
| `--conv-tol FLOAT` | SCF収束許容値（Hartree） | `1e-9` |
| `--grid-level INT` | PySCF数値積分グリッドレベル | `3` |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_dft/` |
| `--engine [gpu\|cpu\|auto]` | バックエンドポリシー: GPU4PySCF優先、CPUのみ、または自動 | `gpu` |
| `--convert-files {True\|False}` | インターフェース一貫性のために受け入れ; `dft` はPDB/GJF出力を生成しない | `True` |
| `--ref-pdb FILE` | 原子数検証とXYZ/GJF入力のリガンド電荷導出を有効にする参照PDBトポロジー（出力変換は行わない） | _None_ |
| `--args-yaml FILE` | YAMLオーバーライド | _None_ |

## 出力
```
out_dir/ (デフォルト: ./result_dft/)
├─ input_geometry.xyz   # PySCFに送信された構造スナップショット
└─ result.yaml          # 収束/エンジンメタデータを含むエネルギー/電荷/スピンサマリー
```
- `result.yaml` には以下が含まれます:
  - `energy`: Hartree/kcal·mol⁻¹、収束フラグ、実行時間、エンジン情報（`gpu4pyscf`/`pyscf(cpu)`、`used_gpu`）
  - `charges`: Mulliken/meta-Löwdin/IAO原子電荷（失敗時は `null`）
  - `spin_densities`: Mulliken/meta-Löwdin/IAOスピン密度（UKSのみ、失敗時は `null`）
- 電荷・多重度（2S）、汎関数/基底、収束設定、出力ディレクトリも要約されます。

## 注意事項
- GPU4PySCFは利用可能な場合は常に使用されます。`--engine auto` はGPU優先→CPUフォールバックを行い、GPU import/runtime エラー時はCPUで再試行します（`--engine cpu` はCPU固定）。
- **Blackwellアーキテクチャ**GPUが検出された場合、現在のGPU4PySCFがサポートされていない可能性があるため警告が出力されます。
- GPU4PySCFのコンパイル済みホイールはBlackwellをサポートしない場合があり、非x86環境ではソースビルドが必要です。該当環境ではCPUバックエンドまたは自身でのビルドを推奨します（参照: https://github.com/pyscf/gpu4pyscf）。
- 密度フィッティングは常にPySCFデフォルトで試行されます（補助基底の推定は実装されていません）。
- YAML入力はトップレベルに `dft` を持つマッピングである必要があります（非マッピングは `load_yaml_dict` でエラー）。
- IAOの電荷/スピン解析が難しい系では失敗することがあり、`result.yaml` の該当列は `null` となり警告が出力されます。


## YAML設定（`--args-yaml`）
YAMLはトップレベルに `dft` を持つマッピングで指定します（任意で `geom` を追加可能）。YAML値はCLIを上書きします。

`dft` キー（括弧内はデフォルト）:
- `func` (`"wb97m-v"`): 交換相関汎関数
- `basis` (`"def2-tzvpd"`): 基底セット名
- `func_basis` (_None_): `FUNC/BASIS` 形式の統合指定（`func`/`basis` を上書き）
- `conv_tol` (`1e-9`): SCF収束閾値（Hartree）
- `max_cycle` (`100`): 最大SCF反復
- `grid_level` (`3`): PySCF `grids.level`
- `verbose` (`0`): PySCF冗長度（0–9）
- `out_dir` (`"./result_dft/"`): 出力ディレクトリ

_汎関数/基底の既定は `wb97m-v/def2-tzvpd` ですがCLIで上書き可能です。電荷/スピンは `.gjf` テンプレートがあればそれを継承し、`-q` が省略され `--ligand-charge` が与えられている場合は `extract.py` の電荷サマリーで総電荷を導出します。明示的な `-q` は常に優先され、`.gjf` 以外で `--ligand-charge` が使えない場合は中断します。多重度は省略時 `1` がデフォルトです。_

```yaml
geom:
  coord_type: cart       # optional geom_loader settings
dft:
  func: wb97m-v         # exchange–correlation functional
  basis: def2-tzvpd     # basis set name (alternatively use func_basis: "FUNC/BASIS")
  conv_tol: 1.0e-09     # SCF convergence tolerance (Hartree)
  max_cycle: 100        # maximum SCF iterations
  grid_level: 3         # PySCF grid level
  verbose: 0            # PySCF verbosity (0-9)
  out_dir: ./result_dft/  # output directory root
```
