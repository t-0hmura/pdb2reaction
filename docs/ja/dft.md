# `dft`

## 概要

> **要約:** GPU4PySCF または CPU PySCF を使用して DFT 一点計算を実行します。デフォルトの汎関数/基底関数は ωB97M-V/def2-TZVPD です。結果にはエネルギーと電子密度解析（Mulliken、meta-Löwdin、IAO 電荷）が含まれます。

`pdb2reaction dft` は PySCF（CPU）または GPU4PySCF（GPU）を使用して DFT 一点計算を実行します。デフォルトの汎関数/基底関数は ωB97M-V/def2-TZVPD です。結果にはエネルギーと電子密度解析（Mulliken、meta-Löwdin、IAO 電荷）が含まれます。

バックエンドは `--engine` で制御します:
- `gpu`（デフォルト）: GPU4PySCF を使用します。**GPU が利用できない場合はエラーになります。**
- `cpu`: CPU PySCF を強制的に使用します。
- `auto`（移植性重視の場合に推奨）: GPU4PySCF を試行し、GPU が利用できない場合は CPU PySCF にフォールバックします。

総エネルギーに加えて、Mulliken、meta-Löwdin、IAO の原子電荷およびスピン密度も報告します。XYZ/GJF 入力では `--ref-pdb` で参照 PDB トポロジーを指定でき、原子数の検証や電荷導出に使われますが、DFT 段階自体は PDB/GJF 出力を生成しません。

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
1. **入力処理** – `geom_loader` でロード可能な任意のファイル（.pdb/.xyz/.trj/…）を受け入れ、座標は `input_geometry.xyz` として再エクスポートされます。XYZ/GJF 入力では `--ref-pdb` が参照 PDB トポロジーを提供し、原子数検証や（`--ligand-charge` 使用時の）電荷導出に使われます。DFT 段階自体は PDB/GJF 出力を生成しません。
2. **設定マージ** – デフォルト → CLI → YAML（`dft` ブロック）。YAMLがCLIより優先されます。電荷/多重度は `.gjf` があればそのメタデータを継承し、`-q` が省略され `--ligand-charge` が与えられている場合は酵素–基質複合体として扱って `extract.py` の電荷サマリーから総電荷を導出します。明示的な `-q` は常に優先され、`.gjf` 以外で `-q` が無く `--ligand-charge` も使えない場合は中断します。多重度は省略時 `1` がデフォルトです。
3. **SCFビルド** – `--func-basis` を汎関数と基底に解析し、密度フィッティングは PySCF のデフォルト設定で自動的に有効化されます。`--engine` でGPU/CPUの優先度を制御します（`gpu` はGPU4PySCF必須、`cpu` はCPU固定、`auto` はGPU→CPUの順）。非局所補正（例: VV10）はバックエンドのデフォルト設定を超える明示的な設定は行いません。
4. **電子密度解析 & 出力** – 収束後（または失敗後）、エネルギー（Hartree/kcal·mol⁻¹）、収束メタデータ、タイミング、バックエンド情報、および原子ごとのMulliken/meta-Löwdin/IAO電荷とスピン密度を要約する `result.yaml` を書き込みます。解析に失敗した列は `null` に設定され、警告が出力されます。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | PySCFに提供される総電荷。`.gjf` テンプレートまたは `--ligand-charge`（PDB 入力または `--ref-pdb` 付きXYZ/GJF）が提供しない限り必須。両方指定時は `-q` が優先 | テンプレート/導出が適用されない限り必須 |
| `--ligand-charge TEXT` | `-q` が省略された場合に使用される総電荷または残基名ごとのマッピング。PDB 入力（または `--ref-pdb` 付きXYZ/GJF）でextract方式の電荷導出を有効化 | _None_ |
| `-m, --multiplicity INT` | スピン多重度（2S+1）。PySCF用に `2S` に変換 | `.gjf` テンプレート値または `1` |
| `--func-basis TEXT` | `FUNC/BASIS` 形式の汎関数/基底ペア | `wb97m-v/def2-tzvpd` |
| `--max-cycle INT` | 最大SCF反復 | `100` |
| `--conv-tol FLOAT` | SCF収束許容値（Hartree） | `1e-9` |
| `--grid-level INT` | PySCF数値積分グリッドレベル | `3` |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_dft/` |
| `--engine [gpu\|cpu\|auto]` | バックエンドポリシー: GPU4PySCF優先、CPUのみ、または自動 | `gpu` |
| `--convert-files {True\|False}` | インターフェースの一貫性のために受け付けるが、`dft` では PDB/GJF 出力は生成されない | `True` |
| `--ref-pdb FILE` | 原子数検証とXYZ/GJF 入力のリガンド電荷導出を有効にする参照 PDB トポロジー（出力変換は行わない） | _None_ |
| `--args-yaml FILE` | YAML 上書き | _None_ |

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
- `--engine gpu`（デフォルト）は GPU4PySCF を必要とし、GPU が利用できない場合は**エラーになります**。GPU リソースが検出されない場合に自動フォールバックさせるには `--engine auto` を使用してください。CPU のみで実行するには `--engine cpu` を指定します。
- **Blackwellアーキテクチャ**GPUが検出された場合、現在のGPU4PySCFがサポートされていない可能性があるため警告が出力されます。
- GPU4PySCFのコンパイル済みホイールはBlackwellをサポートしない場合があり、非x86環境ではソースビルドが必要です。該当環境ではCPUバックエンドまたは自身でのビルドを推奨します（参照: https://github.com/pyscf/gpu4pyscf）。
- 密度フィッティングは常にPySCFデフォルトで試行されます（補助基底の推定は実装されていません）。
- YAML 入力ファイルのルートはマッピングでなければなりません。`dft` セクションは任意です。マッピング以外のルートは `load_yaml_dict` でエラーになります。
- IAO の電荷/スピン解析は難しい系で失敗する場合があり、`result.yaml` の該当列は `null` となり警告が出力されます。


## YAML 設定（`--args-yaml`）
マッピングルートで指定します。`dft` セクション（および任意の `geom`）が存在する場合に適用されます。YAML の値は CLI を上書きします。

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

---

## 関連項目

- [freq](freq.md) — UMAベースの振動解析（DFT精密化の前に行うことが多い）
- [all](all.md) — `--dft True` を使用したエンドツーエンドワークフロー
- [YAML リファレンス](yaml-reference.md) — `dft` の完全な設定オプション
- [用語集](glossary.md) — DFT、SP（一点計算）の定義
