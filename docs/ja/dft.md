# `dft`

## 概要

> **要約:** GPU4PySCF または CPU PySCF を使用して DFT 一点計算を実行します。デフォルトの汎関数/基底関数は ωB97M-V/def2-tzvpd です。結果にはエネルギーと電子密度解析（population analysis: Mulliken、meta-Löwdin、IAO 電荷）が含まれます。

### 要点
- **想定場面:** 小規模な活性部位モデルに対して DFT 一点エネルギー（および電子密度解析）が必要なとき。多くは MLIP で最適化した R/TS/P 構造の精密化に使います。
- **計算手法:** PySCF（CPU）または GPU4PySCF（GPU）。バックエンドは `--engine {gpu|cpu}` で選択します。closed-shell の GPU 経路は GPU4PySCF の低メモリ実装 `rks_lowmem.RKS` をデフォルトで使用し（`--lowmem/--no-lowmem`）、open-shell や CPU、`rks_lowmem` 未対応の旧 `gpu4pyscf` では標準 RKS/UKS に自動フォールバックします。
- **主な出力:** `input_geometry.xyz` と `result.yaml`（Hartree/kcal·mol⁻¹ のエネルギー、収束・実行時間・エンジン情報、Mulliken/meta-Löwdin/IAO 電荷とスピン密度）。
- **デフォルト値:** `--engine gpu`、`--func-basis wb97m-v/def2-tzvpd`、`--max-cycle 100`、`--conv-tol 1e-9`、`--grid-level 3`、`--out-dir ./result_dft/`。
- **次のステップ:** `all --dft --thermo` で DFT エネルギーと MLIP の熱補正を組み合わせる（DFT//MLIP Gibbs）か、ここで得た構造をそのまま機構報告に使います。

`pdb2reaction dft` は PySCF（CPU）または GPU4PySCF（GPU）を使用して DFT 一点計算を実行します。デフォルトの汎関数/基底関数は ωB97M-V/def2-tzvpd です。結果にはエネルギーと電子密度解析（Mulliken、meta-Löwdin、IAO 電荷）が含まれます。

> `--engine`（単体の `dft`）と `--dft-engine`（`pdb2reaction all` から転送する場合）の命名規則は {ref}`ja-engine-vs-dft-engine` を参照してください。

バックエンドは `--engine` で制御します:
- `gpu`（デフォルト）: GPU4PySCF を使用します。**GPU が利用できない場合はエラーになります。** GPU アクセラレーションを保証したい本番計算に最適です。
- `cpu`: CPU PySCF を強制的に使用します。GPU が利用できない場合や、決定的な CPU のみの実行が必要な場合（移植性やデバッグなど）に使用します。

> **前提条件:** DFT 依存パッケージ（PySCF、GPU4PySCF）はデフォルトではインストールされません。`pip install "pdb2reaction[dft]"` でインストールしてください。

総エネルギーに加え、Mulliken、meta-Löwdin、IAO の原子電荷およびスピン密度も報告します。XYZ/GJF 入力では `--ref-pdb` で参照 PDB トポロジーを指定でき、原子数の検証や電荷導出に使用されますが、DFT 段階自体は PDB/GJF 出力を生成しません。

### 実用上限

DFT 一点計算は基底関数のコストとシステムサイズの双方で制約されます。以下の閾値はデフォルト設定（`wb97m-v/def2-tzvpd`、密度フィッティング、グリッドレベル 3）を前提としています。

- **デフォルト基底のコスト:** `def2-tzvpd` はトリプルゼータのディフューズ拡張セットであり、大きな系では計算コストが高くなります。探索的な計算には小さい基底（例: `6-31g**` や `def2-svp`）を検討してください。
- **GPU メモリ、def2-TZVPD:** 16–24 GB GPU ではデフォルトの `def2-tzvpd` で **150 原子以上** の系は OOM になります。代替として `--func-basis 'wb97m-v/def2-svp'` を使用してください。def2-SVP と def2-TZVPD のバリアハイト差は通常 1–3 kcal/mol です。
- **GPU メモリ、小規模活性部位（≲150 原子）:** 厳しい `def2-tzvpd` 設定は十分な VRAM を持つ GPU 上で**小さい**活性部位モデルにのみ適しています。16–24 GB GPU でより大きな系を扱うとこの組み合わせは OOM になります。`def2-svp` に切り替えるか、全系を本番運用するなら外部 DFT プログラム（ORCA, Gaussian）を使用してください。
- **Blackwell アーキテクチャ GPU（RTX 50xx）:** GPU4PySCF は小規模な系（~100 原子）でもメモリ不足エラーが発生する場合があります。これらの GPU では `--engine cpu` または外部 DFT プログラム（ORCA, Gaussian）を使用してください。
- **CPU バックエンド:** `--engine cpu` は活性部位モデル（**≲150 原子**）と小さい基底関数（例: `def2-svp`）に限り実用的で、より大きな系を CPU で計算すると非常に低速になるため、全系計算には外部 DFT プログラムの利用を推奨します。
- **総合的なシステムサイズ上限:** DFT 一点計算は **約 300 原子** までのシステムで実用的です。それ以上のシステムでは計算時間とメモリ使用量が実用範囲を超え、A100 や H200 等の高性能 GPU を搭載した HPC クラスタの利用が必要になる場合があります。酵素系では、DFT 実行前に小さな活性部位モデル（バインディングポケット）を抽出してください。

## 最小例

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 --engine gpu --out-dir ./result_dft
```

## 出力の見方

- `result_dft/input_geometry.xyz`
- `result_dft/result.yaml`
- `result.yaml` 内のエンジン情報（`gpu4pyscf` / `pyscf(cpu)`）

## よくある例

1. 大きい基底と厳しい SCF 条件で実行する。

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 \
 --func-basis 'wb97m-v/def2-tzvpd' --conv-tol 1e-10 --max-cycle 200 \
 --engine gpu --out-dir ./result_dft_tight
```

> **注意:** 上記の `def2-tzvpd` 設定は、十分な VRAM を持つ GPU 上で**小さい**活性部位モデル（≲150 原子）にのみ適しています。サイズ／基底／バックエンドごとの閾値の総覧は上の「実用上限」を参照してください。

2. 移植性重視で CPU バックエンドを強制する。

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 --engine cpu --out-dir ./result_dft_cpu
```

3. `-q` を省略し、リガンド定義から総電荷を導出する。

```bash
pdb2reaction dft -i input.pdb -l 'LIG:0' -m 1 \
 --engine gpu --out-dir ./result_dft_ligand
```

## 使用法
```bash
pdb2reaction dft -i INPUT.{pdb|xyz|gjf|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULTIPLICITY] \
 [--func-basis 'FUNC/BASIS'] \
 [--max-cycle N] [--conv-tol Eh] [--grid-level L] \
 [--out-dir DIR] [--engine gpu|cpu] [--convert-files/--no-convert-files] \
 [--ref-pdb FILE] [--config FILE] [--show-config] [--dry-run]
```

### 例
```bash
# 明示的な汎関数/基底を使用したデフォルトGPU優先ポリシー
pdb2reaction dft -i input.pdb -q 0 -m 1 --func-basis 'wb97m-v/6-31g**'

# より厳密な制御、大きい基底、CPUのみバックエンド
pdb2reaction dft -i input.pdb -q 1 -m 2 \
 --func-basis 'wb97m-v/def2-tzvpd' --max-cycle 150 \
 --conv-tol 1e-9 --engine cpu
```

## ワークフロー
1. **入力処理** – `geom_loader` でロード可能な任意のファイル（.pdb/.xyz/_trj.xyz/…）を受け入れ、座標は `input_geometry.xyz` として再エクスポートされます。XYZ/GJF 入力では `--ref-pdb` が参照 PDB トポロジーを提供し、原子数検証や（`--ligand-charge` 使用時の）電荷導出に使われます。DFT 段階自体は PDB/GJF 出力を生成しません。
2. **SCFビルド** – `--func-basis` を汎関数と基底に解析し、密度フィッティングは PySCF のデフォルト設定で自動的に有効化されます。`--engine` でGPU/CPUを制御します（`gpu` はGPU4PySCF必須でエラー終了、`cpu` はCPU固定）。非局所補正（例: VV10）はバックエンドのデフォルト設定を超える明示的な設定は行いません。
3. **電子密度解析 & 出力** – 収束後（または失敗後）、エネルギー（Hartree/kcal·mol⁻¹）、収束メタデータ、タイミング、バックエンド情報、および原子ごとのMulliken/meta-Löwdin/IAO電荷とスピン密度を要約する `result.yaml` を書き込みます。解析に失敗した列は `null` に設定され、警告が出力されます。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | PySCFに提供される総電荷。`.gjf` テンプレートまたは `--ligand-charge`（PDB 入力または `--ref-pdb` 付きXYZ/GJF）が提供しない限り必須。両方指定時は `-q` が優先 | テンプレート/導出が適用されない限り必須 |
| `-l, --ligand-charge TEXT` | 残基別電荷マッピング（例: `GPP:-3,SAM:1`）。PDB の残基電荷から全系の電荷を自動導出します（手動計算不要）。`-q` 省略時に使用（PDB 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `-m, --multiplicity INT` | スピン多重度（2S+1）。PySCF用に `2S` に変換 | `.gjf` テンプレート値または `1` |
| `--func-basis TEXT` | `FUNC/BASIS` 形式の汎関数/基底ペア | `wb97m-v/def2-tzvpd` |
| `--max-cycle INT` | 最大SCF反復 | `100` |
| `--conv-tol FLOAT` | SCF収束許容値（Hartree） | `1e-9` |
| `--grid-level INT` | PySCF数値積分グリッドレベル | `3` |
| `-o, --out-dir TEXT` | 出力ディレクトリ | `./result_dft/` |
| `--engine [gpu\|cpu]` | SCFバックエンド: gpu (GPU4PySCF) または cpu (PySCF)。`--engine` と `--dft-engine` の命名規則は {ref}`ja-engine-vs-dft-engine` を参照。 | `gpu` |
| `--lowmem/--no-lowmem` | closed-shell の GPU 経路で `gpu4pyscf.dft.rks_lowmem.RKS` を使用（密度フィッティングを使わず、メモリ効率の良い直接 JK を使用）。open-shell や CPU エンジン、`rks_lowmem` 未搭載の旧 `gpu4pyscf` では標準 RKS/UKS に自動フォールバック。 | `True` |
| `--convert-files/--no-convert-files` | **`dft` では no-op。** 他のサブコマンドとのインターフェース整合性のためだけに受け付けられます。`dft` は PDB や GJF を一切出力せず（`input_geometry.xyz` + `result.yaml` のみ）、このフラグの値は無視されます。 | `True` |
| `--ref-pdb FILE` | 原子数検証とXYZ/GJF 入力のリガンド電荷導出を有効にする参照 PDB トポロジー（出力変換は行わない） | _None_ |
| `--config FILE` | 明示的な CLI オプション適用前に読み込むベース YAML | _None_ |
| `--show-config/--no-show-config` | 解決済み設定を表示して実行を継続 | `False` |
| `--out-json/--no-out-json` | `out_dir` に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに設定検証と実行計画表示のみ行う。 | `False` |

## 出力
```
out_dir/ (デフォルト:./result_dft/)
├─ input_geometry.xyz # PySCFに送信された構造スナップショット
├─ result.yaml # 収束/エンジンメタデータを含むエネルギー/電荷/スピンサマリー
```
- `result.yaml` には以下が含まれます:
 - `energy`: Hartree/kcal·mol⁻¹、収束フラグ、実行時間、エンジン情報（`gpu4pyscf`/`pyscf(cpu)`、`used_gpu`）
 - `charges`: Mulliken/meta-Löwdin/IAO原子電荷（失敗時は `null`）
 - `spin_densities`: Mulliken/meta-Löwdin/IAOスピン密度（UKSのみ、失敗時は `null`）
- 電荷・多重度（2S）、汎関数/基底、収束設定、出力ディレクトリも要約されます。

## 終了コード

終了コードは CLI 規約の {ref}`ja-exit-codes` を参照。

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- `--engine gpu`（デフォルト）は GPU4PySCF を必要とし、GPU が利用できない場合は**エラーになります**。CPU のみで実行するには `--engine cpu` を指定します。
- 基底関数のコスト、GPU/CPU メモリ上限、Blackwell GPU、約 300 原子の総合上限については、上の「実用上限」を参照してください。
- GPU4PySCF のコンパイル済みホイールは非 x86 環境では動作しない場合があります。ソースからビルドしてください（参照: https://github.com/pyscf/gpu4pyscf）。
- 密度フィッティングは常に PySCF のデフォルト設定で試行されます（補助基底の推定は未実装）。
- YAML 入力ファイルのルートはマッピングでなければなりません。`dft` セクションは任意です。マッピング以外のルートは `load_yaml_dict` でエラーになります。
- IAO の電荷/スピン解析は難しい系で失敗する場合があり、`result.yaml` の該当列は `null` となり警告が出力されます。


マッピングルートで指定します。`dft` セクション（および任意の `geom`）が存在する場合に適用されます。マージ順は次の通りです。
- defaults
- `--config`
- 明示的に指定した CLI オプション

`dft` キー（括弧内はデフォルト）:
- `func` (`"wb97m-v"`): 交換相関汎関数
- `basis` (`"def2-tzvpd"`): 基底セット名
- `func_basis` (_None_): `FUNC/BASIS` 形式の統合指定（`func`/`basis` を上書き）
- `conv_tol` (`1e-9`): SCF収束閾値（Hartree）
- `max_cycle` (`100`): 最大SCF反復
- `grid_level` (`3`): PySCF `grids.level`
- `verbose` (`0`): PySCF冗長度（0–9）
- `out_dir` (`"./result_dft/"`): 出力ディレクトリ
- `lowmem` (`True`): closed-shell の GPU 経路で `gpu4pyscf.dft.rks_lowmem.RKS` を使用（密度フィッティングをスキップ）。open-shell、CPU エンジン、`rks_lowmem` 非搭載の旧 `gpu4pyscf` では標準 RKS/UKS に自動フォールバック。

_汎関数/基底のデフォルトは `wb97m-v/def2-tzvpd` ですが CLI で上書き可能です。電荷/スピンは `.gjf` テンプレートがあればそれを継承し、`-q` が省略され `--ligand-charge` がある場合は `extract.py` の電荷サマリーから総電荷を導出します。明示的な `-q` は常に最優先です。`.gjf` 以外の入力で `--ligand-charge` もない場合は中断します。多重度は省略時 `1` がデフォルトです。_

```yaml
geom:
 coord_type: cart # optional geom_loader settings
dft:
 func: wb97m-v # exchange–correlation functional
 basis: def2-tzvpd # basis set name (alternatively use func_basis: "FUNC/BASIS")
 conv_tol: 1.0e-09 # SCF convergence tolerance (Hartree)
 max_cycle: 100 # maximum SCF iterations
 grid_level: 3 # PySCF grid level
 verbose: 0 # PySCF verbosity (0-9)
 out_dir: ./result_dft/ # output directory root
```

---

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け

- [freq](freq.md) — MLIPベースの振動解析（DFT精密化の前に行うことが多い）
- [all](all.md) — `--dft` を使用した一気通貫ワークフロー
- [YAML リファレンス](yaml-reference.md) — `dft` の完全な設定オプション
- [用語集](glossary.md) — DFT、SP（一点計算）の定義
