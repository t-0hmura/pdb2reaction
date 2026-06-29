# `dft`

GPU4PySCF または CPU PySCF を使用して DFT 一点計算を実行し、エネルギーと布居解析（population analysis: Mulliken、meta-Löwdin、IAO 電荷）を出力します。デフォルトの汎関数/基底関数は ωB97M-V/def2-tzvpd です。小規模な活性部位モデルの DFT 一点エネルギー（および布居解析）を得たい場面で使用します。多くは、MLIP で最適化した R/TS/P 構造の精密化に用います。バックエンドは `--engine`（デフォルト `gpu`）で選択します。GPU が利用できない場合や移植性・デバッグ目的の実行には `cpu` を使用します。

> `--engine`（単体の `dft`）と `--dft-engine`（`pdb2reaction all` から転送する場合）の命名規則は {ref}`ja-engine-vs-dft-engine` を参照してください。

> **前提条件:** DFT 依存パッケージ（PySCF、GPU4PySCF）はデフォルトではインストールされません。`pip install "pdb2reaction[dft]"` でインストールしてください。

## 実行例

コマンド形式:

```bash
pdb2reaction dft -i INPUT.{pdb|xyz|gjf|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULTIPLICITY] \
 [--func-basis 'FUNC/BASIS'] \
 [--max-cycle N] [--conv-tol Eh] [--grid-level L] \
 [--out-dir DIR] [--engine gpu|cpu] [--convert-files/--no-convert-files] \
 [--ref-pdb FILE] [--config FILE] [--show-config] [--dry-run]
```

基本的な GPU 一点計算。

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 --engine gpu --out-dir ./result_dft
```

大きい基底と厳しい SCF 条件で実行する。

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 \
 --func-basis 'wb97m-v/def2-tzvpd' --conv-tol 1e-10 --max-cycle 200 \
 --engine gpu --out-dir ./result_dft_tight
```

> **注意:** 上記の `def2-tzvpd` 設定は、十分な VRAM を持つ GPU 上で**小さい**活性部位モデル（≲150 原子）にのみ適しています。サイズ／基底／バックエンドごとの閾値の総覧は下記の注意事項を参照してください。

移植性重視で CPU バックエンドを強制する。

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 --engine cpu --out-dir ./result_dft_cpu
```

`-q` を省略し、リガンド定義から総電荷を導出する。

```bash
pdb2reaction dft -i input.pdb -l 'LIG:0' -m 1 \
 --engine gpu --out-dir ./result_dft_ligand
```

`-q` が省略され `--ligand-charge/-l` がある場合、入力は酵素−基質複合体として扱われ、`extract.py` の電荷サマリーから総電荷を計算します。明示的な `-q` は常に最優先です。`.gjf` 以外の入力で `--ligand-charge/-l` もない場合は中断します。

## 処理の流れ

1. **入力処理** – `geom_loader` でロード可能な任意のファイル（.pdb/.xyz/_trj.xyz/…）を受け入れ、座標は `input_geometry.xyz` として再エクスポートされます。XYZ/GJF 入力では `--ref-pdb` が参照 PDB トポロジーを提供し、原子数検証や（`--ligand-charge` 使用時の）電荷導出に使われます。DFT 段階自体は PDB/GJF 出力を生成しません。
2. **SCF ビルド** – `--func-basis` を汎関数と基底に解析します。`--engine` で GPU/CPU を制御します（`gpu` は GPU4PySCF 必須でエラー終了、`cpu` は CPU 固定）。closed-shell + GPU + `--lowmem`（デフォルト）では SCF オブジェクトに `gpu4pyscf.dft.rks_lowmem.RKS` を使用し、メモリ効率の良い直接 JK で密度フィッティングをスキップします。open-shell GPU、CPU、または `--no-lowmem` の経路では密度フィッティングが PySCF のデフォルト設定で自動的に有効化されます。非局所補正（例: VV10）はバックエンドのデフォルトに従い、明示的な上書きは行いません。
3. **布居解析 & 出力** – 収束後（または失敗後）、エネルギー（Hartree/kcal·mol⁻¹）、収束メタデータ、バックエンド情報、および原子ごとの Mulliken/meta-Löwdin/IAO 電荷とスピン密度を要約する `result.yaml` を書き込みます。解析に失敗した項目は `null` に設定され、警告が出力されます。

## 出力

```
out_dir/ (デフォルト:./result_dft/)
├─ input_geometry.xyz # PySCFに送信された構造スナップショット
├─ result.yaml # 収束/エンジンメタデータを含むエネルギー/電荷/スピンサマリー
```

- `result.yaml` には以下が含まれます:
 - `energy`: Hartree/kcal·mol⁻¹、収束フラグ、エンジン情報（`engine`: `gpu4pyscf(rks_lowmem)`/`gpu4pyscf`/`pyscf(cpu)`、`used_gpu`、`used_lowmem`）
 - `charges`: Mulliken/meta-Löwdin/IAO 原子電荷（失敗時は `null`）
 - `spin_densities`: Mulliken/meta-Löwdin/IAO スピン密度（UKS のみ、失敗時は `null`）
- 電荷・多重度（2S）、汎関数/基底、収束設定、出力ディレクトリも要約されます。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル（`.pdb`/`.xyz`/`_trj.xyz`/`.gjf`/…） | 必須 |
| `-q, --charge INT` | PySCF に提供される総電荷。`.gjf` テンプレートまたは `--ligand-charge`（PDB 入力または `--ref-pdb` 付き XYZ/GJF）が提供しない限り必須。両方指定時は `-q` が優先 | テンプレート/導出が適用されない限り必須 |
| `-l, --ligand-charge TEXT` | 単一の整数（例: `-1`）でリガンド総電荷を指定するか、残基別マッピング（例: `GPP:-3,SAM:1`）で PDB 残基電荷から全系の電荷を導出。`-q` 省略時に使用（PDB 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `-m, --multiplicity INT` | スピン多重度（2S+1）。PySCF 用に `2S` に変換 | `.gjf` テンプレート値または `1` |
| `--func-basis TEXT` | `FUNC/BASIS` 形式の汎関数/基底ペア | `wb97m-v/def2-tzvpd` |
| `--max-cycle INT` | 最大 SCF 反復 | `100` |
| `--conv-tol FLOAT` | SCF 収束許容値（Hartree） | `1e-9` |
| `--grid-level INT` | PySCF 数値積分グリッドレベル | `3` |
| `-o, --out-dir TEXT` | 出力ディレクトリ | `./result_dft/` |
| `--engine [gpu\|cpu]` | SCF バックエンド: gpu (GPU4PySCF) または cpu (PySCF)。`--engine` と `--dft-engine` の命名規則は {ref}`ja-engine-vs-dft-engine` を参照 | `gpu` |
| `--lowmem/--no-lowmem` | closed-shell の GPU 経路で `gpu4pyscf.dft.rks_lowmem.RKS` を使用（密度フィッティングを使わず、メモリ効率の良い直接 JK を使用）。open-shell や CPU エンジン、`rks_lowmem` 未搭載の旧 `gpu4pyscf` では標準 RKS/UKS に自動フォールバック | `True` |
| `--convert-files/--no-convert-files` | **`dft` では no-op。** 他のサブコマンドとのインターフェース整合性のためだけに受け付けられます。`dft` は PDB や GJF を一切出力せず（`input_geometry.xyz` + `result.yaml` のみ）、このフラグの値は無視されます | `True` |
| `--ref-pdb FILE` | 原子数検証と XYZ/GJF 入力のリガンド電荷導出を有効にする参照 PDB トポロジー（出力変換は行わない） | _None_ |
| `--config FILE` | 明示的な CLI オプション適用前に読み込むベース YAML | _None_ |
| `--show-config/--no-show-config` | 解決済み設定を表示して実行を継続 | `False` |
| `--out-json/--no-out-json` | `out_dir` に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに設定検証と実行計画表示のみ行う | `False` |

## YAML 設定

マッピングルートで指定します。`dft` セクション（および任意の `geom`）が存在する場合に適用されます。マージ順は次の通りです。

- defaults
- `--config`
- 明示的に指定した CLI オプション

```yaml
geom:
 coord_type: cart # optional geom_loader settings
dft:
 func: wb97m-v # exchange–correlation functional
 basis: def2-tzvpd # basis set name (alternatively use func_basis: "FUNC/BASIS")
 conv_tol: 1.0e-09 # SCF convergence tolerance (Hartree)
 max_cycle: 100 # maximum SCF iterations
 grid_level: 3 # PySCF grid level
 verbose: 0 # PySCF verbose レベル (0-9); CLI -v 2/3 では実行時 PySCF verbose レベル が >=4
 out_dir: ./result_dft/ # output directory root
```

`dft` キー（括弧内はデフォルト）:
- `func` (`"wb97m-v"`): 交換相関汎関数
- `basis` (`"def2-tzvpd"`): 基底セット名
- `func_basis` (_None_): `FUNC/BASIS` 形式の統合指定（`func`/`basis` を上書き）
- `conv_tol` (`1e-9`): SCF 収束閾値（Hartree）
- `max_cycle` (`100`): 最大 SCF 反復
- `grid_level` (`3`): PySCF `grids.level`
- `verbose` (`0`): PySCF verbose レベル（0–9）。デフォルトは quiet。CLI `-v 2/3` では実行時に PySCF verbose レベルが最低 `4` へ上がります。
- `out_dir` (`"./result_dft/"`): 出力ディレクトリ
- `lowmem` (`True`): closed-shell の GPU 経路で `gpu4pyscf.dft.rks_lowmem.RKS` を使用（密度フィッティングをスキップ）。open-shell、CPU エンジン、`rks_lowmem` 非搭載の旧 `gpu4pyscf` では標準 RKS/UKS に自動フォールバック。

## 終了コード

終了コードは CLI 規約の {ref}`ja-exit-codes` を参照。

(ja-notes)=
## 注意事項

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- **デフォルト基底のコスト:** `def2-tzvpd` はトリプルゼータのディフューズ拡張セットであり、大きな系では計算コストが高くなります。探索的な計算には小さい基底（例: `6-31g**` や `def2-svp`）を検討してください。
- **GPU メモリ、def2-TZVPD:** 16–24 GB GPU ではデフォルトの `def2-tzvpd` で **150 原子以上** の系は OOM になります。厳しい設定は十分な VRAM を持つ GPU 上で**小さい**活性部位モデルにのみ適します。代替として `--func-basis 'wb97m-v/def2-svp'` を使用するか（def2-SVP と def2-TZVPD の反応障壁の差は通常 1–3 kcal/mol）、全系を本番計算するなら外部 DFT プログラム（ORCA, Gaussian）を使用してください。
- **Blackwell アーキテクチャ GPU（RTX 50xx）:** GPU4PySCF は小規模な系（~100 原子）でもメモリ不足エラーが発生する場合があります。これらの GPU では `--engine cpu` または外部 DFT プログラム（ORCA, Gaussian）を使用してください。
- **CPU バックエンド:** `--engine cpu` は活性部位モデル（**≲150 原子**）と小さい基底関数（例: `def2-svp`）に限り実用的で、より大きな系を CPU で計算すると非常に低速になるため、全系計算には外部 DFT プログラムの利用を推奨します。
- **系サイズの上限:** DFT 一点計算は **約 300 原子** までの系で実用的です。それ以上の系では計算時間とメモリ使用量が実用範囲を超え、A100 や H200 等の高性能 GPU を搭載した HPC クラスタの利用が必要になる場合があります。酵素系では、DFT 実行前に小さな活性部位モデル（バインディングポケット）を抽出してください。
- **HPC スクラッチ領域:** PySCF / GPU4PySCF は積分や中間ファイルを `$PYSCF_TMPDIR`（未設定なら `$TMPDIR`、最後は `/tmp`）に書き出します。HPC ノードの `/tmp` は容量が小さい、または `tmpfs` (memory-backed) であることが多く、数百原子規模では SCF 途中で枯渇します。`dft` を起動する前に `PYSCF_TMPDIR` をジョブの作業ファイルシステム配下に向けてください（例: `export PYSCF_TMPDIR="$PBS_O_WORKDIR"`）。
- GPU4PySCF のコンパイル済みホイールは非 x86 環境では動作しない場合があります。ソースからビルドしてください（参照: https://github.com/pyscf/gpu4pyscf）。
- 補助基底の推定は未実装です。密度フィッティングの挙動は処理の流れ（SCF ビルド）と `--lowmem` CLI オプションで説明しています。
- YAML 入力ファイルのルートはマッピングでなければなりません。`dft` セクションは任意です。マッピング以外のルートは `load_yaml_dict` でエラーになります。
- IAO の電荷/スピン解析は難しい系で失敗する場合があり、`result.yaml` の該当項目は `null` となり警告が出力されます。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — 一般的な失敗モードの詳細な対処
- [freq](freq.md) — MLIP ベースの振動解析（DFT 精密化の前に行うことが多い）
- [all](all.md) — `--dft` を使用した一気通貫ワークフロー
- [YAML リファレンス](yaml-reference.md) — `dft` の完全な設定オプション
- [用語集](glossary.md) — DFT、SP（一点計算）の定義
