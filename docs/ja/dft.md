# `dft`

GPU4PySCF または CPU PySCF を使用して DFT 一点計算を実行し、エネルギーと布居解析（population analysis: Mulliken、meta-Löwdin、IAO 電荷）を出力します。デフォルトの汎関数/基底関数は ωB97M-V/def2-tzvpd です。小規模な活性部位モデルの DFT 一点エネルギー（および布居解析）を得たい場面で使用します。多くは、MLIP で最適化した R/TS/P 構造上の DFT 一点エネルギー評価に用います。バックエンドは `--engine`（デフォルト `gpu`）で選択します。GPU が利用できない場合や移植性・デバッグ目的の実行には `cpu` を使用します。

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

> **注意:** 上記の `def2-tzvpd` 設定は高costです。普遍的な
> atom-count/VRAM cutoffはないため、代表構造でpilotし、下記の注意事項を
> 参照してください。

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

1. **入力処理** – 共通bridgeがPDB/mmCIFと`geom_loader`対応形式を受け入れ、座標を`input_geometry.xyz`へ再出力します。XYZ/GJF入力では`--ref-pdb`にPDBまたはmmCIF topologyを指定し、原子数検証と電荷導出に使用できます。DFT 段階自体はPDB/CIF/GJF出力を生成しません。
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
- 電荷・多重度・スピン(2S)、汎関数/基底、収束設定、出力ディレクトリも要約されます。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力bridgeが受け入れる構造（`.pdb`/`.cif`/`.mmcif`/`.xyz`/`_trj.xyz`/`.gjf`/…） | 必須 |
| `-q, --charge INT` | PySCF に提供される総電荷。`.gjf` テンプレートまたは `--ligand-charge`（PDB/mmCIF 入力または `--ref-pdb` 付き XYZ/GJF）が提供しない限り必須。両方指定時は `-q` が優先 | テンプレート/導出が適用されない限り必須 |
| `-l, --ligand-charge TEXT` | 単一の整数（例: `-1`）でリガンド総電荷を指定するか、残基別マッピング（例: `GPP:-3,SAM:1`）で PDB/mmCIF 残基電荷から全系の電荷を導出。`-q` 省略時に使用（PDB/mmCIF 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `-m, --multiplicity INT` | スピン多重度（2S+1）。PySCF 用に `2S` に変換 | `.gjf` テンプレート値または `1` |
| `--func-basis TEXT` | `FUNC/BASIS` 形式の汎関数/基底ペア | `wb97m-v/def2-tzvpd` |
| `--max-cycle INT` | 最大 SCF 反復 | `100` |
| `--conv-tol FLOAT` | SCF 収束許容値（Hartree） | `1e-9` |
| `--grid-level INT` | PySCF 数値積分グリッドレベル | `3` |
| `-o, --out-dir TEXT` | 出力ディレクトリ | `./result_dft/` |
| `--engine [gpu\|cpu]` | SCF バックエンド: gpu (GPU4PySCF) または cpu (PySCF)。`--engine` と `--dft-engine` の命名規則は {ref}`ja-engine-vs-dft-engine` を参照 | `gpu` |
| `--lowmem/--no-lowmem` | closed-shell の GPU 経路で `gpu4pyscf.dft.rks_lowmem.RKS` を使用（密度フィッティングを使わず、メモリ効率の良い直接 JK を使用）。open-shell や CPU エンジン、`rks_lowmem` 未搭載の旧 `gpu4pyscf` では標準 RKS/UKS に自動フォールバック | `True` |
| `--convert-files/--no-convert-files` | **`dft` では no-op。** 他のサブコマンドとのインターフェース整合性のためだけに受け付けられます。`dft` は PDB や GJF を一切出力せず（`input_geometry.xyz` + `result.yaml` のみ）、このフラグの値は無視されます | `True` |
| `--ref-pdb FILE` | XYZ/GJF入力の原子数検証とリガンド電荷導出に使う参照PDBまたはmmCIF topology（出力変換なし） | _None_ |
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

全keyとdefaultは [YAMLリファレンス](yaml-reference.md) を参照してください。

## 終了コード

終了コードは CLI 規約の {ref}`ja-exit-codes` を参照。

(ja-notes)=
## 注意事項

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- **system size / basis cost:** `def2-tzvpd` は高コストですが、普遍的な atom-count/VRAM cutoff はありません。basis-function 数、元素、functional、grid、density-fitting path、GPU に依存します。代表構造をpilotし、peak memoryを監視してください。`def2-svp` など小さい基底は安価ですが method 自体が変わるため、基底変更に一律の barrier error を割り当てないでください。
- **新しい GPU architecture:** OOM や unsupported-kernel error は、実メモリ需要だけでなく package/kernel compatibility が原因の場合があります。engine を変更する前に GPU4PySCF/CuPy version と traceback を確認し、全 Blackwell card に同じ既知不具合があると扱わないでください。
- **CPU backend:** `--engine cpu` は対応していますが、実用性は method/system/hardware に依存します。固定の atom-count cutoff ではなく代表 single point を計測してください。
- **HPC scratch:** PySCF / GPU4PySCF は積分や中間fileを `$PYSCF_TMPDIR`（未設定なら `$TMPDIR`、最後は `/tmp`）へ書きます。代表runの実使用量とsite quotaを確認し、必要なら `PYSCF_TMPDIR` をjob filesystem配下へ向けてください（例: `export PYSCF_TMPDIR="$PBS_O_WORKDIR"`）。
- GPU4PySCF のコンパイル済みホイールは非 x86 環境では動作しない場合があります。ソースからビルドしてください（参照: https://github.com/pyscf/gpu4pyscf）。
- 補助基底の推定は未実装です。密度フィッティングの挙動は処理の流れ（SCF ビルド）と `--lowmem` CLI オプションで説明しています。
- YAML 入力ファイルのルートはマッピングでなければなりません。`dft` セクションは任意です。マッピング以外のルートは `load_yaml_dict` でエラーになります。
- IAO の電荷/スピン解析は難しい系で失敗する場合があり、`result.yaml` の該当項目は `null` となり警告が出力されます。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — 一般的な失敗モードの詳細な対処
- [freq](freq.md) — MLIP ベースの振動解析（DFT 一点エネルギー評価の前に行うことが多い）
- [all](all.md) — `--dft` を使用した一気通貫ワークフロー
- [YAML リファレンス](yaml-reference.md) — `dft` の完全な設定オプション
- [用語集](glossary.md) — DFT、SP（一点計算）の定義
