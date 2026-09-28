# `dft`

GPU4PySCF または CPU PySCF を使用して DFT 一点計算を実行し、エネルギーとポピュレーション解析（population analysis: Mulliken、meta-Löwdin、IAO 電荷）を出力します。デフォルトの汎関数/基底関数は ωB97M-V/def2-svp です。小規模な活性部位モデルの DFT 一点エネルギー（およびポピュレーション解析）を得たい場面で使用します。多くは、MLIP で最適化した R/TS/P 構造上の DFT 一点エネルギー評価に用います。バックエンドは `--dft-engine`（デフォルト `gpu`）で選択します。GPU が利用できない場合や移植性・デバッグ目的の実行には `cpu` を使用します。

> **前提条件:** native CUDA 13 GPU4PySCF は `pdb2reaction[dft]`、CUDA 12 site では `pdb2reaction[dft-cuda12]` をインストールします。

> **溶媒:** `--solvent NAME --solvent-model pcm|smd`はPySCF native implicit solventを
> 使用します。MLIP backendのxTB solvent-delta補正とは別経路です。

## 実行例

コマンド形式:

```bash
pdb2reaction dft -i INPUT.{pdb|xyz|gjf|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULTIPLICITY] \
 [--func-basis 'FUNC/BASIS'] \
 [--scf-max-cycles N] [--scf-tol Eh] [--grid-level L] \
 [--out-dir DIR] [--dft-engine gpu|cpu] \
 [--solvent NAME] [--solvent-model pcm|smd] \
 [--ref-pdb FILE] [--config FILE] [--show-config] [--dry-run]
```

基本的な GPU 一点計算。

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 --dft-engine gpu --out-dir ./result_dft
```

大きい基底と厳しい SCF 条件で実行する。

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 \
 --func-basis 'wb97m-v/def2-tzvpd' --scf-tol 1e-10 --scf-max-cycles 200 \
 --dft-engine gpu --out-dir ./result_dft_tight
```

> **注意:** 上記の `def2-tzvpd` 設定は高costです。普遍的な
> atom-count/VRAM cutoffはないため、代表構造でpilotし、下記の注意事項を
> 参照してください。

移植性重視で CPU バックエンドを強制する。

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 --dft-engine cpu --out-dir ./result_dft_cpu
```

`-q` を省略し、リガンド定義から総電荷を導出する。

```bash
pdb2reaction dft -i input.pdb -l 'LIG:0' -m 1 \
 --dft-engine gpu --out-dir ./result_dft_ligand
```

`-q` が省略され `--ligand-charge/-l` がある場合、入力は酵素−基質複合体として扱われ、`extract.py` の電荷サマリーから総電荷を計算します。明示的な `-q` は常に最優先です。どちらの CLI 電荷指定もない場合は YAML `calc.charge`、GJF ヘッダーの順に参照し、電荷が決まらなければ中断します。

## 処理の流れ

1. **入力処理** – 共通bridgeがPDB/mmCIFと`geom_loader`対応形式を受け入れ、座標を`input_geometry.xyz`へ再出力します。XYZ/GJF入力では`--ref-pdb`にPDBまたはmmCIF topologyを指定し、原子数検証と電荷導出に使用できます。DFT 段階自体はPDB/CIF/GJF出力を生成しません。
2. **SCF ビルド** – `--func-basis` を汎関数と基底に解析します。`--dft-engine` で GPU/CPU を制御します。低メモリモードは既定で有効です。PCM/SMDを含むclosed-shell GPU計算は`gpu4pyscf.dft.rks_lowmem.RKS`、open-shell GPUとCPUはDF tensorを保持しない標準direct-JK RKS/UKSを使います。十分なメモリがある場合、`--no-dft-low-memory`でdensity fittingを有効にすると難しいSCFの収束が改善することがあります。CPU thread数とhost RAM上限はscheduler/process制約から自動検出し、`--dft-nprocs`と`--dft-memory`で上書きできます。これらの資源値は記録されますが、科学的checkpoint identityには入りません。
3. **ポピュレーション解析 & 出力** – 収束後（または失敗後）、エネルギー（Hartree/kcal·mol⁻¹）、収束メタデータ、バックエンド情報、および原子ごとの Mulliken/meta-Löwdin/IAO 電荷とスピン密度を要約する `result.yaml` を書き込みます。解析に失敗した項目は `null` に設定され、警告が出力されます。

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
| `-q, --charge INT` | PySCF に提供される総電荷。優先順位は `-q` → `--ligand-charge/-l` による残基電荷の導出 → YAML `calc.charge` → GJF ヘッダー | 他の指定から電荷が決まらなければ必須 |
| `-l, --ligand-charge TEXT` | 単一の整数（例: `-1`）でリガンド総電荷を指定するか、残基別マッピング（例: `GPP:-3,SAM:1`）で PDB/mmCIF 残基電荷から全系の電荷を導出。`-q` 省略時に使用（PDB/mmCIF 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `-m, --multiplicity INT` | スピン多重度（2S+1）。PySCF 用に `2S` に変換 | YAML `calc.spin` → GJF → `1` |
| `--func-basis TEXT` | `FUNC/BASIS` 形式の汎関数/基底ペア | `wb97m-v/def2-svp` |
| `--scf-max-cycles INT` | 最大 SCF 反復 | `100` |
| `--scf-tol FLOAT` | SCF 収束許容値（Hartree） | `1e-9` |
| `--grid-level INT` | PySCF 数値積分グリッドレベル | `3` |
| `-o, --out-dir TEXT` | 出力ディレクトリ | `./result_dft/` |
| `--dft-engine [gpu\|cpu]` | SCF バックエンド: gpu (GPU4PySCF) または cpu (PySCF)。 | `gpu` |
| `--solvent TEXT` | PySCF native implicit-solvent名。`none`で無効。 | `none` |
| `--solvent-model [pcm\|smd]` | PySCF native implicit-solvent model。 | `smd` |
| `--dft-low-memory/--no-dft-low-memory` | PCM/SMDを含むclosed-shell GPU経路で`gpu4pyscf.dft.rks_lowmem.RKS`を使用。open-shell GPUとCPUは標準direct-JK RKS/UKSを使い、`--no-dft-low-memory`でdensity fittingを有効化 | `True` |
| `--dft-nprocs INT` | PySCF/OpenMP の CPU thread 数。省略時は scheduler/affinity/host から自動検出 | `auto` |
| `--dft-memory SIZE` | PySCF host RAM 上限（例: `64GB`、`120000MB`）。GPU VRAM ではありません | `auto` |
| `--ref-pdb FILE` | XYZ/GJF入力の原子数検証とリガンド電荷導出に使う参照PDBまたはmmCIF topology（出力変換なし） | _None_ |
| `--config FILE` | 明示的な CLI オプション適用前に読み込むベース YAML | _None_ |
| `--show-config/--no-show-config` | 読み込んだ YAML ファイルとその最上位の key を表示して実行を継続 | `False` |
| `--out-json/--no-out-json` | `out_dir` に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照 | `False` |
| `--dry-run/--no-dry-run` | 実行せずにオプションと入力を検証する | `False` |

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
 basis: def2-svp # basis set name (alternatively use func_basis: "FUNC/BASIS")
 lowmem: true # direct-JK低メモリmode。falseでdensity fitting
 nprocs: auto # 必要なら正の整数で明示
 memory: auto # 必要なら64GBなどのhost RAM上限
 conv_tol: 1.0e-09 # SCF convergence tolerance (Hartree)
 max_cycle: 100 # maximum SCF iterations
 grid_level: 3 # PySCF grid level
 pyscf: {mf: {level_shift: 0.2}} # 任意の PySCF object attribute
 verbose: 0 # PySCF verbose レベル (0-9); CLI -v 2/3 では実行時 PySCF verbose レベル が >=4
 out_dir: ./result_dft/ # output directory root
```

`dft` subcommand は `dft.pyscf` を、`-b dft` の calculator workflow は同じ PySCF object 名を `calc.dft.pyscf` から読みます。

全keyとdefaultは [YAMLリファレンス](yaml-reference.md) を参照してください。

## 終了コード

終了コードは CLI 規約の {ref}`ja-exit-codes` を参照。

(ja-notes)=
## 注意事項

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- **system size / basis cost:** `def2-tzvpd` は高コストですが、普遍的な atom-count/VRAM cutoff はありません。basis-function 数、元素、functional、grid、density-fitting path、GPU に依存します。代表構造をpilotし、peak memoryを監視してください。`def2-svp` など小さい基底は安価ですが method 自体が変わるため、基底変更に一律の barrier error を割り当てないでください。
- **新しい GPU architecture:** OOM や unsupported-kernel error は、実メモリ需要だけでなく package/kernel compatibility が原因の場合があります。engine を変更する前に GPU4PySCF/CuPy version と traceback を確認し、全 Blackwell card に同じ既知不具合があると扱わないでください。
- **CPU backend:** `--dft-engine cpu` は対応していますが、実用性は method/system/hardware に依存します。固定の atom-count cutoff ではなく代表 single point を計測してください。
- **HPC scratch:** PySCF / GPU4PySCF は積分や中間fileを `$PYSCF_TMPDIR`（未設定なら `$TMPDIR`、最後は `/tmp`）へ書きます。代表runの実使用量とsite quotaを確認し、必要なら `PYSCF_TMPDIR` をjob filesystem配下へ向けてください（例: `export PYSCF_TMPDIR="$PBS_O_WORKDIR"`）。
- GPU4PySCF のコンパイル済みホイールは非 x86 環境では動作しない場合があります。ソースからビルドしてください（参照: https://github.com/pyscf/gpu4pyscf）。
- 補助基底の推定は未実装です。密度フィッティングの挙動は処理の流れ（SCF ビルド）と `--dft-low-memory` CLI オプションで説明しています。
- YAML 入力ファイルのルートはマッピングでなければなりません。`dft` セクションは任意です。マッピング以外のルートは `load_yaml_dict` でエラーになります。
- IAO の電荷/スピン解析は難しい系で失敗する場合があり、`result.yaml` の該当項目は `null` となり警告が出力されます。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — 一般的な失敗モードの詳細な対処
- [freq](freq.md) — MLIP ベースの振動解析（DFT 一点エネルギー評価の前に行うことが多い）
- [all](all.md) — `--dft` を使用した一気通貫ワークフロー
- [YAML リファレンス](yaml-reference.md) — `dft` の完全な設定オプション
- [用語集](glossary.md) — DFT、SP（一点計算）の定義
