# MLIP バックエンド

計算バックエンドの選び方と、バックエンドごとのインストール、モデル名、精度、再現性の設定、Hessian の計算方式をまとめたページです。デフォルトのバックエンドは **UMA**（Meta の Universal Models for Atoms）で、`-b/--backend` で **ORB**、**MACE**、**AIMNet2** も選べます。4 つとも機械学習原子間ポテンシャル（MLIP）です。

## バックエンドごとの特性

バックエンドは、計算を行うどのコマンドでも `-b/--backend` で選びます。

```bash
# UMA（デフォルト）
pdb2reaction opt -i input.pdb -q 0

# ORB
pdb2reaction opt -i input.pdb -q 0 -b orb

# MACE
pdb2reaction opt -i input.pdb -q 0 -b mace

# AIMNet2
pdb2reaction opt -i input.pdb -q 0 -b aimnet2
```

| バックエンド | インストール | モデル名 | `--precision` | 解析 Hessian | 複数ワーカー |
|---------|---------|------------------|------------------|--------------------|------------------|
| `uma` | 同梱（`fairchem-core` ≥ 2.22 は本体の依存）＋ [Hugging Face へのログイン](installation.md#必須) | `uma-s-1p2`（デフォルト）/ `uma-m-1p1` | `fp32` / `fp64` | あり（autograd） | あり |
| `orb` | `pip install "pdb2reaction[orb]"` | `orb_v3_conservative_omol`（エネルギーの勾配から力を求める conservative モデルのみ） | `fp32` / `fp64` | あり（autograd） | なし |
| `mace` | 専用の conda 環境に pdb2reaction を入れてから `pip uninstall -y fairchem-core && pip install 'mace-torch>=0.3.8'`（この環境では UMA は動きません） | `MACE-OMOL-0` | `fp32` / `fp64` | あり | なし |
| `aimnet2` | `pip install "pdb2reaction[aimnet]"` | `aimnet2` | `fp32` のみ | あり | なし |

`--backend-model NAME` は、選んだ `--backend` のモデルを替えます（例：`--backend uma --backend-model uma-m-1p1`）。

実行時には、読み込むバックエンドとモデルが `[backend] Preparing MLIP model (<backend> / <model>)...` の形で表示されます。デフォルトのモデルでは括弧の中が `UMA / UMA-S-1.2 (OMol)`、`ORB / ORB-v3-conservative-OMol`、`MACE / MACE-OMOL-0`、`AIMNet2 / aimnet2` になります。

(ja-precision-by-gpu-class)=
### 精度（precision）

`--precision fp32|fp64` は、どのバックエンドでも MLIP の推論の浮動小数点精度を決めます。値は UMA の `precision`、ORB の `precision`、MACE の `default_dtype` に渡されます。AIMNet2 には精度の設定がなく、`fp32` を指定しても何も変わりません。

`--precision` を指定しないときは、バックエンドごとのデフォルト値を使います。

| バックエンド | デフォルト値 | 理由 |
|---------|------|------|
| `uma` | fp32 | 上流の fairchem の基準の設定です。 |
| `orb` | fp64 | ORB の `float32-high` を使うときは `--precision fp32` を明示します。 |
| `mace` | fp64 | MACE は上流で `default_dtype="float64"` をデフォルトにしています。 |
| `aimnet2` | fp32 | 精度の設定がありません。 |

どの値を選ぶかは目的で決めます。

| 目的 | 推奨 | 理由 |
| --- | --- | --- |
| 通常の計算 | 指定しない | 上のデフォルト値（UMA・AIMNet2 は fp32、ORB・MACE は fp64）のままにします。 |
| 速さを優先するスクリーニング | 必要なときだけ `--precision fp32` | ORB・MACE の精度が下がります（[使用上の注意点](#使用上の注意点)）。 |
| 最終の TS と Hessian | 指定しない。UMA で n_imag ≥ 2 のときは `--precision fp64` と比べる（{ref}`tsopt <ja-wrong-imaginary-mode-count>`） | 精度によらず、`tsopt` の最後の Hessian で n_imag を確かめ、IRC と端点の最適化で TS が狙った R と P をつなぐことを確かめます。 |

fp64 は次のように指定します。

```bash
pdb2reaction tsopt -i ts.pdb -q 0 --precision fp64 ...
pdb2reaction freq  -i opt.pdb -q 0 --precision fp64 ...
pdb2reaction irc   -i ts.pdb -q 0 --precision fp64 ...
```

YAML では次のように書きます。

```yaml
calc:
  precision: fp64
```

## 決定論的実行と再現性

`--deterministic` を付けると、同じソフトウェアと GPU で、同じ入力から同じ結果が得られます。付けないと、GPU では同じ入力の 2 回の計算でも最後の桁が違うことがあります。

`--deterministic` は、PyTorch の決定論的アルゴリズム（`torch.use_deterministic_algorithms`）を有効にし、GPU で決定論的に動く版の無い PyTorch の演算 1 つを置き換えます。DFT の計算、自作の ASE calculator、PyTorch の外の GPU のコードは制御しません。

```bash
pdb2reaction opt -i input.pdb -q 0 --deterministic
pdb2reaction all -i r.pdb p.pdb -q -1 --tsopt --deterministic
```

- プロセス全体に効きます。`all` に付ければ `all` が実行する MLIP の段すべてに効くので、段ごとに付ける必要はありません。
- 計算は遅くなります。決定論的な GPU の演算は処理が遅いので、繰り返しの計算を一致させる必要があるときだけ使ってください。
- 実行する演算に PyTorch の決定論的な版が無いときは、再現しない結果を黙って出さずに、エラーで止まります。

| バックエンド | `--deterministic` |
|---|---|
| `uma` | 対応 |
| `orb` / `mace` | PyTorch の決定論的モードは有効になります。入れた版で 2 回の計算が一致するか確かめてください |
| `aimnet2` | **非対応**：エラーで止まります（[使用上の注意点](#使用上の注意点)） |
| `custom` | 自作の ASE calculator しだいで、このフラグでは保証できません |

## ワーカーと Hessian の計算方式

`--uma-workers N`（デフォルト 1）は UMA の予測器を N 個並列に動かし、`--uma-workers-per-node`（デフォルト 1）はそのうち 1 ノードで動かす数を決めます。どちらのフラグも `opt`、`tsopt`、`freq`、`irc`、`sp`、`all`、`path-opt`、`path-search`、`scan`、`scan2d`、`scan3d` にあります。ORB、MACE、AIMNet2 はこれらを警告を出して無視します。ワーカーを増やすと速くなるかは [HPC 実行例 › ウォールタイム見積り](hpc-example.md#ウォールタイム見積り) にあります。

(ja-hessian-evaluation)=
### Hessian の計算方式

`--hessian-calc-mode` で Hessian の計算方式を選びます。`FiniteDifference`（デフォルト）は力の中心差分を取り、`Analytical` は選んだデバイスで 2 階の自動微分を行います。UMA、ORB、MACE、AIMNet2、DFT は解析 Hessian を計算でき、自作の calculator は `FiniteDifference` だけに対応します。

## xTB 溶媒補正

MLIP のバックエンドでは、`--solvent NAME` で xTB の溶媒和エネルギー `E_xTB(solvent) - E_xTB(vacuum)` と、それに対応する力と Hessian の差を MLIP のポテンシャル面に足します。主な用途は小分子の溶液中の計算です。`--solvent-model` で `alpb`（デフォルト）か `cpcmx` を選びます。`--solvent-xtb-cmd` には追加の引数を含めた xTB のコマンドを渡せます。xTB の SCC（自己無撞着電荷）の反復が収束しないときは、`'xtb --etemp 1000'` のように指定してください。

補正のたびに xTB を 2 回（溶媒中と真空中）実行するので、200〜300 原子くらいが実用の目安です（上限ではありません）。実際の系と計算機で計算時間を測ってください。この補正は名前の付いたバルクの溶媒を表すので、酵素クラスターに使うのはその環境を仮定する根拠があるときだけにしてください。溶液中の障壁とクラスターの障壁を比べるときは、反応する化学種、電荷、多重度、バックエンドとモデル、エネルギーの基準をそろえてください。

## DFT バックエンド

`sp`、`opt`、`tsopt`、`irc`、`freq`、`scan`、`scan2d`、`scan3d`、`path-opt`、`path-search`、`all` は `-b dft --func-basis FUNCTIONAL/BASIS --dft-engine gpu|cpu`（デフォルトは `wb97m-v/def2-svp` と `gpu`）を受け付け、すべてのエネルギーと力を PySCF/GPU4PySCF で計算します。ポピュレーション解析付きの一点計算には、別の `pdb2reaction dft` コマンドがあります。低メモリモード、CPU のスレッド数とホスト RAM、SCF のチェックポイントは [MLIP の TS を DFT で確かめる](dft-backend.md) にあります。

(ja-backends-custom-calculator)=
## カスタムバックエンド — 任意の ASE Calculator を使う（`--calc-file`）

組み込みの MLIP バックエンドのほかに、`--calc-file` で任意の [ASE](https://wiki.fysik.dtu.dk/ase/) Calculator を実行時に使えます。pdb2reaction 本体を変える必要はありません。GFN-xTB（`tblite` か `xtb-python` 経由）、DFTB+、ORCA、Psi4 など、ASE に対応した計算エンジンをつなげます。受け渡しは標準の ASE Calculator の形（エネルギーは eV、力は eV/Å）です。

ASE Calculator を返す `get_calculator` 関数を持つ Python ファイルを書きます。

```python
# my_calc.py（最小の例）
from ase.calculators.emt import EMT

def get_calculator(charge=0, spin=1, device="auto", **kwargs):
    return EMT()
```

`EMT()` を、GFN-xTB の `tblite.ase.TBLite(...)`、DFTB+ の ASE calculator、`ase.calculators.orca.ORCA(...)` など、使いたいエンジンに替えてください。このファイルを各コマンドか `all` に渡すと `custom` バックエンドが選ばれ、`--backend` の指定より優先されます。

```bash
pdb2reaction sp     -i model.xyz --calc-file my_calc.py -q 0 -m 1
pdb2reaction opt    -i model.xyz --calc-file my_calc.py -q 0 -m 1
pdb2reaction tsopt  -i ts.xyz    --calc-file my_calc.py -q 0 -m 1
pdb2reaction freq   -i ts.xyz    --calc-file my_calc.py -q 0 -m 1
pdb2reaction all    -i R.pdb P.pdb -c 'A:LIG' --calc-file my_calc.py -q 0 -m 1
```

- 関数が引数か `**kwargs` で受け取れば、`charge`、`spin`、`device` が渡されるので、総電荷が要るエンジン（xTB など）も設定できます。`spin` は多重度で、`mult`・`multiplicity` の名前でも渡します。関数の名前は `--calc-file-func-name NAME` で変えられ、その名前に Calculator のインスタンスを置いてもかまいません。
- Hessian は力の有限差分で求めるので、`freq` と `tsopt --opt-mode hess` はどのエンジンでも動きます。`--freeze-links`・`--freeze-atoms` による固定もふつうどおり効きます。
- `all`、`sp`、`opt`、`tsopt`、`freq`、`irc`、`scan`・`scan2d`・`scan3d`、`path-opt`、`path-search` で使えます。`all` は、calculator を使うすべての段に同じ関数を渡します。独自の `--backend` 名を持つ、インストールできるバックエンドにするときは [開発者向け](#開発者向け) を見てください。

## Python API

### クイックスタート

```python
import numpy as np
from pdb2reaction.backends.uma import UMACalculator

# 例: 中性一重項の2原子系（GPUが利用可能ならGPU、なければCPU）
calc = UMACalculator(charge=0, spin=1, model="uma-s-1p2", device="auto")

# UMACalculator には Bohr 単位の座標（形状: [n_atoms, 3]）を渡します
coords_bohr = np.array([
 [0.0, 0.0, 0.0],
 [2.2, 0.0, 0.0], # 約 1.16 Å
])

symbols = ["C", "O"]

# 注: これらのメソッドは dict を返すため、適切なキーで値を取り出します
energy_h = calc.get_energy(symbols, coords_bohr)["energy"] # float (Hartree)
forces_h_bohr = calc.get_forces(symbols, coords_bohr)["forces"] # ndarray (Hartree/Bohr)
hessian_h_bohr2 = calc.get_hessian(symbols, coords_bohr)["hessian"] # ndarray (Hartree/Bohr²)
```

- 座標は **Bohr** で与えます。内部で Å に変換して UMA で計算し、エネルギーと微分を Hartree、Hartree/Bohr、Hartree/Bohr² に戻します。
- `device="auto"` は、CUDA が使えれば GPU を、使えなければ CPU を選びます。
- `pysisyphus`（同梱の最適化ライブラリ）の geometry オブジェクトに付けるか、上のように直接呼び出します。

### Calculator ファクトリ

`backends` モジュールには、MLIP の calculator をプログラムから作るファクトリがあります。

```python
from pdb2reaction.backends import create_calculator, create_ase_calculator
```

| 関数 | 説明 |
|----------|-------------|
| `create_calculator(backend="uma", **kwargs)` | pysisyphus 用の MLIP calculator を作ります。受け付ける kwargs はバックエンドごとに違い、選んだバックエンドが受け付けないキーは警告を出して捨てます。ただし UMA 用のキー（`task_name`、`max_neigh` など）は、ほかのバックエンドでは警告なしに捨てます。Python から直接渡す `freeze_atoms` の番号は 0 始まりです（CLI と YAML は 1 始まりの値を変換します）。 |
| `create_ase_calculator(backend="uma", **kwargs)` | ASE 用の calculator を作ります。受け付ける kwargs はバックエンドごとに違い、使えないキーは警告なしに捨てられます。UMA・ORB・MACE の calculator は構造ごとの電荷とスピンを `atoms.info` から読み、AIMNet2 は `charge`・`spin` を作るときの引数で受け取ります。 |

```python
from pdb2reaction.backends import create_calculator

# 解析 Hessian 付きの UMA の calculator
calc = create_calculator(
    backend="uma",
    charge=0,
    spin=1,
    device="auto",
    hessian_calc_mode="Analytical",
)

```

返される calculator は pysisyphus の calculator の形を持ちます。`get_energy`、`get_forces`、`get_hessian` は `(atoms: List[str], coords: np.ndarray)` を受け取り、coords は **Bohr** 単位です。戻り値は `"energy"`（Hartree）、`"forces"`（Hartree/Bohr）、`"hessian"`（Hartree/Bohr²）を持つ dict です。固定原子の力は 0 になり、Hessian は動ける原子のブロック（`return_partial_hessian=True`）か、固定原子の行と列を 0 にした全体の行列です。

(ja-configuration-reference)=
### 設定の一覧

calculator の主な引数です。YAML の `calc` 節のキーと同じです。xTB の溶媒のキーを含む `calc` のすべてのキーは [YAML 設定の一覧 › calc](yaml-reference.md#calc) にあります。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `backend` | MLIP バックエンド | `"uma"` |
| `charge` | 総電荷。YAML に書いたときだけ使い、`-q`・`-l` が優先 | なし |
| `spin` | スピン多重度（2S+1） | `1` |
| `model` | 選んだバックエンドのモデル（`--backend-model`）。UMA のデフォルトのままなら、ORB・MACE・AIMNet2 は `orb_v3_conservative_omol`・`MACE-OMOL-0`・`aimnet2` を使います | `"uma-s-1p2"` |
| `precision` | MLIP の数値精度（`"fp32"` か `"fp64"`）。`"auto"` は UMA・AIMNet2 で fp32、ORB・MACE で fp64 | `"auto"` |
| `task_name` | UMA のバッチに記録するタスクの名前 | `"omol"` |
| `device` | `"cuda"`、`"cpu"`、または自動選択 | `"auto"` |
| `workers` / `workers_per_node` | 並列の UMA 予測器（UMA だけ。ORB・MACE・AIMNet2 は警告を出して無視します） | `1` / `1` |
| `max_neigh`, `radius`, `r_edges` | UMA の近傍の作り方の上書き | `None`, `None`, `False` |
| `freeze_atoms` | 固定する原子の番号。Python API では 0 始まり（CLI と YAML は 1 始まり） | _None_ |
| `hessian_calc_mode` | Hessian の計算方式（`"Analytical"` か `"FiniteDifference"`） | `"FiniteDifference"` |
| `return_partial_hessian` | 全体の行列でなく、動ける原子の Hessian だけを返す | `True` |
| `hessian_double` | Hessian を float64 で組み立てて返す | `True` |
| `out_hess_torch` | Hessian を `torch.Tensor` で返す | `True` |
| `print_timing` | Hessian の計算時間の内訳を表示 | `True` |
| `print_vram` | Hessian の計算中の CUDA VRAM の使用量を表示（UMA だけ） | `True` |

## 開発者向け

### バックエンドディスパッチャのパターン

```python
from pdb2reaction.backends import create_calculator, create_ase_calculator

calc = create_calculator(
    backend="uma",        # one of: "uma", "orb", "mace", "aimnet2", "auto"
    charge=0, spin=1,
    device="cuda", workers=1,
    model="uma-s-1p2",
)
# calc is a pysisyphus-compatible MLIPCalculator.

# ASE-based stages (e.g. DMF path optimization) use the ASE factory:
ase_calc = create_ase_calculator(backend="uma", model="uma-s-1p2", device="cuda")
```

pysisyphus を使う構造と経路の段は `create_calculator(...)` を、DMF（Direct Max Flux）のように ASE を使う段は `create_ase_calculator(...)` を使います。`backend="auto"` は UMA、ORB、MACE、AIMNet2 の順に試し、最初にインポートできたものを使います。YAML の `calc.backend: auto` も同じで、`-b` は `auto` を受け付けません。

### ファイルマップ

| ファイル | 役割 |
|------|------|
| `pdb2reaction/backends/__init__.py` | `BACKEND_REGISTRY` の dict、`create_calculator()`・`create_ase_calculator()` のファクトリ、UMA から順に試す `resolve_backend('auto')` |
| `pdb2reaction/backends/base.py` | `MLIPCalculator(pysisyphus.calculators.Calculator)` の基底クラス：固定原子の扱い、有限差分 Hessian の組み立て、単位の変換、バックエンドの失敗を表すエラーの型 |
| `pdb2reaction/backends/uma.py` | UMA（Meta FAIR の fairchem-core）：autograd の Hessian |
| `pdb2reaction/backends/orb.py` | Orb（Orbital Materials）：precision と compile_model |
| `pdb2reaction/backends/mace.py` | MACE：default_dtype |
| `pdb2reaction/backends/aimnet2.py` | AIMNet2：電荷を入力に取るモデル |
| `pdb2reaction/backends/pyscf_dft.py` | PySCF/GPU4PySCF の DFT/HF の calculator。段の間で SCF の状態を引き継ぎ、同じ座標の結果を再利用し、必要なら SCF のチェックポイントを書きます |

組み込みのバックエンドを独自の `--backend` 名で足すときは、[CONTRIBUTING](https://github.com/t-0hmura/pdb2reaction/blob/main/CONTRIBUTING.md) のレシピ 3.2「Add an MLIP backend」に従ってください。

## 使用上の注意点

- ORB と MACE の `--precision fp32` はスクリーニング専用です。有限差分 Hessian のノイズが増えるので、結果を使う前に n_imag を確かめてください。
- AIMNet2 は `--precision fp64` にも `--deterministic` にも対応せず、どちらもエラーで止まります。AIMNet2 はモデルへの入力を float32 にし、力を PyTorch の決定論的モードの外にある独自の CUDA のコードで計算するので、力はビット単位で再現しません（エネルギーは再現します）。繰り返しの計算を一致させたいときは、UMA、ORB、MACE のいずれかで `--deterministic` を付け、同じ環境で 2 回実行して比べてください。

(ja-workers-analytical-error)=
- UMA で `--uma-workers` を 2 以上にすると、`--hessian-calc-mode Analytical` とは併用できず、エラーで止まります。並列の予測器は autograd のモデルを持たないためです。解析 Hessian には `--uma-workers 1` を、複数のワーカーには `FiniteDifference` を使ってください。
- モデルの精度と Hessian の精度は別の設定です。エネルギーと力は常に float64 で返り、Hessian もデフォルトでは float64 で組み立てます。`calc.hessian_double: false` にすると、モデル本来の dtype（ふつうは float32）で返します。`--precision fp64` のときは Hessian も常に float64 になり、設定ファイルの `hessian_double: false` は警告を出して上書きされます。
- CI のジョブや Python API の `create_calculator` では、環境変数 `PDB2REACTION_STRICT_DETERMINISTIC=1` で `--deterministic` と同じモードになります。

## 関連ドキュメント

- [アーキテクチャ](architecture.md)：ディレクトリの構成と依存の向き
- [HPC 実行例](hpc-example.md)：PBS、Open MPI、Ray で `workers`・`workers_per_node` を複数のノードに広げるテンプレート
- [MLIP の TS を DFT で確かめる](dft-backend.md)：DFT の設定、メモリ、チェックポイント
- [トラブルシューティング](troubleshooting.md)：計算が失敗したとき
- [opt](opt.md)、[path-opt](path-opt.md)、[all](all.md)：選んだバックエンドで動くコマンド
