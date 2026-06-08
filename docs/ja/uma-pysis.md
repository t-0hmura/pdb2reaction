# MLIP 計算機

## 概要
`pdb2reaction` は複数の機械学習原子間ポテンシャル（MLIP）を pysisyphus 向けの計算機バックエンドとしてサポートします。デフォルトバックエンドは **UMA**（Meta の Universal Models for Atoms）ですが、**ORB**、**MACE**、**AIMNet2** も利用可能です。各バックエンドはエネルギー / 力 / ヘシアンを Hartree 単位で返し、GPU/CPU ディスパッチと bohr↔Å 変換を内部処理します。クラスターサイズ（数百原子）の系では、augmented-Hessian 固有値解（RFO / RS-I-RFO）と IRC 伝播が扱う 3N × 3N の Hessian テンソルを on-device に保持しないと、ホスト–デバイス同期の繰り返しで wall-clock が支配されがちです。そのため Hessian は pysisyphus 経由で on-device に保持します（エネルギー / 力 はスカラー / numpy 配列で host へ戻すため、そのコストは無視できます）。これらの計算機は `pdb2reaction` の最適化・経路探索・熱化学・軌跡後処理など幅広く用いられます。

`pdb2reaction` は GPU 加速版の `pysisyphus` fork を同梱しています。RFO / RS-I-RFO 単一構造 optimizer、EulerPC IRC integrator、振動モード（Hessian）対角化が `device="cuda"`（または GPU ホストでの `"auto"`）の場合に CUDA で実行されます。CPU ホストでは upstream の NumPy/SciPy パスへ透過的にフォールバックします。

## クイックスタート
```python
import numpy as np
from pdb2reaction.backends.uma import UMACalculator

# 例: 中性一重項の2原子系（GPUが利用可能ならGPU、なければCPU）
calc = UMACalculator(charge=0, spin=1, model="uma-s-1p1", device="auto")

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

- 座標は **Bohr** で与えます。ラッパー内部で Å に変換し、UMA 計算後に Hartree / Hartree·Bohr⁻¹ / Hartree·Bohr⁻² に戻します。
- `pysisyphus` の geometry オブジェクトにアタッチするか、上記のように直接呼び出せます。

## Python API: Calculator Factory

`backends` モジュールは MLIP 計算機をプログラムから生成するファクトリを提供します:

```python
from pdb2reaction.backends import create_calculator, create_ase_calculator
```

| 関数 | 説明 |
|----------|-------------|
| `create_calculator(backend="uma", **kwargs)` | pysisyphus 互換の MLIP 計算機を生成します。`charge`、`spin`、`model`、`device`、`solvent`、`solvent_model`、`hessian_calc_mode`、`freeze_atoms` などバックエンド固有の kwargs も受け付けます。未知のキーはバックエンドごとに黙って除外されます |
| `create_ase_calculator(backend="uma", **kwargs)` | ASE 互換の MLIP 計算機を生成します（DMF ワークフローや ASE ベースのツールで使用）。kwargs は `create_calculator` と同じです |

### 例

```python
from pdb2reaction.backends import create_calculator

# 解析ヘシアン付き UMA 計算機
calc = create_calculator(
    backend="uma",
    charge=0,
    spin=1,
    device="auto",
    hessian_calc_mode="Analytical",
)

# 暗黙溶媒付き ORB 計算機
calc_orb = create_calculator(
    backend="orb",
    charge=-1,
    spin=1,
    solvent="water",
    solvent_model="alpb",
)
```

返される計算機は pysisyphus の計算機インターフェース（`get_energy`, `get_forces`, `get_hessian`）を実装しています。すべてのメソッドは `(atoms: List[str], coords: np.ndarray)` を受け取り、coords は **Bohr** 単位です。戻り値は `"energy"`（Hartree）、`"forces"`（Hartree/Bohr）、`"hessian"`（Hartree/Bohr²）を含む dict です。

## バックエンド選択

任意のコマンドで `-b/--backend` を指定するか、YAML で `calc.backend` を設定します:

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

| バックエンド | インストール | 解析ヘシアン | マルチワーカー | 備考 |
|---------|---------|-------------------|-------------|-------|
| **UMA** | 同梱 | あり（autograd） | あり | `fairchem-core` による完全機能 |
| **ORB** | `pip install "pdb2reaction[orb]"` | あり（autograd） | なし | orb-models（conservative モデルのみ） |
| **MACE** | `pip uninstall -y fairchem-core && pip install mace-torch` | あり（`calc.get_hessian`） | なし | mace-torch >= 0.3.8 |
| **AIMNet2** | `pip install "pdb2reaction[aimnet]"` | あり（native） | なし | aimnet |

### 暗黙溶媒補正

すべてのバックエンドで xTB ベースの暗黙溶媒補正が `--solvent` により利用可能です:

```bash
pdb2reaction opt -i input.pdb -q 0 --solvent water
pdb2reaction opt -i input.pdb -q 0 -b orb --solvent water --solvent-model cpcmx
```

デルタアプローチによる補正: ΔE = E_xTB(溶媒) - E_xTB(真空) を MLIP エネルギー/力/ヘシアンに加算します。`xtb` が `PATH` 上にインストールされている必要があります。

## 主な特徴

- **MLIP バックエンド** – デフォルトの UMA バックエンドは `fairchem-core` の `pretrained_mlip` ヘルパーで UMA チェックポイントを読み込み、AtomicData バッチに電荷/スピン情報を付与。代替バックエンド（ORB、MACE、AIMNet2）は `-b/--backend` で利用可能。
- **デバイス処理** – `device="auto"` は CUDA があれば GPU、なければ CPU を選択。グラフ構築は選択デバイス上で行い、`workers>1` では並列予測器が転送を管理。
- **凍結原子** – `freeze_atoms` に 1 始まりの原子インデックスを渡すと、凍結原子の力がゼロ化。`return_partial_hessian=True` で凍結自由度を除いたヘシアンを返すか、フル行列で該当行/列をゼロ化できます。
- **精度制御** – エネルギー/力は常に float64。`hessian_double=False` でヘシアンをモデルのネイティブ dtype（通常 float32）で返します。
- **マルチワーカー推論** – `workers>1` で `fairchem-core` の `ParallelMLIPPredictUnit` を起動し、`workers_per_node` をノードごとに指定可能。バッチ処理速度の向上に有効です。

(ja-workers-fd-downgrade)=
### `workers > 1` による暗黙的な FD ダウングレード（UMA バックエンド）

```{warning}
UMA バックエンドを `workers > 1` で使用する場合、`hessian_calc_mode="Analytical"` を明示的に指定していても解析ヘシアンは暗黙的に有限差分へ切り替わります（内部フラグ、ユーザーから設定不可）。**このダウングレード発生時にログマーカーは出力されません。**診断する唯一の方法は、同じ原子数の解析ヘシアン基準ランと比較してヘシアン計算時間が FD 相当に長くなっていることの確認です。このダウングレード規則は `--workers` / `--workers-per-node` を持つすべてのサブコマンド（`opt`, `tsopt`, `freq`, `irc`, `all`, `path-opt`, `path-search`, scan 系）に適用されます。UMA 以外のバックエンド（ORB / MACE / AIMNet2）では `workers` / `workers_per_node` は `_BACKEND_ACCEPTED_KEYS` で暗黙的に除外されるため、このダウングレード規則は該当しません。
```

(ja-hessian-evaluation)=
### ヘシアン評価モード

`hessian_calc_mode="Analytical"` は選択されたデバイス上で 2 階自動微分を行い、`"FiniteDifference"`（デフォルト）は力の中心差分を計算します。複数の推論ワーカーを要求した状態では、`Analytical` を明示指定していても解析ヘシアンは警告なく有限差分へダウングレードされます（上記の注記を参照）。解析ヘシアンが必要な場合は `workers = 1` を指定してください。

## HPC での使用例: PBS + Open MPI + Ray

複数ノードで `workers` / `workers_per_node` を Ray クラスタに分散する PBS + Open MPI のテンプレートは [HPC 実行例のページ](hpc-example.md) を参照してください。モジュール名、conda パス、PBS リソース要求などを環境に合わせて調整すればそのまま使えます。

(ja-configuration-reference)=
## 設定リファレンス
代表的なコンストラクタ引数（右端はデフォルト値）:

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `backend` | MLIP バックエンドエンジン | `"uma"` |
| `charge` | 総電荷 | `0` |
| `spin` | スピン多重度（2S+1） | `1` |
| `model` | UMA モデル名 (`uma-s-1p1`, `uma-m-1p1`) | `"uma-s-1p1"` |
| `task_name` | UMA バッチに記録されるタスクタグ | `"omol"` |
| `device` | `"cuda"` / `"cpu"` / `"auto"` | `"auto"` |
| `workers` / `workers_per_node` | 並列 UMA 予測器（UMA バックエンド限定。ORB / MACE / AIMNet2 では無視されます）。`workers>1` の場合、解析ヘシアンは**警告なく**有限差分へダウングレードされます。`hessian_calc_mode` のデフォルトはそもそも FD のため、`Analytical` を明示的に選んだ場合のみ影響があります | `1` / `1` |
| `max_neigh`, `radius`, `r_edges` | 近傍構築のオプション上書き | `None`, `None`, `False` |
| `freeze_atoms` | 1 始まりの凍結原子インデックス | _None_ |
| `hessian_calc_mode` | `"Analytical"` または `"FiniteDifference"` | `"FiniteDifference"` |
| `return_partial_hessian` | アクティブ自由度のみ返す | `True` |
| `hessian_double` | ヘシアンを float64 で返す | `True` |
| `out_hess_torch` | ヘシアンを `torch.Tensor` で返す | `True` |
| `print_timing` | ヘシアン計算のタイミング内訳を表示 | `True` |
| `print_vram` | ヘシアン計算中の CUDA VRAM 使用量を表示（UMA バックエンド限定） | `True` |
| `solvent` | 暗黙溶媒名（例: `"water"`）または `"none"` | `"none"` |
| `solvent_model` | xTB 溶媒モデル: `"alpb"` または `"cpcmx"` | `"alpb"` |
| `xtb_cmd` | 溶媒補正で使用する xTB 実行コマンド | `"xtb"` |
| `xtb_acc` | 溶媒補正実行時の xTB 精度設定 | `0.2` |


## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [opt](opt.md) -- MLIP バックエンドを使う単一構造最適化
- [path-opt](path-opt.md) -- MLIP バックエンドを使う MEP 最適化
- [all](all.md) -- MLIP を複数段で使う一気通貫ワークフロー
