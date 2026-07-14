# MLIP バックエンド

pdb2reaction は、すべてのワークフローステージ（`opt`、`scan`、`tsopt`、`freq`、`irc`、`path-search`,...）を単一の `MLIPCalculator` アダプタ経由でディスパッチします。本ページでは、バックエンドの選択方法・バックエンドごとの kwargs・新しいバックエンドの追加方法を説明します。

## バックエンドディスパッチャのパターン

```python
from pdb2reaction.backends import create_calculator, create_ase_calculator

calc = create_calculator(
    backend="uma",        # one of: "uma", "orb", "mace", "aimnet2", "auto"
    charge=0, spin=1,
    device="cuda", workers=1,
    model="uma-s-1p2",
)
# calc is a pysisyphus-compatible MLIPCalculator; pass it to any stage runner.

# ASE-based stages (e.g. DMF path optimization) use the ASE factory:
ase_calc = create_ase_calculator(backend="uma", model="uma-s-1p2", device="cuda")
```

`create_calculator(...)` は `**kwargs` を各バックエンドの
`_BACKEND_ACCEPTED_KEYS` セットと照合してフィルタするため、ワークフロー側のコードは
スーパーセットを渡しても問題なく、未知のキーは警告なく破棄されます。同じパターンが
`create_ase_calculator()` にも適用されます。

`backend="auto"` は UMA → Orb → MACE → AIMNet2 の順で解決し、最初にインポートに
成功したものを選択します。

## ファイルマップ

| file | role |
|------|------|
| `pdb2reaction/backends/__init__.py` | `BACKEND_REGISTRY` dict + `create_calculator()` / `create_ase_calculator()` ファクトリ + `resolve_backend('auto')` の UMA 優先フォールバック |
| `pdb2reaction/backends/base.py` | `MLIPCalculator(pysisyphus.calculators.Calculator)` ABC + FD-Hessian 組み立て + 単位変換 + `BackendError(RuntimeError)` |
| `pdb2reaction/backends/uma.py` | UMA (Meta FAIR fairchem-core) — autograd Hessian path |
| `pdb2reaction/backends/orb.py` | Orb (Orbital Materials) — precision / compile_model |
| `pdb2reaction/backends/mace.py` | MACE — default_dtype |
| `pdb2reaction/backends/aimnet2.py` | AIMNet2 — charge-aware（p2r 論文の 5-backend ベンチマークからは除外） |
| `pdb2reaction/backends/solvent.py` | xTB ALPB ラッパー — 任意のベース calculator を溶媒 ΔE 補正でラップ |
| `pdb2reaction/backends/xtb_alpb_correction.py` | スタンドアロンの xTB ALPB 補正（上記ラッパーとは別の API。混同しないこと） |

## バックエンドごとの特性

| backend | install | model identifier | precision option |
|---------|---------|------------------|------------------|
| `uma` | `pip install fairchem-core` + HF auth | `uma-s-1p2` / `uma-s-1p1` | `precision="fp32" \| "fp64"` |
| `orb` | `pip install orb-models` | `orb_v3_conservative_omol` | `precision="float32-high" \|...` |
| `mace` | 専用 conda env で pdb2reaction を入れた後、`pip uninstall -y fairchem-core && pip install 'mace-torch>=0.3.8'`（現行版どうしも `e3nn` 要件が競合） | `MACE-OMOL-0` | `default_dtype="float64"` |
| `aimnet2` | `pip install aimnet` | `aimnet2` | n/a |

### 精度（precision）

`--precision` はバックエンド非依存です。UMA の `precision` /
ORB の `precision` / MACE の `default_dtype` へ自動的にルーティングされます。AIMNet2 では
fp32 は何もしない指定（no-op）として扱われ、fp64 は拒否されます（モデル入力が前段で float32 にキャストされるため）。

`--precision` を指定しない場合、既定値はバックエンドごとに決まります。

| backend | 既定 | 理由 |
|---------|------|------|
| `uma` | fp32 | 上流 fairchem のベースライン。 |
| `orb` | fp64 | ORB の fp32 は縮約された `float32-high`（TF32）matmul であり、その力のノイズが有限差分 Hessian に偽の虚振動を生じさせる。 |
| `mace` | fp64 | MACE は上流で `default_dtype="float64"` を既定とする。 |
| `aimnet2` | fp32 | 精度の切り替えを持たない。 |

OMol で学習された UMA を既定の fp32 から fp64 に切り替えると、TSopt + Hessian に
無視できない影響が生じることがあります。以下で有効化します。

```bash
pdb2reaction tsopt -i ts.pdb -q 0 --precision fp64 ...
pdb2reaction freq  -i opt.pdb -q 0 --precision fp64 ...
pdb2reaction irc   -i ts.pdb -q 0 --precision fp64 ...
```

逆に `--precision fp32` はスループットのために ORB / MACE の精度を明示的に落とします（スクリーニング用途以外では非推奨）。使用する場合は、Hessian のノイズが増えるため、虚振動の本数を確認してください。

`--backend-model NAME` は選択中の `--backend` のモデル変種を上書きします
（例: `--backend uma --backend-model uma-s-1p2`）。未指定ならバックエンドデフォルトのモデルを使用します。

または YAML 設定で:

```yaml
calc:
  precision: fp64
```

`InferenceSettings` API のため `fairchem-core ≥ 2.0` が必要です。

## カスタムバックエンド — 任意の ASE Calculator を使う（`--calc-file`）

組み込みの MLIP バックエンドに加えて、`--calc-file` で任意の
[ASE](https://wiki.fysik.dtu.dk/ase/) Calculator を実行時に指定できます
（pdb2reaction 本体の変更は不要）。これにより GFN-xTB（`tblite` / `xtb-python`
経由）、DFTB+、ORCA、Psi4 など ASE 互換の任意エンジンと結合できます。境界は標準の
ASE Calculator インターフェース（エネルギー eV、力 eV/Å）です。

ASE Calculator を返す `get_calculator` ファクトリを持つ Python ファイルを用意します:

```python
# my_calc.py（最小の例）
from ase.calculators.emt import EMT

def get_calculator(charge=0, spin=1, device="auto", **kwargs):
    return EMT()
```

`EMT()` を使いたいエンジンに差し替えてください — 例えば GFN-xTB なら
`tblite.ase.TBLite(...)`、DFTB+ の ASE calculator、`ase.calculators.orca.ORCA(...)`
など。このファイルを単一ステージのサブコマンドに渡すと、`custom` バックエンドが
選択され `--backend` を上書きします:

```bash
pdb2reaction sp     -i model.xyz --calc-file my_calc.py -q 0 -m 1
pdb2reaction opt    -i model.xyz --calc-file my_calc.py -q 0 -m 1
pdb2reaction tsopt  -i ts.xyz    --calc-file my_calc.py -q 0 -m 1
pdb2reaction freq   -i ts.xyz    --calc-file my_calc.py -q 0 -m 1
```

補足:

- ファクトリには、シグネチャが受け取る場合（または `**kwargs` を宣言している場合）に
  `charge`・`spin`（多重度。`mult` / `multiplicity` でも渡されます）・`device` が
  渡されるため、全電荷が必要なエンジン（xTB など）も設定できます。ファクトリ名を
  変える場合は `--calc-factory NAME`、モジュール直下の Calculator インスタンスも
  受け付けます。
- Hessian は `MLIPCalculator` から継承する有限差分経路を使うため、`freq` や
  `tsopt --opt-mode hess` も任意エンジンで動作します。凍結原子（`--freeze-links` /
  `--freeze-atoms`）も通常どおり尊重されます。
- 単一ステージのサブコマンド（`sp`・`opt`・`tsopt`・`freq`・`irc`・`scan` /
  `scan2d` / `scan3d`・`path-opt`・`path-search`）で利用できます。独自の
  `--backend` 名を持つ恒久的なバックエンドにする場合は、以下のレシピを参照してください。

## バックエンド追加レシピ（5 ステップ）

`--backend xyz` として公開する新しいバックエンド `XYZModel` を追加するには:

1. **`pdb2reaction/backends/xyz.py` を作成**:
   `XYZCalculator(MLIPCalculator)`（pysisyphus パス）と `XYZASECalculator(...)`
   （ASE パス）を実装します。どちらも共通の kwargs である `charge / spin / device /
   freeze_atoms / hessian_calc_mode / return_partial_hessian / hessian_double /
   print_timing / model` と、バックエンド固有の kwargs（`precision`、
   `default_dtype` など）を受け付ける必要があります。
2. **`MLIPCalculator` を継承**（`backends/base.py`）し、サブクラスフックである
   `_compute_energy_forces_ev(elem, coord_ang) -> (energy_eV,
   forces_eV_Ang)` を実装します。バックエンドが解析的 Hessian を公開している場合は、
   必要に応じて `_compute_analytical_hessian_ev(elem,
   coord_ang) -> hessian_eV_Ang2` をオーバーライドします。
   そうでなければ、基底クラスが FD-Hessian の組み立て + 単位変換
   （eV/Å → Hartree/Bohr）を行います。
3. **`BACKEND_REGISTRY` に登録**（`backends/__init__.py`）: 既存のエントリの隣に
   `"xyz": {"module": "pdb2reaction.backends.xyz", "pysis_cls": "XYZCalculator",
   "ase_cls": "XYZASECalculator"}` を追加します。
4. **受け付ける kwargs を宣言**: `_BACKEND_ACCEPTED_KEYS["xyz"]` と
   `_ASE_ACCEPTED_KEYS["xyz"]` にセットを追加します。`resolve_backend('auto')` の
   フォールバック順に参加させたい場合は `"xyz"` を追加します。
5. **ドキュメント化 + smoke**: 本ページのファイルマップ / バックエンドごとの
   表にエントリを追加し、model identifier + インストールコマンドを記載し、
   新しいバックエンドが end-to-end で実行されるよう `tests/smoke/run.sh` に
   `xyz` の行を追加します。

## 関連項目

- [アーキテクチャ](architecture.md) — 6 層のディレクトリマップ + 依存方向。
- [CONTRIBUTING](https://github.com/t-0hmura/pdb2reaction/blob/main/CONTRIBUTING.md) — Recipe 3.2「Add an MLIP backend」、ゲートサイクルの完全な参照付き。
