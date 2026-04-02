# JSON 出力リファレンス

pdb2reaction は AI エージェント、スクリプト、下流ツールからの利用に向けた機械可読 JSON 出力を提供します。

## `--out-json` フラグ

MLIP を使用するすべてのサブコマンドが `--out-json / --no-out-json`（デフォルト: off）に対応しています。
有効にすると、出力ディレクトリに `result.json` が生成されます。

```bash
pdb2reaction opt -i r.pdb -q -1 --out-json --out-dir result_opt
cat result_opt/result.json | python -m json.tool
```

`all` / `path-search` は常に `summary.json` を出力します（`--out-json` 不要）。

## 共通エンベロープ

すべての `result.json` に自動付与されるフィールド:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `command` | string | サブコマンド名（例: `"opt"`） |
| `pdb2reaction_version` | string | パッケージバージョン |
| `elapsed_seconds` | float | 実行時間（秒） |
| `environment` | object | ハードウェア情報（下表参照） |

**`environment`**:

| フィールド | 型 | 例 |
|-----------|------|------|
| `device` | string | `"cuda"` または `"cpu"` |
| `gpu_name` | string | `"NVIDIA GeForce RTX 5080"` |
| `gpu_vram_gb` | float | `16.6` |
| `cuda_version` | string | `"12.9"` |
| `cpu` | string | `"AMD Ryzen 9 7950X 16-Core Processor"` |
| `n_cpus` | int | `32` |
| `ram_gb` | float | `133.7` |

## サブコマンド別スキーマ

### `opt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"converged"` / `"not_converged"` |
| `energy_hartree` | float | 最終エネルギー (Hartree) |
| `n_opt_cycles` | int | 最適化サイクル数 |
| `opt_mode` | string | `"grad"` / `"hess"` |
| `backend` | string | MLIP バックエンド |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `model` | string | MLIP モデル名 |
| `n_atoms` | int | 原子数 |
| `n_freeze_atoms` | int | 凍結原子数 |
| `solvent` | string | 暗黙溶媒 or `"none"` |
| `thresh` | string | 収束閾値プリセット名 |
| `max_cycles` | int | 最大サイクル数 |
| `input_file` | string | 入力ファイル名 |
| `final_max_force` | float | 最終 max gradient (Hartree/Bohr) |
| `final_rms_force` | float | 最終 RMS gradient |
| `final_max_step` | float | 最終 max 変位 (Bohr) |
| `final_rms_step` | float | 最終 RMS 変位 |
| `convergence_thresholds` | object | 収束閾値の数値 |
| `files` | object | 出力ファイルマップ |

### `tsopt`

`opt` と同じフィールドに加え:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_imaginary_modes` | int | 虚振動数 |
| `imaginary_frequencies_cm` | float[] | 虚振動数 (cm$^{-1}$, 負の値) |
| `opt_mode` | string | `"rsirfo"` / `"dimer"` |

`files` には `imaginary_mode_files`（vib ファイルリスト）を含む場合があります。

### `freq`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"completed"` |
| `n_modes` | int | 全基準振動数 |
| `n_imaginary` | int | 虚振動数 |
| `frequencies_cm` | float[] | 全振動数 (cm$^{-1}$) |
| `imaginary_frequencies_cm` | float[] | 負の振動数のみ |
| `thermochemistry` | object\|null | 熱化学データ（下表参照） |
| `backend` | string | MLIP バックエンド |
| `n_atoms` | int | 原子数 |
| `temperature_K` | float | 温度 (K) |
| `pressure_atm` | float | 圧力 (atm) |
| `files` | object | `{"frequencies_txt": "frequencies_cm-1.txt"}` |

**`thermochemistry`** (thermoanalysis 利用不可時は null):

| フィールド | 型 | 単位 |
|-----------|------|------|
| `electronic_energy_ha` | float | Hartree |
| `zpe_correction_ha` | float | Hartree |
| `thermal_correction_energy_ha` | float | Hartree |
| `thermal_correction_enthalpy_ha` | float | Hartree |
| `thermal_correction_free_energy_ha` | float | Hartree |
| `sum_EE_and_ZPE_ha` | float | Hartree |
| `sum_EE_and_thermal_energy_ha` | float | Hartree |
| `sum_EE_and_thermal_enthalpy_ha` | float | Hartree |
| `sum_EE_and_thermal_free_energy_ha` | float | Hartree |
| `E_thermal_cal_per_mol` | float | cal/mol |
| `Cv_cal_per_mol_K` | float | cal/(mol K) |
| `S_cal_per_mol_K` | float | cal/(mol K) |

### `irc`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"completed"` |
| `n_frames_forward` / `backward` / `total` | int | IRC フレーム数 |
| `energy_reactant_hartree` | float | 反応物エネルギー |
| `energy_ts_hartree` | float | TS エネルギー |
| `energy_product_hartree` | float | 生成物エネルギー |
| `forward_converged` / `backward_converged` | bool | IRC 収束判定 |
| `backend` | string | MLIP バックエンド |
| `bond_changes` | object | `{formed: [...], broken: [...]}` |
| `files` | object | 軌跡ファイル (xyz + pdb) |

### `scan`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"completed"` |
| `n_stages` | int | スキャンステージ数 |
| `stages` | object[] | ステージごとのデータ |
| `backend` | string | MLIP バックエンド |
| `files` | object | 出力ファイル |

**`stages[]`**: `n_steps`, `converged`, `pairs_1based`, `energies_hartree`, `final_energy_hartree`, `bond_changes`

### `scan2d` / `scan3d`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_grid_points` | int | グリッド点数 |
| `grid_shape` | int[] | グリッド次元 |
| `pair1`, `pair2` (,`pair3`) | object | `{i, j, low, high}` |
| `min_energy_hartree` | float | 表面最小エネルギー |
| `backend` | string | MLIP バックエンド |
| `files` | object | CSV + プロットファイル |

### `path-opt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `converged` | bool | 収束判定 |
| `mep_mode` | string | `"dmf"` / `"gsm"` |
| `image_energies_hartree` | float[] | 全イメージエネルギー |
| `barrier_kcal` | float | 前方障壁 (kcal/mol) |
| `delta_kcal` | float | 反応エネルギー (kcal/mol) |
| `backend` | string | MLIP バックエンド |
| `files` | object | 軌跡 + HEI ファイル |

### `dft`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `converged` | bool | SCF 収束? |
| `energy_hartree` | float | DFT エネルギー |
| `xc_functional` | string | 汎関数 |
| `basis_set` | string | 基底関数 |
| `used_gpu` | bool | GPU 使用? |
| `charges` | object | `{mulliken, lowdin, iao}` 原子電荷配列 |
| `spin_densities` | object | `{mulliken, lowdin, iao}` スピン密度配列 |
| `n_atoms` | int | 原子数 |
| `grid_level` | int | DFT グリッドレベル |
| `conv_tol` | float | SCF 収束閾値 |

### `extract`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"ok"` |
| `n_atoms_extracted` | int | 抽出後の原子数 |
| `total_charge` | float | 合計電荷 |
| `protein_charge` | float | タンパク質電荷 |
| `ligand_total_charge` | float | リガンド電荷合計 |
| `ion_total_charge` | float | イオン電荷合計 |
| `unknown_residue_charges` | object | `{残基名: 電荷}` |
| `center` | string | 中心残基 |
| `radius` | float | 抽出半径 (angstrom) |

## `summary.json` (`path-search` / `all`)

`all` / `path-search` は `summary.json` を出力します:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"success"` / `"partial"` |
| `segments` | object[] | セグメントごとの障壁、反応エネルギー、結合変化 |
| `energy_diagrams` | object[] | エネルギーダイアグラム（ラベル + kcal/mol） |
| `mlip_backend` | string | モデル名 |
| `environment` | object | ハードウェア情報 |

`all` はさらに以下を含みます:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `rate_limiting_step` | object | 律速段階のセグメント番号と障壁 |
| `overall_reaction_energy_kcal` | float | 全体反応エネルギー |
| `post_segments` | list | セグメントごとの TS/IRC/freq/DFT 結果 |

## 使用例

### Python

```python
import json

with open("result_opt/result.json") as f:
    result = json.load(f)

if result["status"] == "converged":
    print(f"Energy: {result['energy_hartree']:.6f} Hartree")
else:
    print(f"Not converged after {result['n_opt_cycles']} cycles")
```

### jq

```bash
# 収束確認
jq '.status' result.json

# 障壁エネルギー取得
jq '.barrier_kcal' result.json

# 虚振動数の確認
jq '.imaginary_frequencies_cm' result.json

# 自由エネルギー取得
jq '.thermochemistry.sum_EE_and_thermal_free_energy_ha' result.json
```
