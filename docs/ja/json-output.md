# JSON 出力リファレンス

pdb2reaction は、AI エージェント・スクリプト・下流ツールがプログラムから扱える機械可読の JSON 出力を提供します。

## `--out-json` フラグ

MLIP を使用するすべてのサブコマンドが `--out-json / --no-out-json`（デフォルト: off）に対応しています。
有効にすると、出力ディレクトリに `result.json` が生成されます。

```bash
pdb2reaction opt -i r.pdb -q -1 --out-json --out-dir result_opt
cat result_opt/result.json | python -m json.tool
```

`all` / `path-search` は常に `summary.json` を出力します（`--out-json` 不要）。

### `summary.json` ミラー

`write_result_json` は、ステージごとの `result.json` のペイロードをすべて同じディレクトリの `summary.json` にミラーリングします。これにより、エージェントスクリプトはどのサブコマンドでも単一のファイル名（`summary.json`）を読めばよく、隣に書き出される `result.json` も同一のペイロードを持ちます。

## 共通エンベロープ

すべての `result.json` に自動付与されるフィールド:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `schema_version` | string | エンベロープのスキーマバージョン。現在値は `pdb2reaction.core.utils.RESULT_JSON_SCHEMA_VERSION` にあります（このドキュメント中のリテラルではなく、この定数を参照してください）。値が上がった場合は構造変更を意味します。 |
| `command` | string | サブコマンド名（例: `"opt"`） |
| `pdb2reaction_version` | string | パッケージバージョン |
| `status` | string | `success` / `partial` / `error` / `unknown` のいずれか |
| `elapsed_seconds` | float | 実行時間（秒） |
| `environment` | object | ハードウェア情報（下表参照） |

**`environment`**:

| フィールド | 型 | 例 |
|-----------|------|------|
| `device` | string | `"cuda"` または `"cpu"` |
| `gpu_name` | string | `"<gpu model>"` |
| `gpu_vram_gb` | float | `<vram in GB>` |
| `cuda_version` | string | `"<cuda version>"` |
| `cpu` | string | `"<cpu model>"` |
| `n_cpus` | int | `<int>` |
| `ram_gb` | float | `<ram in GB>` |

### エラーエンベロープ（`status == "error"` のとき）

| フィールド | 型 | 説明 |
|-----------|------|------|
| `error` | string | 元の例外の `str(exc)` |
| `error_type` | string | 例外クラス名 |
| `error_class_chain` | list[string] | MRO 全体のクラス名。エージェントがテキストを解析せずに階層を照合できます |
| `error_module` | string | 例外クラスが定義されたモジュール |
| `error_label` | string | 高レベルの CLI ステージラベル |

## エラー処理

ジョブが失敗した場合（クラッシュ、OOM、収束失敗による `sys.exit` など）でも、`"status": "error"` と失敗種別を表す `"error_type"` を含む `result.json` が書き出されます。詳細なトレースバックは `.out` ログファイルを参照してください。失敗判定には `result.json` の不在ではなく `status == "error"` を使用してください。

収束しなかったが完了したジョブでは、`"status": "not_converged"` と最終 force/step 値を含む `result.json` が書き出されるため、AI エージェントはサイクル数を増やして再試行するかどうかをこの情報をもとに判断できます。

## サブコマンド別スキーマ

### `opt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"converged"` / `"not_converged"` |
| `energy_hartree` | float | 最終エネルギー (Hartree) |
| `n_opt_cycles` | int | 最適化サイクル数 |
| `opt_mode` | string | `"grad"` / `"hess"` |
| `backend` | string | MLIP バックエンド (`"uma"`, `"orb"`, `"mace"`, `"aimnet2"`) |
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
| `convergence_thresholds` | object | `{max_force_thresh, rms_force_thresh, max_step_thresh, rms_step_thresh}` (Hartree/Bohr) |
| `files` | object | 出力ファイルマップ |

### `tsopt`

`opt` と同じフィールドに加え:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_imaginary_modes` | int | 虚振動の数 |
| `imaginary_frequencies_cm` | float[] | 虚振動数 (cm⁻¹, 負の値) |
| `opt_mode` | string | `"rsprfo"` (default) / `"rsirfo"` / `"trim"` / `"dimer"` |

`files` には `imaginary_mode_files`（vib ファイルリスト）を含む場合があります。
収束詳細 (force/step) は rsirfo モードで利用可能です。dimer モードも `runner.is_converged` に応じて `status` に `"converged"` / `"not_converged"` を返しますが、`n_opt_cycles` のみを出力し、rsirfo モードが出すサイクルごとの力・ステップ収束の詳細キーは省略されます。

### `freq`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"completed"` |
| `n_modes` | int | 基準振動モードの総数 |
| `n_imaginary` | int | 虚モード数 |
| `frequencies_cm` | float[] | 全振動数 (cm⁻¹) |
| `imaginary_frequencies_cm` | float[] | 負の振動数のみ |
| `thermochemistry` | object\|null | 熱化学データ（下表参照） |
| `backend` | string | MLIP バックエンド |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `model` | string | MLIP モデル名 |
| `n_atoms` | int | 原子数 |
| `n_freeze_atoms` | int | 凍結原子数 |
| `solvent` | string | 暗黙溶媒 or `"none"` |
| `temperature_K` | float | 温度 (K) |
| `pressure_atm` | float | 圧力 (atm) |
| `input_file` | string | 入力ファイル名 |
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
| `n_frames_forward` | int | 前方 IRC フレーム数 |
| `n_frames_backward` | int | 後方 IRC フレーム数 |
| `n_frames_total` | int | 全フレーム数 |
| `energy_reactant_hartree` | float | 反応物エネルギー |
| `energy_ts_hartree` | float | TS エネルギー |
| `energy_product_hartree` | float | 生成物エネルギー |
| `forward_converged` | bool \| null | 前方 IRC 収束? インテグレータがフラグを公開しない場合は `null` |
| `backward_converged` | bool \| null | 後方 IRC 収束? インテグレータがフラグを公開しない場合は `null` |
| `backend` | string | MLIP バックエンド |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `model` | string | MLIP モデル名 |
| `n_freeze_atoms` | int | 凍結原子数 |
| `solvent` | string | 暗黙溶媒 or `"none"` |
| `bond_changes` | object | `{formed: [...], broken: [...]}` の各リストは元素記号付き 1 始まりの原子ペア文字列（例 `"C7-O12"`）。比較が失敗または `finished_first.xyz`/`finished_last.xyz` が存在しない場合はキー自体が省略されます。 |
| `step_length` | float | IRC ステップ長 (Bohr) |
| `max_cycles` | int | 最大 IRC ステップ数 |
| `input_file` | string | 入力ファイル名 |
| `files` | object | 軌跡ファイル (xyz + pdb) |

### `scan`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"completed"` |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `backend` | string | MLIP バックエンド |
| `model` | string | MLIP モデル名 |
| `solvent` | string | 暗黙溶媒 or `"none"` |
| `preopt` | bool | 事前最適化を実行したか |
| `max_step_size_angstrom` | float | 1 ステップ当たりの最大結合長変位 (Å) |
| `n_stages` | int | スキャンステージ数 |
| `stages` | object[] | ステージごとのデータ（下記参照） |
| `files` | object | 出力ファイル |

**`stages[]`**:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `index` | int | 1-based ステージインデックス |
| `n_steps` | int | ステップ数 |
| `converged` | bool | 拘束最適化の収束? |
| `pairs_1based` | list | 原子ペア (1-based) |
| `initial_distances_angstrom` | list | 初期距離 |
| `target_distances_angstrom` | list | 目標距離 |
| `final_energy_hartree` | float | 最終エネルギー |
| `energies_hartree` | float[] | ステップごとのエネルギー |
| `bond_changes` | object | `{"changed": bool \| null, "summary": str}`（`has_bond_change` の自由記述サマリ。比較が走らなかった場合は `null`/`""`）。 |

### `scan2d` / `scan3d`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"completed"` |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `backend` | string | MLIP バックエンド |
| `model` | string | MLIP モデル名 |
| `solvent` | string | 暗黙溶媒 or `"none"` |
| `max_step_size_angstrom` | float | 1 ステップ当たりの最大結合長変位 (Å, `scan2d` のみ) |
| `n_grid_points` | int | グリッド点数 |
| `grid_shape` | int[] | グリッド次元 (`scan3d --csv` 再プロット時には省略) |
| `pair1`, `pair2` (,`pair3`) | object | `{i, j, low, high}` (オプション: `label_i`, `label_j`)。`scan3d` で `--csv` 再プロット時は省略 |
| `min_energy_hartree` | float | 表面最小エネルギー |
| `files` | object | CSV + プロットファイル |

### `path-opt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"converged"` / `"not_converged"` / `"completed"` |
| `converged` | bool | 収束判定 |
| `mep_mode` | string | `"dmf"` / `"gsm"` |
| `backend` | string | MLIP バックエンド |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `model` | string | MLIP モデル名 |
| `solvent` | string | 暗黙溶媒 or `"none"` |
| `image_energies_hartree` | float[] | 全イメージエネルギー |
| `n_images` | int | イメージ数 |
| `hei_index` | int | 最高エネルギーイメージのインデックス |
| `hei_energy_hartree` | float | HEI エネルギー |
| `barrier_kcal` | float | 前方障壁 (kcal/mol) |
| `delta_kcal` | float | 反応エネルギー (kcal/mol) |
| `files` | object | 軌跡 + HEI ファイル |

### `path-search`

`path-search` は `--out-json` フラグを持たず、`summary.json` を出力ディレクトリに**常に**書き出します。共通エンベロープ（`command`, `pdb2reaction_version`, `environment`）に加え、以下を含みます:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"success"` / `"partial"` |
| `n_segments` | int | 再帰 MEP のセグメント数 |
| `segments` | object[] | セグメントごとの `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, `bond_changes`（`{title: [entries]}` dict のリスト。bridge セグメントは `""`）。 |
| `energy_diagrams` | object[] | セグメントごとのラベル付きエネルギープロファイル (kcal/mol) |
| `mlip_backend` | string | バックエンド / モデル名 |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |

`all` がさらに追加するフィールドは下の [`summary.json` セクション](#ja-summary-json-path-search-all) を参照してください。

### `dft`

> **注:** `dft` は SCF が収束した場合 (exit 0) のみ `result.json` を書き出します。SCF が収束しなかった場合は exit code 3 を返し `result.json` は作成されません。SCF 状態は `converged: bool` フィールドと exit code で表現され、`status` フィールドは持ちません。上記の汎用 "not_converged" ステータスは `dft` には適用されません。ただし、想定外の例外（unhandled exception）が発生した場合は、他サブコマンドと同じ標準の "error" エンベロープ（`result.json` + ミラーの `summary.json`）が書き出されます。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `converged` | bool | SCF 収束? |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `n_atoms` | int | 原子数 |
| `grid_level` | int | DFT グリッドレベル |
| `conv_tol` | float | SCF 収束閾値 |
| `max_cycle` | int | 最大 SCF サイクル数 |
| `input_file` | string | 入力ファイル名 |
| `energy_hartree` | float | DFT エネルギー |
| `energy_kcal_per_mol` | float | DFT エネルギー (kcal/mol) |
| `xc_functional` | string | 汎関数 |
| `basis_set` | string | 基底関数 |
| `engine` | string | 実効エンジンラベル (`"gpu4pyscf(rks_lowmem)"` / `"gpu4pyscf"` / `"pyscf(cpu)"`) |
| `used_gpu` | bool | GPU 使用? |
| `used_lowmem` | bool | `gpu4pyscf.dft.rks_lowmem.RKS` を実際に使用したか? (open-shell, CPU, `--no-lowmem`, あるいは `rks_lowmem` 未搭載の旧 `gpu4pyscf` では False) |
| `charges` | object | `{mulliken, lowdin, iao}` 原子電荷配列 |
| `spin_densities` | object | `{mulliken, lowdin, iao}` スピン密度配列 |
| `files` | object | `{"result_yaml": "result.yaml", "input_geometry_xyz": "input_geometry.xyz"}` |

### `extract`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"ok"` |
| `n_atoms_raw` | int | 入力 PDB の原子数 |
| `n_atoms_extracted` | int | 抽出後の原子数 |
| `total_charge` | float | 合計電荷 |
| `protein_charge` | float | タンパク質電荷 |
| `ligand_total_charge` | float | リガンド電荷合計 |
| `ion_total_charge` | float | イオン電荷合計 |
| `ion_charges` | list | `[[名前, 電荷], ...]` |
| `unknown_residue_charges` | object | `{残基名: 電荷}` |
| `n_link_hydrogens` | int | 切断された C/N 結合に追加されたキャップ水素数 |
| `exclude_backbone` | bool | バックボーンを除外したか |
| `include_h2o` | bool | 結晶水を含めたか |
| `ligand_charge_input` | string | ユーザ指定 `--ligand-charge` マッピング |
| `center` | string | 中心残基 |
| `radius` | float | 抽出半径 (angstrom) |
| `input_files` | string[] | 入力 PDB パス |
| `files` | object | 出力 PDB / クラスターファイル |

### `trj2fig`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"ok"` |
| `n_frames` | int | 軌跡フレーム数 |
| `min_energy_hartree` | float | フレーム中の最小エネルギー |
| `max_energy_hartree` | float | フレーム中の最大エネルギー |
| `backend` | string | MLIP バックエンド |
| `files` | object | 出力プロットファイル |

### `energy-diagram`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"ok"` |
| `n_points` | int | エネルギーデータ点数 |
| `files` | object | 出力ダイアグラムファイル |

### `bond-summary`

`--json` 有効時、`bond-summary` は JSON を**標準出力**に出力します（`result.json` ファイルは書き出しません。永続化したい場合は stdout をリダイレクトしてください）。上の MLIP 系サブコマンドが `out_dir` に `result.json` を書き出すのとは**異なる**挙動です:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"ok"`（全ペアが正常に比較できた場合）または `"partial"`（一部のペアが失敗した場合） |
| `comparisons` | object[] | ペアごとの比較（`structure_a` (string), `structure_b` (string), `bonds_formed` (int), `bonds_broken` (int)） |

(ja-summary-json-path-search-all)=
## `summary.json` (`path-search` / `all`)

`all` / `path-search` は、より構造化された `summary.json` を出力します:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"success"` / `"partial"` / `"failed"` (`all`); `"success"` / `"partial"` (`path-search`) |
| `n_segments` | int | セグメント数 |
| `segments` | object[] | セグメントごとの `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, `bond_changes` |
| `energy_diagrams` | object[] | エネルギーダイアグラム（ラベル + kcal/mol） |
| `mlip_backend` | string | モデル名 |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `environment` | object | ハードウェア情報 |

`all` はさらに以下を含みます:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `rate_limiting_step` | object | 律速段階のセグメント番号と障壁 |
| `overall_reaction_energy_kcal` | float | 全体反応エネルギー |
| `post_segments` | list | セグメントごとの TS/IRC/freq/DFT 結果 |
| `key_output_files` | object | 主要出力ファイル一覧 |

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
    print(f"Max force: {result['final_max_force']:.6f}")
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

## 関連項目

- [CLI 規約](cli-conventions.md) — `--out-json` / `--no-out-json` フラグの規約と終了コード
- [YAML リファレンス](yaml-reference.md) — これらのスキーマに現れる設定入力
- [all](all.md), [path-search](path-search.md) — 常に `summary.json` を書き出すサブコマンド
- [opt](opt.md), [tsopt](tsopt.md), [freq](freq.md), [irc](irc.md), [scan](scan.md), [path-opt](path-opt.md), [dft](dft.md), [extract](extract.md) — `--out-json` 指定時にのみ `result.json` を書き出すサブコマンド
