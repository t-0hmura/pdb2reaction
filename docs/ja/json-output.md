# JSON 出力リファレンス

pdb2reaction は、AI エージェント・スクリプト・下流ツールがプログラムから扱える機械可読の JSON 出力を提供します。

## `--out-json` フラグ

MLIP を使用するすべてのサブコマンドが `--out-json / --no-out-json`（デフォルト: off）に対応しています。
有効にすると、出力ディレクトリに `result.json` が生成されます。

```bash
pdb2reaction opt -i r.pdb -q -1 --out-json --out-dir result_opt
cat result_opt/result.json | python -m json.tool
```

`all` / `path-search` は summary writer まで到達すれば `--out-json` なしで
`summary.json` を出力します。早期の CLI/input 検証で失敗した場合は、file が
作られないことがあります。

### `summary.json` ミラー

`write_result_json` は、ステージごとの `result.json` のペイロードをすべて同じディレクトリの `summary.json` にミラーリングします。これにより、エージェントスクリプトはどのサブコマンドでも単一のファイル名（`summary.json`）を読めばよく、隣に書き出される `result.json` も同一のペイロードを持ちます。

## 共通エンベロープ

shared writerが供給するfieldを示します。「任意」と明記したrowはproducerが対応dataを渡した場合だけ存在します:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `schema_version` | string | エンベロープのスキーマバージョン。現在値は `pdb2reaction.core.utils.RESULT_JSON_SCHEMA_VERSION` にあります（このドキュメント中のリテラルではなく、この定数を参照してください）。値が上がった場合は構造変更を意味します。 |
| `command` | string | leaf envelope はサブコマンド名（例: `"opt"`）、aggregate `all` / `path-search` summary は完全な invocation string |
| `pdb2reaction_version` | string | パッケージバージョン |
| `status` | string | commandごとに異なります（下記参照）。`converged` / `not_converged`（opt, tsopt）、`completed`（irc, freq）、`ok` / `partial`（bond-summary）、`success` / `partial` / `failed`（aggregate workflow）、失敗時の `error` などです |
| `elapsed_seconds` | float | 任意の実行時間（秒）。shared writerへ時間を渡さないproducerでは省略 |
| `environment` | object | ハードウェア情報（下表参照） |
| `mlip_backend` | string | 任意。payloadがMLIP stageを表しbackend provenanceを渡した場合に追加 |
| `mlip_model` | string \| null | 任意。backendと分離した正確なmodel/checkpoint名 |

各leaf schemaの`backend` / `model`も維持されますが、command横断の処理では
`mlip_backend` / `mlip_model`を使用してください。

### Rigid projection provenance

`freq`、`irc`、`tsopt`のresultは`rigid_projection` objectを含み、`opt`では`--flatten`実行時に含みます。`freq --dump`は同じobjectを`thermoanalysis.yaml`にも書きます。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `treatment` | string | 選択した`--tr-projection` mode: `"constrained"` / `"legacy-active"` |
| `algorithm` | string | 射影kernelの識別子 |
| `effective_rank` | int | active Hessianから除去した剛体方向の数 |
| `full_rigid_rank` | int | 凍結anchor拘束前の全系剛体basisのrank |
| `frozen_constraint_rank` | int | 凍結anchor拘束によって除かれたrank |
| `active_atom_count` / `frozen_atom_count` | int | active／frozen原子数 |
| `active_atoms` / `frozen_atoms` | int[] | 射影kernelが使用した0始まり原子index |
| `hessian_space` | string | 入力Hessian空間: `"full"` / `"active"` |
| `hessian_source` / `source` | string | Hessian provenance。`freq`/`irc`は`hessian_source`、`opt`/`tsopt`は`source`を使用 |
| `hessian_shape` / `raw_hessian_shape` | int[2] | 入力Hessian shape。`freq`/`irc`は`hessian_shape`、`opt`/`tsopt`は`raw_hessian_shape`を使用 |

`constrained`は凍結anchorを動かさない全系剛体運動だけを除去します。`legacy-active`はisolated-active比較用であり、near-linear／縮退構造に対する旧結果のbitwise replayは保証しません。詳細は[凍結原子](freeze-atoms.md#凍結境界での剛体モード)を参照してください。

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

捕捉されたruntime exceptionでは、`"status": "error"` と `"error_type"` を
含む `result.json` をbest effortで書きます。usage/validation exitやout-dir確定前の
失敗ではJSONが作られない場合があるため、nonzero exit codeまたは期待JSONの欠損も
failure signalです。stderr/job logを確認してください。

収束しなかったが完了したジョブでは、`"status": "not_converged"` と最終 force/step 値を含む `result.json` が書き出されるため、AI エージェントはサイクル数を増やして再試行するかどうかをこの情報をもとに判断できます。

## サブコマンド別スキーマ

### `sp`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"ok"` |
| `stage` | string | `"sp"` |
| `input` | string | 準備後の公開入力path |
| `backend` / `model` | string / string \| null | MLIP provenance（`mlip_backend` / `mlip_model` にもmirror） |
| `charge` / `spin` | int / int | 総電荷とspin多重度 |
| `energy_au` | float | 一点energy (Hartree) |
| `forces_path` | string | `forces.npy` のpath |
| `hessian_path` | string \| null | `hessian.npy` のpath。`--hess` 無指定時はnull |
| `elapsed` | string | 人間可読の経過時間 |

### `opt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"converged"` / `"not_converged"` |
| `energy_hartree` | float | 最終エネルギー (Hartree) |
| `n_opt_cycles` | int | 最適化サイクル数 |
| `opt_mode` | string | `"grad"` / `"hess"` / `"lbfgs"` / `"rfo"` |
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
| `rigid_projection` | object | 任意。`--flatten`実行時に含む。[projection provenance](#rigid-projection-provenance)を参照 |

### `tsopt`

`opt` と同じフィールドに加え:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_imaginary_modes` | int | 虚振動の数 |
| `imaginary_frequencies_cm` | float[] | 虚振動数 (cm⁻¹, 負の値) |
| `opt_mode` | string | `"rsprfo"` (default) / `"rsirfo"` / `"trim"` / `"dimer"` |
| `reference_mode_file` | string\|null | `--ref-mode` で渡した advanced path-mode ファイル。通常は `all` が生成して内部指定します |
| `safeguards` | object | mode-loss rejection、exact saddle check、最終目的モードの index/overlap、高次鞍点でのMEP再anchor flag、停止理由、有界 path-mode restart の診断情報 |
| `rigid_projection` | object | 剛体モードとexact Hessianのprovenance。[projection provenance](#rigid-projection-provenance)を参照 |

`files` には `imaginary_mode_files`（vib ファイルリスト）を含む場合があります。
Cartesian Hessian mode の `status: "converged"` には、最終 exact PHVA で
有意な虚振動がちょうど1個必要です。対応する内部座標 mode では、代わりに
exact optimizer-space Hessian の負の固有値が1個であることを要求します。
高次鞍点と `n_imag=0` 構造は `not_converged` です。収束詳細 (force/step) は
rsirfo モードで利用可能です。
dimer モードも `status` に `"converged"` / `"not_converged"` を返しますが、
`n_opt_cycles` のみを出力し、Hessian mode の収束詳細と `safeguards` は
省略されます。

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
| `rigid_projection` | object | 剛体モードとHessian provenance。`--dump`時は`thermoanalysis.yaml`にも記録 |

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
| `energy_first_hartree` | float | stitched path の最初の端点エネルギー。standalone IRC は化学的identityを割り当てない |
| `energy_ts_hartree` | float | TS エネルギー |
| `energy_last_hartree` | float | stitched path の最後の端点エネルギー。standalone IRC は化学的identityを割り当てない |
| `endpoint_energy_orientation` | string | `"finished_first_to_finished_last"` |
| `energy_reactant_hartree` / `energy_product_hartree` | float | first / last の旧alias。key名から化学的R/P identityを推定しないこと |
| `forward_converged` | bool \| null | 前方 IRC 収束? インテグレータがフラグを公開しない場合は `null` |
| `backward_converged` | bool \| null | 後方 IRC 収束? インテグレータがフラグを公開しない場合は `null` |
| `forward_energy_increased` | bool \| null | 前方の最終stepでenergyが上昇したか |
| `backward_energy_increased` | bool \| null | 後方の最終stepでenergyが上昇したか |
| `backend` | string | MLIP バックエンド |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `model` | string | MLIP モデル名 |
| `never_stop` | bool | energy上昇／平坦化停止を無視したか |
| `n_freeze_atoms` | int | 凍結原子数 |
| `solvent` | string | 暗黙溶媒 or `"none"` |
| `bond_changes` | object | first→last 方向の `{formed: [...], broken: [...]}`。各リストは元素記号付き1始まりの原子ペア文字列（例 `"C7-O12"`）。比較が失敗または `finished_first.xyz`/`finished_last.xyz` が存在しない場合はキー自体が省略されます。 |
| `bond_changes_direction` | string | `bond_changes` がある場合は `"finished_first_to_finished_last"` |
| `step_length` | float | IRC ステップ長 (Bohr) |
| `max_cycles` | int | 最大 IRC ステップ数 |
| `input_file` | string | 入力ファイル名 |
| `files` | object | 軌跡ファイル (xyz + pdb) |
| `rigid_projection` | object | 剛体モードと初期Hessianのprovenance。[projection provenance](#rigid-projection-provenance)を参照 |

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
| `preopt` | bool | 端点pre-optimizationを有効にしたか |
| `reactant_energy_hartree` | float | 最初のimage energy (Hartree) |
| `product_energy_hartree` | float | 最後のimage energy (Hartree) |
| `image_energies_hartree` | float[] | 全イメージエネルギー |
| `n_images` | int | イメージ数 |
| `hei_index` | int | 最高エネルギーイメージのインデックス |
| `hei_energy_hartree` | float | HEI エネルギー |
| `barrier_kcal` | float | 前方障壁 (kcal/mol) |
| `delta_kcal` | float | 反応エネルギー (kcal/mol) |
| `files` | object | 軌跡 + HEI ファイル |

### `path-search`

`path-search` は `--out-json` フラグを持ちません。summary writerまで到達すれば
共通エンベロープ（`command`, `pdb2reaction_version`, `environment`）を持つ
`summary.json` を書き出しますが、早期のCLI/input検証ではfileが作られない場合があります。
追加fieldは以下です:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"success"` / `"partial"` |
| `n_segments` | int | 再帰 MEP のセグメント数 |
| `segments` | object[] | セグメントごとの `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, `bond_changes`（`{title: [entries]}` dict のリスト。bridge セグメントは `""`）。 |
| `energy_diagrams` | object[] | セグメントごとのラベル付きエネルギープロファイル (kcal/mol) |
| `mlip_backend` | string | バックエンド名 |
| `mlip_model` | string \| null | バックエンドと分離して記録するモデル名 |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |

`all` がさらに追加するフィールドは下の [`summary.json` セクション](#ja-summary-json-path-search-all) を参照してください。

### `dft`

> **注:** `dft` はSCF収束時 (exit 0) のみ `status: "converged"` を含む
> `result.json` を書きます。非収束時はexit code 3でJSONを作らないため、exit codeが
> signalです。汎用 `not_converged` statusは `dft` に適用しません。捕捉された想定外の
> exceptionでは標準error envelopeをbest effortで書きます。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"converged"` |
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
| `n_atoms_raw` | int | 選択残基に含まれるbackbone/truncation filter前の原子数（入力全体ではない） |
| `n_atoms_extracted` | int | truncation後に保持した原子数（cap-H追加前） |
| `total_charge` | float | 合計電荷 |
| `protein_charge` | float | タンパク質電荷 |
| `ligand_total_charge` | float | リガンド電荷合計 |
| `ion_total_charge` | float | イオン電荷合計 |
| `ion_charges` | list | `[[名前, 電荷], ...]` |
| `unknown_residue_charges` | object | `{残基名: 電荷}` |
| `n_link_hydrogens` | int | 炭素原子側に残る切断結合へ追加されたキャップ水素数 |
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
| `energy_source` | string | `"trajectory_comment"` または `"mlip_recomputed"` |
| `backend` | string または null | フレームを再計算した場合のみ MLIP backend。comment energy mode では null |
| `charge` / `multiplicity` | int または null | 再計算時の解決済み電子状態。それ以外は null |
| `solvent` / `solvent_model` | string または null | 再計算 calculator の溶媒設定。それ以外は null |
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
| `mlip_backend` | string | バックエンド名 |
| `mlip_model` | string \| null | バックエンドと分離して記録するモデル名 |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `environment` | object | ハードウェア情報 |

`all` はさらに以下を含みます:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `rate_limiting_step` | object | 全reactive segmentを完全に覆う最高method（`DFT//MLIP_Gibbs` > `DFT` > `MLIP_Gibbs` > `MLIP` > `MEP`）での最大障壁。`method` と raw `mep_barrier_kcal` を明記 |
| `overall_reaction_energy_kcal` | float | 全体反応エネルギー |
| `overall_reaction_energy_method` | string | 全体反応energyのmethod (`MEP`, `MLIP`, `MLIP_Gibbs`, `DFT`, `DFT//MLIP_Gibbs`) |
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
