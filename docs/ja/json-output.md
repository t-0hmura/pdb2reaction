# JSON 出力リファレンス

pdb2reaction は、AI エージェント・スクリプト・下流ツールがプログラムから扱える機械可読の JSON 出力を提供します。

## `--out-json` フラグ

MLIP を使用する主なサブコマンドとレポート系サブコマンド（`opt`、`sp`、
`tsopt`、`freq`、`irc`、scan 系、`path-opt`、`dft`、`extract`、
`trj2fig`、`energy-diagram`）は `--out-json / --no-out-json`（デフォルト:
off）に対応しています。有効にすると、正規の `result.json` と、同一内容の
互換ミラー `summary.json` が通常の出力と同じ場所に生成されます。

```bash
pdb2reaction opt -i r.pdb -q -1 --out-json --out-dir result_opt
cat result_opt/result.json | python -m json.tool
```

`all` / `path-search` は、集約結果を書き込む段階まで到達すれば `--out-json` なしで
`summary.json` を出力します。早期の CLI 引数または入力の検証で失敗した場合は、ファイルが
作られないことがあります。

### `summary.json` ミラー

`write_result_json` はステージごとのペイロードを一度だけシリアライズし、
互換ミラー `summary.json` を先に、正規の `result.json` を最後にアトミックに公開
します。正常終了時は両ファイルのバイト列が同一です。I/O エラーによって公開が
中断された場合、書き込み処理はミラーの不一致を隠さず例外を送出します。
MCP の利用側は、割り当てられている場合には現在の `run_id` も検証してください。

## 共通エンベロープ

共通の書き込み処理が供給するフィールドを示します。「任意」と明記した行は、
生成側が対応するデータを渡した場合にだけ存在します。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `schema_version` | string | エンベロープのスキーマバージョン。現在値は `pdb2reaction.core.utils.RESULT_JSON_SCHEMA_VERSION` にあります（このドキュメント中のリテラルではなく、この定数を参照してください）。値が上がった場合は構造変更を意味します。 |
| `command` | string | leaf envelope はサブコマンド名（例: `"opt"`）、aggregate `all` / `path-search` summary は完全な invocation string |
| `pdb2reaction_version` | string | パッケージバージョン |
| `status` | string | commandごとに異なります（下記参照）。`converged` / `not_converged` / `stalled`（opt, tsopt）、`completed`（irc, freq）、`ok` / `partial` / `failed`（bond-summary）、`success` / `partial` / `failed`（`all`）、`success` / `partial`（`path-search`）、失敗時の `error` などです |
| `run_id` | string | 任意。現在の MCP 呼び出し UUID。プロダクト固有の実行環境が有効な場合のみ注入され、producer 側の値と衝突する場合は拒否されます。 |
| `elapsed_seconds` | float | 任意の実行時間（秒）。shared writerへ時間を渡さないproducerでは省略 |
| `environment` | object | ハードウェア情報（下表参照） |
| `mlip_backend` | string \| null | 任意のMLIP backend識別子。plot-only commandがcalculatorを評価していない場合はnull |
| `mlip_model` | string \| null | 任意。backendと分離した正確なmodel/checkpoint名 |
| `mlip_precision` | string \| null | 実効精度の共通token（`fp32` / `fp64`）。dtypeをuser codeが管理するcustom calculatorではnull |

各leaf schemaの`backend` / `model`も維持されますが、command横断の処理では
`mlip_backend` / `mlip_model` / `mlip_precision`を使用してください。

### 実行結果と科学的妥当性

複数段階のワークフローと scan の出力処理は、構成要素を評価できる場合に以下のフィールドを追加します。出力されるフィールドはコマンドによって異なり、各コマンド固有の `status` も互換性のため維持されます。結果を利用できるか判断する際は、`scientific_status` と各 outcome を確認してください。収束を確認できない個別結果は安全側に倒して扱われ、`usable` にはなりません。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `execution_status` | string | 通常は `completed` または `failed`。必須の構成コマンドが実行されたかを示します。 |
| `scientific_status` | string | `success`、`partial`、`failed`。得られた科学的結果が完全かつ利用可能かを示します。 |
| `scientific_status_reasons` | string[] | 利用できない、または欠落した個別結果の理由。正常終了時は省略されます。集約ワークフローの従来の `status_reasons` とは別です。 |
| `expected_item_ids` / `observed_item_ids` | string[] | 集約結果の欠落を検出するための、期待された項目と観測された項目の ID。 |
| `stage_outcomes` | object[] | `stage`、`item_id`、`required`、`executed`、`converged`、`usable`、`reason`、`artifacts` を持つ段階別 outcome。 |
| `point_outcomes` | object[] | `point_id`、`executed`、`converged`、`energy_valid`、`artifact_written`、`seed_eligible`、`reason` を持つ scan 点別 outcome。 |

`run_id` が存在する場合は、現在の呼び出しを識別します。`all` の集約
`summary.json` では、`current_output_paths` と `key_output_files` も、その
呼び出しが記録した成果物だけを示します。再利用した出力ディレクトリに残る
既存ファイルは現在の結果に含まれません。

### Rigid projection provenance

`freq`、`irc`、`tsopt`のresultは`rigid_projection` objectを含み、`opt`では`--flatten`実行時に含みます。`freq --dump`は同じobjectを`thermoanalysis.yaml`にも書きます。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `treatment` | string | 固定の剛体モード処理: `"constrained"` |
| `algorithm` | string | 射影kernelの識別子 |
| `effective_rank` | int | active Hessian から除去した剛体方向の数 |
| `full_rigid_rank` | int | 凍結anchor拘束前の全系剛体basisのrank |
| `frozen_constraint_rank` | int | 凍結anchor拘束によって除かれたrank |
| `svd_rtol` | float | rank 判定に用いる相対 SVD 許容値 |
| `active_atom_count` / `frozen_atom_count` | int | active／frozen原子数 |
| `active_atoms` / `frozen_atoms` | int[] | 射影kernelが使用した0始まり原子index |
| `hessian_space` | string | 入力 Hessian 空間: `"full"` / `"active"` |
| `hessian_source` / `source` | string | Hessian provenance。`freq`/`irc`は`hessian_source`、`opt`/`tsopt`は`source`を使用 |
| `hessian_shape` / `raw_hessian_shape` | int[2] | 入力 Hessian shape。`freq`/`irc`は`hessian_shape`、`opt`/`tsopt`は`raw_hessian_shape`を使用 |

`constrained`は凍結anchorを動かさない全系剛体運動だけを除去します。詳細は[凍結原子](freeze-atoms.md#凍結境界での剛体モード)を参照してください。

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

捕捉された実行時例外では、`"status": "error"` と `"error_type"` を含む
`result.json` を可能な範囲で書き出します。使用法・入力検証による終了や出力先の
確定前に失敗した場合は JSON が作られないことがあるため、0 以外の終了コードや
期待した JSON の欠落も失敗のシグナルです。標準エラー出力またはジョブログを確認してください。

制御された非収束結果の writer まで到達した job は、`"status": "not_converged"` の `result.json` を書きます。optimizer family は該当する最終 force/step または cycle field を含みますが、DFT と Dimer は所有しない field を出力しません。再試行の判断前に command-specific schema を確認してください。

オプティマイザは `"status": "stalled"` を返すこともあります。これは、設定した force/step の収束基準を満たさないまま、設定ウィンドウにわたってエネルギーが減少しなくなった状態（エネルギープラトー）です。stalled は converged とは別の非収束アウトカムであり、`converged` として報告されることは決してありません。`tsopt` では、停滞した探索を繰り返さないよう flatten/再試行ループも停止します。`opt --flatten` の flatten ループは停止しません — このループは残った虚振動方向へ変位してプラトーから脱出するためのものだからです。存在する場合は `stop_reason` にエネルギー範囲・ウィンドウ・満たせなかった基準が記録されます。stalled は（例えば摂動した構造やより厳しいステップ制御で）再試行し得るものであり、`max_cycles` 枯渇や一般的な失敗のエイリアスではありません。

## サブコマンド別スキーマ

### `sp`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"ok"` |
| `stage` | string | `"sp"` |
| `input` | string | 準備後の公開入力path |
| `backend` / `model` | string / string \| null | local MLIP provenance（共通の`mlip_*` fieldにもmirror） |
| `custom_calculator` | string \| null | `--calc-file`時の`filename:factory`。組み込みbackendではnull |
| `charge` / `spin` | int / int | 総電荷とspin多重度 |
| `n_atoms` | int | 原子数 |
| `energy_au` | float | 一点energy (Hartree) |
| `forces_path` | string | `forces.npy` のpath |
| `hessian_path` | string \| null | `hessian.npy` のpath。`--hess` 無指定時はnull |
| `elapsed` | string | 人間可読の経過時間 |

### `opt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"converged"` / `"not_converged"` / `"stalled"`（エネルギープラトー、上記参照） |
| `stop_reason` | string | 非収束停止（stalled/stopped）時のみ出力。エネルギープラトーの範囲・ウィンドウと満たせなかった基準を記録 |
| `energy_hartree` | float | 最終エネルギー (Hartree) |
| `n_opt_cycles` | int | 最適化サイクル数 |
| `opt_mode` | string | `"grad"` / `"hess"` / `"lbfgs"` / `"rfo"` |
| `backend` | string | MLIP バックエンド (`"uma"`, `"orb"`, `"mace"`, `"aimnet2"`、または `--calc-file` 時の `"custom"`) |
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
| `safeguards` | object | exact saddle check、最終目的modeのindex/overlap、停止理由、および明示的に有効化した場合のmode-loss/recovery診断。これらの回復経路はデフォルト無効 |
| `rigid_projection` | object | 剛体モードとexact Hessian のprovenance。[projection provenance](#rigid-projection-provenance)を参照 |

`files` には `imaginary_mode_files`（vib ファイルリスト）を含む場合があります。
Cartesian Hessian mode の `status: "converged"` には、最終 exact PHVA で
有意な虚振動がちょうど1個必要です。対応する内部座標 mode では、代わりに
exact optimizer-space Hessian の負の固有値が1個であることを要求します。
高次鞍点と `n_imag=0` 構造は `not_converged` です。エネルギープラトーによる
`stalled`（上記参照）も同様に `converged` として報告されず、以降の flatten/
再試行を停止します。収束詳細 (force/step) は
rsirfo モードで利用可能です。
dimer モードも `status` に `"converged"` / `"not_converged"` / `"stalled"` を
返しますが、`n_opt_cycles` のみを出力し、Hessian mode の収束詳細と
`safeguards` は省略されます。

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
| `rigid_projection` | object | 剛体モードと Hessian provenance。`--dump`時は`thermoanalysis.yaml`にも記録 |

**`thermochemistry`** (thermoanalysis 利用不可時は null):

| フィールド | 型 | 単位 |
|-----------|------|------|
| `point_group` | string | 自動判定した分子点群 |
| `point_group_source` | string | `"auto"` または保守的なフォールバックを示す `"auto-fallback"` |
| `symmetry_number` | int | 外部回転対称数 |
| `symmetry_number_source` | string | `"auto"` / `"auto-fallback"` / `"config"` / `"override"` |
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
| `forward_energy_increased` | bool \| null | 前方の最終stepで`irc.energy_increase_thresh`（デフォルト`0` Hartree、上昇はすべて対象）を超えてenergyが上昇したか |
| `backward_energy_increased` | bool \| null | 後方の最終stepで`irc.energy_increase_thresh`（デフォルト`0` Hartree、上昇はすべて対象）を超えてenergyが上昇したか |
| `backend` | string | MLIP バックエンド |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `model` | string | MLIP モデル名 |
| `never_stop` | bool | opt-inの物理的端点停止bypass modeを有効にしたか |
| `never_stop_energy_bypasses` | int | 実際にbypassしたenergy上昇／1 step energy変化量停止event数 |
| `n_freeze_atoms` | int | 凍結原子数 |
| `solvent` | string | 暗黙溶媒 or `"none"` |
| `bond_changes` | object | first→last 方向の `{formed: [...], broken: [...]}`。各リストは元素記号付き1始まりの原子ペア文字列（例 `"C7-O12"`）。比較が失敗または `finished_first.xyz`/`finished_last.xyz` が存在しない場合はキー自体が省略されます。 |
| `bond_changes_direction` | string | `bond_changes` がある場合は `"finished_first_to_finished_last"` |
| `step_length` | float | IRC ステップ長 (Bohr) |
| `max_cycles` | int | 最大 IRC ステップ数 |
| `input_file` | string | 入力ファイル名 |
| `files` | object | 軌跡ファイル（XYZと、利用可能なPDB/CIF companion） |
| `rigid_projection` | object | 剛体モードと初期 Hessian のprovenance。[projection provenance](#rigid-projection-provenance)を参照 |

### `scan`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"completed"` |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `backend` | string | MLIPバックエンド |
| `model` | string | MLIPモデル名 |
| `solvent` | string | 暗黙溶媒または`"none"` |
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
| `charge` | int \| null | 系の電荷。plot-onlyの`scan3d --csv`ではnull |
| `spin` | int \| null | スピン多重度。plot-onlyの`scan3d --csv`ではnull |
| `backend` | string \| null | MLIPバックエンド。plot-onlyの`scan3d --csv`ではnull |
| `model` | string \| null | MLIPモデル名。plot-onlyの`scan3d --csv`ではnull |
| `solvent` | string \| null | 暗黙溶媒または`"none"`。importしたenergyにはcalculator provenanceがないため、plot-onlyの`scan3d --csv`ではnull |
| `max_step_size_angstrom` | float | 1 ステップ当たりの最大結合長変位 (Å, `scan2d` のみ) |
| `n_grid_points` | int | `is_preopt=true` を除くグリッド行数 |
| `execution_status` | string | 実行レベルの完了状態 |
| `n_points_attempted` | int | fresh run の試行グリッド点数（事前最適化を除く） |
| `n_points_usable` | int | 科学計算上再利用可能な fresh-run 点数 |
| `point_outcomes` | object[] | 各点の収束・energy・artifact・eligibility |
| `grid_shape` | int[] | グリッド次元 (`scan3d --csv` 再プロット時には省略) |
| `pair1`, `pair2` (,`pair3`) | object | `{i, j, low, high}` (オプション: `label_i`, `label_j`)。`scan3d` で `--csv` 再プロット時は省略 |
| `min_energy_hartree` | float | 表面最小エネルギー |
| `files` | object | CSV + プロットファイル |

outcome count は fresh scan で出力します。plot-only `scan3d --csv` は
試行数を出力せず、入力 CSV の provenance が不完全なら usable count も
省略する場合があります。

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
| `mlip_precision` | string \| null | 実効`fp32` / `fp64` token。custom calculatorではnull |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |

`all` がさらに追加するフィールドは下の [`summary.json` セクション](#ja-summary-json-path-search-all) を参照してください。

### `dft`

> **注:** `dft` は SCF 収束・非収束の両方で `result.json` と
> `summary.json` を書きます。非収束時は `status: "not_converged"`、
> `converged: false` を記録して exit code 3 で終了します。未処理 exception
> では標準 error envelope を書きます。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"converged"` または `"not_converged"` |
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
| `ligand_charge_input` | string または null | ユーザー指定 `--ligand-charge` mapping。省略時は null |
| `center` | string | 中心残基 |
| `radius` | float | 抽出半径 (angstrom) |
| `input_files` | string[] | 入力 PDB / mmCIF パス |
| `files` | object | 出力 PDB / クラスターファイル |

### `trj2fig`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"ok"` |
| `n_frames` | int | 軌跡フレーム数 |
| `min_energy_hartree` | float | フレーム中の最小エネルギー |
| `max_energy_hartree` | float | フレーム中の最大エネルギー |
| `energy_source` | string | `"trajectory_comment"` または `"mlip_recomputed"` |
| `energy_provenance` | string[] | frame ごとの energy source provenance |
| `energy_unit` | string | 保存 energy unit（`hartree`） |
| `backend` | string または null | フレームを再計算した場合のみ MLIP backend。comment energy mode では null |
| `charge` / `multiplicity` | int または null | 再計算時の解決済み電子状態。それ以外は null |
| `solvent` / `solvent_model` | string または null | 再計算 calculator の溶媒設定。それ以外は null |
| `output_files` | string[] | すべての出力パスを順序どおりに保持する正規フィールド。別ディレクトリに同名ファイルがあっても保持される |
| `files` | object | 後方互換用のベース名からパスへの対応表。同じベース名が重複すると一方だけが残る |

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
| `status` | string | `"ok"`（全ペアが正常に比較できた場合）、`"partial"`（一部のペアが失敗した場合）、`"failed"`（正常に比較できたペアが無い場合） |
| `comparisons` | object[] | ペアごとの比較（`structure_a` (string), `structure_b` (string), `bonds_formed` (int), `bonds_broken` (int)） |

(ja-summary-json-path-search-all)=
## `summary.json` (`path-search` / `all`)

`all` / `path-search` は、より構造化された `summary.json` を出力します:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"success"` / `"partial"` / `"failed"` (`all`); `"success"` / `"partial"` (`path-search`) |
| `execution_status` / `scientific_status` | string / string | 実行の完了度と科学的な利用可能性。従来の `status` とは分けて評価します。 |
| `scientific_status_reasons` | string[] | 不完全または利用できない科学的結果の理由。正常終了時は省略されます。 |
| `expected_item_ids` / `observed_item_ids` | string[] | 期待された集約項目と観測された集約項目。 |
| `n_segments` | int | セグメント数 |
| `segments` | object[] | セグメントごとの `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, `bond_changes` |
| `energy_diagrams` | object[] | エネルギーダイアグラム（ラベル + kcal/mol） |
| `mlip_backend` | string | バックエンド名 |
| `mlip_model` | string \| null | バックエンドと分離して記録するモデル名 |
| `mlip_precision` | string \| null | 実効`fp32` / `fp64` token。custom calculatorではnull |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `environment` | object | ハードウェア情報 |
| `references` | object[] | 解決済みworkflowで実際に使った手法の `{method, citation, doi}` record。同じreference setを `summary.log` と最終標準出力の末尾（elapsed time直前）にまとめて出力します。 |

`all` はさらに以下を含みます:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `rate_limiting_step` | object | 互換性のため維持するキー。全 reactive segment を完全に覆う最高 method（`DFT//MLIP_Gibbs` > `DFT` > `MLIP_Gibbs` > `MLIP` > `MEP`）で、各段階の始状態を基準にした局所障壁の最大値。`method` と raw `mep_barrier_kcal` を明記する。microkinetics に基づく律速段階の判定ではない。 |
| `overall_reaction_energy_kcal` | float | 全体反応エネルギー |
| `overall_reaction_energy_method` | string | 全体反応energyのmethod (`MEP`, `MLIP`, `MLIP_Gibbs`, `DFT`, `DFT//MLIP_Gibbs`) |
| `post_segments` | list | セグメントごとの TS/IRC/freq/DFT 結果 |
| `post_segments[].thermo_symmetry` | object | 子 freq が報告した状態別の点群・回転対称 provenance。有効な対称数 provenance を持つ R/TS/P 状態だけを含み、欠けた状態は省略する。どの状態にも有効な provenance が無い場合だけフィールド全体を省略する。 |
| `current_output_paths` | string[] | `--out-dir` からの相対パスを並べたリスト。現在の呼び出しが記録した成果物だけを含みます。 |
| `key_output_files` | object | 現在の呼び出しの出力索引。ルートファイルはファイル名 → 説明、各 `seg_NN` は `{description, files}` で、`files` はそのセグメントディレクトリからの相対パスです。 |

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
