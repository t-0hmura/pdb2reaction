# JSON 出力の一覧

このページでは、`--out-json` で書き出す `result.json` と `summary.json` の欄（key）を、全コマンドに共通の欄とコマンドごとの欄に分けて示します。

## `--out-json` フラグ

`opt`、`sp`、`tsopt`、`freq`、`irc`、`scan`、`scan2d`、`scan3d`、`path-opt`、`dft`、`extract`、`trj2fig`、`energy-diagram` は `--out-json / --no-out-json`（デフォルト: 無効）に対応しています。有効にすると、通常の出力と同じ場所に `result.json` と `summary.json` を書き出します。2 つのファイルの中身は同じなので、`result.json` を読みます。

```bash
pdb2reaction opt -i r.pdb -q -1 --out-json --out-dir result_opt
cat result_opt/result.json | python -m json.tool
```

`result_opt/result.json` を開き、まず `execution_status` と `scientific_status` を読みます。`opt`、`tsopt`、`path-opt` は `optimization_status` も記録します。`all` と `path-search` が `--out-json` を付けなくても書き出す `summary.json` は、これとは[構造が違う別のファイル](#ja-summary-json-path-search-all)です。

## 共通の欄

どの結果ファイルにも下の欄があります。「任意」と書いた欄は、そのコマンドが対応するデータを持つときだけ出ます。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `schema_version` | string | ファイルのスキーマの版。版が上がると構造が変わったことを示します。 |
| `command` | string | 単体のコマンドはサブコマンド名（例: `"opt"`）、`all` / `path-search` の要約はコマンドライン全体を記録します |
| `pdb2reaction_version` | string | パッケージバージョン |
| `execution_status` | string | 実行の完了状況: `completed` / `failed`。 |
| `scientific_status` | string | 結果の利用可否: `success` / `partial` / `failed`。 |
| `run_id` | string | 任意。現在の呼び出しの UUID。MCP サーバーからコマンドを起動したときと、`all` の実行（各段を含む）で書かれます。 |
| `elapsed_seconds` | float | 任意。実行時間（秒）。時間を記録しないコマンドでは省きます |
| `environment` | object | ハードウェア情報（下表参照） |
| `mlip_backend` | string \| null | 任意。MLIP のバックエンドの識別子。描画だけのコマンドが calculator を評価しなかった場合は null |
| `mlip_model` | string \| null | 任意。バックエンドと分けて記録する、正確なモデル/チェックポイント名 |
| `mlip_model_label` | string \| null | 任意。正確な識別子から導出した論文表記用のモデル名 |
| `mlip_task` | string \| null | 任意。複数ドメインのモデルで使ったバックエンドのタスク。正確な識別子は `mlip_model` に保持 |
| `mlip_precision` | string \| null | 実効精度の共通表記（`fp32` / `fp64`）。dtype をユーザーのコードが管理する自作の calculator では null |

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

### 実行の完了と指定した段の完了

すべての結果に `execution_status` と `scientific_status` を出します。複数段階の計算と scan は、各段の結果も下の欄に残します。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `execution_status` | string | `completed` または `failed`。`all` では、IRC の後の端点の最適化が収束しなくても `completed` のままで、エラーで止まると `failed` になります。 |
| `scientific_status` | string | 指定したすべての段階が収束すれば `success`、そうでなければ `partial` か `failed`。`all` では、n_imag ≥ 2 の TS は `partial` になり、n_imag = 0 の TS では IRC の前で止まるので、`success` なら n_imag = 1 です。単体の `tsopt` は n_imag を見ないので、`n_imaginary_modes` を読んでください。 |
| `scientific_status_reasons` | string[] | 利用できない、または欠落した個別結果の理由。正常終了時は省略されます。 |
| `expected_item_ids` / `observed_item_ids` | string[] | 集約結果の欠落を検出するための、期待された項目と観測された項目の ID。 |
| `stage_outcomes` | object[] | `stage`、`item_id`、`required`、`executed`、`converged`、`usable`、`reason`、`artifacts` を持つ段階別の結果。 |
| `point_outcomes` | object[] | `point_id`、`executed`、`converged`、`energy_valid`、`artifact_written`、`seed_eligible`、`reason` を持つ scan の点ごとの結果。 |

### エラー時の欄（`execution_status == "failed"` のとき）

| フィールド | 型 | 説明 |
|-----------|------|------|
| `error` | string | 元の例外の `str(exc)` |
| `error_type` | string | 例外クラス名 |
| `error_class_chain` | list[string] | 例外のクラスとそのすべての親クラスの名前。エージェントはテキストを解析せずに階層を照合できます |
| `error_module` | string | 例外クラスが定義されたモジュール |
| `error_label` | string | CLI の段の名前 |

## エラー処理

出力ディレクトリを用意した後に例外で実行が止まった場合は、`--out-json` が無くても `"execution_status": "failed"` と `"error_type"` を含む `result.json` と `summary.json` を書き出します。それより前の失敗は[使用上の注意点](#使用上の注意点)を参照してください。

最適化が収束せずに終わった場合、`result.json` には `"optimization_status": "not_converged"` が記録されます。オプティマイザの結果には該当する最終の力/ステップやサイクルの欄が入り、DFT と Dimer は持たない欄を出しません。再実行の前に何を変えるかは、{ref}`トラブルシューティングの計算 / 収束の問題 <ja-calculation-convergence-problems>` を参照してください。

オプティマイザは `"optimization_status": "stalled"` を返すこともあります。これは、力/ステップの収束基準を満たさないまま、設定ウィンドウにわたってエネルギーが減少しなくなった状態（エネルギープラトー）です。停滞は未収束の一種で、`converged` にはなりません。止まった理由は `stop_reason` に記録されます。

## サブコマンド別スキーマ

### `sp`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `stage` | string | `"sp"` |
| `input` | string | 計算に使った入力のパス |
| `backend` / `model` | string / string \| null | MLIP のバックエンドとモデル。共通の `mlip_*` 欄にも同じ値が入ります |
| `custom_calculator` | string \| null | `--calc-file` 時の `filename:factory`。組み込みのバックエンドでは null |
| `charge` / `spin` | int / int | 総電荷とスピン多重度 |
| `n_atoms` | int | 原子数 |
| `energy_au` | float | 一点エネルギー (Hartree) |
| `forces_path` | string | `forces.npy` のパス |
| `hessian_path` | string \| null | `hessian.npy` のパス。`--hess` 無指定時は null |
| `elapsed` | string | 読みやすい形の経過時間 |

### `opt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `optimization_status` | string | `"converged"` / `"not_converged"` / `"stalled"`（エネルギープラトー、上記参照） |
| `stop_reason` | string | オプティマイザが収束せずに早く止まったとき（`stalled` か `not_converged`）のみ出力。止まった理由（エネルギープラトーの範囲・ウィンドウと満たせなかった基準など）を記録 |
| `energy_hartree` | float | 最終エネルギー (Hartree) |
| `n_opt_cycles` | int | 最適化サイクル数 |
| `opt_mode` | string | `"grad"` / `"hess"` / `"lbfgs"` / `"rfo"` |
| `backend` | string | calculator のバックエンド（`"uma"`, `"orb"`, `"mace"`, `"aimnet2"`, `"dft"`、または `--calc-file` 時の `"custom"`）。`dft` では `model` が `FUNCTIONAL/BASIS` で、エンジンは別の欄に記録し、MLIP の精度は null です |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `model` | string | MLIP モデル名。`dft` では `FUNCTIONAL/BASIS` |
| `n_atoms` | int | 原子数 |
| `n_freeze_atoms` | int | 固定原子数 |
| `solvent` | string | 陰溶媒または `"none"` |
| `thresh` | string | 収束閾値プリセット名 |
| `max_cycles` | int | 最大サイクル数 |
| `input_file` | string | 入力ファイル名 |
| `final_max_force` | float | 最終 max 勾配 (Hartree/Bohr) |
| `final_rms_force` | float | 最終 RMS 勾配 |
| `final_max_step` | float | 最終 max 変位 (Bohr) |
| `final_rms_step` | float | 最終 RMS 変位 |
| `convergence_thresholds` | object | `{max_force_thresh, rms_force_thresh, max_step_thresh, rms_step_thresh}`。力/ステップの単位は Cartesian で Hartree/Bohr と Bohr、角度成分で Hartree/rad と rad。 |
| `files` | object | 出力ファイルマップ |
| `rigid_projection` | object | 任意。`--flatten` 実行時に含む。[剛体モードの射影の記録](#剛体モードの射影の記録)を参照 |

### `tsopt`

`opt` と同じフィールドに加え:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `optimization_status` | string | 数値オプティマイザの結果: `"converged"` / `"not_converged"` / `"stalled"`。鞍点次数とは独立 |
| `saddle_validation` | string | 終端の厳密な PHVA（部分 Hessian 振動解析）による `"first_order"` / `"higher_order"` / `"no_imaginary"` / `"unavailable"` |
| `hessian_status` | string | `"completed"` / `"failed"` / `"skipped"` / `"unavailable"`。失敗理由は`hessian_error` |
| `reaction_mode_index` | int\|null | `all` が IRC でたどる虚振動の、PHVA の振動数の一覧での 0 始まりの番号。オプティマイザが追ったモードが虚振動ならそのモード、そうでなければ最も低い虚振動です。どちらかは `reaction_mode_source`（`"mep-reference-overlap"` / `"lowest-imaginary"`）に入ります。虚振動が無ければ `null` |
| `n_imaginary_modes` | int\|null | 虚振動の数。PHVA を実行しなかった場合は `null` |
| `n_negative_modes` | int\|null | 大きさを問わない負の振動数の数（閾値以内も含む）。`n_imaginary_modes` と並べて見る診断用の値で、PHVA を実行しなかった場合は `null` |
| `imaginary_frequencies_cm` | float[]\|null | 虚振動数 (cm⁻¹, 負の値)。PHVA 未実行時は `null` |
| `frequency_zero_cutoff_cm` | float | 閾値（cm⁻¹、デフォルト `5.0`）。デフォルトでは ν < −5.00 cm⁻¹ だけを虚振動として数えます |
| `imaginary_mode_criterion` | string | 数え方の規則の名前。値は `"frequency_cutoff_cm"` |
| `imaginary_frequency_threshold_cm` | float | 符号を負にした同じ閾値（デフォルト `-5.0`） |
| `opt_mode` | string | `"rsprfo"`（デフォルト） / `"rsirfo"` / `"trim"` / `"dimer"` |
| `opt_mode_requested` | string | CLI で要求したプリセット（`grad` / `hess` / 明示したアルゴリズム） |
| `optimizer` | string | 実際に使用したオプティマイザのアルゴリズム |
| `reference_mode_file` | string\|null | `--ref-mode` で渡した、経路から求めた反応モードのファイル（上級者向け）。通常は `all` が生成して内部指定します |
| `safeguards` | object | 厳密な鞍点の確認、最後に追ったモードの番号と重なり、止まった理由と、指定したときだけ働くモードの消失・回復の処理の記録。これらの回復の処理はデフォルトでは無効 |
| `rigid_projection` | object | 剛体モードと厳密な Hessian の記録。[剛体モードの射影の記録](#剛体モードの射影の記録)を参照 |

`files` には `imaginary_mode_files`（vib ファイルの一覧）と `hessian_npy`（`--dump-hess` のファイルの絶対パス）が入ることがあります。TS 最適化が成功すると、反応モードの虚振動が 1 つ出ます。このとき `saddle_validation: "first_order"`、`n_imaginary_modes: 1` です。最後の PHVA は、オプティマイザが収束したときと、エネルギープラトーで止まったとき（`stalled`）に実行します。それ以外はスキップとして記録します。`optimization_status` と `saddle_validation` は独立なので、収束した実行が `saddle_validation: "higher_order"` で終わることがあり、これは一次の TS ではありません。`all` が IRC に進む条件は [tsopt の TS の判定](tsopt.md#ts-の判定) を参照してください。Dimer の結果も同じ TS の欄を持ちますが、サイクルごとの力/ステップの詳細と `safeguards` はありません。

### `freq`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_modes` | int | 基準振動モードの総数 |
| `n_imaginary` | int | 虚振動の数（n_imag）。−`freq.zero_cutoff_cm` より小さい振動数の数 |
| `n_negative_modes` | int | 大きさを問わない負の振動数の数（閾値以内も含む） |
| `frequencies_cm` | float[] | 全振動数 (cm⁻¹) |
| `imaginary_frequencies_cm` | float[] | `n_imaginary` に数えた振動数 |
| `thermochemistry` | object\|null | 熱化学データ（下表参照） |
| `backend` | string | MLIP バックエンド |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `model` | string | MLIP モデル名 |
| `n_atoms` | int | 原子数 |
| `n_freeze_atoms` | int | 固定原子数 |
| `solvent` | string | 陰溶媒または `"none"` |
| `temperature_K` | float | 温度 (K) |
| `pressure_atm` | float | 圧力 (atm) |
| `input_file` | string | 入力ファイル名 |
| `files` | object | `{"frequencies_txt": "frequencies_cm-1.txt"}`。`--dump-hess` で書いたときは `hessian_npy`（絶対パス）を含む |
| `rigid_projection` | object | 剛体モードと Hessian の記録。`--dump` 時は `thermoanalysis.yaml` にも記録 |

**`thermochemistry`**:

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

IRC は `execution_status` と `scientific_status` に加え、方向ごとの停止理由と軌跡を残します。`all` は、その後の端点最適化を `endpoint_opt` に記録します。つないだ経路は `finished_first`（前方の端）から TS を通って `finished_last`（後方の端）へ並びます。どちらが R でどちらが P かは並び順では決まりません。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_frames_forward` | int | 前方 IRC フレーム数 |
| `forward_short_branch` / `backward_short_branch` | bool | サイクル上限前に 3 フレーム以内で停止した分岐。診断用のみ |
| `n_frames_backward` | int | 後方 IRC フレーム数 |
| `n_frames_total` | int | 全フレーム数 |
| `energy_first_hartree` | float | `finished_first` のエネルギー |
| `energy_ts_hartree` | float | TS エネルギー |
| `energy_last_hartree` | float | `finished_last` のエネルギー |
| `endpoint_energy_orientation` | string | `"finished_first_to_finished_last"` |
| `forward_requested` / `backward_requested` | bool | 各方向を要求したか |
| `forward_integration_converged` / `backward_integration_converged` | bool \| null | RMS 勾配の停止条件を満たして止まったか。診断専用で、`--never-stop` はこの判定を迂回するため常に `false`。`*_downhill_departure_valid` と合わせて見ると、その分岐が TS からエネルギーの下がる向きに離れ、かつこの判定を満たしたかを確かめられる |
| `forward_downhill_departure_valid` / `backward_downhill_departure_valid` | bool \| null | TS からエネルギーの下がる向きに離れたことを確認できたか |
| `forward_integration_stop_reason` / `backward_integration_stop_reason` | string \| null | 数値の積分が失敗したときだけ入る理由 |
| `forward_energy_increased` | bool \| null | 前方の最終ステップで `irc.energy_increase_thresh`（デフォルト `0` Hartree、上昇はすべて対象）を超えてエネルギーが上昇したか |
| `backward_energy_increased` | bool \| null | 後方の最終ステップで `irc.energy_increase_thresh`（デフォルト `0` Hartree、上昇はすべて対象）を超えてエネルギーが上昇したか |
| `backend` | string | MLIP バックエンド |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `model` | string | MLIP モデル名 |
| `never_stop` | bool | 端点での物理的な停止を迂回するモード（指定したときだけ有効）を使ったか |
| `never_stop_energy_bypasses` | int | エネルギーの上昇、または 1 ステップのエネルギー変化量による停止を、実際に迂回した回数 |
| `n_freeze_atoms` | int | 固定原子数 |
| `solvent` | string | 陰溶媒または `"none"` |
| `bond_changes` | object | first→last 方向の `{formed: [...], broken: [...]}`。各リストは元素記号付き 1 始まりの原子ペア文字列（例 `"C7-O12"`）。比較が失敗または `finished_first.xyz`/`finished_last.xyz` が存在しない場合はキー自体が省略されます。 |
| `bond_changes_direction` | string | `bond_changes` がある場合は `"finished_first_to_finished_last"` |
| `step_length` | float | IRC ステップ長 (Bohr) |
| `max_cycles` | int | 最大 IRC ステップ数 |
| `input_file` | string | 入力ファイル名 |
| `files` | object | 軌跡ファイル（XYZ と、あれば PDB/CIF 形式のファイル） |
| `rigid_projection` | object | 剛体モードと初期 Hessian の記録。[剛体モードの射影の記録](#剛体モードの射影の記録)を参照 |

### `scan`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `scan_opt_mode` | string | 拘束緩和に使用したオプティマイザのプリセット |
| `scan_optimizer` | string | 実効のオプティマイザの種類（`lbfgs` / `rfo`） |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `backend` | string | MLIP バックエンド |
| `model` | string | MLIP モデル名 |
| `solvent` | string | 陰溶媒または `"none"` |
| `preopt` | bool | 事前最適化を実行したか |
| `max_step_size_angstrom` | float | 1 ステップ当たりの最大結合長変位 (Å) |
| `n_stages` | int | スキャンステージ数 |
| `stages` | object[] | ステージごとのデータ（下記参照） |
| `files` | object | 出力ファイル |

**`stages[]`**:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `index` | int | 1 始まりのステージインデックス |
| `n_steps` | int | ステップ数 |
| `converged` | bool | 拘束最適化が収束したか |
| `pairs_1based` | list | 原子ペア (1 始まり) |
| `initial_distances_angstrom` | list | 初期距離 |
| `target_distances_angstrom` | list | 目標距離 |
| `final_energy_hartree` | float | 最終エネルギー |
| `energies_hartree` | float[] | ステップごとのエネルギー |
| `bond_changes` | object | `{"changed": bool \| null, "summary": str}`（自由記述の要約。比較が走らなかった場合は `null`/`""`）。 |

### `scan2d` / `scan3d`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `charge` | int \| null | 系の電荷。描画だけの `scan3d --csv` では null |
| `spin` | int \| null | スピン多重度。描画だけの `scan3d --csv` では null |
| `backend` | string \| null | MLIP バックエンド。描画だけの `scan3d --csv` では null |
| `model` | string \| null | MLIP モデル名。描画だけの `scan3d --csv` では null |
| `solvent` | string \| null | 陰溶媒または `"none"`。読み込んだエネルギーには calculator の記録が無いため、描画だけの `scan3d --csv` では null |
| `max_step_size_angstrom` | float | 1 ステップ当たりの最大結合長変位 (Å, `scan2d` のみ) |
| `n_grid_points` | int | `is_preopt=true` を除くグリッド行数 |
| `execution_status` | string | 実行レベルの完了状態 |
| `n_points_attempted` | int | 新規実行で試行したグリッド点数（事前最適化を除く） |
| `n_points_usable` | int | 新規実行の点のうち、科学計算上再利用可能な点数 |
| `point_outcomes` | object[] | 各点の収束・エネルギー・成果物・再利用可否 |
| `grid_points` | object[] | 各格子点のインデックス、距離、エネルギー、収束、`geometry_file` の対応 |
| `current_output_paths` | string[] | 今回の実行が書いた CSV/HTML/PNG と格子点の構造。前の実行で残ったファイルは入りません |
| `grid_shape` | int[] | グリッド次元 (`scan3d --csv` 再プロット時には省略) |
| `pair1`, `pair2` (,`pair3`) | object | `{i, j, low, high}` (オプション: `label_i`, `label_j`)。`scan3d` で `--csv` 再プロット時は省略 |
| `min_energy_hartree` | float | エネルギー曲面上の最小エネルギー |
| `files` | object | CSV + プロットファイル |

試行した点と使える点の数の欄は、新しく scan したときに出力します。描画だけの `scan3d --csv` は試行数を出力せず、入力 CSV の記録が不完全なら使える点の数も省略する場合があります。

### `path-opt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `optimization_status` | string | `"converged"` / `"not_converged"` / `"completed"` |
| `converged` | bool \| null | 収束判定: エンジン自身の収束シグナルによる `true` / `false`。読み取れない場合は `null`（`optimization_status` は `"completed"` となり、収束を主張しない） |
| `mep_mode` | string | `"dmf"` / `"gsm"` |
| `backend` | string | MLIP バックエンド |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `model` | string | MLIP モデル名 |
| `solvent` | string | 陰溶媒または `"none"` |
| `preopt` | bool | 端点の事前最適化を有効にしたか |
| `reactant_energy_hartree` | float | 最初のイメージのエネルギー (Hartree) |
| `product_energy_hartree` | float | 最後のイメージのエネルギー (Hartree) |
| `image_energies_hartree` | float[] | 全イメージエネルギー |
| `n_images` | int | イメージ数 |
| `hei_index` | int | 最高エネルギーイメージのインデックス |
| `hei_energy_hartree` | float | HEI エネルギー |
| `barrier_kcal` | float | 前方障壁 (kcal/mol) |
| `delta_kcal` | float | 反応エネルギー (kcal/mol) |
| `files` | object | 軌跡 + HEI ファイル |

### `path-search`

`path-search` には `--out-json` フラグがありません。`summary.json` を書き出し、その欄は [`summary.json` (`path-search` / `all`)](#ja-summary-json-path-search-all) に示します。

### `dft`

> **注:** `--out-json` 指定時、`dft` は SCF 収束・非収束の両方で `result.json` と `summary.json` を書き、非収束時は `scientific_status: "failed"`、`converged: false` を記録します。SCF が収束しなかった実行は終了コード 1 で終わります。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `converged` | bool | SCF が収束したか |
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
| `used_gpu` | bool | GPU を使ったか |
| `used_lowmem` | bool | GPU4PySCF の低メモリのソルバーを実際に使ったか（開殻、CPU、`--no-dft-low-memory` では False） |
| `lowmem_requested` | bool | 低メモリモードを要求したか |
| `dft_settings` / `dft_resources` | object | 正規化した計算の設定と、実際に使ったホストの計算資源 |
| `effective_ecp` | string/object \| null | PySCF へ渡した実効 ECP |
| `solvent` / `solvent_model` | string | 実効の組み込み陰溶媒の設定 |
| `charges` | object | `{mulliken, lowdin, iao}` 原子電荷配列 |
| `spin_densities` | object | `{mulliken, lowdin, iao}` スピン密度配列 |
| `files` | object | `{"result_yaml": "result.yaml", "input_geometry_xyz": "input_geometry.xyz"}` |

### `extract`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_atoms_raw` | int | 選択残基に含まれる、主鎖の除外と切断の前の原子数（入力全体ではない） |
| `n_atoms_extracted` | int | 切断後に保持した原子数（キャップ水素の追加前） |
| `total_charge` | float | 合計電荷 |
| `protein_charge` | float | タンパク質電荷 |
| `ligand_total_charge` | float | リガンド電荷合計 |
| `ion_total_charge` | float | イオン電荷合計 |
| `ion_charges` | list | `[[名前, 電荷], ...]` |
| `unknown_residue_charges` | object | `{残基名: 電荷}` |
| `n_link_hydrogens` | int | 炭素原子側に残る切断結合へ追加されたキャップ水素数。モデルの原子数は `n_atoms_extracted` + `n_link_hydrogens` です |
| `exclude_backbone` | bool | 主鎖を除外したか |
| `include_h2o` | bool | 結晶水を含めたか |
| `ligand_charge_input` | string \| null | ユーザーが指定した `--ligand-charge` の対応表。省略時は null |
| `center` | string | 中心残基 |
| `radius` | float | 抽出半径 (angstrom) |
| `input_files` | string[] | 入力 PDB / mmCIF パス |
| `files` | object | 出力 PDB / クラスターファイル |

### `trj2fig`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_frames` | int | 軌跡フレーム数 |
| `min_energy_hartree` | float | フレーム中の最小エネルギー |
| `max_energy_hartree` | float | フレーム中の最大エネルギー |
| `energy_source` | string | `"trajectory_comment"` または `"mlip_recomputed"` |
| `mlip_backend` / `mlip_model` / `mlip_model_label` / `mlip_task` / `mlip_precision` | string \| null | 再計算の設定。軌跡コメントを読む場合はすべて null |
| `energy_provenance` | string[] | フレームごとのエネルギーの出どころ |
| `energy_unit` | string | 保存したエネルギーの単位（`hartree`） |
| `backend` | string \| null | フレームを再計算した場合のみ MLIP のバックエンド。軌跡コメントのエネルギーを読む場合は null |
| `charge` / `multiplicity` | int \| null | エネルギーを再計算したときの電荷と多重度。それ以外は null |
| `solvent` / `solvent_model` | string \| null | 再計算 calculator の溶媒設定。それ以外は null |
| `output_files` | string[] | すべての出力パスを順序どおりに保持する正規フィールド。別ディレクトリに同名ファイルがあっても保持される |
| `files` | object | ベース名からパスへの対応表。同じベース名の出力が 2 つあると一方だけが残るため、`output_files` を使ってください |

### `energy-diagram`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_points` | int | エネルギーデータ点数 |
| `files` | object | 出力ダイアグラムファイル |

### `bond-summary`

`--json` 有効時、`bond-summary` は JSON を**標準出力**に出し、上のサブコマンドと違って `result.json` を書き出しません:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `execution_status` / `scientific_status` | string / string | すべての組を比較できれば `completed` / `success`。比較できない組があると `execution_status` は `failed`、`scientific_status` は `partial` か `failed` |
| `comparisons` | object[] | ペアごとの比較（`structure_a` (string), `structure_b` (string), `bonds_formed` (int), `bonds_broken` (int)） |

### 剛体モードの射影の記録

`freq`、`irc`、`tsopt` の結果は `rigid_projection` object を含み、`opt` では `--flatten` 実行時に含みます。`freq --dump` は同じ object を `thermoanalysis.yaml` にも書きます。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `treatment` | string | 固定の剛体モード処理: `"constrained"` |
| `algorithm` | string | 射影の方法の名前 |
| `effective_rank` | int | 動ける原子の Hessian から除いた剛体方向の数 |
| `full_rigid_rank` | int | 固定原子を考える前の、系全体の剛体運動のランク |
| `frozen_constraint_rank` | int | 固定原子を動かさない条件で除かれたランク |
| `svd_rtol` | float | ランクの判定に用いる相対 SVD 許容値 |
| `active_atom_count` / `frozen_atom_count` | int | 動ける原子と固定原子の数 |
| `active_atoms` / `frozen_atoms` | int[] | 動ける原子と固定原子の 0 始まりのインデックス |
| `hessian_space` | string | 入力 Hessian 空間: `"full"` / `"active"` |
| `hessian_source` / `source` | string | Hessian の出どころ。`freq`/`irc` は `hessian_source` で、`"file"`（`--read-hess`）、`"cache"`（同じ実行の前の段）、`"fresh"`（新規計算）のいずれか。`opt`/`tsopt` は `source` を使用 |
| `hessian_shape` / `raw_hessian_shape` | int[2] | 入力 Hessian の形状。`freq`/`irc` は `hessian_shape`、`opt`/`tsopt` は `raw_hessian_shape` を使用 |
| `near_zero_mode_count` / `near_zero_frequencies_cm` | int / float[] | ±`frequency_zero_cutoff_cm`（デフォルト 5.00 cm⁻¹）以内のモードの数と値。これらのモードは全振動数の一覧にも入ります |

`constrained` は、固定原子を動かさない系全体の剛体運動だけを除きます。詳しくは [freq](freq.md#固定境界での剛体モード) を参照してください。

(ja-summary-json-path-search-all)=
## `summary.json` (`path-search` / `all`)

`all` / `path-search` は、より構造化された `summary.json` を出力します:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `execution_status` / `scientific_status` | string / string | 実行が完了したかと、指定した計算の段がすべて完了したか。 |
| `scientific_status_reasons` | string[] | 要求した結果の欠損・未収束などの理由。正常終了時は省略されます。 |
| `pipeline_stop` | object \| なし | 早期停止時のみ存在。`stage` は `post`（`reason` は `no_segments` / `no_reactive_segment`）、`before_irc`（TSOPT の理由と `segment`・`tsopt_result`）、または `endpoint_opt`（`endpoint_execution_failed` と端点別 `failures`）。`summary.log` では `Pipeline stop` |
| `expected_item_ids` / `observed_item_ids` | string[] | 期待された集約項目と観測された集約項目。 |
| `config` | object | 実効設定。`mep_mode` は GSM/DMF、`ts_opt_mode` / `endpoint_opt_mode` は設定した後処理のプリセットを示す。汎用の `opt_mode*` は実効の CLI の値を記録する。`path_opt_mode` は端点の事前最適化に使う単一構造オプティマイザであり（`preopt` を参照）、MEP の経路アルゴリズムではない。 |
| `scan` | object \| なし | scanを初期経路に使う`all`の、前段scanの状態・各段階の結果・診断。 |
| `n_segments` | int | セグメント数 |
| `search_max_depth` | int | 実効の再帰分割階層上限。`0` は分割無効 |
| `path_optimizers` | string[] | 経路の準備・精密化で実際に使用した単一構造オプティマイザ（`lbfgs`, `rfo`）。`all` ではスキャン・アライメントの実行も含む。`path-opt` の `result.json` にも記録 |
| `preopt_requested` / `preopt_converged` | bool / bool \| null | 端点事前最適化を実行したか、および全端点が収束したか。読み取れない端点があれば `null`。`all` では `preopt_converged` が `scientific_status` の判定に入ります。`--tsopt` で、すべての反応セグメントの TS と両端点の最適化が収束すると、判定から外れます |
| `segments` | object[] | セグメントごとの `index`（1 始まり。`all` はセグメント 1 を `segments/seg_01/` に書きます）、`tag`（名前。番号は `index` と一致するとは限りません）、`kind`（反応セグメントは `"seg"`、共有結合が変わらない[キンク（kink）](path-search.md#処理の仕組みと計算仕様)は `"kink"`、[ブリッジセグメント](path-search.md#処理の仕組みと計算仕様)は `"bridge"`、TS-only モードは `"tsopt"`）、`converged`（そのセグメントを作った最適化がすべて収束したか。収束の信号を読めなければ `null`）、`barrier_kcal`、`delta_kcal`、`bond_changes`（`{title: [entries]}` dict のリスト。ブリッジセグメントは `""`）。`barrier_kcal` は TS 最適化の前の MEP 上の障壁です（キンクでは `null`）。TS-only モードでは TS − R で、R は IRC の両端のうちエネルギーが高いほうです（向きの名前で、化学的な反応の向きではありません。[all の使用上の注意点](all.md#使用上の注意点)を参照）。 |
| `energy_diagrams` | object[] | エネルギーダイアグラム（ラベル + kcal/mol） |
| `mlip_backend` | string | バックエンド名 |
| `mlip_model` | string \| null | バックエンドと分離して記録するモデル名 |
| `mlip_model_label` | string \| null | 論文表記用のモデル名 |
| `mlip_task` | string \| null | 複数ドメインのモデルで使ったバックエンドのタスク |
| `mlip_precision` | string \| null | 実効の `fp32` / `fp64` 表記。自作の calculator では null |
| `charge` | int | 系の電荷 |
| `spin` | int | スピン多重度 |
| `environment` | object | ハードウェア情報 |
| `references` | object[] | 実行で実際に使った手法の `{method, citation, doi}` の記録。同じ文献の一覧を `summary.log` と最後の標準出力の末尾（経過時間の直前）にまとめて出力します。 |

`all` はさらに以下を含みます:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `rate_limiting_step` | object | 反応のあるすべてのセグメントで使える最上位の手法（`DFT//MLIP_Gibbs` > `DFT` > `MLIP_Gibbs` > `MLIP` > `MEP`）で求めた局所障壁の最大値。使った手法 `method` と、MEP から求めた障壁 `mep_barrier_kcal` も含みます。微視的速度論による律速段階の判定ではありません。 |
| `overall_reaction_energy_kcal` | float | 全体反応エネルギー |
| `overall_reaction_energy_method` | string | 全体反応エネルギーを求めた手法（`MEP`、`MLIP`、`MLIP_Gibbs`、`DFT`、`DFT//MLIP_Gibbs`） |
| `post_segments` | list | セグメントごとの TS/IRC/freq/DFT 結果 |
| `post_segments[].tsopt.n_imaginary_modes` / `.imaginary_frequencies_cm` | int / float[] | 最適化した TS の n_imag と虚振動数 (cm⁻¹, 負の値) |
| `post_segments[].mlip` / `.gibbs_mlip` / `.dft` / `.gibbs_dft_mlip` | object | 1 つのレベルでの R・TS・P のエネルギーと、`barrier_kcal`・`delta_kcal`・`energies_kcal`。順に MLIP の電子エネルギー（`--tsopt`）、MLIP の Gibbs エネルギー（`--thermo`）、DFT のエネルギー（`--dft`）、DFT//MLIP の Gibbs エネルギー（`--thermo` と `--dft`）です |
| `post_segments[].tsopt.energy_valid` / `.structure_valid` | bool | `energy_valid`：最終の TS のエネルギーが有限の値。`structure_valid`：最終の TS の構造ファイルがあり、座標が有限の値 |
| `post_segments[].tsopt.n_opt_cycles` / `.max_cycles` | int / int\|null | TS 最適化で実行したサイクル数と設定上限。通常の非収束時にも記録します。 |
| `post_segments[].irc` / `.endpoint_assignment` / `.endpoint_opt` | object | 順に IRC 停止診断、端点の向き付け、端点 OPT の収束記録。`endpoint_opt.reactant` / `.product` に `optimization_status`, `n_opt_cycles`, `max_cycles`, `stop_reason`（存在する場合）を記録します。IRC がどう止まったかと結合変化が合うかは `scientific_status` に入りません。最適化した端点が狙った R と P かは自分で確かめてください。 |
| `post_segments[].thermo_symmetry` | object | freq が R・TS・P ごとに報告した点群と回転対称数とその出どころ。有効な対称数の出どころを持つ R/TS/P 状態だけを含み、欠けた状態は省略する。どの状態にも有効な出どころが無い場合だけフィールド全体を省略する。 |
| `current_output_paths` | string[] | `--out-dir` からの相対パスを並べたリスト。現在の呼び出しが記録した成果物だけを含みます。 |
| `key_output_files` | object | 現在の呼び出しの出力索引。ルートファイルはファイル名 → 説明、各 `seg_NN` は `{description, files}` で、`files` はそのセグメントディレクトリからの相対パスです。 |

## 使用例

### Python

`opt --out-json` が書き出す `result_opt/result.json` を読む例です。

```python
import json

with open("result_opt/result.json") as f:
    result = json.load(f)

status = result.get("optimization_status")
if result["execution_status"] == "failed":
    raise RuntimeError(f"{result['error_type']}: {result['error']}")
elif status == "converged":
    print(f"Energy: {result['energy_hartree']:.6f} Hartree")
elif status in {"not_converged", "stalled"}:
    print(f"Not converged after {result['n_opt_cycles']} cycles")
    print(f"Max force: {result['final_max_force']:.6f}")
else:
    print(f"Status: {status}")
```

### jq

```bash
# 収束確認
jq '{execution_status, scientific_status}' result.json

# path-opt の障壁
jq '.barrier_kcal' result.json

# tsopt の虚振動数
jq '.imaginary_frequencies_cm' result.json

# freq の自由エネルギー
jq '.thermochemistry.sum_EE_and_thermal_free_energy_ha' result.json

# all の各セグメントの障壁（--tsopt の後）
jq '.post_segments[] | {index, barrier_kcal: .mlip.barrier_kcal}' result_all/summary.json
```

## 使用上の注意点

- CLI のオプションや入力を確かめている段階（出力ディレクトリを用意する前）で失敗すると、JSON を書かずに止まることがあります。0 以外の終了コードは失敗として扱い、標準エラー出力かジョブログでメッセージを確認してください。
- `all` と `path-search` は要約の段に着いてから `summary.json` を書くため、早い段階の入力エラーではファイルが残りません。
- コマンドごとの節にも、それぞれの `backend` / `model` 欄があります。複数のコマンドの結果を読むスクリプトでは、どのコマンドでも同じ意味の `mlip_backend` / `mlip_model` / `mlip_precision` を使ってください。
- [IRC がどう止まったか](irc.md#irc-の成否の判定)は `scientific_status` に入りません。単体の `irc` は、積分が終われば `--max-cycles` に達した場合も `success` で、`all` は TS と端点の最適化から判定します。

## 関連ドキュメント

- {ref}`終了コード <ja-exit-codes>` — 終了コードの意味
- [出力ディレクトリのレイアウト](output-layout.md) — `result.json` と `summary.json` を書き出す場所
- [トラブルシューティング](troubleshooting.md) — 失敗・未収束の後に何を変えるか
- [YAML 設定の一覧](yaml-reference.md) — これらのスキーマに現れる設定入力
- [all](all.md), [path-search](path-search.md) — `--out-json` なしで `summary.json` を書き出すサブコマンド
- [opt](opt.md), [sp](sp.md), [tsopt](tsopt.md), [freq](freq.md), [irc](irc.md), [scan](scan.md), [scan2d](scan2d.md), [scan3d](scan3d.md), [path-opt](path-opt.md), [dft](dft.md), [extract](extract.md), [trj2fig](trj2fig.md), [energy-diagram](energy-diagram.md) — `--out-json` 指定時にのみ `result.json` を書き出すサブコマンド
