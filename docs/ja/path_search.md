# `path-search`


## 概要

> **要約:** **2 構造以上**から、GSM（デフォルト）または DMF（`--mep-mode dmf`）で連続的な MEP を構築します。共有結合変化のある領域のみを自動で精密化し、最高エネルギー画像（HEI）を TS 候補として出力します（freq/IRC で検証）。

### 要点
- **想定場面:** R → … → P のように **2 構造以上**を入力として、自動精密化を含めた連続 MEP を構築したい場合に使います。
- **手法:** GSM/DMF セグメントを連鎖し、結合変化が残る区間だけを再帰的に精密化します。
- **主な出力:** `mep_trj.xyz`（主軌跡）、`summary.yaml`（セグメントごとの結果）、必要に応じてプロットやマージ済み PDB。
- **既定値:** `--mep-mode gsm`、`--opt-mode light`（LBFGS）、`--preopt`、`--align`、`--thresh gau`。
- **次にやること:** HEI は **TS 候補**です。単独では TS 検証になりません。続けて [tsopt](tsopt.md) → [freq](freq.md) → [irc](irc.md) を実行してください。

`pdb2reaction path-search` は、反応順に並んだ 2 構造以上を入力として連続的な最小エネルギー経路（MEP）を構築します。共有結合変化が検出される領域のみを選択的に精密化し、解決済みのサブパスを連結して 1 本の軌跡にまとめます。


`--convert-files` が有効（デフォルト）な場合、参照 PDB があれば軌跡の `.pdb` コンパニオンを、Gaussian テンプレートがあれば HEI スナップショットの `.gjf` コンパニオンを生成します。XYZ/GJF 入力では `--ref-pdb` がポケット PDB トポロジーを提供し（XYZ 座標は保持）、`--ref-full-pdb` によりフルテンプレートへのマージが可能です（XYZ/GJF 入力では PDB コンパニオンは生成されません）。

**2 端点だけ**で再帰精密化が不要な場合は、[path-opt](path_opt.md) の方がシンプルです。

## 最小例

```bash
pdb2reaction path-search -i reactant.pdb product.pdb -q 0 -m 1 \
 --out-dir ./result_path_search
```

## 出力の見方

- `result_path_search/mep_trj.xyz`
- `result_path_search/summary.yaml`
- `result_path_search/summary.log`
- `result_path_search/mep_plot.png`（プロット生成時）

## よくある例

1. 中間体を明示して多段の経路を与える。

```bash
pdb2reaction path-search -i R.pdb IM1.pdb IM2.pdb P.pdb -q -1 -m 1 \
 --out-dir ./result_path_search_multi
```

2. テンプレート参照を使って全系マージ出力を有効化する。

```bash
pdb2reaction path-search -i R.pdb IM1.pdb P.pdb -q 0 -m 1 \
 --ref-full-pdb holo_template.pdb --out-dir ./result_path_search_merge
```

3. DMF + minima リファインで探索する。

```bash
pdb2reaction path-search -i reactant.pdb product.pdb -q 0 -m 1 \
 --mep-mode dmf --refine-mode minima --out-dir ./result_path_search_dmf
```

## 使用法

```bash
pdb2reaction path-search -i R.pdb [I.pdb...] P.pdb [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [--multiplicity 2S+1]
 [--workers N] [--workers-per-node N]
 [--mep-mode {gsm|dmf}] [--freeze-links/--no-freeze-links] [--thresh PRESET]
 [--refine-mode {peak|minima}]
 [--max-nodes N] [--max-cycles N] [--climb/--no-climb]
 [--opt-mode light|heavy] [--dump/--no-dump]
 [--out-dir DIR] [--preopt/--no-preopt]
 [--align/--no-align] [--ref-full-pdb FILE...] [--ref-pdb FILE...]
 [--convert-files/--no-convert-files]
 [--show-config/--no-show-config] [--dry-run/--no-dry-run]
```

### 例
- **ポケットのみ**の2つのエンドポイント間のMEP:
 ```bash
 pdb2reaction path-search -i reactant.pdb product.pdb -q 0
 ```
- YAML 上書きとマージされた全系出力を使用した**マルチステップ**探索:
 ```bash
 pdb2reaction path-search \
 -i R.pdb IM1.pdb IM2.pdb P.pdb -q -1 \
 --ref-full-pdb holo_template.pdb --out-dir ./run_ps
 ```


## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 反応順序の2つ以上の構造（反応物 → 生成物）。`-i` を繰り返すか1つのフラグに複数パスを渡す | 必須 |
| `-q, --charge INT` | 総電荷。非`.gjf`入力では `--ligand-charge` の導出が成功しない限り必須。両方指定時は `-q` が優先 | テンプレート/導出が適用されない限り必須 |
| `--ligand-charge TEXT` | `-q` 省略時に使用される総電荷または残基名ごとのマッピング。PDB 入力（または `--ref-pdb` 付き XYZ/GJF）で extract 方式の電荷導出を有効化 | _None_ |
| `--workers`, `--workers-per-node` | UMA予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| `--freeze-links/--no-freeze-links` | PDB ポケット読み込み時、リンク水素の親原子を凍結 | `True` |
| `--max-nodes INT` | MEPセグメントごとの内部ノード | `10` |
| `--max-cycles INT` | 最大MEP最適化サイクル（GSM/DMF） | `300` |
| `--climb/--no-climb` | GSMセグメントのクライミングイメージを有効化（ブリッジは無効） | `True` |
| `--opt-mode TEXT` | HEI±1/kinkノード用の単一構造オプティマイザー（`light`=LBFGS、`heavy`=RFO） | `light` |
| `--mep-mode {gsm\|dmf}` | セグメント生成器: GSM（string）またはDMF（direct flux） | `gsm` |
| `--refine-mode {peak\|minima}` | 精密化シード: `peak` はHEI±1、`minima` はHEIから最寄り局所極小へ外側探索。未指定時はGSMで`peak`、DMFで`minima` | _Auto_ |
| `--dump/--no-dump` | MEP（GSM/DMF）と単一構造軌跡/リスタートをダンプ | `False` |
| `--convert-files/--no-convert-files` | PDB/Gaussian入力のXYZ/TRJ → PDB/GJFコンパニオンを切り替え | `True` |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_path_search/` |
| `--thresh TEXT` | GSMおよびイメージごとの最適化の収束プリセットを上書き（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `gau` |
| `--config FILE` | 明示 CLI 指定より前に適用されるベース YAML | _None_ |
| `--show-config/--no-show-config` | 解決済み設定（YAML レイヤ情報を含む）を表示して実行継続 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画表示のみを行う | `False` |
| `--preopt/--no-preopt` | MEP 探索前に各エンドポイントを事前最適化（推奨） | `True` |
| `--align/--no-align` | 探索前にすべての入力を最初の構造にアライメント | `True` |
| `--ref-full-pdb PATH...` | フルサイズテンプレート PDB（`--align` があれば先頭のみ再利用可） | _None_ |
| `--ref-pdb PATH...` | 入力がXYZ/GJFの場合のポケット参照 PDB（XYZ 座標は保持） | _None_ |

## ワークフロー

1. **ペアごとの初期セグメント（GSM/DMF）** – 各隣接入力（A→B）間で `GrowingString` または DMF を実行し、粗いMEPと最高エネルギーイメージ（HEI）を取得。
2. **HEI周辺の局所緩和** – `refine-mode=peak` なら HEI±1、`refine-mode=minima` なら HEI 近傍の局所極小を、選択した単一構造オプティマイザー（`opt-mode`）で精密化し `End1`/`End2` を得る。
3. **kink vs. 精密化の決定** – `End1` と `End2` 間に共有結合変化がなければ *kink* とみなし、`search.kink_max_nodes` の線形ノードを挿入して個別最適化。結合変化がある場合は **精密化セグメント（GSM/DMF）** を起動。
4. **選択的再帰** – `(A→End1)` と `(End2→B)` の結合変化を `bond` しきい値で比較し、共有結合更新が残るサブ区間のみ再帰的に探索。再帰深度は `search.max_depth` で制限。
5. **スティッチング & ブリッジング** – 解決済みのサブパスを連結し、RMSD ≤ `search.stitch_rmsd_thresh` の重複エンドポイントを除去。RMSDギャップが `search.bridge_rmsd_thresh` を超える場合はブリッジMEPを挿入。境界で結合変化が検出される場合はブリッジではなく新規の再帰セグメントで置換。
6. **アライメント & マージング（オプション）** – `--align`（既定）で事前最適化構造を先頭入力へ剛体アライメントし、`freeze_atoms` を整合。`--ref-full-pdb` を指定するとポケット軌跡をフルサイズPDB テンプレートへマージ（`--align` により先頭テンプレートの再利用が可能）。

結合変化の判定は `bond_changes.compare_structures` を用い、`bond` セクションのしきい値に従います。UMA 計算機は全構造で共有され、効率的に再利用されます。

## 出力
```
out_dir/ (デフォルト:./result_path_search/)
├─ mep_trj.xyz # 主要 MEP 軌跡
├─ mep.pdb # 入力がPDB テンプレートで変換が有効な場合のPDB コンパニオン
├─ mep_w_ref.pdb # マージされた全系MEP（参照 PDB/テンプレートが必要）
├─ mep_w_ref_seg_XX.pdb # 共有結合変化がある場合のマージされたセグメントごとのパス
├─ summary.yaml # すべての再帰セグメントの障壁と分類サマリー
├─ summary.log # 人間が読めるサマリー
├─ mep_plot.png # ΔEプロファイル（kcal/mol、反応物基準）
├─ energy_diagram_MEP.png # MEP状態エネルギーダイアグラムの静的エクスポート
└─ seg_000_*/ # セグメントごとの GSM/DMF ダンプ、HEI スナップショット、kink/精密化の診断情報
```


- コンソールには確定済みの設定ブロック（`geom`, `calc`, `gs`, `opt`, `sopt.*`, `bond`, `search`）が出力されます。

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- 入力は2つ以上が必須。満たさない場合は `click.BadParameter` が発生します。
- `--ref-full-pdb` は1回の指定で複数ファイルを続けて渡せます。`--align` が有効な場合、マージでは先頭テンプレートのみが再利用されます。
- UMA 計算機は全構造で共有され、効率化されます。
- `--dump` が有効な場合、MEP（GSM/DMF）と単一構造最適化の軌跡が出力されます。リスタート YAML は YAML で `dump_restart` を有効にした場合のみ書き出されます。

マージ順は **defaults < config < 明示指定 CLI < override** です。

YAML ルートはマッピングでなければなりません。共通セクションは [YAML リファレンス](yaml_reference.md) を再利用します: `geom`/`calc` は単一構造設定を反映し（PDBでは `--freeze-links` が `geom.freeze_atoms` にマージ）、`opt` は `path-opt`（[path_opt.md](path_opt.md)）に記載の StringOptimizer 設定を継承します。

`gs`（Growing String）は `pdb2reaction.path_opt.GS_KW` の既定値を継承し、`max_nodes`（セグメント内部ノード）、クライミング設定（`climb`, `climb_rms`, `climb_fixed`）、再パラメータ化（`reparam_every_full`, `reparam_check`）を上書きできます。

`sopt` は HEI±1 と kink ノードに使う単一構造オプティマイザーで、`lbfgs` と `rfo` に分かれます。各サブセクションは [YAML リファレンス](yaml_reference.md) と同じキーを持ちますが、デフォルトは `out_dir:./result_path_search/`、`dump: False` です。

`bond` は UMA ベースの結合変化検出パラメータで、[scan](scan.md) の bond セクションと共通の `device`, `bond_factor`, `margin_fraction`, `delta_fraction` を持ちます。


`dmf` は `--mep-mode dmf` 選択時に適用される Direct Max Flux + (C)FB-ENM の設定です。既定値は `DMF_KW` を踏襲し、実行ごとに上書きできます。

### YAML例（デフォルト値）
```yaml
geom:
 coord_type: cart # coordinate type: cartesian vs dlc internals
 freeze_atoms: [] # 0-based frozen atoms merged with CLI/link detection
calc:
 charge: 0 # total charge (CLI/template override)
 spin: 1 # spin multiplicity 2S+1
 model: uma-s-1p1 # UMA model tag
 task_name: omol # UMA task name
 device: auto # UMA device selection
 max_neigh: null # maximum neighbors for graph construction
 radius: null # cutoff radius for neighbor search
 r_edges: false # store radial edges
 out_hess_torch: true # request torch-form Hessian
 freeze_atoms: null # calculator-level frozen atoms
 hessian_calc_mode: FiniteDifference # Hessian mode selection
 return_partial_hessian: false # full Hessian (avoids shape mismatches)
gs:
 fix_first: true # keep the first endpoint fixed during optimization
 fix_last: true # keep the last endpoint fixed during optimization
 max_nodes: 10 # maximum string nodes
 perp_thresh: 0.005 # perpendicular displacement threshold
 reparam_check: rms # reparametrization check metric
 reparam_every: 1 # reparametrization stride
 reparam_every_full: 1 # full reparametrization stride
 param: equi # parametrization scheme
 max_micro_cycles: 10 # micro-iteration limit
 reset_dlc: true # rebuild delocalized coordinates each step
 climb: true # enable climbing image
 climb_rms: 0.0005 # climbing RMS threshold
 climb_lanczos: true # Lanczos refinement for climbing
 climb_lanczos_rms: 0.0005 # Lanczos RMS threshold
 climb_fixed: false # keep climbing image fixed
 scheduler: null # optional scheduler backend
opt:
 type: string # optimizer type label
 stop_in_when_full: 300 # early stop threshold when string is full
 align: false # alignment toggle (kept off)
 scale_step: global # step scaling mode
 max_cycles: 300 # maximum optimization cycles
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 reparam_thresh: 0.0 # reparametrization threshold
 coord_diff_thresh: 0.0 # coordinate difference threshold
 out_dir:./result_path_search/ # output directory
 print_every: 10 # logging stride
dmf:
 max_cycles: 300 # DMF/IPOPT の最大反復数
 correlated: true # correlated DMF propagation
 sequential: true # sequential DMF execution
 fbenm_only_endpoints: false # run FB-ENM beyond endpoints
 fbenm_options:
 delta_scale: 0.2 # FB-ENM displacement scaling
 bond_scale: 1.25 # bond cutoff scaling
 fix_planes: true # enforce planar constraints
 two_hop_mode: sparse # neighbor traversal strategy
 cfbenm_options:
 bond_scale: 1.25 # CFB-ENM bond cutoff scaling
 corr0_scale: 1.1 # Correlation scale for corr0
 corr1_scale: 1.5 # Correlation scale for corr1
 corr2_scale: 1.6 # Correlation scale for corr2
 eps: 0.05 # Correlation epsilon
 pivotal: true # Pivotal residue handling
 single: true # Single-atom pivots
 remove_fourmembered: true # Prune four-membered rings
 two_hop_mode: sparse # Neighbor traversal strategy
 dmf_options:
 remove_rotation_and_translation: false # Keep rigid-body motions
 mass_weighted: false # Toggle mass weighting
 parallel: false # Enable parallel DMF
 eps_vel: 0.01 # Velocity tolerance
 eps_rot: 0.01 # Rotational tolerance
 beta: 10.0 # Beta parameter for DMF
 update_teval: false # Update transition evaluation
 k_fix: 300.0 # Harmonic constant for restraints
sopt:
 lbfgs:
 thresh: gau # LBFGS convergence preset
 max_cycles: 10000 # iteration limit
 print_every: 100 # logging stride
 min_step_norm: 1.0e-08 # minimum accepted step norm
 assert_min_step: true # assert when steps stagnate
 rms_force: null # explicit RMS force target
 rms_force_only: false # rely only on RMS force convergence
 max_force_only: false # rely only on max force convergence
 force_only: false # skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # RMS threshold when targeting geometry
 overachieve_factor: 0.0 # tighten thresholds
 check_eigval_structure: false # validate Hessian eigenstructure
 line_search: true # enable line search
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 prefix: "" # filename prefix
 out_dir:./result_path_search/ # output directory
 keep_last: 7 # history size for LBFGS buffers
 beta: 1.0 # initial damping beta
 gamma_mult: false # multiplicative gamma update toggle
 max_step: 0.3 # maximum step length
 control_step: true # control step length adaptively
 double_damp: true # double damping safeguard
 mu_reg: null # regularization strength
 max_mu_reg_adaptions: 10 # cap on mu adaptations
 rfo:
 thresh: gau # RFOptimizer convergence preset
 max_cycles: 10000 # iteration cap
 print_every: 100 # logging stride
 min_step_norm: 1.0e-08 # minimum accepted step norm
 assert_min_step: true # assert when steps stagnate
 rms_force: null # explicit RMS force target
 rms_force_only: false # rely only on RMS force convergence
 max_force_only: false # rely only on max force convergence
 force_only: false # skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # RMS threshold when targeting geometry
 overachieve_factor: 0.0 # tighten thresholds
 check_eigval_structure: false # validate Hessian eigenstructure
 line_search: true # enable line search
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 prefix: "" # filename prefix
 out_dir:./result_path_search/ # output directory
 trust_radius: 0.1 # trust-region radius
 trust_update: true # enable trust-region updates
 trust_min: 0.0 # minimum trust radius
 trust_max: 0.1 # maximum trust radius
 max_energy_incr: null # allowed energy increase per step
 hessian_update: bfgs # Hessian update scheme
 hessian_init: calc # Hessian initialization source
 hessian_recalc: 200 # rebuild Hessian every N steps
 hessian_recalc_adapt: null # adaptive Hessian rebuild factor
 small_eigval_thresh: 1.0e-08 # eigenvalue threshold for stability
 alpha0: 1.0 # initial micro step
 max_micro_cycles: 50 # micro-iteration limit
 rfo_overlaps: false # enable RFO overlaps
 gediis: false # enable GEDIIS
 gdiis: true # enable GDIIS
 gdiis_thresh: 0.0025 # GDIIS acceptance threshold
 gediis_thresh: 0.01 # GEDIIS acceptance threshold
 gdiis_test_direction: true # test descent direction before DIIS
 adapt_step_func: true # adaptive step scaling toggle
bond:
 device: cuda # UMA device for bond analysis
 bond_factor: 1.2 # covalent-radius scaling
 margin_fraction: 0.05 # tolerance margin for comparisons
 delta_fraction: 0.05 # minimum relative change to flag bonds
search:
 max_depth: 10 # recursion depth limit
 stitch_rmsd_thresh: 0.0001 # RMSD threshold for stitching segments
 bridge_rmsd_thresh: 0.0001 # RMSD threshold for bridging nodes
 max_nodes_segment: 10 # max nodes per segment
 max_nodes_bridge: 5 # max nodes per bridge
 kink_max_nodes: 3 # max nodes for kink optimizations
 max_seq_kink: 2 # max sequential kinks
 refine_mode: null # optional refinement strategy (auto-chooses when null)
```

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け

- [path-opt](path_opt.md) — 単一パスMEP最適化（再帰的精密化なし）
- [tsopt](tsopt.md) — HEIを遷移状態として最適化
- [extract](extract.md) — path-search入力用のポケットPDBを生成
- [all](all.md) — 内部でpath-searchを呼び出すエンドツーエンドワークフロー
- [YAML リファレンス](yaml_reference.md) — `gs`、`dmf`、`bond`、`search` の完全な設定オプション
- [用語集](glossary.md) — MEP、GSM、DMF、HEIの定義
