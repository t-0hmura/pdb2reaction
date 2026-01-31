# `path-search`


## 概要
反応座標に沿って順序付けられた**2つ以上**の構造にわたって連続的な最小エネルギー経路（MEP）を構築します。`path-search` はGSM**またはDMF**セグメントを連鎖させ、共有結合変化のある領域のみを選択的に精密化し、（オプションで）PDBポケットをフルサイズテンプレートにマージします。`--mep-mode` でどちらを選んでも同じ再帰ワークフローが動作し、**GSMがデフォルト**です。デフォルトの `--opt-mode` は **light**（LBFGS）です。RFOを使用する場合は `--opt-mode heavy` を指定してください。フォーマット対応の変換は、PDB参照がある場合に軌跡を `.pdb` に、Gaussianテンプレートがある場合にXYZスナップショット（例: HEI）を `.gjf` にミラーします（`--convert-files` 有効時、デフォルト）。XYZ/GJF入力では `--ref-pdb` がポケットPDBトポロジーを提供しXYZ座標を保持するため、`--ref-full-pdb` が与えられればフルテンプレートマージが可能になります（XYZ/GJF入力ではPDBコンパニオン自体は生成されません）。

## 使用法

```bash
pdb2reaction path-search -i R.pdb [I.pdb ...] P.pdb [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [--multiplicity 2S+1]
                         [--workers N] [--workers-per-node N]
                         [--mep-mode {gsm|dmf}] [--freeze-links {True\|False}] [--thresh PRESET]
                         [--refine-mode {peak|minima}]
                         [--max-nodes N] [--max-cycles N] [--climb {True\|False}]
                         [--opt-mode light|heavy] [--dump {True\|False}]
                         [--out-dir DIR] [--preopt {True\|False}]
                         [--align {True\|False}] [--ref-full-pdb FILE ...] [--ref-pdb FILE ...]
                         [--convert-files {True\|False}]
                         [--args-yaml FILE]
```

### 例
- **ポケットのみ**の2つのエンドポイント間のMEP:
  ```bash
  pdb2reaction path-search -i reactant.pdb product.pdb -q 0
  ```
- YAMLオーバーライドとマージされた全系出力を使用した**マルチステップ**探索:
  ```bash
  pdb2reaction path-search \
      -i R.pdb IM1.pdb IM2.pdb P.pdb -q -1 \
      --args-yaml params.yaml --ref-full-pdb holo_template.pdb --out-dir ./run_ps
  ```


## CLIオプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 反応順序の2つ以上の構造（反応物 → 生成物）。`-i` を繰り返すか1つのフラグに複数パスを渡す | 必須 |
| `-q, --charge INT` | 総電荷。非`.gjf`入力では `--ligand-charge` の導出が成功しない限り必須。両方指定時は `-q` が優先 | テンプレート/導出が適用されない限り必須 |
| `--ligand-charge TEXT` | `-q` が省略された場合に使用する総電荷または残基名マッピング。PDB入力ではextract方式の電荷導出を有効化 | _None_ |
| `--workers`, `--workers-per-node` | UMA予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| `--freeze-links {True\|False}` | PDBポケット読み込み時、リンク水素の親原子を凍結 | `True` |
| `--max-nodes INT` | MEPセグメントごとの内部ノード | `10` |
| `--max-cycles INT` | 最大MEP最適化サイクル（GSM/DMF） | `300` |
| `--climb {True\|False}` | GSMセグメントのクライミングイメージを有効化（ブリッジは無効） | `True` |
| `--opt-mode TEXT` | HEI±1/kinkノード用の単一構造オプティマイザー（`light`=LBFGS、`heavy`=RFO） | `light` |
| `--mep-mode {gsm\|dmf}` | セグメント生成器: GSM（string）またはDMF（direct flux） | `gsm` |
| `--refine-mode {peak\|minima}` | 精密化シード: `peak` はHEI±1、`minima` はHEIから最寄り局所極小へ外側探索。未指定時はGSMで`peak`、DMFで`minima` | _Auto_ |
| `--dump {True\|False}` | MEP（GSM/DMF）と単一構造軌跡/リスタートをダンプ | `False` |
| `--convert-files {True\|False}` | PDB/Gaussian入力のXYZ/TRJ → PDB/GJFコンパニオンを切り替え | `True` |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_path_search/` |
| `--thresh TEXT` | GSMおよびイメージごとの最適化の収束プリセットを上書き（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `gau` |
| `--args-yaml FILE` | YAMLオーバーライド（下記参照） | _None_ |
| `--preopt {True\|False}` | MEP探索前に各エンドポイントを事前最適化（推奨） | `True` |
| `--align {True\|False}` | 探索前にすべての入力を最初の構造にアライメント | `True` |
| `--ref-full-pdb PATH...` | フルサイズテンプレートPDB（`--align` があれば先頭のみ再利用可） | _None_ |
| `--ref-pdb PATH...` | 入力がXYZ/GJFの場合のポケット参照PDB（XYZ座標は保持） | _None_ |

## ワークフロー

1. **ペアごとの初期セグメント（GSM/DMF）** – 各隣接入力（A→B）間で `GrowingString` または DMF を実行し、粗いMEPと最高エネルギーイメージ（HEI）を取得。
2. **HEI周辺の局所緩和** – `refine-mode=peak` なら HEI±1、`refine-mode=minima` なら HEI 近傍の局所極小を、選択した単一構造オプティマイザー（`opt-mode`）で精密化し `End1`/`End2` を得る。
3. **kink vs. 精密化の決定** – `End1` と `End2` 間に共有結合変化がなければ *kink* とみなし、`search.kink_max_nodes` の線形ノードを挿入して個別最適化。結合変化がある場合は **精密化セグメント（GSM/DMF）** を起動。
4. **選択的再帰** – `(A→End1)` と `(End2→B)` の結合変化を `bond` しきい値で比較し、共有結合更新が残るサブ区間のみ再帰的に探索。再帰深度は `search.max_depth` で制限。
5. **スティッチング & ブリッジング** – 解決済みのサブパスを連結し、RMSD ≤ `search.stitch_rmsd_thresh` の重複エンドポイントを除去。RMSDギャップが `search.bridge_rmsd_thresh` を超える場合はブリッジMEPを挿入。境界で結合変化が検出される場合はブリッジではなく新規の再帰セグメントで置換。
6. **アライメント & マージング（オプション）** – `--align`（既定）で事前最適化構造を先頭入力へ剛体アライメントし、`freeze_atoms` を整合。`--ref-full-pdb` を指定するとポケット軌跡をフルサイズPDBテンプレートへマージ（`--align` により先頭テンプレートの再利用が可能）。

結合変化の判定は `bond_changes.compare_structures` を用い、`bond` セクションのしきい値に従います。UMA計算機は全構造で共有され、効率的に再利用されます。

## 出力
```
out_dir/ (デフォルト: ./result_path_search/)
├─ mep.trj                  # プライマリMEP軌跡
├─ mep.pdb                  # 入力がPDBテンプレートで変換が有効な場合のPDBコンパニオン
├─ mep_w_ref.pdb            # マージされた全系MEP（参照PDB/テンプレートが必要）
├─ mep_w_ref_seg_XX.pdb     # 共有結合変化がある場合のマージされたセグメントごとのパス
├─ summary.yaml             # すべての再帰セグメントの障壁と分類サマリー
├─ mep_plot.png             # ΔEプロファイル（kcal/mol、反応物基準）
├─ energy_diagram_MEP.png   # MEP状態エネルギーダイアグラムの静的エクスポート
└─ seg_000_*/               # セグメントごとのGSM/DMFダンプ、HEIスナップショット
```


- 設定解決済みの `geom`, `calc`, `gs`, `opt`, `sopt.*`, `bond`, `search` ブロックはコンソールに出力されます。

## 注意事項
- 入力は2つ以上が必須。満たさない場合は `click.BadParameter` が発生します。
- `--ref-full-pdb` は1回の指定で複数ファイルを続けて渡せます。`--align` が有効な場合、マージでは先頭テンプレートのみが再利用されます。
- UMA計算機は全構造で共有され、効率化されます。
- `--dump` が有効な場合、MEP（GSM/DMF）と単一構造最適化の軌跡および再開用YAMLが出力されます。
- 電荷/スピンは `.gjf` テンプレートがあればそれを継承します。`-q` が省略され `--ligand-charge` が与えられている場合、入力は酵素–基質複合体として扱われ、PDB入力では `extract.py` の電荷サマリーで総電荷が導出されます。明示的な `-q` は常に優先されます。`.gjf` 以外で `--ligand-charge` が使えない場合は実行が中断され、多重度は省略時に `1` がデフォルトです。

## YAML設定（`--args-yaml`）
YAMLルートはマッピングである必要があります。YAML値はCLIを上書きします。共通セクションは [YAMLリファレンス](yaml-reference.md) を再利用します: `geom`/`calc` は単一構造設定を反映し（PDBでは `--freeze-links` が `geom.freeze_atoms` にマージ）、`opt` は `path_opt` に記載のStringOptimizer設定を継承します。

`gs`（Growing String）は `pdb2reaction.path_opt.GS_KW` の既定値を継承し、`max_nodes`（セグメント内部ノード）、クライミング設定（`climb`, `climb_rms`, `climb_fixed`）、再パラメータ化（`reparam_every_full`, `reparam_check`）を上書きできます。

`sopt` は HEI±1 と kink ノードに使う単一構造オプティマイザーで、`lbfgs` と `rfo` に分かれます。各サブセクションは [YAMLリファレンス](yaml-reference.md) と同じキーを持ちますが、デフォルトは `out_dir: ./result_path_search/`、`dump: False` です。

`bond` は UMA ベースの結合変化検出で、{ref}`scan <section-bond>` と共通の `device`, `bond_factor`, `margin_fraction`, `delta_fraction` を持ちます。

`search` は再帰ロジックを制御します: `max_depth`, `stitch_rmsd_thresh`, `bridge_rmsd_thresh`, `max_nodes_segment`, `max_nodes_bridge`, `kink_max_nodes`, `max_seq_kink`, `refine_mode`（`null` の場合は GSM→`peak`、DMF→`minima` を自動選択）。旧 `rmsd_align` フラグは互換性のため保持されますが無視されます。

`dmf` は `--mep-mode dmf` 選択時に適用される Direct Max Flux + (C)FB-ENM の設定です。既定値は `DMF_KW` を踏襲し、実行ごとに上書きできます。

### YAML例（デフォルト値）
```yaml
geom:
  coord_type: cart           # coordinate type: cartesian vs dlc internals
  freeze_atoms: []           # 0-based frozen atoms merged with CLI/link detection
calc:
  charge: 0                  # total charge (CLI/template override)
  spin: 1                    # spin multiplicity 2S+1
  model: uma-s-1p1           # UMA model tag
  task_name: omol            # UMA task name
  device: auto               # UMA device selection
  max_neigh: null            # maximum neighbors for graph construction
  radius: null               # cutoff radius for neighbor search
  r_edges: false             # store radial edges
  out_hess_torch: true       # request torch-form Hessian
  freeze_atoms: null         # calculator-level frozen atoms
  hessian_calc_mode: FiniteDifference   # Hessian mode selection
  return_partial_hessian: false         # full Hessian (avoids shape mismatches)
gs:
  fix_first: true            # keep the first endpoint fixed during optimization
  fix_last: true             # keep the last endpoint fixed during optimization
  max_nodes: 10              # maximum string nodes
  perp_thresh: 0.005         # perpendicular displacement threshold
  reparam_check: rms         # reparametrization check metric
  reparam_every: 1           # reparametrization stride
  reparam_every_full: 1      # full reparametrization stride
  param: equi                # parametrization scheme
  max_micro_cycles: 10       # micro-iteration limit
  reset_dlc: true            # rebuild delocalized coordinates each step
  climb: true                # enable climbing image
  climb_rms: 0.0005          # climbing RMS threshold
  climb_lanczos: true        # Lanczos refinement for climbing
  climb_lanczos_rms: 0.0005  # Lanczos RMS threshold
  climb_fixed: false         # keep climbing image fixed
  scheduler: null            # optional scheduler backend
opt:
  type: string               # optimizer type label
  stop_in_when_full: 300     # early stop threshold when string is full
  align: false               # alignment toggle (kept off)
  scale_step: global         # step scaling mode
  max_cycles: 300            # maximum optimization cycles
  dump: false                # dump trajectory/restart data
  dump_restart: false        # dump restart checkpoints
  reparam_thresh: 0.0        # reparametrization threshold
  coord_diff_thresh: 0.0     # coordinate difference threshold
  out_dir: ./result_path_search/   # output directory
  print_every: 10            # logging stride
dmf:
  correlated: true           # correlated DMF propagation
  sequential: true           # sequential DMF execution
  fbenm_only_endpoints: false   # run FB-ENM beyond endpoints
  fbenm_options:
    delta_scale: 0.2         # FB-ENM displacement scaling
    bond_scale: 1.25         # bond cutoff scaling
    fix_planes: true         # enforce planar constraints
    two_hop_mode: sparse     # neighbor traversal strategy
  cfbenm_options:
    bond_scale: 1.25         # CFB-ENM bond cutoff scaling
    corr0_scale: 1.1         # Correlation scale for corr0
    corr1_scale: 1.5         # Correlation scale for corr1
    corr2_scale: 1.6         # Correlation scale for corr2
    eps: 0.05                # Correlation epsilon
    pivotal: true            # Pivotal residue handling
    single: true             # Single-atom pivots
    remove_fourmembered: true  # Prune four-membered rings
    two_hop_mode: sparse     # Neighbor traversal strategy
  dmf_options:
    remove_rotation_and_translation: false  # Keep rigid-body motions
    mass_weighted: false               # Toggle mass weighting
    parallel: false                    # Enable parallel DMF
    eps_vel: 0.01                      # Velocity tolerance
    eps_rot: 0.01                      # Rotational tolerance
    beta: 10.0                         # Beta parameter for DMF
    update_teval: false                # Update transition evaluation
  k_fix: 300.0                         # Harmonic constant for restraints
sopt:
  lbfgs:
    thresh: gau                # LBFGS convergence preset
    max_cycles: 10000          # iteration limit
    print_every: 100           # logging stride
    min_step_norm: 1.0e-08     # minimum accepted step norm
    assert_min_step: true      # assert when steps stagnate
    rms_force: null            # explicit RMS force target
    rms_force_only: false      # rely only on RMS force convergence
    max_force_only: false      # rely only on max force convergence
    force_only: false          # skip displacement checks
    converge_to_geom_rms_thresh: 0.05   # RMS threshold when targeting geometry
    overachieve_factor: 0.0    # tighten thresholds
    check_eigval_structure: false   # validate Hessian eigenstructure
    line_search: true          # enable line search
    dump: false                # dump trajectory/restart data
    dump_restart: false        # dump restart checkpoints
    prefix: ""                 # filename prefix
    out_dir: ./result_path_search/   # output directory
    keep_last: 7               # history size for LBFGS buffers
    beta: 1.0                  # initial damping beta
    gamma_mult: false          # multiplicative gamma update toggle
    max_step: 0.3              # maximum step length
    control_step: true         # control step length adaptively
    double_damp: true          # double damping safeguard
    mu_reg: null               # regularization strength
    max_mu_reg_adaptions: 10   # cap on mu adaptations
  rfo:
    thresh: gau                # RFOptimizer convergence preset
    max_cycles: 10000          # iteration cap
    print_every: 100           # logging stride
    min_step_norm: 1.0e-08     # minimum accepted step norm
    assert_min_step: true      # assert when steps stagnate
    rms_force: null            # explicit RMS force target
    rms_force_only: false      # rely only on RMS force convergence
    max_force_only: false      # rely only on max force convergence
    force_only: false          # skip displacement checks
    converge_to_geom_rms_thresh: 0.05   # RMS threshold when targeting geometry
    overachieve_factor: 0.0    # tighten thresholds
    check_eigval_structure: false   # validate Hessian eigenstructure
    line_search: true          # enable line search
    dump: false                # dump trajectory/restart data
    dump_restart: false        # dump restart checkpoints
    prefix: ""                 # filename prefix
    out_dir: ./result_path_search/   # output directory
    trust_radius: 0.1          # trust-region radius
    trust_update: true         # enable trust-region updates
    trust_min: 0.0             # minimum trust radius
    trust_max: 0.1             # maximum trust radius
    max_energy_incr: null      # allowed energy increase per step
    hessian_update: bfgs       # Hessian update scheme
    hessian_init: calc         # Hessian initialization source
    hessian_recalc: 200        # rebuild Hessian every N steps
    hessian_recalc_adapt: null # adaptive Hessian rebuild factor
    small_eigval_thresh: 1.0e-08   # eigenvalue threshold for stability
    alpha0: 1.0                # initial micro step
    max_micro_cycles: 50       # micro-iteration limit
    rfo_overlaps: false        # enable RFO overlaps
    gediis: false              # enable GEDIIS
    gdiis: true                # enable GDIIS
    gdiis_thresh: 0.0025       # GDIIS acceptance threshold
    gediis_thresh: 0.01        # GEDIIS acceptance threshold
    gdiis_test_direction: true # test descent direction before DIIS
    adapt_step_func: true      # adaptive step scaling toggle
bond:
  device: cuda                # UMA device for bond analysis
  bond_factor: 1.2            # covalent-radius scaling
  margin_fraction: 0.05       # tolerance margin for comparisons
  delta_fraction: 0.05        # minimum relative change to flag bonds
search:
  max_depth: 10               # recursion depth limit
  stitch_rmsd_thresh: 0.0001  # RMSD threshold for stitching segments
  bridge_rmsd_thresh: 0.0001  # RMSD threshold for bridging nodes
  rmsd_align: true            # legacy alignment flag (ignored)
  max_nodes_segment: 10       # max nodes per segment
  max_nodes_bridge: 5         # max nodes per bridge
  kink_max_nodes: 3           # max nodes for kink optimizations
  max_seq_kink: 2             # max sequential kinks
  refine_mode: null           # optional refinement strategy (auto-chooses when null)
```