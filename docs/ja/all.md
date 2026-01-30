# `all` サブコマンド

## 概要
`pdb2reaction all` は、**すべてのパイプラインステージ**を一度に実行するコマンドです: ポケット抽出 → オプションの段階的UMAスキャン → 再帰的MEP探索（`path_search`、GSM/DMF） → 全系マージ → オプションのTS最適化 + IRC（`tsopt`） → オプションの振動解析（`freq`） → オプションの一点DFT（`dft`）。このコマンドは複数構造アンサンブルを受け入れ、単一構造スキャンを順序付けられた中間体に変換し、TSOPTのみのポケットワークフローにフォールバックできます。すべての下流ツールは単一のCLIサーフェスを共有するため、1つの呼び出しから長い反応キャンペーンを調整できます。すべてのステージにわたるフォーマット対応XYZ/TRJ → PDB/GJF変換は、共有の `--convert-files {True\|False}` フラグ（デフォルトで有効）によって制御されます。

主要モード:
- **エンドツーエンドアンサンブル** – 反応順序で2つ以上のPDB/GJF/XYZファイルと基質定義を入力; コマンドはポケットを抽出し、GSM/DMF MEP探索を実行し、親PDBにマージし、オプションで反応セグメントごとにTSOPT/freq/DFTを実行
- **単一構造 + 段階的スキャン** – 1つ以上の `--scan-lists` と共に1つの構造を提供; 抽出されたポケットでのUMAスキャンがMEPエンドポイントとなる中間体を生成
- 単一の `--scan-lists` リテラルは1ステージスキャンを実行; 複数のリテラルは順次ステージを実行（1つのフラグの後に複数の値として提供、フラグの繰り返しは不可）
- **TSOPTのみポケット精密化** – 1つの入力構造を提供し、`--scan-lists` を省略し、`--tsopt True` を有効化; `-c/--center` が指定されている場合はポケットを抽出し、その単一システムでTS最適化 + IRC（オプションでfreq/DFT）のみを実行

## 使用法
```bash
pdb2reaction all -i INPUT1 [INPUT2 ...] -c SUBSTRATE [options]
```

### 例
```bash
# 明示的なリガンド電荷と後処理を伴う複数構造アンサンブル
pdb2reaction all -i reactant.pdb product.pdb -c 'GPP,MMT' \
    --ligand-charge 'GPP:-3,MMT:-1' --mult 1 --freeze-links True \
    --max-nodes 10 --max-cycles 100 --climb True --opt-mode light \
    --out-dir result_all_${date} --tsopt True --thermo True --dft True

# 単一構造段階的スキャン + GSM/DMF + TSOPT/freq/DFT
pdb2reaction all -i single.pdb -c '308,309' \
    --scan-lists '[("TYR,285,CA","MMT,309,C10",2.20),("TYR,285,CB","MMT,309,C11",1.80)]' \
    --opt-mode heavy --tsopt True --thermo True --dft True

# TSOPTのみワークフロー（経路探索なし）
pdb2reaction all -i reactant.pdb -c 'GPP,MMT' \
    --ligand-charge 'GPP:-3,MMT:-1' --tsopt True --thermo True --dft True
```

## ワークフロー
1. **活性部位ポケット抽出**（`-c/--center` が提供された場合）
   - 基質はPDBパス、残基ID（`123,124` または `A:123,B:456`）、または残基名（`GPP,MMT`）で指定可能
   - 抽出オプション: `--radius`、`--radius-het2het`、`--include-H2O`、`--exclude-backbone`、`--add-linkH`、`--selected_resn`、`--verbose`
   - 入力ごとのポケットPDBは `<out-dir>/pockets/` に保存。複数構造が提供された場合、ポケットは残基選択ごとに統合
   - **最初のポケットの総電荷**がスキャン/MEP/TSOPTに伝播

2. **オプションの段階的スキャン（単一入力のみ）**
   - 各 `--scan-lists` 引数はUMAスキャンステージを記述する `(i,j,target_Å)` タプルのPythonライクなリスト
   - 単一リテラルは1ステージスキャンを実行; 複数リテラルは**順次**実行
   - スキャンは電荷/スピン、`--freeze-links`、UMA最適化プリセット（`--opt-mode`）、`--args-yaml`、`--preopt` を継承。`--dump` はこのコマンドで明示指定された場合のみスキャンへ転送され、それ以外は scan 側のデフォルト（`False`）を使用
   - `--scan-out-dir`、`--scan-one-based`、`--scan-max-step-size`、`--scan-bias-k`、`--scan-relax-max-cycles`、`--scan-preopt`、`--scan-endopt` などの上書きフラグが利用可能
   - ステージエンドポイント（`stage_XX/result.pdb`）が後続MEPステップに供給される順序付き中間体となる

3. **ポケットでのMEP探索（再帰的GSM/DMF）**
   - 抽出されたポケット（または抽出がスキップされた場合は元の全構造）を使用してデフォルトで `path_search` を実行
   - `--refine-path False` で再帰的精密化なしのシングルパス `path-opt` GSM/DMFチェーンに切り替え

4. **ポケットを全系にマージ**
   - 参照PDBテンプレートが存在する場合、マージされた `mep_w_ref*.pdb` およびセグメントごとの `mep_w_ref_seg_XX.pdb` ファイルが `<out-dir>/path_search/` に出力

5. **オプションのセグメントごとの後処理**
   - `--tsopt True`: 各HEIポケットでTS最適化を実行、EulerPC IRCで追跡し、セグメントエネルギーダイアグラムを出力
   - `--thermo True`: (R, TS, P) で `freq` を呼び出し振動/熱化学データとUMA Gibbsダイアグラムを取得
   - `--dft True`: (R, TS, P) で一点DFTを起動しDFTダイアグラムを構築。`--thermo True` と組み合わせるとDFT//UMA Gibbsダイアグラムも生成
   - 共有の上書きには `--opt-mode`、`--opt-mode-post`（TSOPT/IRC後最適化のプリセット上書き）、`--flatten-imag-mode`、`--hessian-calc-mode`、`--tsopt-max-cycles`、`--tsopt-out-dir`、`--freq-*`、`--dft-*`、`--dft-engine`（GPU優先）などが含まれる
   - VRAMが十分な場合は `--hessian-calc-mode` を `Analytical` に設定することを強く推奨

6. **TSOPTのみモード**（単一入力、`--tsopt True`、`--scan-lists` なし）
   - MEP/マージステージをスキップ。ポケット（または抽出がスキップされた場合は完全入力）で `tsopt` を実行し、EulerPC IRCを実行
   - 高エネルギー側のIRC終端を反応物 (R) として識別し、同じ種類のエネルギーダイアグラムとオプションの freq/DFT 出力を生成


### 電荷とスピンの優先順位
- 抽出あり: ポケット電荷 = 抽出器のモデル#1総電荷（`--ligand-charge` があればそれを使用）。スピンは `--mult`（デフォルト1）。
- 抽出なし: 明示的な `-q/--charge` が最優先。省略され、`--ligand-charge` が与えられている場合は**全複合体を酵素–基質アダクトとして扱い**、先頭入力がPDBなら `extract.py` の電荷サマリーで総電荷を導出。導出に失敗しても `--ligand-charge` が数値なら総電荷として採用。それ以外は `.gjf` の電荷かデフォルト0。スピンの優先順位は `--mult` → `.gjf` → 1。

### 入力要件
- 抽出有効（`-c/--center`）: 残基同定のため入力は **PDB** が必須。
- 抽出なし: **PDB/XYZ/GJF** を使用可能。
- 複数構造実行は 2 つ以上の構造が必要。


## CLIオプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 反応順序の2つ以上の完全構造（`--scan-lists` または `--tsopt True` を使用する場合のみ単一入力可） | 必須 |
| `-c, --center TEXT` | 基質指定（PDBパス、残基ID `123,124` / `A:123,B:456`、または残基名 `GPP,MMT`） | 抽出に必須 |
| `--out-dir PATH` | トップレベル出力ディレクトリ | `./result_all/` |
| `-r, --radius FLOAT` | ポケット包含カットオフ（Å） | `2.6` |
| `--radius-het2het FLOAT` | 独立したヘテロ–ヘテロカットオフ（Å） | `0.0` |
| `--include-H2O BOOLEAN` | 水分子を含める（`False` で HOH/WAT/TIP3/SOL を除外） | `True` |
| `--exclude-backbone BOOLEAN` | 非基質アミノ酸の主鎖原子を除去 | `True` |
| `--add-linkH BOOLEAN` | 切断結合にリンク水素を付加（炭素のみ） | `True` |
| `--selected_resn TEXT` | 強制包含残基（カンマ/空白区切り、チェーン/挿入コード可） | `""` |
| `--verbose BOOLEAN` | 抽出器のINFOログを有効化 | `True` |
| `--ligand-charge TEXT` | 未知残基の残基別マッピングまたは総電荷（推奨）。`-q` 省略時にPDB入力でextract方式の電荷導出を行い、数値指定は総電荷のフォールバックとしても利用される | _None_ |
| `-q, --charge INT` | 総電荷を強制上書き（extractor/`.gjf`/`--ligand-charge` を上書き、警告を出力） | _None_ |
| `--workers`, `--workers-per-node` | UMA予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --mult INT` | すべての下流ステップに転送されるスピン多重度 | `1` |
| `--freeze-links BOOLEAN` | ポケットPDBでリンクHの親を凍結（scan/tsopt/freq でも再利用） | `True` |
| `--max-nodes INT` | MEP内部ノード数（GSM string / DMF images） | `10` |
| `--max-cycles INT` | MEP最大最適化サイクル（GSM/DMF） | `300` |
| `--climb BOOLEAN` | 各ペアの最初のセグメントでTSクライミングを有効化 | `True` |
| `--opt-mode [light\|heavy]` | scan/tsopt/path_search で共有する最適化プリセット（light → LBFGS/Dimer、heavy → RFO/RSIRFO） | `light` |
| `--opt-mode-post [light\|heavy]` | TSOPT/IRC後最適化のプリセット上書き。未指定なら `--opt-mode` が明示されていればそれを適用、そうでなければ TSOPT は `heavy` | _None_ |
| `--dump BOOLEAN` | MEP(GSM/DMF)軌跡を出力。`path_search`/`path-opt` へ常に転送。scan/tsopt へは明示指定時のみ。freq はこのラッパーでは `dump=True` が既定で `thermoanalysis.yaml` を書くため、`--dump False` で明示的に無効化 | `False` |
| `--convert-files {True\|False}` | テンプレート利用可能時のXYZ/TRJ → PDB/GJFコンパニオンのグローバルトグル | `True` |
| `--thresh TEXT` | MEPおよび単一構造最適化の収束プリセット（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `gau` |
| `--thresh-post TEXT` | IRC後エンドポイント最適化の収束プリセット | `baker` |
| `--args-yaml FILE` | `path_search`/`scan`/`tsopt`/`freq`/`dft` へそのまま転送されるYAML | _None_ |
| `--preopt BOOLEAN` | MEP前にポケット端点を事前最適化（scan の preopt 既定にもなる） | `True` |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | tsopt/freq に転送される共有UMAヘシアンエンジン | `FiniteDifference` |
| `--tsopt BOOLEAN` | セグメントごとのTS最適化 + IRC を実行、またはTSOPTのみモードを有効化 | `False` |
| `--thermo BOOLEAN` | R/TS/Pで振動解析を実行しUMA Gibbsダイアグラムを構築 | `False` |
| `--dft BOOLEAN` | R/TS/Pで一点DFTを実行しDFT図を構築 | `False` |
| `--dft-engine [gpu\|cpu\|auto]` | DFTステージのバックエンド（`auto` はGPU優先でCPUへフォールバック） | `gpu` |
| `--tsopt-max-cycles INT` | 各TS精密化で `tsopt --max-cycles` を上書き | `10000` |
| `--tsopt-out-dir PATH` | tsopt出力サブディレクトリの上書き（相対パスは `<out-dir>` 配下に解決） | _None_ |
| `--flatten-imag-mode {True\|False}` | 余分な虚数モードのフラット化（light: dimer ループ、heavy: RSIRFO後処理） | `False` |
| `--freq-out-dir PATH` | freq 出力ディレクトリの上書き | _None_ |
| `--freq-max-write INT` | `freq --max-write` 上書き | `10` |
| `--freq-amplitude-ang FLOAT` | `freq --amplitude-ang` 上書き（Å） | `0.8` |
| `--freq-n-frames INT` | `freq --n-frames` 上書き | `20` |
| `--freq-sort [value\|abs]` | freq のソート方法上書き | `value` |
| `--freq-temperature FLOAT` | 熱化学温度 (K) 上書き | `298.15` |
| `--freq-pressure FLOAT` | 熱化学圧力 (atm) 上書き | `1.0` |
| `--dft-out-dir PATH` | DFT出力ディレクトリ上書き | _None_ |
| `--dft-func-basis TEXT` | `dft --func-basis` 上書き | `wb97m-v/def2-tzvpd` |
| `--dft-max-cycle INT` | `dft --max-cycle` 上書き | `100` |
| `--dft-conv-tol FLOAT` | `dft --conv-tol` 上書き | `1e-9` |
| `--dft-grid-level INT` | `dft --grid-level` 上書き | `3` |
| `--scan-lists TEXT...` | 抽出ポケットでの段階的スキャンを記述するPythonライクなリスト（単一入力のみ）。各要素は `(i,j,target_Å)`。単一リテラルは1ステージ、複数リテラルは順次ステージ。`i`/`j` は整数インデックスまたは PDB セレクタで内部的に再マップ | _None_ |
| `--scan-out-dir PATH` | scan 出力ディレクトリ上書き（相対パスはデフォルト親配下） | _None_ |
| `--scan-one-based BOOLEAN` | scan のインデックスを1始まり/0始まりに強制（省略時は scan デフォルト=1始まり） | `True` |
| `--scan-max-step-size FLOAT` | scan の `--max-step-size` 上書き（Å） | `0.20` |
| `--scan-bias-k FLOAT` | 調和バイアス強度 `k` 上書き（eV/Å²） | `100` |
| `--scan-relax-max-cycles INT` | scan の緩和サイクル上限上書き | `10000` |
| `--scan-preopt BOOLEAN` | scan 事前最適化の上書き（省略時は `--preopt` に追従） | `True` |
| `--scan-endopt BOOLEAN` | scan ステージ終端最適化の上書き | `True` |

## 出力
```
out_dir/ (デフォルト: ./result_all/)
├─ summary.log               # クイック検査用フォーマット済みサマリー
├─ summary.yaml              # YAMLバージョンサマリー
├─ pockets/                  # 抽出実行時の入力ごとのポケットPDB
├─ scan/                     # 段階的ポケットスキャン結果（--scan-lists提供時）
├─ path_search/              # MEP結果: 軌跡、マージPDB、ダイアグラム
├─ path_search/post_seg_XX/  # 後処理出力（TS最適化、IRC、freq、DFT）
└─ tsopt_single/             # TSOPTのみ出力とIRCエンドポイント
```


- コンソールにはポケット電荷の解決結果、YAML内容、スキャン段数、MEP進行状況（GSM/DMF）、各ステージ時間のサマリーが出力されます。

### `summary.log` の読み方
ログは番号付きセクションで構成されます:
- **[1] グローバルMEP概要** – 画像/セグメント数、MEP軌跡プロットのパス、MEP全体のエネルギーダイアグラム。
- **[2] セグメント別MEPサマリー（UMAパス）** – セグメントごとの障壁（`ΔE‡`）、反応エネルギー（`ΔE`）、結合変化サマリー。
- **[3] セグメント別後処理（TSOPT / Thermo / DFT）** – TS虚数振動数チェック、IRC出力、UMA/熱化学/DFTのエネルギーテーブル。
- **[4] エネルギーダイアグラム（概要）** – MEP/UMA/Gibbs/DFT 系の図表と、任意の横断サマリー表。
- **[5] 出力ディレクトリ構造** – 生成ファイルを注釈付きでまとめたツリー。

### `summary.yaml` の読み方
YAMLは機械可読サマリーです。代表的なトップレベルキーは以下です:
- `out_dir`, `n_images`, `n_segments` – 実行メタデータと総数。
- `segments` – `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, `bond_changes` を含むセグメント配列。
- `energy_diagrams`（任意） – `labels`, `energies_kcal`, `energies_au`, `ylabel`, `image` などを含む図表データ。

`summary.yaml` には `summary.log` にある整形テーブルやファイルツリーは含まれません。

## 注意事項
- 形式電荷が推定できない場合は `--ligand-charge`（数値または残基別マッピング）を必ず指定し、scan/MEP/TSOPT/DFTへ正しい総電荷を伝播させてください。
- マージ用の参照PDBテンプレートは元の入力から自動導出されます。`path_search` の `--ref-full-pdb` はこのラッパーでは意図的に隠されています。
- 収束プリセット: `--thresh` の既定は `gau`、`--thresh-post` の既定は `baker`。
- 抽出半径: `--radius` または `--radius-het2het` に `0` を渡すと、内部で `0.001 Å` にクランプされます。
- ダイアグラムのエネルギーは反応物（最初の状態）基準の kcal/mol で報告されます。
- `-c/--center` を省略すると抽出をスキップして全構造を MEP/tsopt/freq/DFT に渡しますが、単一構造実行には `--scan-lists` か `--tsopt True` が引き続き必要です。
- `--args-yaml` で全計算器を単一設定から制御できます。YAMLはCLIを上書きします。

## YAML設定（`--args-yaml`）
同じYAMLファイルが**すべての**呼び出されるサブコマンドにそのまま転送されます。各ツールは独自のドキュメントに記載されているセクションを読み取ります:

- [`path_search`](path_search.md): `geom`, `calc`, `gs`, `opt`, `sopt`, `bond`, `search`
- [`scan`](scan.md): `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`, `bond`
- [`tsopt`](tsopt.md): `geom`, `calc`, `opt`, `hessian_dimer`, `rsirfo`
- [`freq`](freq.md): `geom`, `calc`, `freq`, `thermo`
- [`dft`](dft.md): `dft`

すべてのYAMLオプションの完全なリファレンスについては、[YAML設定リファレンス](yaml-reference.md) を参照してください。

YAMLルートには必要なセクションだけを含めてください。`geom`/`calc`/`opt` などの共通セクションは複数モジュールで共有され、`freq` や `dft` などのモジュール固有ブロックは対応するステージのみで適用されます。CLIとYAMLの両方が指定された場合、YAMLが優先されます。

共有・固有セクションを組み合わせた例:

```yaml
geom:
  coord_type: cart                     # coordinate type: cartesian vs dlc internals
calc:
  model: uma-s-1p1                     # UMA model tag
  hessian_calc_mode: FiniteDifference  # Hessian mode selection
gs:
  max_nodes: 12                        # maximum string nodes
freq:
  max_write: 8                         # maximum modes written
dft:
  grid_level: 6                        # PySCF grid level
```

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
```
