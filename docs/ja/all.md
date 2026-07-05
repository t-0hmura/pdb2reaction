# `all`

`pdb2reaction all` は、抽出から解析までの一連の処理を **まとめて実行する最上位コマンド** です。`extract` → `scan` / `path-search` → `tsopt` → `irc` / `freq` / `dft` を手動で連結する代わりに、構造から検証済みの機構までを 1 コマンドで得られます。1 つ以上の PDB 入力から活性部位モデル（バインディングポケット）を抽出し、（任意の）段階的スキャンを行い、MEP 探索（デフォルトで単一パス `path-opt`、`--refine-path True` で再帰的 `path-search` に切り替え）を実行したうえで、必要に応じて TS 最適化・IRC・振動解析・DFT 一点計算まで連結します。MLIP バックエンドはデフォルトで UMA を使用しますが、`-b/--backend` オプションで ORB・MACE・AIMNet2 も選択できます。

`all` は与える入力に応じて次の 3 つのモードのいずれかで動作します。

- **複数構造 MEP** — 反応順に並べた 2 構造以上（PDB/GJF/XYZ）と基質定義を与え、活性部位モデル抽出 → GSM/DMF MEP 探索 → 全系テンプレートへのマージまで一括実行する場合。
- **単一構造 + 段階的スキャン** — 1 つの構造に対して `-s/--scan-lists` を 1 つ以上与え、（段階的）スキャンで得た中間体列を MEP の端点として用いる場合。
- **TSOPT のみ** — 1 つの入力構造に `--scan-lists` を省略して `--tsopt` を指定し、MEP/マージをスキップして TS 最適化 + IRC（必要に応じて freq / DFT）だけ実行する場合。

```{important}
`--tsopt` **なし**の `all` ワークフローは **TS 候補**（MEP 探索の最高エネルギー画像 / HEI）を出力します。`--tsopt` を追加すると、これらを虚振動数チェックで検証済みの最適化 TS 構造に精密化し、続いて IRC でエンドポイントを検証します。結果（虚振動数の本数と端点の接続性）は機構解釈の前に必ず目視で確認してください。
```

## 実行例

[`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples) ディレクトリに GPP C6-メチル基転移酵素 BezA（[Tsutsumi et al., *Angew. Chem. Int. Ed.* 2022, 61, e202111217](https://doi.org/10.1002/anie.202111217)）の完全な `all` ワークフロースクリプト（MEP およびスキャンパイプライン）があります。

コマンド形式:

```bash
pdb2reaction all -i INPUT1 [INPUT2 ...] -c SUBSTRATE [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] [options]
```

TS 最適化・IRC・熱化学・DFT まで一括実行する複数構造 MEP:

```bash
# TS 最適化・IRC・熱化学・DFT まで一括実行する複数構造 MEP
pdb2reaction all -i 1.R.pdb 3.P.pdb -c "SAM,GPP,MG" -l "SAM:1,GPP:-3" \
 --tsopt --thermo --dft --out-dir ./result_mep
```

単一構造 + 段階的スキャン（2 ステージ）:

```bash
# 単一構造 + 段階的スキャン（2 ステージ）
pdb2reaction all -i 1.R.pdb -c "SAM,GPP,MG" -l "SAM:1,GPP:-3" \
 -s '[("CS1 SAM 320","GPP 321 C7",1.60)]' '[("GPP 321 H11","GLU 186 OE2",0.90)]' \
 --tsopt --thermo --out-dir ./result_scan
```

TSOPT のみワークフロー（経路探索なし）:

```bash
# TSOPT のみワークフロー（経路探索なし）
pdb2reaction all -i TS_candidate.pdb -c 'SAM,GPP,MG' \
 -l 'SAM:1,GPP:-3' --tsopt --thermo --dft
```

テンプレートがある場合の XYZ/TRJ → PDB/GJF 変換（付随ファイルの生成）は、全ステージ共通の `--convert-files/--no-convert-files`（デフォルト: `True`）で制御できます。

ヘルプ出力は `pdb2reaction all --help` で主要オプションを、`pdb2reaction all --help-advanced` で全オプションを確認できます。

## 処理の流れ

```text
全系入力 (PDB/XYZ/GJF)
 │
 ├─ (任意) 活性部位モデル抽出 [`extract`](extract.md) ← --center/-c を使う場合は PDB が必要
 │ ↓
 │ 活性部位モデル/クラスターモデル (PDB)
 │ │
 │ ├─ (任意) 段階的スキャン [`scan`](scan.md) ← 単一構造ワークフロー
 │ │ ↓
 │ │ 順序付けられた中間体
 │ │ ↓
 │ └─ MEP 探索 [`path-opt`](path-opt.md) または [`path-search`](path-search.md)
 │ ↓
 │ MEP 経路 (mep_trj.xyz) + エネルギーダイアグラム
 │ ↓
 └─ (任意) TS 最適化 + IRC [`tsopt`](tsopt.md) → [`irc`](irc.md)
 └─ (任意) 熱化学 [`freq`](freq.md)
 └─ (任意) DFT 一点計算 [`dft`](dft.md)
```

`all` は次のステージを順に実行します。各ステージはサブコマンドとして単独でも実行できます。

0. **事前チェック**（自動）
 - `all` は他の処理の前に `add-elem-info` と `fix-altloc` を自動実行します。`add-elem-info` は元素記号が欠落している PDB 入力（列 77–78 の元素欄が空）に対してのみ実行され、欠落した元素記号を補完します。`fix-altloc` は代替配座（altLoc）を含む PDB 入力に対してのみ実行され、代替配座を解消します。個別のサブコマンド（`extract`、`opt` など）を使用する場合は、必要に応じてこれらを手動で実行してください。

1. **活性部位モデル抽出**（`-c/--center` が指定された場合）
 - 基質は PDB パス、残基 ID（`123,124` または `A:123,B:456`）、または残基名（`GPP,SAM`）で指定可能
 - 抽出オプション: `--radius`、`--radius-het2het`、`--include-h2o`、`--exclude-backbone`、`--add-linkh`、`--selected-resn`、`--verbose`
 - 入力ごとの活性部位モデル PDB は `<out-dir>/_work/models/` に保存。複数構造が提供された場合、活性部位モデルは残基選択ごとに統合
 - **最初の活性部位モデルの総電荷**がスキャン/MEP/TSOPT に伝播

2. **オプションの段階的スキャン（単一入力のみ）**
 - 各 `--scan-lists` 引数は MLIP スキャンステージを記述する `(i,j,target_Å)` タプルの Python ライクなリスト。原子インデックスは元の入力順序（1 始まり）を参照し、活性部位モデル順序に自動変換されます。PDB 入力の場合、`i`/`j` は整数インデックスまたは `'TYR,285,CA'` のようなセレクタ文字列を使用可能です。セレクタはスペース/カンマ/スラッシュ/バッククォート/バックスラッシュ（` ` `,` `/` `` ` `` `\`）を区切り文字として受け付け、トークン順序は任意です（フォールバックは resname, resseq, atom と仮定）。
 - 単一リテラルは 1 ステージスキャンを実行し、複数リテラルは**順次**実行されるため、ステージ 2 はステージ 1 の結果から開始されます。複数リテラルは 1 つの `-s/--scan-lists` に並べて指定します（例: `-s '[(…)]' '[(…)]'`）。
 - ステージエンドポイント（`stage_XX/result.pdb`）が、後続 MEP ステップへ渡される順序付き中間体となる

3. **活性部位モデルでの MEP 探索（デフォルトで単一パス `path-opt`、`--refine-path True` で再帰的 `path-search`）**
 - デフォルトでは、単一パス `path-opt`（GSM/DMF）を実行します。エンジン生出力は `<out-dir>/_work/path_opt/` に書かれ、マージ済み成果物（`mep.pdb`、`mep_trj.xyz`、`energy_diagram_MEP.png`）はルート直下へ配置します。
 - `--refine-path True` を指定すると、再帰的 `path-search` に切り替わり、多段階反応を自動検出して各素反応の詳細な MEP を構築します（エンジン生出力は `<out-dir>/_work/path_search/`）。複雑な多段階反応では手動での試行錯誤が必要な場合があります。
 - 複数入力 PDB の場合、参照マージ用の全系テンプレートが MEP エンジン（デフォルトは `path-opt`、`--refine-path True` 時は `path-search`）に自動的に渡されます（全系マージ自体は `--refine-path True` 時のみ実行）。単一構造スキャンの場合は、元の全系 PDB テンプレートが全ステージで再利用されます。

4. **活性部位モデルを全系にマージ**（`--refine-path True` 使用時のみ）
 - `--refine-path True` を指定し、かつ参照 PDB テンプレートがある場合、マージ済みの `mep_w_ref.pdb` が `<out-dir>/` 直下へ配置し、セグメントごとの `mep_w_ref_seg_NN.pdb` は `<out-dir>/_work/path_search/` に残ります。デフォルトの `path-opt` モードでは全系マージは行われません。

5. **オプションのセグメントごとの後処理**（反応セグメントのみ — 結合変化のあるセグメント。ブリッジセグメントはスキップ）
 - `--tsopt`: 各 HEI 活性部位モデルで TS 最適化を実行し、EulerPC IRC で追跡した後、IRC エンドポイントを `--thresh-post`（デフォルト `baker`）で再最適化します。エンドポイント最適化の作業ディレクトリは完了後に自動削除されます。
 - `--thermo`: (R, TS, P) で `freq` を呼び出し、振動/熱化学データと MLIP Gibbs ダイアグラムを取得
 - `--dft`: (R, TS, P) で DFT 一点計算を実行し、DFT ダイアグラムを構築。`--thermo` と組み合わせると DFT//MLIP Gibbs ダイアグラムも生成
  - 共有の上書きオプション: `--opt-mode`、`--opt-mode-post`（TSOPT/IRC 後最適化のプリセット上書き）、`--flatten/--no-flatten`、`--hessian-calc-mode`、`--tsopt-max-cycles`、`--tsopt-out-dir`、`--freq-*`、`--dft-*`、`--dft-engine`（GPU 優先）など
 - Hessian評価モードの詳細は {ref}`ja-hessian-evaluation` を参照してください。

6. **TSOPT のみモード**（単一入力、`--tsopt`、`--scan-lists` なし）
 - MEP/マージステージをスキップし、活性部位モデル（または抽出がスキップされた場合は全入力構造）で `tsopt` → EulerPC IRC を実行し、高エネルギー側の IRC 終端を反応物 (R) として識別したうえで、エネルギーダイアグラム一式とオプションの freq/DFT 出力を生成します。

## 出力

ツリーは 3 つのゾーンで構成されます: **ルート直下の成果物**、**`segments/seg_NN/` 配下のセグメント別成果物**、**`_work/` 配下のパイプライン作業領域**（必要な結果を取り出したあとは `rm -rf` で削除して構いません）。

```text
out_dir/ (デフォルト:./result_all/)
├─ summary.log                  # 結果要約（ルート直下に生成）
├─ summary.json                 # JSON 結果
├─ mep.pdb                      # マージ済み MEP 経路（エンジンから配置）
├─ mep_w_ref.pdb               # 全系テンプレートへマージした MEP（複数入力時）
├─ mep_trj.xyz                 # MEP 全体軌道
├─ energy_diagram_MEP.png      # 全セグメントの MEP 障壁
├─ energy_diagram_*.png        # 集約後処理ダイアグラム（UMA / Gibbs / DFT、--tsopt 等で生成）
├─ segments/                    # 反応セグメント別の成果物（ブリッジセグメントはスキップ）
│  └─ seg_NN/                   # 2 桁インデックス、例: seg_01, seg_02
│     ├─ reactant.{pdb,xyz,gjf}   #   正規 R/TS/P（出力形式は入力形式と一致）
│     ├─ ts.{pdb,xyz,gjf}
│     ├─ product.{pdb,xyz,gjf}
│     ├─ ts/                    # TS 最適化出力と振動解析（--tsopt）
│     ├─ irc/                   # IRC 軌道とプロット（--tsopt）
│     ├─ freq/{R,TS,P}/         # frequencies_cm-1.txt + thermoanalysis.yaml（--thermo）
│     └─ dft/                   # DFT 一点計算結果（--dft）
└─ _work/                       # パイプライン作業領域（削除可）
   ├─ models/                   # 抽出実行時の活性部位モデル PDB（model_<input_stem>.pdb）
   ├─ scan/                     # 段階的スキャン結果（--scan-lists 提供時）
   ├─ add_elem_info/            # 前処理: 元素記号補完
   ├─ fix_altloc/               # 前処理: altLoc 解決
   └─ path_opt/                 # MEP エンジン生出力（--refine-path True 時は path_search/）
```

**TSOPT のみモード**（単一入力 + `--tsopt`、`--scan-lists` なし）では MEP ステージが無く、最適化済み R/TS/P と `ts/`・`irc/`・`freq/`・`dft/` は `segments/seg_01/` 直下に生成され、MEP 作業ディレクトリ（`_work/path_opt/`）は存在しません。

```{note}
**正規構造は `segments/seg_NN/reactant.*`・`ts.*`・`product.*`** です — 機構を報告する際はこれらを引用してください。同じ `seg_NN/` 内の `ts/`・`irc/`・`freq/`・`dft/` サブディレクトリは各ステージの作業ファイル（例: `ts/vib/imag_*_trj.xyz`、`irc/*_trj.xyz`）を保持し、特定ステージのデバッグに使います。`_work/path_opt/` 配下の MEP エンジン生出力は作業領域であり、必要な成果物（`mep.pdb`、`mep_trj.xyz`、`energy_diagram_MEP.png`）は既にルートへ配置済みです。
```

`-v 2` ではコンソールに活性部位モデルの電荷解決結果、YAML 設定、スキャンステージ、MEP（GSM/DMF）の進行状況、各ステージの所要時間が出力されます。詳細は {ref}`ja-verbosity-levels` を参照してください。

### プロットファイルの命名規則

エネルギーダイアグラムファイルは手法とスコープに基づいて命名されます:

| ファイル名 | 生成タイミング | 内容 |
|---|---|---|
| `energy_diagram_MEP.png` | path-opt/path-search 完了時 | 全セグメント MEP 障壁（生の GSM/DMF 値） |
| `energy_diagram_UMA.png` | セグメントごとの tsopt+IRC 完了時 | R→TS→P（MLIP エネルギー） |
| `energy_diagram_G_UMA.png` | セグメントごとの thermo 完了時 | R→TS→P（MLIP ギブズ自由エネルギー） |
| `energy_diagram_DFT.png` | セグメントごとの DFT 完了時 | R→TS→P（DFT エネルギー） |
| `energy_diagram_G_DFT_plus_UMA.png` | セグメントごとの DFT+thermo 完了時 | R→TS→P（DFT エネルギー + MLIP 熱補正） |
| `energy_diagram_UMA_all.png` | 全セグメント集約時 | 全セグメント統合（MLIP） |
| `energy_diagram_G_UMA_all.png` | 全セグメント + thermo | 全セグメント統合（MLIP ギブズ） |
| `energy_diagram_DFT_all.png` | 全セグメント + DFT | 全セグメント統合（DFT） |
| `energy_diagram_G_DFT_plus_UMA_all.png` | 全セグメント + DFT + thermo | 全セグメント統合（DFT//MLIP ギブズ） |
| `irc_plot.png`（`segments/seg_NN/irc/` 配下） | セグメント IRC 完了 | セグメントごとの IRC プロファイル（MLIP エネルギー） |
| `irc_plot_all.png` | 全セグメント集約 | 全セグメントの IRC プロファイル連結 |

### `summary.log` の読み方
ログは番号付きセクションで構成されます:
- **[1] グローバル MEP 概要** – イメージ/セグメント数、MEP 軌跡プロットのパス、MEP 全体のエネルギーダイアグラム。
- **[2] セグメント別 MEP サマリー（MLIP パス）** – セグメントごとの障壁（`ΔE‡`）、反応エネルギー（`ΔE`）、結合変化サマリー。
- **[3] セグメント別後処理（TSOPT / Thermo / DFT）** – TS 虚振動数チェック、IRC 出力、MLIP/熱化学/DFT のエネルギーテーブル。
- **[4] エネルギーダイアグラム（概要）** – MEP/MLIP/Gibbs/DFT 系の図表と、任意の横断サマリー表。
- **[5] 出力ディレクトリ構造** – 生成ファイルを注釈付きでまとめたツリー。

### `summary.json` の読み方
JSON 結果の代表的なトップレベルキーは以下のとおりです。
- `out_dir`, `n_images`, `n_segments` – 実行メタデータと総数。
- `segments` – `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, `bond_changes` を含むセグメント配列。
- `energy_diagrams`（任意） – `labels`, `energies_kcal`, `energies_au`, `ylabel`, `image` などを含む図表データ。

`summary.json` には `summary.log` にある整形テーブルやファイルツリーは含まれません。

## CLI オプション

入力の要件:

- 抽出有効（`-c/--center`）: 残基同定のため入力は **PDB** が必須。
- 抽出なし: **PDB/XYZ/GJF** を使用可能。
- 複数構造実行は 2 つ以上の構造が必要。完全な入力ファイル要件（水素、元素列、原子順序の一致）は [CLI 規約](cli-conventions.md) を参照してください。

電荷の解決順序の詳細は {ref}`CLI 規約: 電荷の指定 <ja-charge-specification>` を参照してください。`all` コマンドでは、活性部位モデル抽出（`-c` 指定時）による電荷導出が追加の優先度レイヤーとして機能します。スピンの解決順序は `--multiplicity`（CLI）→ `.gjf` テンプレート → デフォルト（1）。非標準の基質には `--ligand-charge/-l` を必ず指定し、scan/MEP/TSOPT/DFT へ正しい総電荷を伝播させてください。

### 入出力オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 反応順序の 2 つ以上の完全構造（`--scan-lists` または `--tsopt` のみ単一入力可） | 必須 |
| `--ref-pdb FILE` | `-i` で XYZ 入力を使用する場合のトポロジー参照 PDB | _None_ |
| `-o, --out-dir PATH` | トップレベル出力ディレクトリ | `./result_all/` |
| `--convert-files/--no-convert-files` | XYZ/TRJ → 対応する PDB/GJF の全体切替 | `True` |
| `--dump/--no-dump` | MEP(GSM/DMF)軌跡を出力。`path-search`/`path-opt` には常時転送され、`scan`/`tsopt` には明示指定時のみ転送。`freq` はデフォルトで dump=True なので `--no-dump` で無効化 | `False` |
| `--config FILE` | 先に適用するベース YAML | _None_ |
| `--show-config/--no-show-config` | 実行前に解決済み設定を表示 | `False` |
| `--dry-run/--no-dry-run` | 実行せず検証と計画表示のみ行う | `False` |

### 電荷・スピンオプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-l, --ligand-charge TEXT` | 総電荷または残基別マッピング（`-q` 省略時に使用、推奨）。PDB 入力（または `--ref-pdb` 付き XYZ/GJF）で extract と同じ全系電荷導出を起動します | _None_ |
| `-q, --charge INT` | 総電荷を強制上書き（`--ligand-charge` より優先） | _None_ |
| `-m, --multiplicity INT` | 全下流ステップへ転送されるスピン多重度 | `1` |

### 活性部位モデル抽出オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-c, --center TEXT` | 基質指定: PDB パス、残基 ID（`123,124` / `A:123,B:456`）、または残基名（`GPP,SAM`）。省略すると抽出をスキップし、全構造を下流へ渡す | 抽出に必須 |
| `-r, --radius FLOAT` | 活性部位モデル包含カットオフ（Å） | `2.6` |
| `--radius-het2het FLOAT` | ヘテロ–ヘテロカットオフ（Å）。`0` を渡すと空の選択を避けるため内部で `0.001 Å` に自動補正されます（単体の `extract` と同じ挙動） | `0.0` |
| `--include-h2o/--no-include-h2o` | 水分子を含める（HOH/WAT/TIP3/SOL） | `True` |
| `--exclude-backbone/--no-exclude-backbone` | 非基質アミノ酸の主鎖原子を除去 | `False` |
| `--add-linkh/--no-add-linkh` | 切断結合にキャップ水素を付加 | `True` |
| `--selected-resn TEXT` | 強制包含残基。**名前とは裏腹にこのフラグは残基 ID（コロン区切り整数、オプションでチェーン/挿入コード付き、例 `A:123A`）を受け付け、3 文字残基名は受け付けません。** 残基名ベースの選択には `-c/--center 'GPP,SAM'` を使用してください | `""` |
| `--modified-residue TEXT` | 修飾アミノ酸残基名をカンマ区切りで指定（任意で電荷付き）。主鎖切断と電荷計算にアミノ酸として扱う。例: `HD1,HD2,HD3` または `HD1:0,SEP:-2` | `""` |
| `--freeze-links/--no-freeze-links` | 活性部位モデル PDB でキャップ H の親を凍結 | `True` |

(ja-mep-search-options)=
### MEP 探索オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--mep-mode [gsm\|dmf]` | MEP 探索アルゴリズム: GSM（Growing String Method）または DMF（Direct Max Flux） | `gsm` |
| `--max-nodes INT` | MEP 内部ノード数。**GSM:** 総イメージ数 = `max_nodes + 2`（端点 2 つは固定）。**DMF:** チェーン上の *可動* イメージ数（端点の暗黙的拡張なし） | `20` |
| `--max-cycles INT` | MEP 最大最適化サイクル | `300` |
| `--climb/--no-climb` | 標準 GSM セグメントでクライミングイメージを有効化（ブリッジセグメントは常に無効） | `True` |
| `--opt-mode [grad\|hess]` | ワークフロープリセット（`grad` → L-BFGS/Dimer、`hess` → RFO/RSPRFO）。コマンド個別実行では `opt --opt-mode grad|hess`、`tsopt --opt-mode grad|hess` を推奨。トークンのマッピングはスコープ依存で、`all` の pre-opt デフォルト（`grad`）と `tsopt` のデフォルト（`hess`）は一致しません。詳細は {ref}`ja-opt-mode-semantics` を参照してください | `grad` |
| `--thresh TEXT` | 収束プリセット（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `gau` |
| `--preopt/--no-preopt` | MEP 前に活性部位モデル端点を事前最適化。単体の `scan`、`scan2d`、`scan3d` では `--preopt` のデフォルトは `False`（`--preopt` を渡すと有効化） | `True` |
| `--refine-path BOOL` | `True` の場合は再帰的 `path-search`、`False`（デフォルト）の場合は `path-opt` を連結して再帰的精密化なしで実行。 | `False` |

### MLIP 計算機オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--workers`, `--workers-per-node` | MLIP 予測器の並列度（workers > 1 で解析Hessian無効; UMA バックエンドのみ）。診断上の注意は {ref}`ja-workers-fd-downgrade` を参照 | `1`, `1` |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | 共有 MLIP Hessianエンジン | `FiniteDifference` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--solvent TEXT` | xTB 補正用の暗黙溶媒名（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |

### 後処理オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--tsopt/--no-tsopt` | セグメントごとの TS 最適化+ IRC を実行 | `False` |
| `--thermo/--no-thermo` | R/TS/P で振動解析を実行 | `False` |
| `--dft/--no-dft` | R/TS/P で DFT 一点計算を実行 | `False` |
| `--opt-mode-post [grad\|hess]` | TSOPT/IRC 後最適化のプリセット上書き（`grad` → Dimer/L-BFGS、`hess` → RSPRFO/RFO） | `hess` |
| `--thresh-post TEXT` | IRC 後エンドポイント最適化の収束プリセット（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `baker` |
| `--flatten/--no-flatten` | 余分な虚振動モードのフラット化 | `False` |

```{warning}
`--dft` による DFT 一点計算（PySCF/GPU4PySCF）は、約 300 原子を超えるモデルでは計算コストが非常に大きくなります。そのような系では、A100 や H200 等の高性能 GPU を搭載した HPC クラスタの利用が必要になる場合があります。
```

TSOPT の最適化モードは、`--opt-mode-post`（指定時）→ `--opt-mode`（明示指定時のみ）→ TSOPT のデフォルト（`hess` → `rsprfo`）の順で決まります。

例: `--opt-mode grad --opt-mode-post hess` は、経路最適化に L-BFGS、TS 精密化に RS-P-RFO を使用します。

### TSOPT 上書き

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--tsopt-max-cycles INT` | `tsopt --max-cycles` 上書き | `10000` |
| `--tsopt-out-dir PATH` | tsopt 出力サブディレクトリ | _None_ |

### Freq 上書き

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--freq-out-dir PATH` | freq 出力ディレクトリ上書き | _None_ |
| `--freq-max-write INT` | 最大モード出力数 | `10` |
| `--freq-amplitude-ang FLOAT` | モードアニメーション振幅（Å） | `0.8` |
| `--freq-n-frames INT` | モードアニメーションフレーム数 | `20` |
| `--freq-sort [value\|abs]` | モードソート方法 | `value` |
| `--freq-temperature FLOAT` | 熱化学温度（K） | `298.15` |
| `--freq-pressure FLOAT` | 熱化学圧力（atm） | `1.0` |

### DFT 上書き

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--dft-engine [gpu\|cpu]` | DFT バックエンド: gpu (GPU4PySCF) または cpu (PySCF)。`all` ラッパーではプレフィックス付きで `--dft-engine` と名付けられていますが、単体の `dft` サブコマンドでは同じオプションが `--engine` という名前になります | `gpu` |
| `--dft-out-dir PATH` | DFT 出力ディレクトリ上書き | _None_ |
| `--dft-func-basis TEXT` | 汎関数/基底関数ペア | `wb97m-v/def2-tzvpd` |
| `--dft-max-cycle INT` | 最大 SCF サイクル | `100` |
| `--dft-conv-tol FLOAT` | SCF 収束閾値 | `1e-9` |
| `--dft-grid-level INT` | PySCF グリッドレベル | `3` |

(ja-scan-options-single-input-runs)=
### スキャンオプション（単一入力）

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-s, --scan-lists TEXT...` | 段階的スキャン: `(i,j,target_Å)` タプル | _None_ |
| `--scan-out-dir PATH` | scan 出力ディレクトリ上書き | _None_ |
| `--scan-one-based/--no-scan-one-based` | 1 始まり/0 始まりインデックス | _None_ |
| `--scan-max-step-size FLOAT` | 最大ステップサイズ（Å） | `0.20` |
| `--scan-bias-k FLOAT` | 調和バイアス強度（eV·Å⁻²） | `300` |
| `--scan-relax-max-cycles INT` | 緩和サイクル上限 | `10000` |
| `--scan-preopt/--no-scan-preopt` | scan の事前最適化トグルを上書き | _None_ |
| `--scan-endopt/--no-scan-endopt` | scan のステージ終端最適化トグルを上書き | _None_ |

## YAML 設定

`all` は YAML の多層指定をサポートします:

- `--config FILE`: ベース設定。

適用順序:

`defaults < config < CLI`

解決後の YAML は呼び出されるすべてのサブコマンドに転送されます。各ツールが読み取るセクションは以下のとおりです:

| サブコマンド | YAML セクション |
|------------|-----------------|
| [`path-opt`](path-opt.md) | `geom`, `calc`, `gs`, `dmf`, `stopt`, `opt` |
| [`path-search`](path-search.md) | `geom`, `calc`, `gs`, `stopt`, `opt`, `bond`, `search` |
| [`scan`](scan.md) | `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`, `bond` |
| [`tsopt`](tsopt.md) | `geom`, `calc`, `opt`, `hessian_dimer`, `rsirfo` |
| [`freq`](freq.md) | `geom`, `calc`, `freq`, `thermo` |
| [`dft`](dft.md) | `dft` |
| [`irc`](irc.md) | `geom`, `calc`, `irc` |

**最小例:**
```yaml
calc:
 model: uma-s-1p1 # uma-s-1p1 | uma-m-1p1
 hessian_calc_mode: Analytical # VRAM に余裕がある場合推奨
gs:
 max_nodes: 12
 climb: true
dft:
 grid_level: 6
```

すべての YAML オプションの完全なリファレンスについては、**[YAML 設定リファレンス](yaml-reference.md)** を参照してください。

## 注記

```{tip}
大規模な活性部位モデルでは、単一構造スキャンワークフロー（`--scan-lists`）の方が、複数構造 MEP ワークフローよりも信頼性の高い反応障壁を得やすい傾向があります。複数の PDB を入力すると、反応座標と無関係な領域の構造差異が蓄積し、障壁を過大評価する可能性があります。スキャンワークフローは単一構造から出発して関連する座標のみを駆動するため、無関係な構造ノイズを最小化できます。この影響はモデルサイズが大きくなるほど顕著になります。
```

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。
- 形式電荷を推定できない場合は `--ligand-charge`（数値または残基別マッピング）を必ず指定し、scan/MEP/TSOPT/DFT へ正しい総電荷を伝播させてください。
- マージ用の参照 PDB テンプレートは元の入力から自動的に導出されます。`path-search` の `--ref-full-pdb` はこのラッパーでは意図的に非公開です。
- 収束プリセット: `--thresh` のデフォルトは `gau`、`--thresh-post` のデフォルトは `baker`。
- 抽出半径: `--radius` または `--radius-het2het` に `0` を渡すと、内部で `0.001 Å` にクランプされます。
- エネルギーダイアグラムは反応物（最初の状態）基準の kcal/mol で表示されます。
- `-c/--center` を省略すると抽出をスキップし、全構造をそのまま MEP/tsopt/freq/DFT に渡します。ただし単一構造実行では `--scan-lists` か `--tsopt` が必要です。

## 関連項目

- [インストール](installation.md) — セットアップと依存関係
- [はじめに](getting-started.md) — 初回実行、ワークフロー概要、主要概念
- [extract](extract.md) — 単独の活性部位モデル抽出（`all` が内部で呼び出し）
- [scan](scan.md) — 単独の段階的距離スキャン
- [path-opt](path-opt.md) — 単一パス MEP 最適化（GSM/DMF）
- [path-search](path-search.md) — 再帰的 MEP 探索（`all` が内部で呼び出し）
- [tsopt](tsopt.md) — 単独の TS 最適化
- [irc](irc.md) — 単独の IRC 計算
- [freq](freq.md) — 単独の振動解析
- [dft](dft.md) — 単独の DFT 計算
- [典型エラー別レシピ](recipes-common-errors.md) — 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと対処法
- [YAML リファレンス](yaml-reference.md) — 全 YAML 設定オプション
- [用語集](glossary.md) — MEP、TS、IRC、GSM、DMF の定義
