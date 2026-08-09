# `all`

`pdb2reaction all` は一連の処理を **まとめて実行する最上位コマンド** です。`-c` を指定した場合だけ活性部位モデル（バインディングポケット）を抽出し、省略時は入力構造全体を使います。入力に応じて段階的スキャンまたは MEP 探索（デフォルトは単一パス `path-opt`、`--refine-path` で再帰的 `path-search`）を行い、必要に応じて TS 最適化・IRC・振動解析・DFT 一点計算まで連結します。MLIP バックエンドはデフォルトで UMA を使用しますが、`-b/--backend` で ORB・MACE・AIMNet2 も選択できます。

`all` は与える入力に応じて次の 3 つのモードのいずれかで動作します。

- **複数構造 MEP** — 反応順に並べた 2 構造以上（PDB/mmCIF/GJF/XYZ）を与える場合。`-c` があれば活性部位モデルを抽出し、省略時は入力全体を使って GSM/DMF MEP 探索を行います。
- **単一構造 + 段階的スキャン** — 1 つの構造に対して `-s/--scan-lists` を 1 つ以上与え、（段階的）スキャンで得た中間体列を MEP の端点として用いる場合。
- **TSOPT のみ** — 1 つの入力構造に `--scan-lists` を省略して `--tsopt` を指定し、MEP/マージをスキップして TS 最適化 + IRC（必要に応じて freq / DFT）だけ実行する場合。高エネルギー側の IRC 端点を反応物として提示します。

```{note}
TSOPT のみモードの反応物/生成物ラベルは**エネルギー順に基づく表示上の慣例**であり、化学的に確定した反応方向ではありません。高エネルギー側の IRC 端点を反応物として提示します（エネルギーが厳密に等しい場合は左側の端点を反応物とする決定的な規則）。R/P ラベル、`reactant_irc`/`product_irc` のファイル名、障壁・ΔE はこの慣例のもとで計算されます。機械可読サマリーの `endpoint_assignment`（`policy = "higher_energy_endpoint_as_reactant"`、`chemical_direction_known = false`）にこの方針が明示されるので、ラベルだけから化学的方向を読み取らず、このフィールドを参照してください。中立的な端点名への変更は将来のメジャースキーマに委ねます。
```

```{important}
`--tsopt` **なし**の `all` ワークフローは **TS 候補**（MEP 探索の最高エネルギー画像 / HEI）を出力します。`--tsopt` を追加すると、これらを虚振動数チェックで検証済みの最適化 TS 構造に精密化します。IRC は optimizer が収束した一次鞍点（`n_imag = 1`）を報告した場合だけ開始し、それ以外では downstream 後処理の前に停止します。結果（虚振動数の本数と端点の接続性）は機構解釈の前に必ず目視で確認してください。
```

## 実行例

[`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples) ディレクトリに GPP C6-メチル基転移酵素 BezA（[Tsutsumi et al., *Angew. Chem. Int. Ed.* 2022, 61, e202111217](https://doi.org/10.1002/anie.202111217)）の完全な `all` ワークフロースクリプト（MEP およびスキャンパイプライン）があります。

コマンド形式:

```bash
pdb2reaction all -i INPUT1 [INPUT2 ...] [-c SUBSTRATE] [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] [options]
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

テンプレートがある場合の XYZ/TRJ → PDB/CIF/GJF 変換（付随ファイルの生成）は、全ステージ共通の `--convert-files/--no-convert-files`（デフォルト: `True`）で制御できます。mmCIF／oversized-PDB入力では、内部連携用PDBに加えて元IDを復元したCIFを生成します。

ヘルプ出力は `pdb2reaction all --help` で主要オプションを、`pdb2reaction all --help-advanced` で全オプションを確認できます。

## 処理の流れ

```text
全系入力 (PDB/mmCIF/XYZ/GJF)
 │
 ├─ (任意) 活性部位モデル抽出 [`extract`](extract.md) ← --center/-c は PDB/mmCIF
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

`all` は次のステージを順に実行します。名前の付いた計算ステージはサブコマンドとして単独でも実行できますが、全系へのマージなど `all` 内部だけの処理もあります。

0. **構造ブリッジと事前チェック**（自動）
 - mmCIF、PDB固定幅を超える構造、altLocを含むPDBは、安全に再採番した内部PDBへ一度だけ正規化します。altLocは全geometry workflow共通の入力ブリッジが残基単位で一貫して選択するため、通常は事前の`fix-altloc`は不要です。通常PDBの元素欄（列77–78）が空の場合だけ`all`が`add-elem-info`を実行します。別途cleaned PDBが必要な場合に限りstandalone `fix-altloc`を使用してください。

1. **活性部位モデル抽出**（`-c/--center` が指定された場合）
 - 基質は PDB/mmCIF パス、残基ID/名、`A:SAM`、または `A:SAM:123` で指定可能
 - 抽出オプション: `--radius`、`--radius-het2het`、`--include-h2o`、`--exclude-backbone`、`--add-linkh`、`--selected-resn`、`--verbose`
 - 入力ごとの内部PDBは `_work/models/` に保存し、mmCIF/oversized-PDB入力では元IDを復元したCIFも生成
 - **最初の活性部位モデルの総電荷**がスキャン/MEP/TSOPT に伝播

2. **オプションの段階的スキャン（単一入力のみ）**
 - 各 `--scan-lists` 引数は MLIP スキャンステージを記述する `(i,j,target_Å)` タプルの Python ライクなリスト。原子インデックスは元の入力順序を参照し、デフォルトでは 1 始まりです（`--no-scan-one-based` を指定すると 0 始まりとして読みます）。いずれの場合も活性部位モデル順序に自動変換されます。3-field selector（例: `'TYR,285,CA'`）はtoken順を問いません。残基名や番号が重複するときは、位置固定の`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`（例: `A:SAM:320:C1`）でchainを明示します。
 - 単一リテラルは 1 ステージスキャンを実行し、複数リテラルは**順次**実行されるため、ステージ 2 はステージ 1 の結果から開始されます。複数リテラルは 1 つの `-s/--scan-lists` に並べて指定します（例: `-s '[(…)]' '[(…)]'`）。
 - ステージエンドポイント（`stage_XX/result.pdb`）が、後続 MEP ステップへ渡される順序付き中間体となる

3. **活性部位モデルでの MEP 探索（デフォルトで単一パス `path-opt`、`--refine-path` で再帰的 `path-search`）**
 - デフォルトでは、単一パス `path-opt`（GSM/DMF）を実行します。エンジン生出力は `<out-dir>/_work/path_opt/` に書かれ、マージ済み成果物（`mep.pdb`、`mep_trj.xyz`、`energy_diagram_MEP.png`）はルート直下へ配置します。
 - `--refine-path` を指定すると、再帰的 `path-search` に切り替わり、結合変化に基づく多段階反応の候補セグメントを構築します。この分割だけで素反応が確定するわけではなく、TS／虚振動／IRC の検証が必要です。粗い MEP から得た HEI で TSOPT が失敗する場合の精密化に有効です。一方、悪い／ノイズの多い path を不要な複数 segment へ分割して計算時間を大幅に増やすことがあるため、意図せぬ cost 増大を避けてデフォルト OFF です（エンジン生出力は `<out-dir>/_work/path_search/`）。
 - 複数入力 PDB の場合、参照マージ用の全系テンプレートが MEP エンジン（デフォルトは `path-opt`、`--refine-path` 時は `path-search`）に自動的に渡されます（全系マージ自体は `--refine-path` 時のみ実行）。単一構造スキャンの場合は、元の全系 PDB テンプレートが全ステージで再利用されます。

4. **活性部位モデルを全系にマージ**（`--refine-path` 使用時のみ）
 - `--refine-path` と参照テンプレートがある場合、`mep_w_ref.pdb` を生成し、bridge入力では `mep_w_ref.cif` も生成します。デフォルトの `path-opt` では全系マージを行いません。

5. **オプションのセグメントごとの後処理**（反応セグメントのみ — 結合変化のあるセグメント。ブリッジセグメントはスキップ）
 - `--tsopt`: 各 HEI 活性部位モデルで TS 最適化を実行します。機械可読な exact-Hessian 結果が `status=converged` かつ `n_imag=1` を報告しない場合は IRC 前に停止します。検証済み TS は EulerPC IRC で追跡し、IRC エンドポイントを `--thresh-post`（デフォルト `baker`）で再最適化します。Hessian TS 最適化には MEP energy-upwinding Cartesian接線を反応参照モードとしてデフォルトで渡します（energyを読めない旧trajectoryだけは正規化secantの二等分線へfallback）。`--no-tsopt-from-mep-tan` では、TSOPT が初期構造の Hessian を計算し、その振動モードから初期モードを選びます。エンドポイント最適化の作業ディレクトリは完了後に自動削除されます。エンドポイント RFO の上り坂拒否はデフォルトで無効で、`--reject-uphill` によりエンドポイント再最適化についてのみ有効化できます。
 - `--thermo`: (R, TS, P) で `freq` を呼び出し、振動/熱化学データと MLIP Gibbs ダイアグラムを取得
 - `--dft`: (R, TS, P) で DFT 一点計算を実行し、DFT ダイアグラムを構築。`--thermo` と組み合わせると DFT//MLIP Gibbs ダイアグラムも生成
  - 共有の上書きオプション: `--opt-mode`、`--opt-mode-post`（TSOPT/IRC 後最適化のプリセット上書き）、`--flatten/--no-flatten`、`--hessian-calc-mode`、`--tsopt-max-cycles`、`--tsopt-out-dir`、`--freq-*`、`--dft-*`、`--dft-engine`（GPU 優先）など。Cartesian PHVA の剛体モードは、凍結anchorを尊重する constrained 処理に固定されています。
 - Hessian 評価モードの詳細は {ref}`ja-hessian-evaluation` を参照してください。

6. **TSOPT のみモード**（単一入力、`--tsopt`、`--scan-lists` なし）
 - MEP/マージステージをスキップし、活性部位モデル（または抽出がスキップされた場合は全入力構造）で `tsopt` → EulerPC IRC を実行し、高エネルギー側の IRC 終端を反応物 (R) として識別したうえで、エネルギーダイアグラム一式とオプションの freq/DFT 出力を生成します。

## 出力

ツリーは 3 つのゾーンで構成されます: **ルート直下の成果物**、**`segments/seg_NN/` 配下のセグメント別成果物**、**`_work/` 配下のパイプライン作業領域**（必要な結果を取り出したあとは `rm -rf` で削除して構いません）。

```text
out_dir/ (デフォルト:./result_all/)
├─ summary.log                  # 結果要約（ルート直下に生成）
├─ summary.json                 # JSON 結果
├─ mep.pdb                      # マージ済み MEP 経路（エンジンから配置）
├─ mep.cif                      # mmCIF/oversized-PDB入力時。元IDを復元
├─ mep_w_ref.pdb               # --refine-path + 参照テンプレート使用時のみ
├─ mep_w_ref.cif               # --refine-path + bridge参照テンプレート使用時のみ
├─ mep_trj.xyz                 # MEP 全体軌道
├─ energy_diagram_MEP.png      # 全セグメントの MEP 障壁
├─ energy_diagram_*.png        # 集約後処理ダイアグラム（MLIP / Gibbs / DFT、--tsopt 等で生成）
├─ segments/                    # 反応セグメント別の成果物（ブリッジセグメントはスキップ）
│  └─ seg_NN/                   # 2 桁インデックス、例: seg_01, seg_02
│     ├─ reactant.{pdb,cif,xyz,gjf} # 正規R/TS/P。bridge入力はPDB+CIF
│     ├─ ts.{pdb,cif,xyz,gjf}
│     ├─ product.{pdb,cif,xyz,gjf}
│     ├─ ts/                    # TS 最適化出力と振動解析（--tsopt）
│     ├─ irc/                   # IRC 軌道とプロット（--tsopt）
│     ├─ freq/{R,TS,P}/         # frequencies_cm-1.txt + thermoanalysis.yaml（--thermo）
│     └─ dft/                   # DFT 一点計算結果（--dft）
└─ _work/                       # パイプライン作業領域（削除可）
   ├─ models/                   # 抽出実行時の活性部位モデル PDB（model_<input_stem>.pdb）
   ├─ scan/                     # 段階的スキャン結果（--scan-lists 提供時）
   ├─ add_elem_info/            # 前処理: 元素記号補完
   └─ path_opt/                 # MEP エンジン生出力（--refine-path 時は path_search/）
```

**TSOPT のみモード**（単一入力 + `--tsopt`、`--scan-lists` なし）では MEP ステージが無く、最適化済み R/TS/P と `ts/`・`irc/`・`freq/`・`dft/` は `segments/seg_01/` 直下に生成され、MEP 作業ディレクトリ（`_work/path_opt/`）は存在しません。

```{note}
**正規構造は `segments/seg_NN/reactant.*`・`ts.*`・`product.*`** です — 機構を報告する際はこれらを引用してください。同じ `seg_NN/` 内の `ts/`・`irc/`・`freq/`・`dft/` サブディレクトリは各ステージの作業ファイル（例: `ts/vib/imag_*_trj.xyz`、`irc/*_trj.xyz`）を保持し、特定ステージのデバッグに使います。`_work/path_opt/` 配下の MEP エンジン生出力は作業領域であり、必要な成果物（`mep.pdb`、bridge入力時の`mep.cif`、`mep_trj.xyz`、`energy_diagram_MEP.png`）は既にルートへ配置済みです。
```

`-v 2` ではコンソールに活性部位モデルの電荷解決結果、YAML 設定、スキャンステージ、MEP（GSM/DMF）の進行状況、各ステージの所要時間が出力されます。詳細は {ref}`ja-verbosity-levels` を参照してください。

### プロットファイルの命名規則

エネルギーダイアグラムファイルは手法とスコープに基づいて命名されます:

| ファイル名 | 生成タイミング | 内容 |
|---|---|---|
| `energy_diagram_MEP.png` | path-opt/path-search 完了時 | 全セグメント MEP 障壁（生の GSM/DMF 値） |
| `energy_diagram_MLIP.png` | セグメントごとの tsopt+IRC 完了時 | R→TS→P（MLIP エネルギー） |
| `energy_diagram_G_MLIP.png` | セグメントごとの thermo 完了時 | R→TS→P（MLIP ギブズ自由エネルギー） |
| `energy_diagram_DFT.png` | セグメントごとの DFT 完了時 | R→TS→P（DFT エネルギー） |
| `energy_diagram_G_DFT_plus_MLIP.png` | セグメントごとの DFT+thermo 完了時 | R→TS→P（DFT エネルギー + MLIP 熱補正） |
| `energy_diagram_MLIP_all.png` | 全セグメント集約時 | 全セグメント統合（MLIP） |
| `energy_diagram_G_MLIP_all.png` | 全セグメント + thermo | 全セグメント統合（MLIP ギブズ） |
| `energy_diagram_DFT_all.png` | 全セグメント + DFT | 全セグメント統合（DFT） |
| `energy_diagram_G_DFT_plus_MLIP_all.png` | 全セグメント + DFT + thermo | 全セグメント統合（DFT//MLIP ギブズ） |
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

- 抽出有効（`-c/--center`）: **PDB/mmCIF** を使用可能。
- 抽出なし: **PDB/mmCIF/XYZ/GJF** を使用可能。
- 複数構造実行は 2 つ以上の構造が必要。完全な入力ファイル要件（水素、元素列、原子順序の一致）は [CLI 規約](cli-conventions.md) を参照してください。

電荷の解決順序の詳細は {ref}`CLI 規約: 電荷の指定 <ja-charge-specification>` を参照してください。`all` コマンドでは、活性部位モデル抽出（`-c` 指定時）による電荷導出が追加の優先度レイヤーとして機能します。スピンの解決順序は `--multiplicity`（CLI）→ `.gjf` テンプレート → デフォルト（1）。非標準の基質には `--ligand-charge/-l` を必ず指定し、scan/MEP/TSOPT/DFT へ正しい総電荷を伝播させてください。

### 入出力オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 反応順序の 2 つ以上の完全構造（`--scan-lists` または `--tsopt` のみ単一入力可） | 必須 |
| `--ref-pdb FILE` | XYZ/GJF入力用の参照 PDB/mmCIF topology | _None_ |
| `-o, --out-dir PATH` | トップレベル出力ディレクトリ | `./result_all/` |
| `--convert-files/--no-convert-files` | XYZ/TRJ → 対応する PDB/CIF/GJF companion の全体切替 | `True` |
| `--dump/--no-dump` | MEP（GSM/DMF）軌跡を出力。親で明示したトグルは `path-search`/`path-opt` と `scan`/`tsopt` に転送され、省略時は各子コマンドの YAML/デフォルトを使用。`--thermo` 使用時は内部入力の `thermoanalysis.yaml` を常に保持し、`--no-dump` でもこの channel は抑止しません | `False` |
| `--config FILE` | 先に適用するベース YAML | _None_ |
| `--show-config/--no-show-config` | 実行前に解決済み設定を表示 | `False` |
| `--dry-run/--no-dry-run` | 設定を検証して計画を表示。`--center` 指定時は一時ディレクトリで extract を実行し、導出電荷と電子数 parity を検証する。scan/MEP/TSOPT/freq/DFT は実行せず、永続出力も作らない | `False` |

### 電荷・スピンオプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-l, --ligand-charge TEXT` | 総電荷または残基別マッピング（`-q` 省略時に使用、推奨）。PDB/mmCIF 入力（または `--ref-pdb` 付き XYZ/GJF）で extract と同じ全系電荷導出を起動します | _None_ |
| `-q, --charge INT` | 明示した総電荷を最優先。不一致時は警告して`-q`を使用し、省略時は抽出／workflowの自動導出を使用 | _None_ |
| `-m, --multiplicity INT` | 全下流ステップへ転送されるスピン多重度 | `1` |

### 活性部位モデル抽出オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-c, --center TEXT` | PDB/mmCIFパス、残基ID/名、`CHAIN:RESNAME`、`CHAIN:RESNAME:RESSEQ` | 抽出に必須 |
| `-r, --radius FLOAT` | 活性部位モデル包含カットオフ（Å） | `2.6` |
| `--radius-het2het FLOAT` | ヘテロ–ヘテロカットオフ（Å）。`0` を渡すと空の選択を避けるため内部で `0.001 Å` に自動補正されます（単体の `extract` と同じ挙動） | `0.0` |
| `--include-h2o/--no-include-h2o` | 水分子を含める（HOH/WAT/TIP3/SOL） | `True` |
| `--exclude-backbone/--no-exclude-backbone` | 非基質アミノ酸の主鎖原子を除去 | `False` |
| `--add-linkh/--no-add-linkh` | 切断結合にキャップ水素を付加 | `True` |
| `--selected-resn TEXT` | `--center` と同じID/名前/chain付きselectorで残基を強制包含 | `""` |
| `--modified-residue TEXT` | アミノ酸として扱う残基名をカンマ区切りで指定。`NAME:charge` はこの抽出中の公称電荷を追加または上書きし、電荷を省略した `NAME` は 0 になります | `""` |
| `--freeze-links/--no-freeze-links` | 活性部位モデル PDB でキャップ H の親を凍結 | `True` |
| `--freeze-atoms TEXT` | 全ステージで凍結する 1-based 原子番号（カンマ区切り）。番号は抽出後モデル基準。`--freeze-links` と YAML `geom.freeze_atoms` にマージされる | _None_ |

組み込みのアミノ酸名は、力場で正規化された Amber/CHARMM の命名を前提とします。
raw PDB CCD との名前衝突は自動判別しないため、`--modified-residue NAME:charge`
で意図する公称電荷を明示してください。

(ja-mep-search-options)=
### MEP 探索オプション

```{note}
`pdb2reaction all` では `--max-cycles` を指定せず、各ステージ固有の
デフォルトを使用してください。`--max-cycles` は `opt`、`tsopt`、
`path-opt` などの単発サブコマンドを直接実行するときだけ指定します。
```

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--mep-mode [gsm\|dmf]` | MEP 探索アルゴリズム: GSM（Growing String Method）または DMF（Direct Max Flux） | `gsm` |
| `--max-nodes INT` | GSM/DMF segment ごとの可動内部イメージ数。両エンジンとも端点2つを保持するため、総イメージ数は `max_nodes + 2` | `20` |
| `--max-cycles INT` | MEP 最大最適化サイクル | `300` |
| `--climb/--no-climb` | 標準 GSM セグメントでクライミングイメージを有効化（ブリッジセグメントは常に無効） | `True` |
| `--opt-mode [grad\|hess]` | ワークフロープリセット（`grad` → L-BFGS/Dimer、`hess` → RFO/RSPRFO）。コマンド個別実行では `opt --opt-mode grad|hess`、`tsopt --opt-mode grad|hess` を推奨。トークンのマッピングはスコープ依存で、`all` の pre-opt デフォルト（`grad`）と `tsopt` のデフォルト（`hess`）は一致しません。詳細は {ref}`ja-opt-mode-semantics` を参照してください | `grad` |
| `--thresh TEXT` | 単一構造最適化と scan 緩和の収束プリセット（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `gau` |
| `--thresh-gsm TEXT` | MEP 段の GSM ストリング最適化の収束プリセット（`--thresh` と同じプリセット群） | `gau_loose` |
| `--thresh-dmf TEXT` | DMF MEP 段の IPOPT dual-infeasibility 許容値。`tight`(0.04)、`middle`(0.10)、`loose`(0.20) または正の float。Gaussian プリセットではない | `tight` |
| `--preopt/--no-preopt` | MEP 前に活性部位モデル端点を事前最適化。単体の `scan`、`scan2d`、`scan3d` では `--preopt` のデフォルトは `False`（`--preopt` を渡すと有効化） | `True` |
| `--refine-path / --no-refine-path` | 再帰的 `path-search` を有効化 / デフォルトの単一パス `path-opt` を使用 | 無効 |

### MLIP 計算機オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--workers`, `--workers-per-node` | UMA 予測器の並列度。`workers > 1` と明示的な解析 Hessian は併用できないため、`workers = 1` または有限差分を使用。診断上の注意は {ref}`ja-workers-analytical-error` を参照 | `1`, `1` |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | 共有 MLIP Hessian エンジン | `FiniteDifference` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--solvent TEXT` | xTB 補正用の暗黙溶媒名（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |

### 後処理オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--tsopt/--no-tsopt` | セグメントごとの TS 最適化+ IRC を実行 | `False` |
| `--tsopt-from-mep-tan/--no-tsopt-from-mep-tan` | HEI の MEP 接線から初期 TS root を選ぶ。OFF では初期構造の Hessian 振動モードから選ぶ | `True` |
| `--thermo/--no-thermo` | R/TS/P で振動解析を実行（`--tsopt` が必要） | `False` |
| `--dft/--no-dft` | R/TS/P で DFT 一点計算を実行（`--tsopt` が必要） | `False` |
| `--opt-mode-post [grad\|hess]` | TSOPT/IRC 後最適化のプリセット上書き（`grad` → Dimer/L-BFGS、`hess` → RSPRFO/RFO） | `hess` |
| `--thresh-post TEXT` | IRC 後エンドポイント最適化の収束プリセット（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `baker` |
| `--flatten/--no-flatten` | 余分な虚振動モードのフラット化 | `False` |
| `--reject-uphill/--no-reject-uphill` | IRC 後の**エンドポイント再最適化のみ**で RFO の上り坂ステップ拒否を明示的に有効化（許容値 `1e-4` Hartree、低エネルギー形状へロールバックして trust radius を縮小）。TS 最適化では拒否を常に無効化し、経路探索には影響しない。emergency floor 到達時は、保持したエンドポイントを通常の収束条件で最終確認 | `False` |
| `--irc-step-size FLOAT` | IRC のEulerPC最大step（Bohr）を上書き。数frameですぐ止まる場合は`0.05`など小さい値で再試行 | IRC デフォルト`0.10` |
| `--irc-never-stop/--no-irc-never-stop` | IRCのgradient・energy停止条件を無視し、各branchを最大cycleまで追跡。数値／integration失敗や外部中断では停止 | `False` |

```{warning}
`--dft` のcost/memoryはbasis-function数、元素、functional、grid、engine、
hardwareに依存し、atom countだけではcutoffを決められません。代表stateをpilot実行し、
peak memoryを測定してresourceを選択してください。
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
| `--scan-one-based/--no-scan-one-based` | `--scan-lists` の原子インデックスの読み方: `True` = 1 始まり、`False` = 0 始まり | _None_（1 始まり） |
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
| [`path-opt`](path-opt.md) | `geom`, `calc`, `gs`, `dmf`, `stopt`, `opt`, `lbfgs`, `rfo` |
| [`path-search`](path-search.md) | `geom`, `calc`, `gs`, `dmf`, `stopt`, `opt`, `lbfgs`, `rfo`, `bond`, `search` |
| [`scan`](scan.md) | `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`, `bond` |
| [`tsopt`](tsopt.md) | `geom`, `calc`, `opt`, `hessian_dimer`, `rsirfo` |
| [`freq`](freq.md) | `geom`, `calc`, `freq`, `thermo` |
| [`dft`](dft.md) | `dft` |
| [`irc`](irc.md) | `geom`, `calc`, `irc` |

**最小例:**
```yaml
calc:
 model: uma-s-1p2 # uma-s-1p2 | uma-m-1p1
 hessian_calc_mode: FiniteDifference # デフォルト。Analytical は対象環境で検証して選択
gs:
 max_nodes: 12
 climb: true
dft:
 grid_level: 6
```

すべての YAML オプションの完全なリファレンスについては、**[YAML 設定リファレンス](yaml-reference.md)** を参照してください。

## 注記

独立に準備した全系構造どうしでは、反応座標以外の構造差が得られる障壁に影響することがあります。構造を確認し、モデル化した系に対して選択した path workflow を検証してください。

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。
- 形式電荷を推定できない場合は `--ligand-charge`（数値または残基別マッピング）を必ず指定し、scan/MEP/TSOPT/DFT へ正しい総電荷を伝播させてください。
- マージ用の参照 PDB テンプレートは元の入力から自動的に導出されます。`path-search` の `--ref-full-pdb` はこのラッパーでは意図的に非公開です。
- 収束プリセット: `--thresh` のデフォルトは `gau`、`--thresh-post` のデフォルトは `baker`、MEP 段は `--thresh-gsm`（デフォルト `gau_loose`）と `--thresh-dmf`（デフォルト `tight`）が担当。
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
