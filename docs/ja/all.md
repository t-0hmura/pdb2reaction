# `all`

## 概要

`pdb2reaction all` は、抽出から解析までの一連の処理を **まとめて実行する最上位コマンド** です。典型的なフローは次のとおりです。

活性部位モデル（バインディングポケット）抽出 →（任意）段階的スキャン → MEP 探索（デフォルトで再帰的 `path-search`）→（任意）TS 最適化（`tsopt`）+ IRC →（任意）振動解析・熱化学（`freq`）→（任意）DFT 一点計算（`dft`）。`--refine-path False` を指定すると単一パス `path-opt`（GSM/DMF）に切り替わります。

MLIP バックエンドはデフォルトで UMA を使用しますが、`-b/--backend` オプションで ORB・MACE・AIMNet2 も選択可能です。

```{important}
`--tsopt` **なし**の `all` ワークフローは **TS 候補**（MEP 探索の最高エネルギー画像 / HEI）を出力します。`--tsopt` を追加すると、これらを虚振動数チェックで検証済みの最適化 TS 構造に精密化し、続いて IRC でエンドポイントを検証します。結果の虚振動モードと端点極小は必ず目視で確認してください。
```

### 要点
- **想定場面:** PDB からの end-to-end 一括実行（活性部位モデル抽出 → MEP → TS 最適化 → IRC → 熱化学 → DFT）。
- **手法:** 入力とフラグに応じて 3 つのモードを切り替え（複数構造 MEP / 単一構造 + 段階的スキャン / TSOPT のみ）。
- **主な出力:** `summary.log`、`summary.json`、`path_search/mep.pdb`（`--refine-path False` 時は `path_opt/`）。`--tsopt`/`--thermo`/`--dft` 有効時はセグメント別 `seg_XX/` と後処理ワークスペース `path_search/post_seg_XX/` も生成。
- **デフォルト値:** バックエンド `uma`、`--mep-mode gsm`、`--opt-mode grad`、`--refine-path True`、`--preopt True`、`--thresh gau`、`--thresh-post baker`。`--tsopt`/`--thermo`/`--dft` はいずれも `False`。
- **次にやること:** `--tsopt` なしでは結果は **TS 候補**（HEI）です。`--tsopt`（虚振動数チェック）+ IRC で検証し、必要に応じて `--thermo`・`--dft` を追加してください。

## ワークフローの全体像

多くのケースでは次の流れになります。

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

各ステージはサブコマンドとして単独実行できます。

主なモードは 3 つあります。

- **end-to-end（複数構造）** — 反応順に並べた 2 構造以上（PDB/GJF/XYZ）と基質定義を与えます。`all` が活性部位モデル抽出 → GSM/DMF による MEP 探索 → 全系テンプレートへのマージまで実行し、必要に応じてセグメントごとに TSOPT / freq / DFT を行います。
- **単一構造 + 段階的スキャン** — 1 つの構造に対して `-s/--scan-lists` を 1 つ以上与えます。（段階的）スキャンで得られた中間体列を MEP の端点として用います。

  - `-s/--scan-lists` を 1 つだけ渡すと 1 ステージになります。
  - 複数ステージは 1 つの `-s/--scan-lists` に複数リテラルを並べて指定します（例: `-s '[(…)]' '[(…)]'`）。

- **TSOPT のみ（活性部位モデル TS 最適化）** — 1 つの入力構造に対し、`--scan-lists` を省略して `--tsopt` を指定します。`-c/--center` がある場合は活性部位モデルを抽出し、その系で TS 最適化+ IRC（必要に応じて freq / DFT）のみ実行します。

```{tip}
大規模な活性部位モデルでは、単一構造スキャンワークフロー（`--scan-lists`）の方が、複数構造 MEP ワークフローよりも信頼性の高い反応障壁を得やすい傾向があります。複数の PDB を入力すると、反応座標と無関係な領域の構造差異が蓄積し、障壁を過大評価する可能性があります。スキャンワークフローは単一構造から出発して関連する座標のみを駆動するため、無関係な構造ノイズを最小化できます。この影響はモデルサイズが大きくなるほど顕著になります。
```

> **実行例:** [`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples) ディレクトリに GPP C6-メチル基転移酵素 BezA（[Tsutsumi et al., *Angew. Chem. Int. Ed.* 2022, 61, e202111217](https://doi.org/10.1002/anie.202111217)）の完全な `all` ワークフロースクリプト（MEP およびスキャンパイプライン）があります。

## 最小例

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c "SAM,GPP,MG" -l "SAM:1,GPP:-3" \
 --out-dir ./result_all
```

## 出力の見方

- `result_all/summary.log`
- `result_all/summary.json`
- `result_all/path_search/mep.pdb`（`--refine-path False` 使用時は `result_all/path_opt/`）

## よくある例

1. TS 最適化・IRC・熱化学・DFT まで一括実行する。

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c "SAM,GPP,MG" -l "SAM:1,GPP:-3" \
 --tsopt --thermo --dft --out-dir ./result_mep
```

2. 単一構造 + 段階的スキャンを実行する。

```bash
pdb2reaction all -i 1.R.pdb -c "SAM,GPP,MG" -l "SAM:1,GPP:-3" \
 -s '[("CS1 SAM 320","GPP 321 C7",1.60)]' '[("GPP 321 H11","GLU 186 OE2",0.90)]' \
 --tsopt --thermo --out-dir ./result_scan
```

テンプレートがある場合の XYZ/TRJ → PDB/GJF 変換（付随ファイルの生成）は、全ステージ共通の `--convert-files/--no-convert-files`（デフォルト: `True`）で制御できます。


## 使用法
```bash
pdb2reaction all -i INPUT1 [INPUT2 ...] -c SUBSTRATE [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] [options]
```

ヘルプ出力は `pdb2reaction all --help` で主要オプションを、`pdb2reaction all --help-advanced` で全オプションを確認できます。

### 例
```bash
# 明示的なリガンド電荷と後処理を伴う複数構造 MEP
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' \
 -l 'SAM:1,GPP:-3' --multiplicity 1 --freeze-links \
 --max-nodes 10 --max-cycles 100 --climb --opt-mode grad \
 --out-dir ./result_mep --tsopt --thermo --dft

# 単一構造段階的スキャン + GSM/DMF + TSOPT/freq/DFT
pdb2reaction all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 -s '[("CS1 SAM 320","GPP 321 C7",1.60)]' '[("GPP 321 H11","GLU 186 OE2",0.90)]' \
 --opt-mode hess --tsopt --thermo --dft

# TSOPT のみワークフロー（経路探索なし）
pdb2reaction all -i TS_candidate.pdb -c 'SAM,GPP,MG' \
 -l 'SAM:1,GPP:-3' --tsopt --thermo --dft
```

## ワークフロー

0. **事前チェック**（自動）
 - `all` は全 PDB 入力に対し、他の処理の前に `add-elem-info`（PDB 列 77–78 の元素記号を補完）と `fix-altloc`（代替配座の解消）を自動実行します。個別のサブコマンド（`extract`、`opt` など）を使用する場合は、必要に応じてこれらを手動で実行してください。

1. **活性部位モデル抽出**（`-c/--center` が指定された場合）
 - 基質は PDB パス、残基 ID（`123,124` または `A:123,B:456`）、または残基名（`GPP,SAM`）で指定可能
 - 抽出オプション: `--radius`、`--radius-het2het`、`--include-h2o`、`--exclude-backbone`、`--add-linkh`、`--selected-resn`、`--verbose`
 - 入力ごとの活性部位モデル PDB は `<out-dir>/models/` に保存。複数構造が提供された場合、活性部位モデルは残基選択ごとに統合
 - **最初の活性部位モデルの総電荷**がスキャン/MEP/TSOPT に伝播

2. **オプションの段階的スキャン（単一入力のみ）**
 - 各 `--scan-lists` 引数は MLIP スキャンステージを記述する `(i,j,target_Å)` タプルの Python ライクなリスト。原子インデックスは元の入力順序（1始まり）を参照し、活性部位モデル順序に自動変換されます。PDB 入力の場合、`i`/`j` は整数インデックスまたは `'TYR,285,CA'` のようなセレクタ文字列を使用可能です。セレクタはスペース/カンマ/スラッシュ/バッククォート/バックスラッシュ（` ` `,` `/` `` ` `` `\`）を区切り文字として受け付け、トークン順序は任意です（フォールバックは resname, resseq, atom と仮定）。
 - 単一リテラルは 1 ステージスキャンを実行し、複数リテラルは**順次**実行されるため、ステージ 2 はステージ 1 の結果から開始されます。複数リテラルは 1 つの `-s/--scan-lists` に並べて指定します（例: `-s '[(…)]' '[(…)]'`）。
 - ステージエンドポイント（`stage_XX/result.pdb`）が、後続 MEP ステップへ渡される順序付き中間体となる

3. **活性部位モデルでの MEP 探索（再帰的 `path-search`）**
 - デフォルトでは、再帰的 `path-search` を実行し、多段階反応を自動検出して各素反応の詳細な MEP を構築します（出力は `<out-dir>/path_search/`）。複雑な多段階反応では手動での試行錯誤が必要な場合があります。
 - `--refine-path False` を指定すると、単一パス `path-opt` GSM/DMF に切り替わります（出力は `<out-dir>/path_opt/`）。
 - 複数入力 PDB の場合、全系テンプレートが参照マージ用に `path-search` に自動的に渡されます。単一構造スキャンの場合は、元の全系 PDB テンプレートが全ステージで再利用されます。

4. **活性部位モデルを全系にマージ**（デフォルトの `--refine-path` 使用時）
 - `--refine-path` が True（デフォルト）かつ参照 PDB テンプレートがある場合、マージ済みの `mep_w_ref*.pdb` とセグメントごとの `mep_w_ref_seg_XX.pdb` が `<out-dir>/path_search/` に出力されます。`--refine-path False`（`path-opt` モード）では全系マージは行われません。

5. **オプションのセグメントごとの後処理**（反応セグメントのみ — 結合変化のあるセグメント。ブリッジセグメントはスキップ）
 - `--tsopt`: 各 HEI 活性部位モデルで TS 最適化を実行し、EulerPC IRC で追跡した後、IRC エンドポイントを `--thresh-post`（デフォルト `baker`）で再最適化します。エンドポイント最適化の作業ディレクトリは完了後に自動削除されます。
 - `--thermo`: (R, TS, P) で `freq` を呼び出し、振動/熱化学データと MLIP Gibbs ダイアグラムを取得
 - `--dft`: (R, TS, P) で DFT 一点計算を実行し、DFT ダイアグラムを構築。`--thermo` と組み合わせると DFT//MLIP Gibbs ダイアグラムも生成
  - 共有の上書きオプション: `--opt-mode`、`--opt-mode-post`（TSOPT/IRC 後最適化のプリセット上書き）、`--flatten/--no-flatten`、`--hessian-calc-mode`、`--tsopt-max-cycles`、`--tsopt-out-dir`、`--freq-*`、`--dft-*`、`--dft-engine`（GPU 優先）など
 - ヘシアン評価モードの詳細は {ref}`ja-hessian-evaluation` を参照してください。

6. **TSOPT のみモード**（単一入力、`--tsopt`、`--scan-lists` なし）
 - MEP/マージステージをスキップし、活性部位モデル（または抽出がスキップされた場合は全入力構造）で `tsopt`→ EulerPC IRC を実行
 - 高エネルギー側の IRC 終端を反応物 (R) として識別し、エネルギーダイアグラム一式とオプションの freq/DFT 出力を生成


### 電荷とスピンの優先順位

電荷の解決順序の詳細は {ref}`CLI 規約: 電荷の指定 <ja-charge-specification>` を参照してください。`all` コマンドでは、活性部位モデル抽出（`-c` 指定時）による電荷導出が追加の優先度レイヤーとして機能します。

**スピンの解決順序:** `--multiplicity`（CLI）→ `.gjf` テンプレート → デフォルト（1）

> **ヒント:** 非標準の基質には `--ligand-charge` を必ず指定し、正しい電荷伝播を確保してください。

### 入力要件
- 抽出有効（`-c/--center`）: 残基同定のため入力は **PDB** が必須。
- 抽出なし: **PDB/XYZ/GJF** を使用可能。
- 複数構造実行は 2 つ以上の構造が必要。


## CLI オプション

> **注記:** 表示されているデフォルト値は、オプション未指定時に使用される内部デフォルトです。

### 入出力オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 反応順序の2つ以上の完全構造（`--scan-lists` または `--tsopt` のみ単一入力可） | 必須 |
| `--ref-pdb FILE` | `-i` で XYZ 入力を使用する場合のトポロジー参照 PDB | _None_ |
| `-o, --out-dir PATH` | トップレベル出力ディレクトリ | `./result_all/` |
| `--convert-files/--no-convert-files` | XYZ/TRJ → PDB/GJFコンパニオンの全体切替 | `True` |
| `--dump/--no-dump` | MEP(GSM/DMF)軌跡を出力。`path-search`/`path-opt` には常時転送され、`scan`/`tsopt` には明示指定時のみ転送。`freq` はデフォルトで dump=True なので `--no-dump` で無効化。 | `False` |
| `--config FILE` | 先に適用するベース YAML | _None_ |
| `--show-config/--no-show-config` | 実行前に解決済み設定を表示 | `False` |
| `--dry-run/--no-dry-run` | 実行せず検証と計画表示のみ行う。 | `False` |
| `--resume/--no-resume` | `--out-dir` から前回の実行を再開。出力ファイルが既に存在する完了済みステージはスキップされる。 | `False` |

### 電荷・スピンオプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-l, --ligand-charge TEXT` | 残基別電荷マッピング（例: `GPP:-3,SAM:1`、推奨）。PDB の残基電荷から全系の電荷を自動導出します（手動計算不要）。`-q` 省略時に使用（PDB 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `-q, --charge INT` | 総電荷を強制上書き（`--ligand-charge` より優先） | _None_ |
| `-m, --multiplicity INT` | 全下流ステップへ転送されるスピン多重度 | `1` |

### 活性部位モデル抽出オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-c, --center TEXT` | 基質指定（PDBパス、残基ID、または残基名） | 抽出に必須 |
| `-r, --radius FLOAT` | 活性部位モデル包含カットオフ（Å） | `2.6` |
| `--radius-het2het FLOAT` | ヘテロ–ヘテロカットオフ（Å）。`0` を渡すと空の選択を避けるため内部で `0.001 Å` に自動補正されます（単体の `extract` と同じ挙動）。 | `0.0` |
| `--include-h2o/--no-include-h2o` | 水分子を含める（HOH/WAT/TIP3/SOL） | `True` |
| `--exclude-backbone/--no-exclude-backbone` | 非基質アミノ酸の主鎖原子を除去 | `False` |
| `--add-linkh/--no-add-linkh` | 切断結合にリンク水素を付加 | `True` |
| `--selected-resn TEXT` | 強制包含残基。**名前とは裏腹にこのフラグは残基 ID（コロン区切り整数、オプションでチェーン/挿入コード付き、例 `A:123A`）を受け付け、3 文字残基名は受け付けません。** 残基名ベースの選択には `-c/--center 'GPP,SAM'` を使用してください。 | `""` |
| `--modified-residue TEXT` | 修飾アミノ酸残基名をカンマ区切りで指定（任意で電荷付き）。主鎖切断と電荷計算にアミノ酸として扱う。例: `HD1,HD2,HD3` または `HD1:0,SEP:-2`。 | `""` |
| `--freeze-links/--no-freeze-links` | 活性部位モデルPDBでリンクHの親を凍結 | `True` |
| `--verbose/--no-verbose` | 抽出器のINFOログを有効化 | `True` |

(mep-search-options)=
### MEP 探索オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--mep-mode [gsm\|dmf]` | MEP 探索アルゴリズム: GSM（Growing String Method）または DMF（Direct Max Flux） | `gsm` |
| `--max-nodes INT` | MEP 内部ノード数。**GSM:** 総イメージ数 = `max_nodes + 2`（端点 2 つは固定）。**DMF:** チェーン上の *可動* イメージ数（端点の暗黙的拡張なし）。 | `20` |
| `--max-cycles INT` | MEP最大最適化サイクル | `300` |
| `--climb/--no-climb` | 最初のセグメントでTSクライミングを有効化 | `True` |
| `--opt-mode [grad\|hess]` | ワークフロープリセット（`grad` → LBFGS/Dimer、`hess` → RFO/RSIRFO）。コマンド個別実行では `opt --opt-mode grad|hess`、`tsopt --opt-mode grad|hess` を推奨。トークンのマッピングはスコープ依存で、`all` の pre-opt デフォルト（`grad`）と `tsopt` のデフォルト（`hess`）は一致しません。詳細は {ref}`ja-opt-mode-semantics` を参照してください。 | `grad` |
| `--thresh TEXT` | 収束プリセット（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `gau` |
| `--preopt/--no-preopt` | MEP前に活性部位モデル端点を事前最適化。**注意:** `all` はここで子サブコマンドのデフォルトを上書きします。単体の `path-search`、`path-opt`、`scan`、`scan2d`、`scan3d` では `--preopt` のデフォルトは `False` です。 | `True` |
| `--refine-path BOOL` | `True`（デフォルト）の場合は再帰的 `path-search`、`False` の場合は `path-opt` を連結して再帰的精密化なしで実行。無効化するには `--refine-path False` を渡してください（`click.BOOL` パラメータのため `--no-refine-path` 形式は自動生成されません）。 | `True` |

### MLIP 計算機オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--workers`, `--workers-per-node` | MLIP 予測器の並列度（workers > 1 で解析ヘシアン無効; UMA バックエンドのみ）。診断上の注意は {ref}`ja-workers-fd-downgrade` を参照。 | `1`, `1` |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | 共有 MLIP ヘシアンエンジン | `FiniteDifference` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--solvent TEXT` | xTB 補正用の暗黙溶媒名（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |

### 後処理オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--tsopt/--no-tsopt` | セグメントごとの TS 最適化+ IRC を実行 | `False` |
| `--thermo/--no-thermo` | R/TS/Pで振動解析を実行 | `False` |
| `--dft/--no-dft` | R/TS/PでDFT一点計算を実行 | `False` |
| `--opt-mode-post [grad\|hess]` | TSOPT/IRC後最適化のプリセット上書き（`grad` → Dimer/LBFGS、`hess` → RSIRFO/RFO） | `hess` |
| `--thresh-post TEXT` | IRC後エンドポイント最適化の収束プリセット（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `baker` |
| `--flatten/--no-flatten` | 余分な虚振動モードのフラット化 | `False` |

```{warning}
`--dft` による DFT 一点計算（PySCF/GPU4PySCF）は、約 300 原子を超えるモデルでは計算コストが非常に大きくなります。そのような系では、A100 や H200 等の高性能 GPU を搭載した HPC クラスタの利用が必要になる場合があります。
```

TSOPT の最適化モードは、`--opt-mode-post`（指定時）→ `--opt-mode`（明示指定時のみ）→ TSOPT のデフォルト（`hess` → `rsirfo`）の順で決まります。

例: `--opt-mode grad --opt-mode-post hess` は、経路最適化に LBFGS、TS 精密化に RS-I-RFO を使用します。

### TSOPT 上書き

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--tsopt-max-cycles INT` | `tsopt --max-cycles` 上書き | `10000` |
| `--tsopt-out-dir PATH` | tsopt出力サブディレクトリ | _None_ |

### Freq 上書き

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--freq-out-dir PATH` | freq出力ディレクトリ上書き | _None_ |
| `--freq-max-write INT` | 最大モード出力数 | `10` |
| `--freq-amplitude-ang FLOAT` | モードアニメーション振幅（Å） | `0.8` |
| `--freq-n-frames INT` | モードアニメーションフレーム数 | `20` |
| `--freq-sort [value\|abs]` | モードソート方法 | `value` |
| `--freq-temperature FLOAT` | 熱化学温度（K） | `298.15` |
| `--freq-pressure FLOAT` | 熱化学圧力（atm） | `1.0` |

### DFT 上書き

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--dft-engine [gpu\|cpu]` | DFTバックエンド: gpu (GPU4PySCF) または cpu (PySCF)。`all` ラッパーではプレフィックス付きで `--dft-engine` と名付けられていますが、単体の `dft` サブコマンドでは同じオプションが `--engine` という名前になります。 | `gpu` |
| `--dft-out-dir PATH` | DFT出力ディレクトリ上書き | _None_ |
| `--dft-func-basis TEXT` | 汎関数/基底関数ペア | `wb97m-v/def2-tzvpd` |
| `--dft-max-cycle INT` | 最大SCFサイクル | `100` |
| `--dft-conv-tol FLOAT` | SCF収束閾値 | `1e-9` |
| `--dft-grid-level INT` | PySCFグリッドレベル | `3` |

(ja-scan-options-single-input-runs)=
### スキャンオプション（単一入力）

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-s, --scan-lists TEXT...` | 段階的スキャン: `(i,j,target_Å)` タプル | _None_ |
| `--scan-out-dir PATH` | scan出力ディレクトリ上書き | _None_ |
| `--scan-one-based/--no-scan-one-based` | 1始まり/0始まりインデックス | _None_ |
| `--scan-max-step-size FLOAT` | 最大ステップサイズ（Å） | `0.20` |
| `--scan-bias-k FLOAT` | 調和バイアス強度（eV·Å⁻²） | `300` |
| `--scan-relax-max-cycles INT` | 緩和サイクル上限 | `10000` |
| `--scan-preopt/--no-scan-preopt` | scan事前最適化 | _None_ |
| `--scan-endopt/--no-scan-endopt` | scanステージ終端最適化 | _None_ |

## 出力
```text
out_dir/ (デフォルト:./result_all/)
├─ summary.log                  # 結果要約
├─ summary.json                 # JSON 結果
├─ models/                      # 抽出実行時の活性部位モデル PDB
│  └─ model_<input_stem>.pdb   #   -i 入力ごとに 1 つ
├─ scan/                        # 段階的スキャン結果（--scan-lists 提供時）
├─ seg_XX/                      # セグメント XX の TS 最適化・IRC 後の構造
│  ├─ reactant.{pdb,xyz,gjf}   #   出力形式は入力形式と一致
│  ├─ ts.{pdb,xyz,gjf}
│  └─ product.{pdb,xyz,gjf}
├─ path_search/                 # MEP 結果（再帰的 path-search、デフォルト）; --refine-path False 時は path_opt/
│  └─ post_seg_XX/              # 後処理: TS 最適化、IRC、freq、DFT
│     ├─ structures/            # IRC 端点の最適化済み R/TS/P 構造
│     ├─ irc/                   # IRC 軌道とプロット
│     ├─ ts/                    # TS 最適化出力と振動解析
│     ├─ freq/{R,TS,P}/         # 状態別 frequencies_cm-1.txt + thermoanalysis.yaml
│     └─ dft/                   # DFT 一点計算結果（--dft 有効時）
└─ tsopt_single/                # TSOPT のみ出力（IRC エンドポイント）
   ├─ structures/               # R/TS/P の作業コピー（正規版は top-level seg_NN/ 配下）
   ├─ ts/                       # TS 最適化出力と振動解析
   ├─ irc/                      # IRC 軌道とプロット
   ├─ freq/{R,TS,P}/            # 状態別 frequencies_cm-1.txt + thermoanalysis.yaml
   └─ dft/                      # DFT 一点計算結果（--dft 有効時）
```

```{note}
**`seg_XX/` と `path_search/post_seg_XX/` の違い。** 2 つのセグメント別ツリーは重複ではなく、それぞれ異なる役割を持ちます:

- **`seg_XX/`**（`out_dir` 直下）は *セグメントごとの統合結果* であり、反応セグメント `XX` の最終的な reactant/TS/product 構造を集約したものです。後処理パイプラインが完了した後に書き出され、機構を報告する際に引用すべき正規の `reactant.{pdb,xyz,gjf}`、`ts.{pdb,xyz,gjf}`、`product.{pdb,xyz,gjf}` を保持します。
- **`path_search/post_seg_XX/`** は *セグメントごとの後処理ワークスペース* です。各ステージの中間生成物（`ts/`、`irc/`、`structures/`、`freq/`、`dft/`）を保存し、特定のステージをデバッグしたい場合（例えば `ts/vib/imag_*_trj.xyz` や `irc/*_trj.xyz` を確認したい場合）に参照します。反応セグメントのみ生成され、結合変化のないブリッジセグメントはスキップされます。

`--refine-path False` を指定した場合、ワークスペースは `path_opt/post_seg_XX/` 配下に移動します。
```

- コンソールには電荷解決結果、YAML 設定、MEP 進行状況、各ステージの所要時間が出力されます。

### エネルギーダイアグラムの命名規則

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

### `summary.log` の読み方
ログは番号付きセクションで構成されます:
- **[1] グローバル MEP 概要** – イメージ/セグメント数、MEP 軌跡プロットのパス、MEP 全体のエネルギーダイアグラム。
- **[2] セグメント別MEPサマリー（MLIPパス）** – セグメントごとの障壁（`ΔE‡`）、反応エネルギー（`ΔE`）、結合変化サマリー。
- **[3] セグメント別後処理（TSOPT / Thermo / DFT）** – TS 虚振動数チェック、IRC 出力、MLIP/熱化学/DFT のエネルギーテーブル。
- **[4] エネルギーダイアグラム（概要）** – MEP/MLIP/Gibbs/DFT 系の図表と、任意の横断サマリー表。
- **[5] 出力ディレクトリ構造** – 生成ファイルを注釈付きでまとめたツリー。

### `summary.json` の読み方
JSON 結果の代表的なトップレベルキーは以下のとおりです。
- `out_dir`, `n_images`, `n_segments` – 実行メタデータと総数。
- `segments` – `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, `bond_changes` を含むセグメント配列。
- `energy_diagrams`（任意） – `labels`, `energies_kcal`, `energies_au`, `ylabel`, `image` などを含む図表データ。

`summary.json` には `summary.log` にある整形テーブルやファイルツリーは含まれません。

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- 形式電荷を推定できない場合は `--ligand-charge`（数値または残基別マッピング）を必ず指定し、scan/MEP/TSOPT/DFT へ正しい総電荷を伝播させてください。
- マージ用の参照 PDB テンプレートは元の入力から自動的に導出されます。`path-search` の `--ref-full-pdb` はこのラッパーでは意図的に非公開です。
- 収束プリセット: `--thresh` のデフォルトは `gau`、`--thresh-post` のデフォルトは `baker`。
- 抽出半径: `--radius` または `--radius-het2het` に `0` を渡すと、内部で `0.001 Å` にクランプされます。
- エネルギーダイアグラムは反応物（最初の状態）基準の kcal/mol で表示されます。
- `-c/--center` を省略すると抽出をスキップし、全構造をそのまま MEP/tsopt/freq/DFT に渡します。ただし単一構造実行では `--scan-lists` か `--tsopt` が必要です。
- **`--resume`**: 同じコマンドに `--resume` を付けて再実行すると、出力ファイルが既に存在するステージをスキップします。各ステージはセンチネルファイルで判定されます（MEP は `summary.json`、TSOPT/IRC は `final_geometry.*` + `finished_irc_trj.xyz`、freq/DFT は `R/`+`TS/`+`P/` ディレクトリ）。resume 時に抽出がスキップされた場合は `-q/--charge` または `--ligand-charge` を明示的に指定してください。**センチネル破損時の注意:** `--resume` はセンチネルファイルの *存在* のみを確認し、整合性は検証しません。ステージが書き込み途中でキル（SIGKILL、OOM、クラスタのプリエンプション）され、センチネルファイルが中途半端な状態でディスクに残った場合でも、`--resume` はそのステージを完了済みとみなします。部分書き込みが疑われる場合は、再開前に該当ステージの出力ディレクトリ（例: `path_search/post_seg_XX/ts/`）を削除してください。


`all` は YAML の多層指定をサポートします:

- `--config FILE`: ベース設定。

適用順序:

`defaults < config < CLI`

解決後の YAML は呼び出されるすべてのサブコマンドに転送されます。各ツールが読み取るセクションは以下のとおりです:

| サブコマンド | YAML セクション |
|------------|-----------------|
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
- [YAML リファレンス](yaml-reference.md) — 全YAML設定オプション
- [用語集](glossary.md) — MEP、TS、IRC、GSM、DMFの定義
