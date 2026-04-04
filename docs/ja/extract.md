# `extract`

## 概要

> **要約:** タンパク質–リガンド PDB からクラスターモデル（活性部位モデル（バインディングポケット））を抽出します。基質は `-c` で指定します（残基名、残基 ID、または PDB パス）。切断された結合はリンク水素でキャップされます。非標準残基の電荷には `--ligand-charge` を使用してください。

### 要点
- **入力:** 1 つ以上の複合体 PDB（アンサンブル対応。原子順序の整合性が必要）。
- **基質指定（`-c`）:** 残基 ID（例: `A:123A`）、残基名（例: `GPP,SAM`）、または基質 PDB（先頭入力と座標一致が必要）。
- **選択ロジック:** `--radius` を基本に、ヘテロ–ヘテロ近接（`--radius-het2het`）やジスルフィド/PRO 隣接などのセーフガードを併用。
- **切断/キャッピング:** 主鎖/側鎖のトリミングと、リンク水素（`--add-linkh` がデフォルト）。
- **電荷:** 未知残基は 0 がデフォルト。必要に応じて `--ligand-charge`（総電荷または残基名マップ）で指定。

`pdb2reaction extract` は基質近傍の残基を選択して活性部位モデル（クラスターモデル）を生成し、規則に従ってトリミングしたうえで、必要に応じて切断結合をリンク水素でキャップします。単一構造だけでなく、複数 PDB を入力するアンサンブル処理にも対応しています。

残基/原子の命名が非標準で残基分類や電荷サマリーに影響がある場合は、下部の付録（PDB 命名規則と内部参照リスト）を参照してください。

## 使用法
```bash
pdb2reaction extract -i COMPLEX.pdb [COMPLEX2.pdb...]
 -c SUBSTRATE_SPEC
 [-o POCKET.pdb [POCKET2.pdb...]]
 [--radius Å] [--radius-het2het Å]
 [--include-h2o/--no-include-h2o]
 [--exclude-backbone/--no-exclude-backbone]
 [--add-linkh/--no-add-linkh]
 [--selected-resn LIST]
 [--modified-residue LIST]
 [-l, --ligand-charge MAP_OR_NUMBER]
 [--verbose/--no-verbose]
```

### 例
```bash
# 最小（ID基準の基質）+ 明示的な総リガンド電荷
pdb2reaction extract -i complex.pdb -c '123' -o model.pdb -l -3

# PDBとして提供される基質; 残基名ごとの電荷マッピング
pdb2reaction extract -i complex.pdb -c substrate.pdb -o model.pdb -l 'GPP:-3,SAM:1'

# 名前基準の基質選択（すべてのマッチを含む）
pdb2reaction extract -i complex.pdb -c 'GPP,SAM' -o model.pdb -l 'GPP:-3,SAM:1'

# ヘテロ-ヘテロ近接を有効にした複数構造から単一のマルチMODEL出力
pdb2reaction extract -i complex1.pdb complex2.pdb -c 'GPP,SAM' \
 -o model_multi.pdb --radius-het2het 2.6 -l 'GPP:-3,SAM:1'

# ヘテロ-ヘテロ近接を有効にした複数構造から複数出力
pdb2reaction extract -i complex1.pdb complex2.pdb -c 'GPP,SAM' \
 -o model1.pdb model2.pdb --radius-het2het 2.6 -l 'GPP:-3,SAM:1'
```

## ワークフロー
### 残基包含
- `-c/--center` からの基質残基を常に含める
- **標準カットオフ（`--radius`、デフォルト2.6 Å）:**
 - `--no-exclude-backbone` の場合、カットオフ内の任意の原子が残基を対象にする
 - `--exclude-backbone` の場合、アミノ酸残基は**非主鎖**原子（N/H*/CA/HA*/C/O以外）で基質に接触する必要がある
- **独立したヘテロ-ヘテロカットオフ（`--radius-het2het`）:** 基質ヘテロ原子（非C/H）がタンパク質ヘテロ原子の指定Å以内にある場合に残基を追加
- **水処理:** HOH/WAT/H2O/DOD/TIP/TIP3/SOLはデフォルトで含まれる（`--include-h2o`）
- **強制包含:** `--selected-resn` はチェーン/挿入コード付きIDを受け入れる（例: `A:123A`）
- **近傍セーフガード:**
 - `--no-exclude-backbone` で主鎖原子が基質に接触した場合、ペプチド隣接のN/C側残基（C–N ≤ 1.9 Å）を自動的に含める。末端はN/H*またはC/O/OXTのキャップを保持。
 - ジスルフィド結合（SG–SG ≤ 2.5 Å）は両方のCysを包含。
 - 非末端PRO残基は常にN側隣接残基を含め、主鎖除去後もCAを保持します。`--exclude-backbone` の場合は隣接残基のC/O/OXTを残し、ペプチド結合を維持。

### 切断/キャッピング
- 孤立残基は側鎖原子のみを保持; アミノ酸主鎖原子（N, CA, C, O, OXT + N/CA水素）はPRO/HYP保護を除いて除去
- 連続ペプチドストレッチは内部主鎖原子を保持; 末端キャップ（N/H*またはC/O/OXT）のみ除去
- TERを認識し、チェーン切断を跨ぐキャッピングは行わない
- `--exclude-backbone` の場合、**非基質**アミノ酸の主鎖原子を除去（PRO/HYP保護とPRO近傍保持は適用）
- 非アミノ酸残基は主鎖様原子名を持つ原子を失わない

### リンク水素（`--add-linkh`）
- 切断された結合ベクトルに沿って1.09 Åで炭素のみのリンク水素を追加（CB–CA、CA–N、CA–C; PRO/HYPはCA–Cのみ）
- `TER` の後に残基 `LKH`（チェーン `L`）の連続した `HETATM` レコードとして `HL` という名前で挿入されます。シリアル番号は本体ブロックからの連番です
- マルチ構造モードでは全モデルで同じ結合にキャップを付け、座標はモデルごとに保持されます

### 電荷サマリー（`--ligand-charge`）
- アミノ酸と一般的なイオンは内部辞書から電荷を取得; 水はゼロ
- 未知残基は `--ligand-charge` が総電荷または残基名ごとのマッピング（例: `GPP:-3,SAM:1`）を提供しない限りデフォルトで0。総電荷が与えられた場合は未知基質残基に配分され、未知基質が無い場合は未知残基全体に配分されます。
- verboseモードが有効な場合、モデル#1の電荷サマリー（タンパク質/リガンド/イオン/総計）がログに記録されます。


### マルチ構造アンサンブル
- 複数の入力 PDB を受け付けます（先頭/末尾で原子順序の一致を検証）。各構造は独立に処理され、選択残基の**和集合**を全モデルに適用することで出力の一貫性を保ちます。
- 出力ポリシー:
 - `-o` なし & 複数入力 → 構造ごとに `model_<original_basename>.pdb`。
 - `-o` を1つだけ指定 → 単一のマルチMODEL PDB。
 - 入力数と同数の `-o` を指定 → 入力ごとに個別PDB。
- 診断ログにモデルごとの全原子数/保持原子数と残基 ID を出力します。

### 基質指定（`-c/--center`）
- **PDB パス**: 座標が先頭入力と完全一致（許容誤差 1e-3 Å）。残基IDは他構造へ伝播。
- **残基ID**: `'123,124'`, `'A:123,B:456'`, `'123A'`, `'A:123A'`（挿入コード対応）。
- **残基名**: カンマ区切り（大文字小文字は無視）。同名残基が複数ある場合は**すべて**含め、警告を出力。

## CLI オプション

> **注記:** 表示されているデフォルト値は、オプション未指定時に使用される内部デフォルトです。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 1つ以上のタンパク質-リガンドPDB ファイル（同一の原子順序が必要） | 必須 |
| `-c, --center SPEC` | 基質指定（PDB パス、残基ID、または残基名） | 必須 |
| `-o, --output PATH...` | 活性部位モデル PDB 出力。1パス ⇒ マルチMODEL、Nパス ⇒ 入力ごと。複数入力で `-o` 1つの場合は単一のマルチMODEL PDB を生成。N個の `-o` が N個の入力と一致する場合は N個の個別 PDB を生成。 | 自動（`model.pdb` または `model_<input>.pdb`） |
| `-r, --radius FLOAT` | 包含のための原子-原子距離カットオフ（Å） | `2.6` |
| `--radius-het2het FLOAT` | 独立したヘテロ-ヘテロカットオフ（Å、非C/H） | `0.0`（0 の場合は内部で 0.001 Å） |
| `--include-h2o/--no-include-h2o` | HOH/WAT/H2O/DOD/TIP/TIP3/SOL水を含める | `True` |
| `--exclude-backbone/--no-exclude-backbone` | 非基質アミノ酸の主鎖原子を除去 | `False` |
| `--add-linkh/--no-add-linkh` | 切断された結合に1.09 Åで炭素のみのリンク水素を追加 | `True` |
| `--selected-resn TEXT` | 強制包含残基（オプションのチェーン/挿入コード付きID） | `""` |
| `--modified-residue TEXT` | 修飾アミノ酸残基名をカンマ区切りで指定（任意で電荷付き）。主鎖切断と電荷計算にアミノ酸として扱う。例: `HD1,HD2,HD3` または `HD1:0,SEP:-2`。電荷省略時は 0。 | `""` |
| `-l, --ligand-charge TEXT` | 総電荷または残基名ごとのマッピング（例: `GPP:-3,SAM:1`） | _None_ |
| `-v, --verbose/--no-verbose` | INFOレベルログを出力（`True`）または警告のみ（`False`） | `True` |

## 出力
```text
<output>.pdb # TERレコード後にオプションのリンク水素を含む活性部位モデル PDB
 # 単一入力 → デフォルトでmodel.pdb
 # -oなしの複数入力 → 構造ごとにmodel_<original_basename>.pdb
 # 複数入力で1つの-oパス → 単一のマルチMODEL PDB
 # 出力ディレクトリは自動作成されません。事前に存在を確認してください
```
- verboseモードが有効な場合、モデル#1の電荷サマリー（タンパク質/リガンド/イオン/総計）がログに記録
- 出力ディレクトリは自動作成されません。必要に応じて事前に作成してください。
- API利用（`extract_api`）では `{"outputs": [...], "counts": [...], "charge_summary": {...}}` を返します。


## MCPB 等で生成された非標準残基を含む系

Amber の `MCPB.py`（Metal Center Parameter Builder）等で金属配位残基のパラメータを生成した場合、金属配位アミノ酸に非標準の残基名（`HD1`, `HE1`, `CM1`, `AP1` 等）が割り当てられます。これらは `extract` の内部辞書 `AMINO_ACIDS` に含まれないため、**主鎖原子の切断・リンク水素の付加が正しく行われません**。

このような系では、`extract` の実行時に以下のような警告が表示されます:

```
[extract] WARNING: Residue HD1 83 may be an amino acid (has N, CA, C, O)
but is not recognized as a standard residue name.
Backbone truncation was not applied.
Consider preparing the active site model manually.
```

### `--modified-residue` オプション

`--modified-residue` を使用すると、非標準の残基名をアミノ酸として登録でき、主鎖切断と電荷割り当てが自動的に適用されます。修飾アミノ酸残基で非標準の3文字コードを持つもの（リン酸化セリン、メチル化残基、特殊な名前の D-アミノ酸、MCPB でリネームされた金属配位残基など）に有用です。

```bash
# HD1, HD2, HD3 をアミノ酸として扱う（電荷はデフォルトで 0）
pdb2reaction extract -i complex.pdb -c 'SUB' -o model.pdb \
  --modified-residue 'HD1,HD2,HD3'

# 各修飾残基に明示的な電荷を指定
pdb2reaction extract -i complex.pdb -c 'SUB' -o model.pdb \
  --modified-residue 'HD1:0,SEP:-2'
```

```{important}
`--modified-residue` で対応できない場合は、**活性部位モデルを手動で構築する**ことを推奨します。
手動構築の手順:

1. 活性部位周辺の残基を選定し、切断箇所を決定する
2. 切断された共有結合の親原子（残る側の原子）に、リンク水素を付加する
3. リンク水素は残基名 `LKH`（チェーン `L`）、原子名 `HL` で記述する
4. 結合方向に沿って **1.09 Å** の位置に配置する
```

## 付録: PDB 命名規則と内部参照リスト

この付録は、`extract` が **非標準の残基名/原子名** により残基分類や電荷割り当てを誤る場合の調査用です。入力が標準的な PDB 命名に従っている場合、多くのユーザーは読み飛ばして問題ありません。

```{important}
`extract` が正しく動作するためには、**入力 PDBの残基名と原子名が標準的なPDB命名規則に準拠している必要があります**。このツールはアミノ酸、イオン、水分子、主鎖原子を認識するために内部辞書を使用しています。非標準の命名を使用すると、残基が誤って分類されたり、電荷が正しく割り当てられない可能性があります。
```

以下の内部定数が認識される名前を定義しています：

### `AMINO_ACIDS`

残基名を公称整数電荷にマッピングする辞書です。この辞書に含まれるかどうかで、残基が主鎖処理、切断、電荷計算においてアミノ酸として扱われるかが決まります。

**標準20アミノ酸**（生理的pHでの電荷）：
- 中性: `ALA`, `ASN`, `CYS`, `GLN`, `GLY`, `HIS`, `ILE`, `LEU`, `MET`, `PHE`, `PRO`, `SER`, `THR`, `TRP`, `TYR`, `VAL`
- 正電荷 (+1): `ARG`, `LYS`
- 負電荷 (−1): `ASP`, `GLU`

**追加の標準残基：**
- `SEC`（セレノシステイン、0）、`PYL`（ピロリシン、+1）

**プロトン化/互変異性体**（Amber/CHARMM形式）：
- `HIP`（+1、完全プロトン化His）、`HID`（0、Nδプロトン化His）、`HIE`（0、Nεプロトン化His）
- `ASH`（0、中性Asp）、`GLH`（0、中性Glu）、`LYN`（0、中性Lys）、`ARN`（0、中性Arg）
- `TYM`（−1、脱プロトン化Tyrフェノラート）

**リン酸化残基：**
- 二価陰イオン（−2）: `SEP`, `TPO`, `PTR`
- 一価陰イオン（−1）: `S1P`, `T1P`, `Y1P`
- リン酸化His（phosaa19SB）: `H1D`（0）、`H2D`（−1）、`H1E`（0）、`H2E`（−1）

**システイン変異体：**
- `CYX`（0、ジスルフィド）、`CSO`（0、スルフェン酸）、`CSD`（−1、スルフィン酸）、`CSX`（0、汎用誘導体）
- `OCS`（−1、システイン酸）、`CYM`（−1、脱プロトン化Cys）

**リシン変異体/カルボキシル化：**
- `MLY`（+1）、`LLP`（+1）、`KCX`（−1、Nz-カルボン酸）

**D-アミノ酸**（19残基）：
- `DAL`, `DAR`, `DSG`, `DAS`, `DCY`, `DGN`, `DGL`, `DHI`, `DIL`, `DLE`, `DLY`, `MED`, `DPN`, `DPR`, `DSN`, `DTH`, `DTR`, `DTY`, `DVA`

**その他の修飾残基：**
- `CGU`（−2、γ-カルボキシグルタミン酸）、`CGA`（−1）、`PCA`（0、ピログルタミン酸）、`MSE`（0、セレノメチオニン）、`OMT`（0、メチオニンスルホン）、`HYP`（0、ヒドロキシプロリン）
- その他: `ASA`, `CIR`, `FOR`, `MVA`, `IIL`, `AIB`, `HTN`, `SAR`, `NMC`, `PFF`, `NFA`, `ALY`, `AZF`, `CNX`, `CYF`

**N末端変異体**（接頭辞 `N`）: `NALA`（+1）、`NARG`（+2）、`NASP`（0）、`NGLU`（0）、`NLYS`（+2）など、および `ACE`（0）、`NTER`（+1、汎用）

**C末端変異体**（接頭辞 `C`）: `CALA`（−1）、`CARG`（0）、`CASP`（−2）、`CGLU`（−2）、`CLYS`（0）など、および `NHE`（0）、`NME`（0）、`CTER`（−1、汎用）

### `BACKBONE_ATOMS`

アミノ酸の主鎖原子と見なされる原子名のセットです。`--exclude-backbone` の場合、非基質残基からどの原子を除去するかを決定するために使用されます：

```
N, C, O, CA, OXT, H, H1, H2, H3, HN, HA, HA2, HA3
```

### `ION`

イオン残基名をその形式電荷にマッピングする辞書です。認識されたイオンは電荷サマリーで自動的に正しい電荷が割り当てられます。

| 電荷 | 残基名 |
|------|--------|
| +1 | `LI`, `NA`, `K`, `RB`, `CS`, `TL`, `AG`, `CU1`, `Ag`, `K+`, `NA+`, `NH4`, `H3O+` |
| +2 | `MG`, `CA`, `SR`, `BA`, `MN`, `FE2`, `CO`, `NI`, `CU`, `ZN`, `CD`, `HG`, `PB`, `BE`, `PD`, `PT`, `SN`, `RA`, `YB2`, `V2+` |
| +3 | `FE`, `AU3`, `AL`, `GA`, `IN`, `CE`, `CR`, `DY`, `EU`, `EU3`, `ER`, `GD3`, `LA`, `LU`, `ND`, `PR`, `SM`, `TB`, `TM`, `Y`, `PU` |
| +4 | `U4+`, `TH`, `HF`, `ZR` |
| −1 | `F`, `CL`, `BR`, `I`, `CL-`, `IOD` |

### `WATER_RES`

水分子として認識される残基名のセットです。水はデフォルトで含まれ（`--include-h2o`）、電荷はゼロが割り当てられます：

```
HOH, WAT, H2O, DOD, TIP, TIP3, SOL
```

---

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- `--radius` のデフォルトは 2.6 Å。`0` を指定すると空選択を避けるため内部で 0.001 Å にクランプされます。`--radius-het2het` もデフォルトでは無効（`0.0`）で、`0` 指定時は同様に 0.001 Å にクランプされます。
- `--no-include-h2o` で水を除外できます。
- 主鎖トリミングとキャッピングはチェーン切断（TER）や PRO/HYP 保護を尊重します。非アミノ酸残基は主鎖様の原子名を保持します。
- リンク水素は炭素切断のみで挿入され、アンサンブルモードでは同一結合パターンを全モデルで再利用します。
- INFO ログに残基選択、切断数、電荷内訳の要約が出力されます。

## リンク水素と凍結原子

pdb2reaction が活性部位モデルを抽出する際、切断された結合は**リンク水素**でキャップされます。デフォルト（`--freeze-links`）では、リンク水素の親原子が最適化や経路探索中に凍結され、境界での非物理的な再配置を防ぎます。

- **力**: 凍結原子の力はゼロ化されます。
- **ヘシアン**: 凍結自由度は除去（`return_partial_hessian: true`）またはフル行列でゼロ化されます。
- **振動解析**: 凍結原子がある場合、`freq` は自動的に部分ヘシアン振動解析（PHVA: Partial Hessian Vibrational Analysis）を行い、活性ブロックのみを対角化します。

凍結原子は `geom.freeze_atoms` YAML キー（1 始まりインデックス）で手動設定も可能です。CLI で検出されたリンク原子は YAML 指定の原子とマージされます。

---

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け

- [all](all.md) — `-c/--center` で内部的にextractを呼び出す一気通貫ワークフロー
- [path-search](path-search.md) — 抽出された活性部位モデルでのMEP 探索
- [scan](scan.md) — 抽出された活性部位モデルでの段階的スキャン
- [add-elem-info](add-elem-info.md) — 抽出前に欠落したPDB元素カラムを修正
- [トラブルシューティング](troubleshooting.md) — よくある抽出エラー
- [用語集](glossary.md) — 活性部位モデル、クラスターモデル、リンク水素の定義
