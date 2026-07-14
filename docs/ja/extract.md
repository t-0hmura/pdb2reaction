# `extract`

タンパク質–リガンド PDB/mmCIF から活性部位クラスターモデルを切り出します。`-c/--center` は残基名、残基ID、chain付き残基名（`A:SAM`）、chain+残基名+番号（`A:SAM:123`）、または PDB/mmCIF パスを受け付けます。mmCIF および PDB 固定幅の上限に達する構造は内部で安全なPDB IDへ再割当てし、元のIDを復元したCIFも出力します。

## 実行例

コマンド形式:

```bash
pdb2reaction extract -i COMPLEX.pdb [COMPLEX2.pdb ...] \
 -c SUBSTRATE_SPEC \
 [-o MODEL.pdb [MODEL2.pdb ...]] \
 [--radius Å] [--radius-het2het Å] \
 [--include-h2o/--no-include-h2o] \
 [--exclude-backbone/--no-exclude-backbone] \
 [--add-linkh/--no-add-linkh] \
 [--selected-resn LIST] \
 [--modified-residue LIST] \
 [-l, --ligand-charge MAP_OR_NUMBER] \
 [--out-json/--no-out-json] \
 [-v LEVEL]
```

最小（ID基準の基質）+ 明示的な総リガンド電荷:

```bash
# 最小（ID基準の基質）+ 明示的な総リガンド電荷
pdb2reaction extract -i complex.pdb -c '123' -o model.pdb -l -3
```

PDB として提供される基質。残基名ごとの電荷マッピング（その他は 0）:

```bash
# PDB として提供される基質。残基名ごとの電荷マッピング（その他は 0）
pdb2reaction extract -i complex.pdb -c substrate.pdb -o model.pdb -l 'GPP:-3,SAM:1'
```

名前基準の基質選択（すべてのマッチを含む。WARNING ログ出力）:

```bash
# 名前基準の基質選択（すべてのマッチを含む。WARNING ログ出力）
pdb2reaction extract -i complex.pdb -c 'GPP,SAM' -o model.pdb -l 'GPP:-3,SAM:1'
```

ヘテロ-ヘテロ近接を有効にした複数構造から単一のマルチ MODEL 出力:

```bash
# ヘテロ-ヘテロ近接を有効にした複数構造から単一のマルチMODEL出力
pdb2reaction extract -i complex1.pdb -i complex2.pdb -c 'GPP,SAM' \
 -o model_multi.pdb --radius-het2het 2.6 -l 'GPP:-3,SAM:1'
# 複数出力にする場合は -o model1.pdb -o model2.pdb を指定
```

## 処理の流れ

### 残基包含

- `-c/--center` からの基質残基を常に含める
- **標準カットオフ（`--radius`、デフォルト 2.6 Å）:**
 - `--no-exclude-backbone` の場合、カットオフ内の任意の原子が残基を対象にする
 - `--exclude-backbone` の場合、アミノ酸残基は**非主鎖**原子（N/H*/CA/HA*/C/O/OXT 以外）で基質に接触する必要がある。非アミノ酸残基は任意の原子で接触判定される。
- **独立したヘテロ-ヘテロカットオフ（`--radius-het2het`）:** 基質ヘテロ原子（非 C/H）がタンパク質ヘテロ原子の指定距離（Å）以内にある場合に残基を追加。`--exclude-backbone` 有効時はタンパク質側原子も非主鎖でなければならない。
- **水処理:** HOH/WAT/H2O/DOD/TIP/TIP3/SOL はデフォルトで含まれる（`--include-h2o`）
- **強制包含:** `--selected-resn` は `--center` と同じ残基ID・残基名・chain付きselectorを受け入れます。詳細は {ref}`ja-selected-resn-takes-ids` を参照。
- **近傍セーフガード:**
 - `--no-exclude-backbone` で主鎖原子が基質に接触した場合、ペプチド隣接の N/C 側残基（C–N ≤ 1.9 Å）を自動的に含める。末端は N/H*または C/O/OXT のキャップを保持。
 - ジスルフィド結合（SG–SG ≤ 2.5 Å）は両方の Cys を包含。
 - 非末端 PRO 残基は常に N 側隣接残基を含め、主鎖除去後も CA を保持します。`--exclude-backbone` の場合は隣接残基の C/O/OXT を残し、ペプチド結合を維持。

### 切断/キャッピング

- 孤立残基は側鎖原子のみを保持; アミノ酸主鎖原子（N, CA, C, O, OXT + N/CA 水素）は PRO/HYP 保護を除いて除去
- 連続ペプチドストレッチは内部主鎖原子を保持; 末端キャップ（N/H*または C/O/OXT）のみ除去
- TER を認識し、チェーン切断を跨ぐキャッピングは行わない
- `--exclude-backbone` の場合、**非基質**アミノ酸の主鎖原子を除去（PRO/HYP 保護と PRO 近傍保持は適用）
- 非アミノ酸残基は主鎖様の原子名（N/CA/HA/H/H1/H2/H3）を持つ原子を失わない

### キャップ水素（`--add-linkh`）

- 切断された結合ベクトル（CB–CA、CA–N、CA–C; PRO/HYP は CA–C のみ）に沿って 1.09 Å のキャップ水素を炭素境界にのみ付加（非炭素境界はキャップしない）
- `TER` の後に残基 `LKH`（チェーン `L`）の連続した `HETATM` レコードとして `HL` という名前で挿入されます。シリアル番号は本体ブロックからの連番です
- マルチ構造モードでは全モデルで同じ結合にキャップを付け、座標はモデルごとに保持されます

### クラスターモデルを手作業で構築／監査する場合

切断結合は単なる距離cutoffではなく、化学的なmodeling判断として扱います。

- タンパク質主鎖断片を残す場合は、両端の主鎖末端が一貫してCα（PDB原子名
  `CA`）になるよう残基範囲を選び、末端原子価をcapで満たします。
- 側鎖・ligand・cofactor境界は、可能な限り非極性の**C–C単結合**
  （典型的には`CA–CB`またはさらに外側の脂肪族C–C）で切断します。距離cutoffを
  跨いだだけでpeptide C–N、極性C–O/C–N、芳香族／共役結合、S–S結合、
  金属配位結合を切らず、結合相手を含めるか境界を移してください。
- cluster側の各境界原子が意図した原子価とcap 1個を持つことを確認し、production
  最適化ではcap親を`--freeze-links`（default）で凍結します。
- R/IM/Pは同一原子・同一順序にします。各状態を独立に再抽出せず、1構造で決めた
  選択／cap patternを全状態へ適用してください。
- truncation後の総電荷とmultiplicityを再計算し、MEP/Hessian前に境界を目視確認します。

自動extractorは一般的なproteinを扱いますが、共有結合cofactor、修飾残基、metal site、
または上記規則に反する境界は手作業で修正してください。

### 電荷サマリー（`--ligand-charge/-l`）

- アミノ酸と一般的なイオンは内部辞書から電荷を取得; 水はゼロ
- 未知残基は `--ligand-charge` が総電荷または残基名ごとのマッピング（例: `GPP:-3,SAM:1`）を提供しない限りデフォルトで 0。総電荷が与えられた場合は未知基質残基に配分され、未知基質が無い場合は未知残基全体に配分されます。

### マルチ構造アンサンブル

- 複数の入力 PDB/mmCIF を受け付けます（全原子のidentity/orderを検証）。各構造は独立に処理され、選択残基の**和集合**を全モデルに適用します。
- 出力ポリシー:
 - `-o` なし & 複数入力 → 構造ごとに `model_<original_basename>.pdb`。
 - `-o` を 1 つだけ指定 → 単一のマルチ MODEL PDB。
 - 入力数と同数の `-o` を指定 → 入力ごとに個別 PDB。
- 診断ログにモデルごとの全原子数/保持原子数と残基 ID を出力します。

## 出力

```text
<output>.pdb # TERレコード後にオプションのキャップ水素を含む活性部位モデル PDB
<output>.cif # mmCIF/oversized-PDB入力時。元のchain/residue IDを復元
 # 単一入力 → デフォルトでmodel.pdb
 # -oなしの複数入力 → 構造ごとにmodel_<original_basename>.pdb
 # 複数入力で1つの-oパス → 単一のマルチMODEL PDB
 # 親出力directoryは自動作成されます
```

- verbose モードが有効な場合、モデル#1 の電荷サマリー（タンパク質/リガンド/イオン/総計）がログに記録されます。
- API 利用（`extract_api`）では `{"outputs": [...], "counts": [...], "charge_summary": {...}, "n_link_hydrogens": N}` を返します。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 1つ以上のタンパク質–リガンド PDB/mmCIF（同一の原子identity/orderが必要） | 必須 |
| `-c, --center SPEC` | PDB/mmCIFパス、残基ID/名、`CHAIN:RESNAME`、`CHAIN:RESNAME:RESSEQ` | 必須 |
| `-o, --output PATH...` | 活性部位モデル PDB 出力。1 パス ⇒ マルチ MODEL、N パス ⇒ 入力ごと。複数入力で `-o` 1 つの場合は単一のマルチ MODEL PDB を生成。N 個の `-o` が N 個の入力と一致する場合は N 個の個別 PDB を生成 | 自動（`model.pdb` または `model_<input>.pdb`） |
| `-r, --radius FLOAT` | 包含のための原子-原子距離カットオフ（Å、0 の場合は内部で 0.001 Å） | `2.6` |
| `--radius-het2het FLOAT` | 独立したヘテロ-ヘテロカットオフ（Å、非 C/H） | `0.0`（0 の場合は内部で 0.001 Å） |
| `--include-h2o/--no-include-h2o` | HOH/WAT/H2O/DOD/TIP/TIP3/SOL 水を含める | `True` |
| `--exclude-backbone/--no-exclude-backbone` | 非基質アミノ酸の主鎖原子を除去 | `False` |
| `--add-linkh/--no-add-linkh` | 切断された結合に 1.09 Å のキャップ水素を炭素境界にのみ付加（非炭素境界はキャップしない） | `True` |
| `--selected-resn TEXT` | `--center` と同じselectorで残基を強制包含 | `""` |
| `--modified-residue TEXT` | 修飾アミノ酸残基名をカンマ区切りで指定（任意で各残基に電荷付き）。主鎖切断と電荷計算でアミノ酸として扱います。例: `HD1,HD2,HD3` または `HD1:0,SEP:-2`。残基ごとに `:charge` 接尾辞を省略した場合、その残基の電荷は `0` になります（例: `HD1,HD2:-1` では `HD1` が電荷 0、`HD2` が電荷 −1）。フラグ全体のデフォルトは空文字列（無効） | `""` |
| `-l, --ligand-charge TEXT` | 総電荷または残基名ごとのマッピング（例: `GPP:-3,SAM:1`） | _None_ |
| `--out-json/--no-out-json` | 抽出された PDB(s) の隣に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照 | `False` |

### 基質指定（`-c/--center`）

- **PDB/mmCIF パス**: 座標が先頭入力と完全一致（許容誤差 1e-3 Å）。残基 ID は他構造へ伝播。
- **残基 ID**: `'123,124'`, `'A:123,B:456'`, `'123A'`, `'A:123A'`（挿入コード対応）。
- **残基名**: カンマ区切り（大文字小文字は無視）。同名残基が複数ある場合は**すべて**含め、警告を出力。
- **chain + 残基名**: `A:SAM` はchain A内のSAMをすべて選択し、`A:SAM:123` は1残基に限定。

## 注記

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。
- 対象系ごとに抽出半径の収束性と化学的完全性を検証してください。`-r` を大きくすると環境を多く含められますが、精度が単調に改善する保証はなく計算costも増えます。
- INFO ログに残基選択、切断数、電荷内訳の要約が出力されます。

## MCPB 等で生成された非標準残基を含む系

Amber の `MCPB.py`（Metal Center Parameter Builder）等で金属配位残基のパラメータを生成した場合、金属配位アミノ酸に非標準の残基名（`HD1`, `HE1`, `CM1`, `AP1` 等）が割り当てられます。これらは `extract` の内部辞書 `AMINO_ACIDS` に含まれないため、**主鎖原子の切断・キャップ水素の付加が正しく行われません**。

このような系では、`extract` の実行時に以下のような警告が表示されます:

```
[extract] WARNING: Residue HD1 83 may be an amino acid (has N, CA, C, O)
but is not recognized as a standard residue name.
Backbone truncation was not applied.
Consider preparing the active site model manually.
```

### `--modified-residue` オプション

`--modified-residue` を使用すると、非標準の残基名をアミノ酸として登録でき、主鎖切断と電荷割り当てが自動的に適用されます。修飾アミノ酸残基で非標準の 3 文字コードを持つもの（リン酸化セリン、メチル化残基、特殊な名前の D-アミノ酸、MCPB でリネームされた金属配位残基など）に有用です。

```bash
# HD1, HD2, HD3 をアミノ酸として扱う（電荷はデフォルトで 0）
pdb2reaction extract -i complex.pdb -c 'SUB' -o model.pdb \
  --modified-residue 'HD1,HD2,HD3'

# 各修飾残基に明示的な電荷を指定
pdb2reaction extract -i complex.pdb -c 'SUB' -o model.pdb \
  --modified-residue 'HD1:0,SEP:-2'
```

```{important}
`--modified-residue` で対応できない場合は、**活性部位モデルを手動で構築**してください。
手動構築の手順:

1. 活性部位周辺の残基を選定し、切断箇所を決定する
2. 切断された共有結合の親原子（残る側の原子）に、キャップ水素を付加する
3. キャップ水素は残基名 `LKH`（チェーン `L`）、原子名 `HL` で記述する
4. 結合方向に沿って **1.09 Å** の位置に配置する
```

## 付録: PDB 命名規則と内部参照リスト

この付録は、`extract` が **非標準の残基名/原子名** により残基分類や電荷割り当てを誤る場合の調査用です。

```{important}
`extract` が正しく動作するためには、**入力 PDBの残基名と原子名が標準的なPDB命名規則に準拠している必要があります**。このツールはアミノ酸、イオン、水分子、主鎖原子を認識するために内部辞書を使用しています。非標準の命名を使用すると、残基の誤分類や電荷の誤割り当てが発生します。
```

以下の内部定数が認識される名前を定義しています：

### `AMINO_ACIDS`

残基名を公称整数電荷にマッピングする辞書です。この辞書に含まれるかどうかで、残基が主鎖処理、切断、電荷計算においてアミノ酸として扱われるかが決まります。

**標準 20 アミノ酸**（生理的 pH での電荷）：
- 中性: `ALA`, `ASN`, `CYS`, `GLN`, `GLY`, `HIS`, `ILE`, `LEU`, `MET`, `PHE`, `PRO`, `SER`, `THR`, `TRP`, `TYR`, `VAL`
- 正電荷 (+1): `ARG`, `LYS`
- 負電荷 (−1): `ASP`, `GLU`

**追加の標準残基：**
- `SEC`（セレノシステイン、0）、`PYL`（ピロリシン、+1）

**プロトン化/互変異性体**（Amber/CHARMM 形式）：
- `HIP`（+1、完全プロトン化 His）、`HID`（0、Nδプロトン化 His）、`HIE`（0、Nεプロトン化 His）
- `ASH`（0、中性 Asp）、`GLH`（0、中性 Glu）、`LYN`（0、中性 Lys）、`ARN`（0、中性 Arg）
- `TYM`（−1、脱プロトン化 Tyr フェノラート）

**リン酸化残基：**
- 二価陰イオン（−2）: `SEP`, `TPO`, `PTR`
- 一価陰イオン（−1）: `S1P`, `T1P`, `Y1P`
- リン酸化 His（phosaa19SB）: `H1D`（0）、`H2D`（−1）、`H1E`（0）、`H2E`（−1）

**システイン変異体：**
- `CYX`（0、ジスルフィド）、`CSO`（0、スルフェン酸）、`CSD`（−1、スルフィン酸）、`CSX`（0、汎用誘導体）
- `OCS`（−1、システイン酸）、`CYM`（−1、脱プロトン化 Cys）

**リシン変異体/カルボキシル化：**
- `MLY`（+1）、`LLP`（+1）、`KCX`（−1、Nz-カルボン酸）

**D-アミノ酸**（19 残基）：
- `DAL`, `DAR`, `DSG`, `DAS`, `DCY`, `DGN`, `DGL`, `DHI`, `DIL`, `DLE`, `DLY`, `MED`, `DPN`, `DPR`, `DSN`, `DTH`, `DTR`, `DTY`, `DVA`

**その他の修飾残基：**
- `CGU`（−2、γ-カルボキシグルタミン酸）、`CGA`（−1）、`PCA`（0、ピログルタミン酸）、`MSE`（0、セレノメチオニン）、`OMT`（0、メチオニンスルホン）、`HYP`（0、ヒドロキシプロリン）
- その他: `ASA`, `CIR`, `FOR`, `MVA`, `IIL`, `AIB`, `HTN`, `SAR`, `NMC`, `PFF`, `NFA`, `ALY`, `AZF`, `CNX`, `CYF`

**N 末端変異体**（接頭辞 `N`）: `NALA`（+1）、`NARG`（+2）、`NASP`（0）、`NGLU`（0）、`NLYS`（+2）など、および `ACE`（0）、`NTER`（+1、汎用）

**C 末端変異体**（接頭辞 `C`）: `CALA`（−1）、`CARG`（0）、`CASP`（−2）、`CGLU`（−2）、`CLYS`（0）など、および `NHE`（0）、`NME`（0）、`CTER`（−1、汎用）

### `BACKBONE_ATOMS`

アミノ酸の主鎖原子と見なされる原子名のセットです。`--exclude-backbone` の場合、非基質残基からどの原子を除去するかを決定するために使用されます：

```
N, C, O, CA, OXT, H, H1, H2, H3, HN, HA, HA2, HA3
```

### `ION`

イオン残基名をその形式電荷にマッピングする辞書です。認識されたイオンは電荷サマリーで自動的に正しい電荷が割り当てられます。

| 電荷 | 残基名 |
|------|--------|
| +1 | `LI`, `NA`, `K`, `RB`, `CS`, `TL`, `AG`, `CU1`, `K+`, `NA+`, `NH4`, `H3O+`, `HE+`, `HZ+` |
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

(ja-link-hydrogen-and-frozen-atoms)=
## キャップ水素と凍結原子

pdb2reaction が活性部位モデルを抽出する際、切断された結合は**キャップ水素**でキャップされます。デフォルト（`--freeze-links`）では、キャップ水素の親原子が最適化や経路探索中に凍結され、境界での非物理的な再配置を防ぎます。

- **力**: 凍結原子の力はゼロ化されます。
- **Hessian**: 凍結自由度は除去（`return_partial_hessian: true`）またはフル行列でゼロ化されます。
- **振動解析**: 凍結原子がある場合、`freq` は自動的に部分Hessian振動解析（PHVA: Partial Hessian Vibrational Analysis）を行い、活性ブロックのみを対角化します。

凍結原子は `geom.freeze_atoms` YAML キー（1 始まりインデックス）で手動設定も可能です。CLI で検出されたキャップ原子は YAML 指定の原子とマージされます。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け

- [all](all.md) — `-c/--center` で内部的に extract を呼び出す一気通貫ワークフロー
- [path-search](path-search.md) — 抽出された活性部位モデルでの MEP 探索
- [scan](scan.md) — 抽出された活性部位モデルでの段階的スキャン
- [add-elem-info](add-elem-info.md) — 抽出前に欠落した PDB 元素カラムを修正
- [トラブルシューティング](troubleshooting.md) — よくある抽出エラー
- [用語集](glossary.md) — 活性部位モデル、クラスターモデル、キャップ水素の定義
