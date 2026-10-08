# `extract`（活性部位モデルの切り出し）

`extract` は、タンパク質–リガンドの PDB/mmCIF から基質の周りの残基を切り出し、切った結合をキャップ水素で埋めて、できたクラスターモデルの電荷を数えます。

---

## 主な用途

* **クラスターモデルを作る**: `all`・`opt`・`tsopt` などで計算する活性部位モデルを作ります。
* **複数の状態を同じ境界で切る**: 原子の並びが同じ反応物と生成物を 1 回で渡すと、どのモデルも同じ残基・同じ境界になります。
* **非標準残基を扱う**: MCPB.py などが付けた残基名を、`--modified-residue` でアミノ酸として登録します。

モデルの大きさの決め方は [クラスターモデルの組み方](model-setup.md) を見てください。

---

## 基本的な実行例

入力には、[すべての水素原子](getting-started.md#入力構造に関する重要事項)と 77–78 列の元素記号が要ります。

### 1. 残基 ID と総電荷で選ぶ

基質を chain:残基名:番号 で、その総電荷を 1 つの数で渡します。

```bash
pdb2reaction extract -i complex.pdb -c 'A:GPP:301' -o model.pdb -l -3 --out-json
```

成功すると終了コード 0 で終わり、端末に `[extract] Atoms after truncation: N` と `[extract] Link-H to add: M` が出ます。モデルの原子数は N + M で、そのうち M 個がキャップ水素です。電荷は `[extract] Total active site model charge` の行に出ます。`model.pdb` をビューアで開き、反応に関わる残基が入っているかを確かめてください。次のコマンドには、`model.pdb` とこの総電荷を `-q` で渡します。

### 2. 基質を PDB ファイルで渡す

基質の PDB ファイルを中心にし、残基名ごとに電荷を渡します。

```bash
pdb2reaction extract -i complex.pdb -c substrate.pdb -o model.pdb -l 'GPP:-3,SAM:1'
```

基質のファイルの座標は、複合体の座標と 0.001 Å 以内で一致している必要があります。

### 3. 残基名で選ぶ

残基名を並べると、その名前の残基がすべて中心になります。

```bash
pdb2reaction extract -i complex.pdb -c 'GPP,SAM' -o model.pdb -l 'GPP:-3,SAM:1'
```

### 4. 複数の構造を 1 回で切る

反応物と生成物を 1 つの `-i` の後に並べると、両方が同じ残基・同じキャップになり、1 つのマルチ MODEL の PDB に書き出されます。

```bash
pdb2reaction extract -i complex_R.pdb complex_P.pdb -c 'A:GPP:301,A:SAM:302' \
    -o model_multi.pdb -l 'GPP:-3,SAM:1'
```

入力ごとに別のファイルにするときは、`-o model_R.pdb -o model_P.pdb` を渡します。

(ja-extract-modified-residue)=
### 5. 非標準残基（`--modified-residue`）

Amber の MCPB.py などは、金属に配位する残基に非標準の名前（`HD1`、`HE1`、`CM1`、`AP1`）を付けます。`extract` はこの名前を知らないので、主鎖を切らず、キャップ水素も付けず、次の警告を出します。

```text
[extract] WARNING: Residue HD1 83 may be an amino acid (has N, CA, C, O) but is not recognized as a standard residue name. Backbone truncation was not applied. Consider preparing the active site model manually.
```

この名前をアミノ酸として登録します。`NAME:charge` で電荷を決め、電荷を書かない `NAME` は 0 になります。

```bash
pdb2reaction extract -i complex.pdb -c 'A:SUB:301' -o model.pdb \
    --modified-residue 'HD1,HE1'
```

`NAME:charge` は組み込みの電荷もこの実行に限って上書きします（例: `LYS:0`）。[付録](#アミノ酸)にある名前は、リン酸化残基や D-アミノ酸も含めて組み込み済みなので、登録は要りません。`--modified-residue` で足りないときは、{ref}`モデルを手で組みます <ja-model-setup-manual>`。

---

## 処理の仕組みと計算仕様

1. **中心**: `-c` に基質・補因子・金属を並べます。各項目は残基の指定で、いちばん具体的な `A:TYR:44`（chain:残基名:番号）から、`A:SAM`、`SAM` のような名前、番号、基質の PDB/mmCIF ファイルまで使えます。`--selected-resn` は同じ形で残基を足し、そこからは距離での探索を始めません。
2. **隣の残基**: 中心の原子から `-r`（既定 2.6 Å）以内に原子がある残基を入れます。水は `--no-include-h2o` を付けない限り数え、`--exclude-backbone` ではアミノ酸の主鎖の原子による接触を数えません。さらに、選んだシステインの S–S 結合の相手（S–S ≤ 2.5 Å）、選んだプロリンの N 側の隣、`--exclude-backbone` でないときは、主鎖の原子が中心に触れたアミノ酸とペプチド結合でつながった両隣の残基を足します。
3. **主鎖の切断**: つながったアミノ酸の並びは内側の主鎖を残し、両端が CA で終わるように切ります。前後の残基がモデルに入らなかった残基は CB で切り、側鎖だけを残します。`-c` のアミノ酸も同じ規則で切りますが、ペプチド結合でつながった隣の残基が `-r` 以内に入ってモデルに加わるので、全部の原子が残ります。プロリンは環を残します。`--exclude-backbone` では、ペプチド結合でつながった `-c` のアミノ酸どうしの間を除き、アミノ酸の主鎖の原子をすべて除きます。水とアミノ酸でない残基は切りません。
4. **キャップ水素**: 切断で CA か CB の結合相手が無くなった所（CB–CA、CA–N、CA–C。プロリンは CA–C だけ）に、その炭素から元の結合の向きに 1.09 Å の位置へ水素を置きます。キャップ水素は `TER` の後に、残基 `LKH`・chain `L` の `HETATM` 原子 `HL` として書かれます。
5. **電荷**: アミノ酸とイオンは組み込みの表から、水は 0、ほかの残基は `-l` で渡さない限り 0 とします。

### 電荷の内訳

`-l` には `'GPP:-3,SAM:1'` のような残基名ごとの電荷か、1 つの数を渡します。未知の残基とは、付録でアミノ酸・イオン・水のどれにも挙がっていない残基です。数を渡すと、`-c` の中の未知の残基に均等に割り、`-c` に未知の残基が無ければ、すべての未知の残基に割ります。例 1 では −3 がすべて GPP に入ります。残基名で渡したときは、書かなかった未知の残基は 0 です。端末には、タンパク質・リガンド・イオンの電荷に続いて `Total active site model charge` が出ます。入力が複数のときは、最初の入力の内訳です。

### 複数の構造

入力が複数のときは、構造ごとに残基を選び、その和集合をすべての構造に当てるので、どのモデルも同じ原子・同じキャップになります。座標はモデルごとのものです。端末には、モデルごとに `[extract:multi] Atoms after truncation (model k): N`、全体で 1 回 `[extract:multi] link-H targets common across models: M` が出ます。

(ja-link-hydrogen-and-frozen-atoms)=
### キャップ水素と凍結原子

`opt`・`tsopt`・`freq`・`irc`・`path-opt`・`path-search`・`scan`・`scan2d`・`scan3d`・`all` は、既定の `--freeze-links` でキャップ水素の親原子を固定し、構造最適化や経路探索の間も境界の形を保ちます。`sp` は固定しません。

* **力**: 固定した原子の力を 0 にします。
* **Hessian**: 固定した原子を Hessian から外します。
* **振動解析**: 固定した原子があると、`freq` は動ける原子で PHVA（部分 Hessian 振動解析）を行います。

`--freeze-atoms` と YAML の `geom.freeze_atoms`（1 始まり）で原子を足せ、どの指定も合わせて使われます。{ref}`原子の固定と距離の拘束 <ja-freeze-atoms-and-restraints>` を見てください。

---

## 主な出力ファイル

```text
./
├─ model.pdb     # クラスターモデル。キャップ水素は TER の後
├─ model.cif     # mmCIF の入力か、PDB の桁に収まらない PDB の入力のとき
├─ result.json   # --out-json のとき。最初の出力ファイルと同じディレクトリ
└─ summary.json  # result.json の写し。result.json を読む（--out-json のとき）
```

| 入力 | `-o` | 出力 |
| --- | --- | --- |
| 1 つ | なし | `model.pdb` |
| 複数 | なし | 入力ごとに `model_<入力の名前>.pdb` |
| 複数 | 1 つ | マルチ MODEL の PDB 1 つ |
| 複数 | 入力と同じ数 | 入力ごとに PDB 1 つ |

`-o` がこれ以外の数だとエラーで止まります。出力先の親ディレクトリは自動で作られます。`result.json` には原子数・電荷・使った設定が入ります。各欄は [JSON 出力リファレンス](json-output.md) にあります。

---

## 主な CLI オプション

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス（複数可） | （必須） | タンパク質–リガンドの PDB/mmCIF。1 つの `-i` の後に並べても、`-i` を繰り返してもよい。原子が同じ順に並んでいること |
| `-c, --center` | 文字列 | （必須） | 中心の残基か、基質の PDB/mmCIF ファイル（例: `'A:TYR:44,A:SAM:301'`） |
| `-o, --output` | パス（複数可） | 上の表 | 出力する PDB のパス |
| `-r, --radius` | 浮動小数点数 | `2.6` | 中心の原子からの距離のしきい値（Å）。`0` では距離で隣の残基を足さない（[使用上の注意点](#使用上の注意点)） |
| `--radius-het2het` | 浮動小数点数 | `0`（無効） | C・H 以外の原子どうしの 2 つ目のしきい値（Å） |
| `--selected-resn` | 文字列 | `""` | 距離で探さずに足す残基。`-c` と同じ形 |
| `--include-h2o/--no-include-h2o` | フラグ | `True` | 水（HOH、WAT、H2O、DOD、TIP、TIP3、SOL）を入れる |
| `--exclude-backbone/--no-exclude-backbone` | フラグ | `False` | アミノ酸から主鎖の原子を除く。ペプチド結合でつながった `-c` のアミノ酸どうしの間は残す |
| `--add-linkh/--no-add-linkh` | フラグ | `True` | 切断で CA か CB の相手が無くなった所にキャップ水素を付ける |
| `--modified-residue` | 文字列 | `""` | アミノ酸として扱う残基名。`NAME` か `NAME:charge` |
| `-l, --ligand-charge` | 文字列 | `None` | 未知の残基（リガンド）の電荷の合計か、残基名ごとの電荷（例: `'GPP:-3,SAM:1'`） |
| `--out-json/--no-out-json` | フラグ | `False` | `result.json` と `summary.json` を書き出す |

全オプションの一覧は [自動生成 CLI リファレンス](../reference/commands/extract.md) を参照してください。

---

## 使用上の注意点

* **`-r 0`** では距離で隣の残基を足さず、`-c` と `--selected-resn` の残基に、手順 2 で足す S–S 結合の相手とプロリンの N 側の隣を加えたモデルになります。
* **モデルの大きさ**: {ref}`モデルを広げて <ja-model-setup-larger>`も結果が変わらないことを系ごとに確かめてください。`-r` を大きくすると計算は重くなり、精度が上がるとは限りません。
* **名前はすべての chain に当たる**: `TYR` のような名前は、どの chain の TYR もすべて選び、複数あれば警告を出します。
* **`TYR:44` は chain TYR と読まれる**: 2 つの欄では最初の欄が必ず chain で、2 つ目は番号か名前なので、`A:TYR:44` と書いてください。chain の欄が空の PDB では、名前か番号だけを使います。
* **1 つの list に 1 つの形**: `'SAM,44'` のように名前と番号を混ぜた list はエラーで止まります。
* **キャップ水素は CA と CB だけ**: ほかの切断にはキャップ水素が付きません。そのような切断が非金属の原子どうしの結合にあると、`extract` は警告を出してモデルをそのまま書き、`all` は計算の前に止まります。警告に並ぶ結合・キャップ・電荷を確かめるか、モデルを手で組んでください。
* **どの入力も同じ原子**: 原子の数や並びが違う入力は `[multi] Atom count mismatch` か `[multi] Atom order mismatch` で止まります。
* **元素の欄**: 元素の欄が空だと `extract` は `Element symbols are missing in '…'` で止まるので、先に [`add-elem-info`](add-elem-info.md) を実行してください。
* **altLoc（別位置の配座）**: `extract` は残基ごとに 1 つの配座を残します。規則は {ref}`mmCIF と大きな構造 <ja-mmcif-input>` にあります。
* **組み込みの残基名**は Amber/CHARMM の命名です。PDB の残基が別の化合物と同じ名前を持つときは、`--modified-residue NAME:charge` で意図する電荷を渡してください。

---

## 関連ドキュメント

* [クラスターモデルの組み方](model-setup.md) — モデルを削る・広げる、原子を固定する
* [all](all.md) — 一括のワークフロー。`-c` で `extract` を実行する
* [path-search](path-search.md) — 切り出したモデルでの最小エネルギー経路（MEP）の探索
* [scan](scan.md) — 切り出したモデルでの段階的なスキャン
* [add-elem-info](add-elem-info.md) — 切り出しの前に元素の欄を埋める
* [共通オプションと残基・原子の指定](cli-conventions.md) — 残基の指定と電荷
* [トラブルシューティング](troubleshooting.md) — 切り出しのエラー
* [用語集](glossary.md) — 活性部位モデル、クラスターモデル、キャップ水素

## 付録: PDB 命名規則と参照リスト

この付録は、非標準の残基名・原子名のために `extract` が残基の分類や電荷を誤るときに使います。標準の PDB の名前なら読み飛ばしてかまいません。

```{important}
`extract` は、アミノ酸・イオン・水・主鎖の原子を PDB の残基名と原子名で見分けます。入力は標準の PDB 化学成分の名前に従う必要があります。
```

### アミノ酸

アミノ酸として扱う残基名と、その公称電荷です。これらの残基だけが、主鎖の切断・キャップ水素・アミノ酸の電荷の対象になります。

**標準 20 アミノ酸**（生理的 pH での電荷）：
- 中性: `ALA`, `ASN`, `CYS`, `GLN`, `GLY`, `HIS`, `ILE`, `LEU`, `MET`, `PHE`, `PRO`, `SER`, `THR`, `TRP`, `TYR`, `VAL`
- 正電荷 (+1): `ARG`, `LYS`
- 負電荷 (−1): `ASP`, `GLU`

**プロトン化/互変異性体**（Amber/CHARMM 形式）：
- `HIP`（+1、完全プロトン化 His）、`HID`（0、Nδプロトン化 His）、`HIE`（0、Nεプロトン化 His）
- `ASH`（0、中性 Asp）、`GLH`（0、中性 Glu）、`LYN`（0、中性 Lys）、`ARN`（0、中性 Arg）
- `TYM`（−1、脱プロトン化 Tyr フェノラート）

**リン酸化残基：**
- 二価陰イオン（−2）: `SEP`, `TPO`, `PTR`
- 一価陰イオン（−1）: `S1P`, `T1P`, `Y1P`
- リン酸化 His（phosaa19SB）: `H1D`（0）、`H2D`（−1）、`H1E`（0）、`H2E`（−1）

**システイン変異体：**
- `CYX`（0、ジスルフィド）、`CSD`（−1、スルフィン酸）
- `OCS`（−1、システイン酸）、`CYM`（−1、脱プロトン化 Cys）

**リシン変異体/カルボキシル化：**
- `MLY`（+1）、`KCX`（−1、Nz-カルボン酸）

**D-アミノ酸**（19 残基）：
- `DAL`, `DAR`, `DSG`, `DAS`, `DCY`, `DGN`, `DGL`, `DHI`, `DIL`, `DLE`, `DLY`, `MED`, `DPN`, `DPR`, `DSN`, `DTH`, `DTR`, `DTY`, `DVA`

**その他の修飾残基：**
- `CGU`（−2、γ-カルボキシグルタミン酸）、`CGA`（−1）、`PCA`（0、ピログルタミン酸）、`MSE`（0、セレノメチオニン）、`OMT`（0、メチオニンスルホン）、`HYP`（0、ヒドロキシプロリン）
- その他（いずれも 0）: `ASA`, `CIR`, `FOR`, `MVA`, `IIL`, `AIB`, `HTN`, `SAR`, `NMC`, `PFF`, `NFA`, `ALY`, `AZF`, `CNX`, `CYF`

**N 末端変異体**（接頭辞 `N`）: `NALA`（+1）、`NARG`（+2）、`NASP`（0）、`NGLU`（0）、`NLYS`（+2）など、および `ACE`（0）、`NTER`（+1、汎用）

**C 末端変異体**（接頭辞 `C`）: `CALA`（−1）、`CARG`（0）、`CASP`（−2）、`CGLU`（−2）、`CLYS`（0）など、および `NHE`（0）、`NME`（0）、`CTER`（−1、汎用）

接頭辞 `N`・`C` の Amber の名前は標準の残基として読み（`NALA` → `ALA`）、末端の電荷は、モデルに N 末端の H1〜H3（プロリンは H2 と H3）か OXT が残るときだけ数えます。

### 主鎖の原子

アミノ酸の主鎖として扱う原子名です。`--exclude-backbone` では、ペプチド結合でつながった `-c` のアミノ酸どうしの間を除き、これらの原子を除きます。

```
N, C, O, CA, OXT, H, H1, H2, H3, HN, HA, HA2, HA3
```

### イオン

イオンとして扱う残基名と、その形式電荷です。

| 電荷 | 残基名 |
|------|--------|
| +1 | `LI`, `NA`, `K`, `RB`, `CS`, `TL`, `AG`, `CU1`, `K+`, `NA+`, `NH4`, `H3O+`, `H3O`, `HE+`, `HZ+` |
| +2 | `MG`, `CA`, `SR`, `BA`, `MN`, `FE2`, `CO`, `NI`, `CU`, `ZN`, `CD`, `HG`, `PB`, `BE`, `PD`, `PT`, `SN`, `RA`, `YB2`, `V2+` |
| +3 | `FE`, `AU3`, `AL`, `GA`, `IN`, `CE`, `CR`, `DY`, `EU`, `EU3`, `ER`, `GD3`, `LA`, `LU`, `ND`, `PR`, `SM`, `TB`, `TM`, `Y`, `PU` |
| +4 | `U4+`, `TH`, `HF`, `ZR` |
| −1 | `F`, `CL`, `BR`, `I`, `CL-`, `IOD` |

### 水

水として扱う残基名です。既定で入り（`--include-h2o`）、電荷は 0 です。

```
HOH, WAT, H2O, DOD, TIP, TIP3, SOL
```
