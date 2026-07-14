# mmCIF と大規模構造

`pdb2reaction` は PDB を受け付ける座標workflowで `.cif` / `.mmcif` も
受け付けます。複数文字chain ID、4桁を超える残基番号、5桁を超える
atom serial、または10,000残基以上のmodelではmmCIFを使用してください。

## 内部変換

計算coreはPDBのままです。入力時に、最初のcoordinate modelを読み、
残基単位でaltLocを解決し、一時的な1文字chainと1–9999の残基番号へ
再割当てします。元のauth chain、残基番号、挿入code、残基名、原子名、
occupancy、B-factor、元素、formal chargeは保持します。

10,000残基以上、99,999原子以上、hybrid-36、または認識可能な固定幅
overflowを持つPDBにも同じbridgeを使用します。内部PDBのIDは実装詳細で、
解析・報告には使用しないでください。

## 出力

`--convert-files` が有効（デフォルト）の場合、座標を生成するworkflowは
stage間で必要なPDBを保持し、元IDを復元したCIF companionも生成します。
`--no-convert-files` では companion 変換を行いません。

```text
final_geometry.pdb  # pipeline stage間で使う正規化表現
final_geometry.cif  # 元のchain/residue identityを復元
```

multi-frame trajectoryは `_atom_site.pdbx_PDB_model_num` を使います。
source CIFの結晶学的refinement categoryは複製せず、計算座標に必要な
atom-site tableを出力します。

## selector

```bash
# LONG_CHAIN内のSAMをすべて選択
pdb2reaction extract -i complex.cif -c 'LONG_CHAIN:SAM' -o model.pdb

# 残基番号10001のSAMだけを選択
pdb2reaction extract -i complex.cif -c 'LONG_CHAIN:SAM:10001' -o model.pdb

# chainと残基番号だけでも指定可能
pdb2reaction extract -i complex.cif -c 'LONG_CHAIN:10001' -o model.pdb
```

`CHAIN:RESNAME` はそのchain内の全matchを選び、複数match時に警告します。
1残基だけを意図する場合は `:RESSEQ` を追加します。scan原子selectorは
`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` の4fieldでchainを指定できます。挿入codeが必要な場合は
残基番号に続けます（例: `A:SAM:12B:C1`）。chain IDは大文字・小文字を
区別します（`A` と `a` は別chainとして扱われます）。

## 制約

- 内部bridgeは最大619,938残基。
- R/IM/P入力のatom identity/order一致は引き続き必須。format変換はatom
  mappingを推定しません。
- 入力は最初のcoordinate modelのみ計算します。
- `fix-altloc` / `add-elem-info` はPDB utilityです。mmCIFのaltLocと
  `type_symbol` は変換時に処理します。

[extract](extract.md) と [CLI規約](cli-conventions.md) も参照してください。
