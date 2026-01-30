# `all` サブコマンド

## 概要
`pdb2reaction all` は、**すべてのパイプラインステージ**を統括する傘コマンドです: ポケット抽出 → オプションの段階的UMAスキャン → 再帰的MEP探索（`path_search`、GSM/DMF） → 全系マージ → オプションのTS最適化 + IRC（`tsopt`） → オプションの振動解析（`freq`） → オプションの一点DFT（`dft`）。このコマンドは複数構造アンサンブルを受け入れ、単一構造スキャンを順序付けられた中間体に変換し、TSOPTのみのポケットワークフローにフォールバックできます。すべての下流ツールは単一のCLIサーフェスを共有するため、1つの呼び出しから長い反応キャンペーンを調整できます。すべてのステージにわたるフォーマット対応XYZ/TRJ → PDB/GJF変換は、共有の `--convert-files {True\|False}` フラグ（デフォルトで有効）によって制御されます。

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
   - VRAMが十分な場合は `--hessian-calc-mode` を `Analytical` に設定することを強く推奨

6. **TSOPTのみモード**（単一入力、`--tsopt True`、`--scan-lists` なし）
   - MEP/マージステージをスキップ。ポケット（または抽出がスキップされた場合は完全入力）で `tsopt` を実行し、EulerPC IRCを実行

## CLIオプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 反応順序の2つ以上の完全構造（`--scan-lists` または `--tsopt True` を使用する場合のみ単一入力可） | 必須 |
| `-c, --center TEXT` | 基質指定（PDBパス、残基ID、または残基名） | 抽出に必須 |
| `--out-dir PATH` | トップレベル出力ディレクトリ | `./result_all/` |
| `-r, --radius FLOAT` | ポケット包含カットオフ（Å） | `2.6` |
| `--ligand-charge TEXT` | 総電荷または未知残基の残基ごとのマッピング（推奨） | _None_ |
| `-q, --charge INT` | 総電荷を強制上書き | _None_ |
| `-m, --mult INT` | すべての下流ステップに転送されるスピン多重度 | `1` |
| `--freeze-links BOOLEAN` | ポケットPDBでリンク親を凍結 | `True` |
| `--max-nodes INT` | セグメントごとのMEP内部ノード | `10` |
| `--max-cycles INT` | MEP最大最適化サイクル | `300` |
| `--climb BOOLEAN` | 各ペアの最初のセグメントでTSクライミングを有効化 | `True` |
| `--opt-mode [light\|heavy]` | オプティマイザープリセット（light → LBFGS/Dimer、heavy → RFO/RSIRFO） | `light` |
| `--dump BOOLEAN` | MEP軌跡をダンプ | `False` |
| `--convert-files {True\|False}` | テンプレート利用可能時のXYZ/TRJ → PDB/GJFコンパニオンのグローバルトグル | `True` |
| `--thresh TEXT` | MEPおよび単一構造最適化の収束プリセット | `gau` |
| `--args-yaml FILE` | 各サブコマンドに転送されるYAML | _None_ |
| `--tsopt BOOLEAN` | 反応セグメントごとにTS最適化 + IRCを実行、またはTSOPTのみモードを有効化 | `False` |
| `--thermo BOOLEAN` | R/TS/Pで振動解析を実行しUMA Gibbsダイアグラムを構築 | `False` |
| `--dft BOOLEAN` | R/TS/Pで一点DFTを実行 | `False` |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | tsoptとfreqに転送される共有UMAヘシアンエンジン | `FiniteDifference` |
| `--scan-lists TEXT...` | 抽出されたポケットでの段階的スキャンを記述するPythonライクなリスト（単一入力実行のみ） | _None_ |

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

## YAML設定（`--args-yaml`）
同じYAMLファイルが**すべての**呼び出されるサブコマンドにそのまま転送されます。各ツールは独自のドキュメントに記載されているセクションを読み取ります:

- [`path_search`](path_search.md): `geom`, `calc`, `gs`, `opt`, `sopt`, `bond`, `search`
- [`scan`](scan.md): `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`, `bond`
- [`tsopt`](tsopt.md): `geom`, `calc`, `opt`, `hessian_dimer`, `rsirfo`
- [`freq`](freq.md): `geom`, `calc`, `freq`, `thermo`
- [`dft`](dft.md): `dft`

すべてのYAMLオプションの完全なリファレンスについては、[YAML設定リファレンス](yaml-reference.md) を参照してください。
