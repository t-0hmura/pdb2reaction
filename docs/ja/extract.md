# `extract`

## 概要
タンパク質-基質複合体から結合ポケット（活性部位）を自動抽出します。このツールは化学的に適切な残基選択（距離カットオフ + ジスルフィド結合、PRO隣接などのヒューリスティクス）を適用し、側鎖/主鎖セグメントを切断し、オプションでリンク水素を追加し、単一構造またはアンサンブルを処理できます。

## 使用法
```bash
pdb2reaction extract -i COMPLEX.pdb [COMPLEX2.pdb ...]
                     -c SUBSTRATE_SPEC
                     [-o POCKET.pdb [POCKET2.pdb ...]]
                     [--radius Å] [--radius-het2het Å]
                     [--include-H2O {True\|False}]
                     [--exclude-backbone {True\|False}]
                     [--add-linkH {True\|False}]
                     [--selected-resn LIST]
                     [--ligand-charge MAP_OR_NUMBER]
                     [--verbose {True\|False}]
```

### 例
```bash
# 最小（ID基準の基質）+ 明示的な総リガンド電荷
pdb2reaction extract -i complex.pdb -c '123' -o pocket.pdb --ligand-charge -3

# PDBとして提供される基質; 残基名ごとの電荷マッピング
pdb2reaction extract -i complex.pdb -c substrate.pdb -o pocket.pdb --ligand-charge 'GPP:-3,SAM:1'

# 名前基準の基質選択（すべてのマッチを含む）
pdb2reaction extract -i complex.pdb -c 'GPP,SAM' -o pocket.pdb --ligand-charge 'GPP:-3,SAM:1'

# ヘテロ-ヘテロ近接を有効にした複数構造から単一のマルチMODEL出力
pdb2reaction extract -i complex1.pdb complex2.pdb -c 'GPP,SAM' -o pocket_multi.pdb --radius-het2het 2.6 --ligand-charge 'GPP:-3,SAM:1'

# ヘテロ-ヘテロ近接を有効にした複数構造から複数出力
pdb2reaction extract -i complex1.pdb complex2.pdb -c 'GPP,SAM' -o pocket1.pdb pocket2.pdb --radius-het2het 2.6 --ligand-charge 'GPP:-3,SAM:1'
```

## ワークフロー
### 残基包含
- `-c/--center` からの基質残基を常に含める
- **標準カットオフ（`--radius`、デフォルト2.6 Å）:**
  - `--exclude-backbone False` の場合、カットオフ内の任意の原子が残基を対象にする
  - `--exclude-backbone True` の場合、アミノ酸残基は**非主鎖**原子（N/H*/CA/HA*/C/O以外）で基質に接触する必要がある
- **独立したヘテロ-ヘテロカットオフ（`--radius-het2het`）:** 基質ヘテロ原子（非C/H）がタンパク質ヘテロ原子の指定Å以内にある場合に残基を追加
- **水処理:** HOH/WAT/H2O/DOD/TIP/TIP3/SOLはデフォルトで含まれる（`--include-H2O True`）
- **強制包含:** `--selected-resn` はチェーン/挿入コード付きIDを受け入れる（例: `A:123A`）
- **近傍セーフガード:**
  - `--exclude-backbone False` で主鎖原子が基質に接触した場合、ペプチド隣接のN/C側残基（C–N ≤ 1.9 Å）を自動的に含める。末端はN/H*またはC/O/OXTのキャップを保持。
  - ジスルフィド結合（SG–SG ≤ 2.5 Å）は両方のCysを包含。
  - 非末端PRO残基は常にN側隣接残基を含め、主鎖除去後もCAを保持します。`--exclude-backbone True` の場合は隣接残基のC/O/OXTを残し、ペプチド結合を維持。

### 切断/キャッピング
- 孤立残基は側鎖原子のみを保持; アミノ酸主鎖原子（N, CA, C, O, OXT + N/CA水素）はPRO/HYP保護を除いて除去
- 連続ペプチドストレッチは内部主鎖原子を保持; 末端キャップ（N/H*またはC/O/OXT）のみ除去
- TERを認識し、チェーン切断を跨ぐキャッピングは行わない
- `--exclude-backbone True` の場合、**非基質**アミノ酸の主鎖原子を除去（PRO/HYP保護とPRO近傍保持は適用）
- 非アミノ酸残基は主鎖様原子名を持つ原子を失わない

### リンク水素（`--add-linkH True`）
- 切断された結合ベクトルに沿って1.09 Åで炭素のみのリンク水素を追加（CB–CA、CA–N、CA–C; PRO/HYPはCA–Cのみ）
- `TER` の後に残基 `LKH`（チェーン `L`）の連続した `HETATM` レコードとして `HL` という名前で挿入し、シリアル番号は本体の続きになります
- マルチ構造モードでは全モデルで同じ結合にキャップを付け、座標はモデルごとに保持されます

### 電荷サマリー（`--ligand-charge`）
- アミノ酸と一般的なイオンは内部辞書から電荷を取得; 水はゼロ
- 未知残基は `--ligand-charge` が総電荷または残基名ごとのマッピング（例: `GPP:-3,SAM:1`）を提供しない限りデフォルトで0。総電荷が与えられた場合は未知基質残基に配分され、未知基質が無い場合は未知残基全体に配分されます。
- verboseモードが有効な場合、モデル#1の電荷サマリー（タンパク質/リガンド/イオン/総計）がログに記録されます。


### マルチ構造アンサンブル
- 複数の入力PDBを受け付けます（先頭/末尾で原子順序が一致することを検証）。各構造は独立に処理され、選択残基の**和集合**が全モデルに適用されるため、出力の一貫性が保たれます。
- 出力ポリシー:
  - `-o` なし & 複数入力 → 構造ごとに `pocket_<original_basename>.pdb`。
  - `-o` を1つだけ指定 → 単一のマルチMODEL PDB。
  - 入力数と同数の `-o` を指定 → 入力ごとに個別PDB。
- 診断ログにモデルごとの生/保持原子数と残基IDを出力します。

### 基質指定（`-c/--center`）
- **PDBパス**: 座標が先頭入力と完全一致（許容誤差 1e-3 Å）。残基IDは他構造へ伝播。
- **残基ID**: `'123,124'`, `'A:123,B:456'`, `'123A'`, `'A:123A'`（挿入コード対応）。
- **残基名**: カンマ区切り（大文字小文字は無視）。同名残基が複数ある場合は**すべて**含め、警告を出力。

## CLIオプション

> **注記:** 表示されているデフォルト値は、オプション未指定時に使用される内部デフォルトです。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 1つ以上のタンパク質-リガンドPDBファイル（同一の原子順序が必要） | 必須 |
| `-c, --center SPEC` | 基質指定（PDBパス、残基ID、または残基名） | 必須 |
| `-o, --output PATH...` | ポケットPDB出力。1パス ⇒ マルチMODEL、Nパス ⇒ 入力ごと | 自動（`pocket.pdb` または `pocket_<input>.pdb`） |
| `-r, --radius FLOAT` | 包含のための原子-原子距離カットオフ（Å） | `2.6` |
| `--radius-het2het FLOAT` | 独立したヘテロ-ヘテロカットオフ（Å、非C/H） | `0.0` |
| `--include-H2O {True\|False}` | HOH/WAT/H2O/DOD/TIP/TIP3/SOL水を含める | `True` |
| `--exclude-backbone {True\|False}` | 非基質アミノ酸の主鎖原子を除去 | `True` |
| `--add-linkH {True\|False}` | 切断された結合に1.09 Åで炭素のみのリンク水素を追加 | `True` |
| `--selected-resn TEXT` | 強制包含残基（オプションのチェーン/挿入コード付きID） | `""` |
| `--ligand-charge TEXT` | 総電荷または残基名ごとのマッピング（例: `GPP:-3,SAM:1`） | _None_ |
| `-v, --verbose` | INFOレベルログを出力 | `true` |

## 出力
```
<output>.pdb  # TERレコード後にオプションのリンク水素を含むポケットPDB
               # 単一入力 → デフォルトでpocket.pdb
               # -oなしの複数入力 → 構造ごとにpocket_<original_basename>.pdb
               # 複数入力で1つの-oパス → 単一のマルチMODEL PDB
```
- verboseモードが有効な場合、モデル#1の電荷サマリー（タンパク質/リガンド/イオン/総計）がログに記録
- 出力ディレクトリは自動作成されません。必要に応じて事前に作成してください。
- API利用（`extract_api`）では `{'outputs': [...], 'counts': [...], 'charge_summary': {...}}` を返します。


## 注意事項
- `--radius` のデフォルトは 2.6 Å。`0` を指定すると空選択を避けるため内部で 0.001 Å に調整されます。`--radius-het2het` もデフォルトは無効で、`0` 指定時は 0.001 Å に調整されます。
- `--include-H2O False` で水を除外できます。
- 主鎖トリミングとキャッピングはチェーン切断や PRO/HYP 保護を尊重し、非アミノ酸残基は主鎖様の原子名を保持します。
- リンク水素は炭素切断のみで挿入され、アンサンブルモードでは同一結合パターンを全モデルで再利用します。
- INFOログに残基選択、切断数、電荷内訳が要約されます。