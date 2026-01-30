# `extract` サブコマンド

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

### 切断/キャッピング
- 孤立残基は側鎖原子のみを保持; アミノ酸主鎖原子（N, CA, C, O, OXT + N/CA水素）はPRO/HYP保護を除いて除去
- 連続ペプチドストレッチは内部主鎖原子を保持; 末端キャップ（N/H*またはC/O/OXT）のみ除去
- 非アミノ酸残基は主鎖様原子名を持つ原子を失わない

### リンク水素（`--add-linkH True`）
- 切断された結合ベクトルに沿って1.09 Åで炭素のみのリンク水素を追加（CB–CA、CA–N、CA–C）
- `TER` の後に残基 `LKH`（チェーン `L`）の連続した `HETATM` レコードとして `HL` という名前で挿入

### 電荷サマリー（`--ligand-charge`）
- アミノ酸と一般的なイオンは内部辞書から電荷を取得; 水はゼロ
- 未知残基は `--ligand-charge` が総電荷または残基名ごとのマッピング（例: `GPP:-3,SAM:1`）を提供しない限りデフォルトで0

## CLIオプション
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
