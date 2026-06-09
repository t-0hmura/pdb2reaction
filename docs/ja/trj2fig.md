# `trj2fig`

XYZ 軌跡のコメント行に書かれたエネルギーを読み取り（または MLIP バックエンドで再計算し）、静的・インタラクティブ図と CSV として出力します。`opt` / `scan` / `path-opt` / `path-search` / `irc` などが出力した軌跡のエネルギーを可視化するときに使います。デフォルトでは各フレームのコメント行に含まれる Hartree エネルギーを読み取り、kcal/mol または Hartree に変換して、必要に応じて基準フレームに対する相対値（ΔE）にします。基準は最初のフレームまたは手動指定（`-r`）で、`--reverse-x` 使用時は最後のフレームが基準になります。`-q/--charge` と `-m/--multiplicity` を与えると、コメント行ではなく MLIP バックエンド（デフォルト: UMA、`-b/--backend` で選択可能）で各フレームのエネルギーを再計算します。

## 実行例

コマンド形式:

```bash
pdb2reaction trj2fig -i TRAJECTORY.xyz [-o OUTPUTS...] [-q CHARGE] [-m MULT] [options]
```

```bash
# デフォルトPNG、最初のフレームを基準としたΔE
pdb2reaction trj2fig -i traj.xyz

# CSV + SVG、フレーム5を基準、Hartree表記
pdb2reaction trj2fig -i traj.xyz -o energy.csv energy.svg -r 5 --unit hartree

# 複数形式 + x軸を反転（基準は最後のフレーム）
pdb2reaction trj2fig -i traj.xyz --reverse-x -o energy.png energy.html energy.pdf

# UMAで全フレームのエネルギーを再計算
pdb2reaction trj2fig -i traj.xyz -q 0 -m 1 -o energy.png
```

## 処理の流れ

1. XYZ 軌跡を解析し、各フレームのコメント行からエネルギーを読み取ります。明示的な `E=` / `Energy:` トークンがあればそれを優先し、なければ行内の数値トークンを採用します（`1.5e-3` などの科学表記に対応）。キーなしで数値トークンが複数ある場合は最後の値を採用し、警告を出力します。`-q/-m` がある場合は MLIP バックエンドで再計算した Hartree エネルギーを使用します。エネルギーが取得できない場合は実行を中断します。
2. 参照指定を正規化:
 - `init` → フレーム `0`（`--reverse-x` が有効な場合は最後のフレーム）
 - `None`/`none`/`null` → 絶対エネルギー（参照なし）
 - 整数リテラル → 0 始まりのフレーム番号
3. エネルギーを kcal/mol または Hartree に変換し、参照指定がある場合は参照値を差し引いて ΔE を作成します。
4. Plotly 図を作成（強調ティック、スプライン補間、マーカー、タイトルなし）し、指定された拡張子へ書き出します。
5. 必要に応じて `frame`, `energy_hartree` と、指定単位での ΔE/E 列を含む CSV を出力します。

## 出力
```
<output>.[png|jpg|jpeg|html|svg|pdf] # 指定された拡張子すべてを出力（デフォルト: energy.png）
<output>.csv # CSV指定時のみ
```
- `-o` や位置引数が無い場合は `energy.png` をカレントディレクトリに書き出します。
- CSV には `frame`, `energy_hartree` と、参照がある場合は `delta_kcal`/`delta_hartree`、参照がない場合は `energy_kcal`/`energy_hartree` が含まれます。
- パース失敗や未対応拡張子がある場合、コンソールに診断メッセージが出力されます。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | コメント行にエネルギーを含む XYZ 軌跡 | 必須 |
| `-o, --out PATH` | 複数指定可能な出力ファイル名（`.png/.jpg/.jpeg/.html/.svg/.pdf/.csv`） | `energy.png` |
| _extra arguments_ | オプション後に続く位置引数（`-o` リストにマージ） | _None_ |
| `--unit {kcal,hartree}` | 出力単位 | `kcal` |
| `-r, --reference TEXT` | 参照指定（`init`、`None`、または 0 始まり整数） | `init` |
| `-q, --charge INT` | 総電荷。指定時は MLIP バックエンドでエネルギーを再計算 | _None_ |
| `-m, --multiplicity INT` | スピン多重度（2S+1）。指定時は MLIP バックエンドでエネルギーを再計算 | _None_ |
| `--reverse-x/--no-reverse-x` | x 軸を反転して最後のフレームを左端に表示し、`init` の参照を最後のフレームに変更 | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | エネルギー再計算用 MLIP バックエンド | `uma` |
| `--solvent TEXT` | xTB 暗黙溶媒（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |
| `--out-json/--no-out-json` | 出力の隣に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照 | `False` |

## 注記
- エネルギーが解析できないコメント行はエラーになります（`E=` / `Energy:` トークンと素の数値の優先順位は処理の流れの手順 1 を参照）。
- 未対応の拡張子がある場合は実行が中断されます。`.png` は Plotly の `scale=2` で高解像度出力されます。
- `--reverse-x` は軸の向きと `-r init` の解釈の両方に影響します。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [path-search](path-search.md) -- MEP 軌跡のプロファイリング
- [irc](irc.md) -- IRC 軌跡のプロファイリング
- [all](all.md) -- 一気通貫ワークフロー
