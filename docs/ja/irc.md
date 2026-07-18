# `irc`

遷移状態（TS）から反応物・生成物方向へ EulerPC（Euler Predictor-Corrector）ベースの固有反応座標（IRC）積分を実行します。`tsopt` で最適化・検証済みの TS 構造を出発点に経路を追跡し、端点接続性（R ↔ TS ↔ P）を確認します。デフォルトで前方・後方の両方向を実行します。`--no-backward`（または `--no-forward`）で一方向のみをたどります。Hessian のデフォルトは有限差分です。解析 autograd の速度とメモリ量は backend・model・系サイズ・precision・GPU に依存するため、対象環境で検証した場合にだけ明示的に選択してください。mmCIF入力は内部PDBで計算し、出力CIFに元IDを復元します。XYZ/GJF入力では`--ref-pdb`にPDBまたはmmCIF topologyを指定できます。一般的な手順は `tsopt` → `irc` です。

## 実行例

コマンド形式:

```bash
pdb2reaction irc -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [--workers N] [--workers-per-node N] [-m 2S+1] \
 [--max-cycles N] [--step-size Δs] [--never-stop/--no-never-stop] [--root k] \
 [--forward/--no-forward] [--backward/--no-backward] \
 [--freeze-links/--no-freeze-links] \
 [--out-dir DIR] [--config FILE] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE] \
 [--hessian-calc-mode Analytical|FiniteDifference] \
 [--show-config] [--dry-run]
```

両方向の基本実行:

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc
```

順方向のみ、有限差分Hessian、大きいステップサイズ:

```bash
# 順方向のみ、有限差分Hessian、大きいステップサイズ
pdb2reaction irc -i ts.pdb -q 0 -m 1 --no-backward \
 --step-size 0.2 --hessian-calc-mode FiniteDifference --out-dir ./irc_fd/
```

解析Hessianを明示的に使用:

```bash
# UMAで解析Hessianを明示する場合はworkersを1にする
pdb2reaction irc -i ts.pdb -q 0 -m 1 \
 --hessian-calc-mode Analytical --out-dir ./result_irc_analytical
```

分岐が1〜2 stepですぐ停止する場合は、まずEulerPCの最大stepを小さくします。
再試行の目安は`0.05` Bohrです。

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --step-size 0.05 \
 --out-dir ./result_irc_small_step
```

小さな数値的な山や平坦な肩を越えて追跡する必要がある場合は
`--never-stop`を追加します。このopt-inモードはenergy上昇／変化量による
停止だけを無視し、gradient/integratorの収束と`--max-cycles`は維持します。

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --step-size 0.05 --never-stop \
 --max-cycles 250 --out-dir ./result_irc_continue
```

## 処理の流れ

1. **入力準備** – 共通bridgeがPDB/mmCIFと`geom_loader`対応形式を受け入れます。参照topologyがある場合は軌跡をPDBへ変換し、bridge入力では元IDを復元したCIFも生成します。`--freeze-links`はcap水素の親原子を凍結します。
2. **EulerPC 積分** – EulerPC 予測子-修正子積分器が遷移状態から IRC 経路をたどります。`--forward`/`--backward` フラグに従って順方向および/または逆方向の分岐を実行します。初期モード選択前に、デフォルトの`constrained`が凍結anchorを動かさない全系並進・回転だけを除去します。各ステップでは、質量加重の最急降下方向に沿う Euler 予測子（勾配は現在のHessianを用いた 2 次の Taylor 展開で近似）を適用し、続いて距離加重補間（DWI）面上で修正 Bulirsch–Stoer 修正子を適用します。
3. **軌跡出力** – 完了済み、順方向、逆方向のIRC軌跡をXYZへ書きます。参照topologyと`--convert-files`があればPDB、bridge入力ならCIF companionも生成します。

## 出力

```
out_dir/ (デフォルト:./result_irc/)
├─ <prefix>finished_irc_trj.xyz   # 完全な IRC 軌跡
├─ <prefix>finished_irc.pdb       # 参照 PDB が利用可能な場合の軌跡に対応する PDB（変換有効時）
├─ <prefix>finished_irc.cif       # bridge入力。元IDを復元
├─ <prefix>forward_irc_trj.xyz    # 順方向分岐が実行された場合
├─ <prefix>forward_irc.pdb        # 順方向分岐に対応する PDB（同条件）
├─ <prefix>forward_irc.cif        # bridge入力
├─ <prefix>backward_irc_trj.xyz   # 逆方向分岐が実行された場合
├─ <prefix>backward_irc.pdb       # 逆方向分岐に対応する PDB（同条件）
└─ <prefix>backward_irc.cif       # bridge入力
```

`irc.prefix`が空でない場合、EulerPCはファイル名との間に`_`を1つ補います。たとえば
`prefix: trial`は`trial_finished_irc_trj.xyz`を生成し、`result.json.files`にも
正規化後の名前を記録します。

コンソールには確定済みの `geom`/`calc`/`irc` 設定と実行時間の要約が表示されます。

## CLI オプション

以下の表は説明が必要なオプションを扱います。全フラグの一覧は生成された [コマンドリファレンス](../reference/commands/index.md) にあります。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる遷移状態構造 | 必須 |
| `-q, --charge INT` | 総電荷。明示的な `-q` は YAML `calc.charge` と `--ligand-charge/-l` より優先。省略時は YAML、残基からの導出、`.gjf` metadata の順に解決 | YAML/テンプレート/導出がない限り必須 |
| `-l, --ligand-charge TEXT` | 単一の整数（例: `-1`）でリガンド総電荷を指定するか、残基別マッピング（例: `GPP:-3,SAM:1`）で PDB/mmCIF 残基電荷から全系の電荷を導出。`-q` 省略時に使用（PDB/mmCIF 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `--workers INT` | UMA 予測器の並列度。`workers > 1` と明示的な解析 Hessian は併用できないため、`workers = 1` または有限差分を使用。{ref}`ja-workers-analytical-error` を参照 | `1` |
| `--workers-per-node INT` | ノードあたりのワーカー数。並列予測器に渡されます | `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1）。明示的な `-m` は YAML `calc.spin` より優先し、省略時は YAML、`.gjf`、`1` の順に解決 | YAML/`.gjf`/`1` |
| `--max-cycles INT` | 最大 IRC ステップ。明示値は YAML `irc.max_cycles` より優先 | `125` |
| `--step-size FLOAT` | ステップ長（Bohr、非質量加重デカルト座標）。明示値は YAML `irc.step_length` より優先 | `0.10` |
| `--never-stop/--no-never-stop` | 一時的なenergy上昇／平坦化による停止を無視し、小さな肩を越えて追跡します。gradient/integratorの収束と`max_cycles`上限は無効化しません | `False` |
| `--root INT` | 射影Hessianの固有値を**昇順**（最も負の値を先頭）に並べたときの**0 始まり**のインデックス。初期 IRC 変位に使用するモードを指定します。虚振動が 1 個だけの妥当な TS では `--root 0`（唯一の負の固有値）のままにしてください。`--root 1`、`--root 2` などは、活性な虚モードがより負のスプリアス（疑似）モードよりも上位にランクされていることが分かっている場合にのみ使用します。明示値は YAML `irc.root` より優先 | `0` |
| `--forward/--no-forward` | 順方向分岐を実行。明示 toggle は YAML `irc.forward` より優先 | `True` |
| `--backward/--no-backward` | 逆方向分岐を実行。明示 toggle は YAML `irc.backward` より優先 | `True` |
| `--freeze-links/--no-freeze-links` | PDB/mmCIF トポロジー用、キャップ H 親を凍結（`geom.freeze_atoms` にマージ）。詳細は [extract](extract.md) を参照 | `True` |
| `--freeze-atoms TEXT` | 凍結する原子の 1 始まりインデックスをカンマ区切りで明示的に指定（例: `'1,3,5'`）。`--freeze-links` と併用可、任意の入力形式に適用 | _None_ |
| `--tr-projection [constrained\|legacy-active]` | 初期frozen／partial Hessianの剛体モード処理。`constrained`は凍結anchorを尊重し、`legacy-active`はisolated-active比較専用 | `constrained` |
| `-o, --out-dir TEXT` | 出力ディレクトリ。明示値は YAML `irc.out_dir` より優先 | `./result_irc/` |
| `--convert-files/--no-convert-files` | 参照PDB/mmCIF topologyがある場合の XYZ/TRJ → PDB/CIF を切り替え | `True` |
| `--ref-pdb FILE` | XYZ/GJF入力に使用する参照PDBまたはmmCIF topology | _None_ |
| `--hessian-calc-mode CHOICE` | MLIP Hessianモード。明示値は YAML `calc.hessian_calc_mode` より優先 | `FiniteDifference` |
| `--config FILE` | 明示 CLI 適用前に読み込むベース YAML | _None_ |
| `--show-config/--no-show-config` | 解決済み YAML レイヤー/設定を表示して続行 | `False` |
| `--out-json/--no-out-json` | `out_dir` に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照 | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--solvent TEXT` | xTB 暗黙溶媒（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画のみ表示 | `False` |

## YAML 設定

設定の優先順位は {ref}`CLI 規約: 設定の優先順位 <ja-configuration-precedence>` を参照してください。

`geom`、`calc`、`irc` の各セクションは [YAML リファレンス](yaml-reference.md) の正規定義から変更ありません: [`geom`](yaml-reference.md#geom)、[`calc`](yaml-reference.md#calc)、[`irc`](yaml-reference.md#irc-section) を参照してください。`--freeze-links` は PDB/mmCIF トポロジーで `geom.freeze_atoms` を拡張し、`--hessian-calc-mode` と CLI の charge/spin 値はマージ済み `calc` ブロックを補完します。

**`irc` 固有の強制上書き**（YAML/CLI マージ後に YAML 値を無視して適用）:

```yaml
geom:
 coord_type: cart # irc では cart に強制（YAML 値は無視）
calc:
 return_partial_hessian: true # irc では true に強制（partial Hessian、active-DOF 処理）
```

## 終了コード

終了コードは CLI 規約の {ref}`ja-exit-codes` を参照。

## 注意事項

- MLIP バックエンド（デフォルト: UMA）は IRC 全体で再利用されます。`step_length` を大きくし過ぎると EulerPC が不安定になることがあります。ほぼ直ちに停止する分岐は、まず小さい `--step-size`（例: `0.05`）で再試行してください。
- `--never-stop` は意図せぬ追跡と計算時間増大を避けるためデフォルトOFFです。小さな山／肩には有効ですが、最寄りの化学的basinを通過する可能性もあるため軌跡と端点接続を確認してください。
- `--freeze-links` が有効な場合、キャップ水素の親原子が自動的に凍結されます（{ref}`キャップ水素と凍結原子 <ja-link-hydrogen-and-frozen-atoms>` を参照）。
- `result.json["rigid_projection"]`にtreatment、effective rank、初期Hessianのsourceとshapeを記録します。詳細は[凍結原子](freeze-atoms.md#凍結境界での剛体モード)を参照してください。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — よくある失敗モードの詳細な対処
- [tsopt](tsopt.md) — IRC 実行前に TS を最適化
- [freq](freq.md) — 完全な振動解析と熱化学補正
- [opt](opt.md) — IRC 端点を真の極小に最適化
- [all](all.md) — tsopt 後に IRC を実行する一気通貫ワークフロー
- [YAML リファレンス](yaml-reference.md) — `irc` の完全な設定オプション
- [用語集](glossary.md) — IRC（固有反応座標）の定義
