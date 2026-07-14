# `freq`

MLIP バックエンド（デフォルト: UMA、`-b/--backend` で ORB ・ MACE ・ AIMNet2 も選択可能）を用いて振動数と熱化学量（ZPE、ギブズ自由エネルギーなど）を計算します。完全な振動解析（極小点に虚振動数がないこと、TS にちょうど 1 つあること等の確認）が必要な場合や、これらの熱化学補正が必要な場合に使用します。デフォルトは有限差分です。`--hessian-calc-mode Analytical` は変位幅による誤差を避けられますが、速度は backend/model/系に依存し、通常はより多くの accelerator memory を使います。対象環境で速度・メモリ・結果を検証して選択してください。虚振動数は負の値で表示されます。

## 実行例

明示的な電荷とスピンでの最小実行:

```bash
# 明示的な電荷とスピンでの最小実行
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 --out-dir ./result_freq
```

freeze-links + 熱化学ダンプを有効化して実行する:

```bash
# freeze-links + 熱化学ダンプを有効化して実行する
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 --freeze-links --dump --out-dir ./result_freq_phva
```

メモリ量と実行時間を検証した上で解析的 Hessian を明示指定する:

```bash
# 解析的 Hessian を明示指定
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 \
 --hessian-calc-mode Analytical --out-dir ./result_freq_analytical
```

## 処理の流れ

- **構造の読み込みと凍結処理**: 構造は共通structure bridgeを経て `pysisyphus.helpers.geom_loader` で読み込まれます。PDB/mmCIF トポロジーでは `--freeze-links` によりキャップ水素を検出して親原子を凍結し、その結果を `geom.freeze_atoms` にマージします。マージされたインデックスはログに表示され、MLIP バックエンドと PHVA に伝播されます。
- **MLIP バックエンド**: `--hessian-calc-mode` で解析的または有限差分Hessianを選択します。MLIP バックエンドは原子が凍結されている場合、部分（活性）Hessianブロックを返すことがあります。Hessian評価モードの詳細は {ref}`ja-hessian-evaluation` を参照してください。
- **PHVAと剛体モード**: 凍結原子がある場合、固有値解析はactive部分空間内で行います。デフォルトの`constrained`は、すべての凍結anchorを動かさない全系剛体運動だけを除去するため、通常の複数anchorクラスターモデルではeffective rankは通常0です。3N×3N Hessianとactive block Hessianの両方に対応します。詳細は[凍結原子](freeze-atoms.md#凍結境界での剛体モード)を参照してください。
- **モードのエクスポート**: `--max-write` でアニメーション化するモード数を制限できます。モードは値順 (`value`) でソートされ、`--sort abs` を指定すると絶対値順になります。すべての入力に `_trj.xyz` を出力し、`--convert-files` 有効時はtopology入力に`.pdb`、mmCIF／oversized-PDB入力に元IDを復元した`.cif`も出力します。
- **熱化学**: `thermoanalysis` がインストールされている場合、QRRHO に準じたサマリー（EE、ZPE、E/H/G 補正、熱容量、エントロピー）が PHVA 振動数に基づいて出力されます。CLI の圧力（atm）は内部で Pa に変換されます。`--dump` を指定すると `thermoanalysis.yaml` も書き込まれます。
- **性能**: GPU メモリ使用量を抑えるため、Hessianは 1 つだけ保持します。

## 出力

```
out_dir/ (デフォルト:./result_freq/)
├─ mode_XXXX_±freqcm-1_trj.xyz # モードごとのアニメーション
├─ mode_XXXX_±freqcm-1.pdb # PDB テンプレートが存在し変換が有効な場合のみ
├─ mode_XXXX_±freqcm-1.cif # mmCIF/oversized-PDB入力
├─ frequencies_cm-1.txt # 選択したソート順での全振動数リスト
└─ thermoanalysis.yaml # thermoanalysisが利用可能で--dumpを有効にした場合
```
- コンソールには `geom`/`calc`/`freq`/熱化学の設定要約ブロックが出力されます。

出力の確認ポイント:

- `result_freq/frequencies_cm-1.txt`
- `result_freq/mode_*_trj.xyz`
- `result_freq/mode_*.pdb`（PDB 入力かつ変換有効時）
- `result_freq/mode_*.cif`（mmCIF／oversized-PDB入力かつ変換有効時）

## 終了コード

終了コードは CLI 規約の {ref}`ja-exit-codes` を参照。

## CLI オプション

以下の表は説明が必要なオプションを扱います。全フラグの一覧は生成された [コマンドリファレンス](../reference/commands/index.md) にあります（ここで手作業で複製しないでください）。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力bridgeが受け入れる構造（`.pdb` / `.cif` / `.mmcif` / `.xyz` / `.trj` / ...） | 必須 |
| `-q, --charge INT` | 総電荷。省略時は `--ligand-charge` から導出可能。明示的な `-q` は導出値より優先される | `.gjf` テンプレートまたは `--ligand-charge` が提供しない限り必須 |
| `-l, --ligand-charge TEXT` | 単一の整数（例: `-1`）でリガンド総電荷を指定するか、残基別マッピング（例: `GPP:-3,SAM:1`）で PDB/mmCIF 残基電荷から全系の電荷を導出。`-q` 省略時に使用（PDB/mmCIF 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `--workers INT` | UMA 予測器の並列度。`workers > 1` と明示的な解析 Hessian は併用できないため、`workers = 1` または有限差分を使用。{ref}`ja-workers-analytical-error` を参照 | `1` |
| `--workers-per-node INT` | ノードあたりのワーカー数。並列予測器に渡されます | `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| `--freeze-links/--no-freeze-links` | PDB/mmCIF 入力（または `--ref-pdb` 付き XYZ/GJF）。キャップ水素の親を凍結し `geom.freeze_atoms` にマージ。キャップ水素の詳細は [extract](extract.md) を参照 | `True` |
| `--freeze-atoms TEXT` | 凍結する原子の 1 始まりインデックスをカンマ区切りで明示的に指定（例: `'1,3,5'`）。`--freeze-links` と併用可、任意の入力形式に適用 | _None_ |
| `--tr-projection [constrained\|legacy-active]` | PHVAの剛体モード処理。`constrained`は凍結anchorを尊重し、`legacy-active`はisolated-active比較用 | `constrained` |
| `--max-write INT` | エクスポートするモード数 | `10` |
| `--amplitude-ang FLOAT` | モードアニメーション振幅（Å） | `0.8` |
| `--n-frames INT` | モードアニメーションのフレーム数 | `20` |
| `--sort CHOICE` | モード順序: `value`（cm⁻¹）または `abs` | `value` |
| `-o, --out-dir TEXT` | 出力ディレクトリ | `./result_freq/` |
| `--temperature FLOAT` | 熱化学計算の温度（K） | `298.15` |
| `--pressure FLOAT` | 熱化学計算の圧力（atm）。CLI では `--pressure` ですが、対応する YAML キー（`thermo:` 配下）は `pressure_atm`（単位接尾辞付き）です。いずれも atm で指定し、内部で Pa に変換されます | `1.0` |
| `--dump/--no-dump` | `thermoanalysis.yaml` を書き込みます。単体 `freq` のデフォルトは OFF です。`pdb2reaction all --thermo` はこのファイルを内部入力として使用するため常に保持し、`all --no-dump` でも抑止しません（任意の scan/MEP/TS 軌跡には `--no-dump` が適用されます） | `False` |
| `--hessian-calc-mode CHOICE` | MLIP Hessianモード（`Analytical` または `FiniteDifference`） | `FiniteDifference` |
| `--convert-files/--no-convert-files` | topology入力に PDB companion、mmCIF／oversized-PDB bridge入力に元IDを復元したCIF companionを生成（GJFは出力しない） | `True` |
| `--ref-pdb FILE` | XYZ/GJF入力に使用する参照PDBまたはmmCIF topology（XYZ座標は保持） | _None_ |
| `--config FILE` | 明示 CLI 適用前に読み込むベース YAML | _None_ |
| `--show-config/--no-show-config` | 解決済み YAML レイヤー/設定を表示して続行 | `False` |
| `--out-json/--no-out-json` | `out_dir` に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照 | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--solvent TEXT` | xTB 暗黙溶媒（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画のみ表示 | `False` |

## YAML 設定

`geom`、`calc`、`freq`、`thermo` の各セクションは [YAML リファレンス](yaml-reference.md) の正規定義から変更ありません: [`geom`](yaml-reference.md#geom)、[`calc`](yaml-reference.md#calc)、[`freq`](yaml-reference.md#freq-section)、[`thermo`](yaml-reference.md#thermo) を参照してください。`freq` ではデフォルトで `calc.return_partial_hessian = true`（PHVA）が自動設定されます（YAML で上書き可）。

`freq` 固有のデフォルトとして異なるのは出力ディレクトリのみです:

```yaml
freq:
 out_dir: ./result_freq/ # freq のデフォルト
```

## 注記

- `tsopt` には虚振動数チェックが内蔵されているため、別途 `freq` を実行するのは主に熱化学量の取得や振動モードの詳細検討のためです。
- 収束した一次の鞍点（TS）では虚振動数が **ちょうど 1 つ**（検出カットオフは `hessian_dimer.neg_freq_thresh_cm`、デフォルト 5 cm⁻¹）になることが期待されます。
- 虚振動数モードは負の振動数として報告されます。`freq` は検出された虚振動数の個数を表示し、`--dump` で詳細を出力します。
- 全原子を凍結した構造にはactiveな振動DOFがないため、明示的なエラーで停止します。
- `--hessian-calc-mode` は **デフォルト < config < 明示 CLI** の優先順位で解決されます。CLI で明示的に指定した値は config YAML の `calc.hessian_calc_mode` より優先されます。
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

## 関連項目

- [tsopt](tsopt.md) — 遷移状態の最適化（内部で虚振動数チェック済み）。続けて IRC で端点の接続性を確認
- [irc](irc.md) — TS からの IRC（端点での freq と組み合わせることが多い）
- [dft](dft.md) — より高精度なエネルギー評価のための DFT 一点計算
- [all](all.md) — `--thermo` を含む一気通貫ワークフロー
- [YAML リファレンス](yaml-reference.md) — `freq` と `thermo` の設定オプション一覧
- [用語集](glossary.md) — ZPE、ギブズ自由エネルギー、エンタルピー、エントロピーの定義
- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — 典型的な失敗モードの詳細な修正方法
