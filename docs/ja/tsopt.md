# `tsopt`

`pdb2reaction tsopt` は、遷移状態（TS）*候補*を一次鞍点に最適化します。虚振動数チェックを内蔵しています。候補には `path-opt` / `path-search` の最高エネルギー像（HEI: highest-energy image）、または自前の構造を使えます。

optimizer は `--opt-mode` で選びます。ほとんどの系では `--opt-mode hess`（デフォルトの **RS-P-RFO**: Restricted-Step Partitioned Rational Function Optimization、Banerjee）を使ってください。完全 Hessian を用いるため、一般的により堅牢です。RS-P-RFO で収束しない場合や完全 Hessian の再計算コストが過大な場合は、`--opt-mode grad`（**Hessian-Guided Dimer**）に切り替えます。候補に複数の虚振動数があり余分なモードの除去が必要な場合は、`--flatten`（デフォルト無効）を有効化します。

`tsopt` は、YAML 上書き後も RFO 系および Dimer optimizer の
`reject_uphill` を常に `false` に固定します。鞍点探索では反応モードに
沿った物理エネルギー上昇を許す必要があるためです。
`--reject-uphill/--no-reject-uphill` は最小値最適化（`opt` と `all` の
IRC 後エンドポイント再最適化）だけに適用されます。

収束後、`tsopt` は最終的な Hessian 計算と虚振動数チェックを自動で行います。妥当な TS では虚振動数が **ちょうど 1 つ** です。別途の [`freq`](freq.md) は、完全な振動解析や熱化学補正が必要な場合にのみ実行します。端点の接続性は必ず [`irc`](irc.md) で確認してください。

TS 初期構造がまず必要な場合は、2 端点なら [path-opt](path-opt.md)、2 構造以上なら [path-search](path-search.md) を実行し、得られた HEI を `tsopt` → `irc` の順で最適化・検証してください。mmCIF入力は内部PDBへ変換され、成果物には元IDを復元したCIFも生成されます。XYZ/GJF入力では`--ref-pdb`にPDBまたはmmCIF topologyを指定できます。

`--ref-mode` は通常の単独 `tsopt` に必要なoptionではなく、主に `all` 内部の MEP→TS handoffです。`all` は MEP energyを使った標準のupwinding Cartesian 3N接線を生成してデフォルトで渡し、energyを読めない旧trajectoryだけは正規化secantの二等分線へfallbackします。`all --no-tsopt-from-mep-tan` では TSOPT が初期構造の Hessian を計算し、その振動モードから初期 root を選びます。外部経路から同じ原子順の非ゼロCartesian 3N方向を用意したexpert standalone runだけが`--ref-mode PATH`を直接指定します。

接線は初期Hessian rootを選び、modeが回転した後もoverlapで追跡するために使います。失敗した探索を別の探索へ自動変換する機能ではありません。デフォルトでは一時的なmode-lossによるtrial棄却、quasi-Newton固有値構造gate、自動saddle recovery、自動変位multistartを実行しません。終端のexact PHVAが合否を決め、`n_imag = 0`または`n_imag > 1`は`not_converged`です。

`--flatten`は余剰虚振動を除くための独立した明示optionです。余分な負方向は除去できますが、欠けた反応modeは生成できません。

> **命名規則の注意:** CLI は `grad|dimer`（= Dimer）、`hess|rsprfo`（= RS-P-RFO、デフォルト）、および `rsirfo`（= RS-I-RFO）/ `trim`（= TRIM）を受け付けます。YAML ではトップレベルの `hessian_dimer:`（Dimer）ブロック、または `rsirfo:` ブロック（RS-P-RFO・RS-I-RFO・TRIM が共用）を直接指定してください。

## TS 候補を得る 2 つの経路

`tsopt` は既に手元にある候補を精密化します。その候補を*構築する*には補完的な 2 つの方法があり、手元の情報に合わせて選びます。

| 経路 | サブコマンド | 使う場面 | 動作 |
| --- | --- | --- | --- |
| (a) MEP / 経路探索 | [`path-search`](path-search.md) | 両端点（反応物**および**生成物）があり、TS を自動でブラケットしたい | 再帰的な最小エネルギー経路探索（GSM / DMF）と結合変化検出。多段階経路を自動分割し、各反応区間を精密化し、区間ごとの最高エネルギー像（`hei_seg_NN.xyz`）を返す |
| (b) 距離拘束スキャン | [`scan`](scan.md) | 反応物のみがある、または特定の反応距離を直接駆動したい | 調和距離拘束 `E = ½k(r − target)²` で各反応距離を完全緩和しながら駆動し、系を TS 候補まで押し上げる |

`opt --restraint` フラグはありません。`opt` は `--dist-freeze`（調和拘束、強さは `--bias-k`）で距離を拘束しますが駆動はせず、距離を駆動する積み上げ経路は `scan`（`--preopt` / `--endopt` で駆動経路まわりの端点を緩和できる）です。いずれの経路で得た候補も `tsopt → freq → irc` に渡して最適化・検証します。

## 実行例

PDB 候補のデフォルト RS-P-RFO 最適化:

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --out-dir ./result_tsopt
```

dimer モード + 解析的 Hessian（VRAM に余裕がある場合）:

```bash
# VRAM に余裕がある場合に dimer モード + 解析的Hessianで実行する
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
 --opt-mode grad --hessian-calc-mode Analytical --out-dir ./result_tsopt_grad
```

RS-P-RFO モードを YAML 上書きと併用:

```bash
# RS-P-RFO モードを YAML 上書きと併用する
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
 --opt-mode hess --config tsopt.yaml --out-dir ./result_tsopt_hess
```

RS-P-RFO モードで flatten を有効化:

```bash
# RS-P-RFO モードで flatten を有効化して実行する
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
 --opt-mode hess --flatten --out-dir ./result_tsopt_flatten
```

最適化軌跡を保存して確認したい場合は `--dump` を追加します。

## 処理の流れ

- **電荷/スピン解決**: 電荷は標準の優先順位チェーンで解決されます。詳細は {ref}`CLI 規約: 電荷の指定 <ja-charge-specification>` を参照してください。
- **構造ロードと freeze-links**: 構造は `pysisyphus.helpers.geom_loader` で読み込まれます。`--freeze-links` が有効な場合、キャップ水素の親原子は自動的に凍結されます（{ref}`キャップ水素と凍結原子 <ja-link-hydrogen-and-frozen-atoms>` を参照）。
- **MLIP Hessian（デフォルト: UMA）**: `--hessian-calc-mode` で解析的 Hessian と有限差分 Hessian を切り替えます。いずれも活性（PHVA）部分空間を考慮します。凍結原子が存在する場合、MLIP バックエンドは活性ブロックのみを返すことがあります。Hessian 評価モードの詳細は {ref}`ja-hessian-evaluation` を参照してください。
- **Dimer モード詳細**:
 - Hessian Guided Dimer 段階は、active部分空間のexact Hessian を周期的に評価してダイマー方向を更新します。剛体モード処理は constrained に固定され、凍結anchorを動かさない全系剛体運動だけを除去します。保存・回転・試行する全方向で凍結Cartesian成分をゼロに保ち、中心外のforce評価でも凍結座標を中心imageと厳密に一致させます。`root == 0` のときは最小固有対に `torch.lobpcg` を優先し、失敗時は `torch.linalg.eigh` にフォールバックします。
 - `--flatten` が有効な場合、フラット化ループはΔx とΔg を用い、Bofill（SR1/MS ↔ PSB ブレンド; `hessian_dimer.flatten_loop_bofill` で切替）で活性 Hessian を更新します。各ループは虚振動数モード推定 → 1 回フラット化 → ダイマー方向再更新 → dimer+L-BFGS マイクロ区間 → （任意で）Bofill 更新を実行します。虚振動数モードが 1 つになったら最終的な正確な Hessian で振動解析を行います。
 - `root != 0` の場合は初期ダイマー方向のみその root を使用し、以降の更新は最も負のモード（`root = 0`）に従います。
- **RS-I-RFO モード**: RS-I-RFO を実行し、任意の Hessian 参照や R+S 分割セーフガード、マイクロサイクル制御は `rsirfo` セクションで設定します。`--flatten` が有効で収束後も虚振動数モードが複数残る場合、追加モードをフラット化して RS-I-RFO を再実行し、虚振動数モードが 1 つになるか上限に達するまで繰り返します。
- **モード出力と変換**: 絶対値が設定した閾値（デフォルト 5 cm⁻¹）未満の虚振動数モードは無視し、それ以外を `vib/imag_*_trj.xyz` に書き出します。変換が有効な場合、PDB入力はPDB companion、mmCIF／oversized-PDB入力はPDBと元IDを復元したCIFを出力します。Gaussian templateでは最終構造のみ`.gjf`を生成します。

## 出力

実行結果は `result.json`、`final_geometry.*` の最終構造、`vib/imag_*` モード（妥当な TS ではちょうど 1 つ）から検証します。

- `result_tsopt/final_geometry.pdb`（または `final_geometry.xyz`）
- `result_tsopt/vib/imag_*_trj.xyz`
- `result_tsopt/vib/imag_*.pdb`（PDB 入力の場合）

```text
out_dir/ (デフォルト:./result_tsopt/)
├─ final_geometry.xyz # 常に書き込み
├─ final_geometry.pdb # 入力がPDBの場合（変換有効時）
├─ final_geometry.cif # mmCIF/oversized-PDB入力（変換有効時）
├─ final_geometry.gjf # 入力がGaussianの場合（変換有効時）
├─ optimization_all_trj.xyz # --dumpがTrueのときの dimer モードダンプ
├─ optimization_all.pdb # PDB 入力の dimer モードに対応する PDB（変換有効時、--dump）
├─ optimization_all.cif # bridge入力の元ID復元CIF
├─ optimization_trj.xyz # --dump時の RS-P-RFO/RS-I-RFO/TRIM 軌跡
├─ optimization.pdb # rsirfo モードに対応する PDB（変換有効時、--dump）
├─ optimization.cif # bridge入力の元ID復元CIF
├─ vib/
│ ├─ imag_±XXXX.Xcm-1_trj.xyz
│ ├─ imag_±XXXX.Xcm-1.pdb
│ └─ imag_±XXXX.Xcm-1.cif # bridge入力
└─.dimer_mode.dat # dimer モード方向シード
```

終了コードは CLI 規約の {ref}`ja-exit-codes` を参照。

## CLI オプション

コマンド形式:

```bash
pdb2reaction tsopt -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m 2S+1] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [--opt-mode grad|hess|dimer|rsirfo|trim|rsprfo] [--flatten/--no-flatten] \
 [--freeze-links/--no-freeze-links] [--max-cycles N] [--thresh PRESET] \
 [--hessian-calc-mode Analytical|FiniteDifference] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

`pdb2reaction tsopt --help` でコアオプション、`pdb2reaction tsopt --help-advanced` で全オプションを表示します。入力ファイルの完全な要件（水素、元素列、原子順序の整合性、電荷指定）は [CLI 規約](cli-conventions.md) を参照してください。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| **入力と電荷** | | |
| `-i, --input PATH` | 単一geometry（`.pdb` / `.cif` / `.mmcif` / `.xyz` / `.gjf`）。trajectory は使用する1 frameを `.xyz` へ抽出してから指定 | 必須 |
| `-q, --charge INT` | 総電荷。`.gjf` テンプレートまたは `--ligand-charge`（PDB/mmCIF 入力または `--ref-pdb` 付き XYZ/GJF）が提供しない限り必須。両方指定時は `-q` が優先 | テンプレート/導出が適用されない限り必須 |
| `-l, --ligand-charge TEXT` | 単一の整数（例: `-1`）でリガンド総電荷を指定するか、残基別マッピング（例: `GPP:-3,SAM:1`）で PDB/mmCIF 残基電荷から全系の電荷を導出。`-q` 省略時に使用（PDB/mmCIF 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| `--ref-pdb FILE` | XYZ/GJF入力に使用する参照PDBまたはmmCIF topology | _None_ |
| **バックエンドと計算** | | |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--workers INT` | UMA 予測器の並列度。`workers > 1` と明示的な解析 Hessian は併用できないため、`workers = 1` または有限差分を使用。{ref}`ja-workers-analytical-error` を参照 | `1` |
| `--workers-per-node INT` | ノードあたりのワーカー数。並列予測器に渡されます | `1` |
| `--hessian-calc-mode CHOICE` | MLIP Hessian モード（`Analytical` または `FiniteDifference`） | `FiniteDifference` |
| `--solvent TEXT` | xTB 暗黙溶媒（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |
| **活性領域の凍結** | | |
| `--freeze-links/--no-freeze-links` | PDB/mmCIF 入力（または `--ref-pdb` 付き XYZ/GJF）。キャップ水素の親を凍結（`geom.freeze_atoms` にマージ）。キャップ水素の詳細は [extract](extract.md) を参照 | `True` |
| `--freeze-atoms TEXT` | 凍結する原子の 1 始まりインデックスをカンマ区切りで明示的に指定（例: `'1,3,5'`）。`--freeze-links` と併用可、任意の入力形式に適用 | _None_ |
| **TS optimizer とモード** | | |
| `--opt-mode TEXT` | TS optimizer プリセット（Choice: `grad` / `hess` / `dimer` / `rsirfo` / `trim` / `rsprfo`）。`grad`/`dimer` → Hessian-Guided Dimer; `hess`/`rsprfo` → RS-P-RFO（Banerjee、デフォルト、non-microiter）; `rsirfo` → RS-I-RFO; `trim` → TRIM（Helgaker、non-microiter）。サブコマンド別の対応表（`opt` は L-BFGS/RFO、`tsopt` は Dimer/RS-P-RFO）は {ref}`ja-opt-mode-semantics` を参照 | `hess` |
| `--ref-mode PATH` | advanced/internal MEP handoff用のCartesian 3N方向（空白区切りtextまたは`.npy`）。`all`がデフォルトで指定し、通常の単独runでは省略します。外部経路を使うexpert runではroot選択とoverlap追跡に使用します | _None_ |
| `--flatten/--no-flatten` | Dimer と RS-P-RFO / RS-I-RFO / TRIM Hessian family の余剰虚振動モード flatten を有効化。`--ref-mode` は保持する負モードを特定するが、それ自体では flatten を有効化しない | `False` |
| `--coord-type TEXT` | 最適化座標系（`cart` / `redund` / `dlc` / `tric`）。`cart` がデフォルトです。`dlc` は条件付けを変えますが、どちらも一律に高速・堅牢ではないため問題のseedで比較してください。Hessian 系`tsopt`は4種類すべて、`path-opt` / `path-search`は`cart` / `dlc`のみ受け付けます | `cart` |
| `--precision [fp32\|fp64]` | MLIP バックエンド精度。バックエンド固有のキー（UMA `precision` / ORB `precision` / MACE `default_dtype`。`aimnet2`: `fp32` は no-op、`fp64` は拒否）へ振り分け。対象系で対応精度を比較してください。{ref}`再現性: GPU クラスによる精度の選択 <ja-precision-by-gpu-class>` を参照 | バックエンドデフォルト (uma `fp32`、orb・mace `fp64`) |
| **閾値とサイクル** | | |
| `--thresh TEXT` | 収束プリセットの上書き（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`、`never`） | `baker` |
| `--max-cycles INT` | `opt.max_cycles` に渡されるマクロサイクル上限 | `10000` |
| **出力と設定** | | |
| `-o, --out-dir TEXT` | 出力ディレクトリ | `./result_tsopt/` |
| `--convert-files/--no-convert-files` | PDB/mmCIF/Gaussian入力用の XYZ/TRJ → PDB/CIF/GJF 出力を切り替え | `True` |
| `--dump/--no-dump` | 軌跡をダンプ | `False` |
| `--out-json/--no-out-json` | `out_dir` に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照 | `False` |
| `--config FILE` | 明示 CLI オプションより前に適用するベース YAML 設定ファイル | _None_ |
| `--show-config/--no-show-config` | 解決後の設定レイヤーを表示して実行を継続 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに入力/設定を検証し、実行計画を表示 | `False` |

(ja-flatten-precedence-caveat)=
### `--flatten` 優先順位の注意

```{note}
**`--flatten` はデフォルトで無効です（優先順位の注意）。** `defaults.py` では `flatten_max_iter: 50` が定義されていますが、CLI は YAML 適用前の初期値を `0` にします。実効値は以下のとおりです:

- CLI `--flatten` **未指定** → YAML で `hessian_dimer.flatten_max_iter` を**明示的に指定**しない限り `flatten_max_iter = 0`（余剰モード除去ループ無効）。フラット化のカウンタは Dimer・RS-I-RFO のどちらの経路でも `hessian_dimer` ブロックからのみ読み取られます。`defaults.py` の値 50 は**無視**されます。
- CLI `--flatten` 指定 → YAML / `defaults.py` の値が有効（デフォルト `flatten_max_iter = 50`）。引き続き YAML で上書きできます。
- CLI `--no-flatten` 指定 → YAML より優先して `flatten_max_iter = 0`。

TS 候補に複数の虚振動数がある場合は、`--flatten` を追加して余分なモードの除去ループを有効にしてください。

pathのHEIからTS最適化がなお失敗する場合は、原因に応じて2通りを試します。

1. 余分な虚振動が残る場合は`--flatten`を追加します。
2. `all` workflowでは`--refine-path`で再帰的`path-search`を実行し、
   TSOPT前のHEIを精密化します。

2番目は意図的にデフォルトOFFです。悪い／ノイズの多いpathを不要な素反応segmentへ
分割し、MEP・TSOPT・IRC・freqの計算量を何倍にもする可能性があります。まず未精密化
MEPを確認し、粗いHEIが原因と判断できる場合に有効化してください。
```

(ja-wrong-imaginary-mode-count)=
### 最適化後に虚振動数の本数が誤っている場合

真の一次鞍点は虚振動数を**ちょうど 1 つ**だけ持ち、そのモードは反応座標に沿って変位します。`tsopt` が代わりに偽の 2 本目の小さい虚振動数を報告したり、支配的な反応モードが無い場合は、以下のレバーを段階的に強めます。これらは補完的なので併用できます。

| レバー | フラグ | 効果 |
| --- | --- | --- |
| 精度を比較する | `--precision fp32|fp64` | 数値挙動はバックエンド・モデル・対象系に依存します。AIMNet2はfp64を受け付けず、どちらの設定も真の負曲率は除去しません |
| 内部座標 | `--coord-type dlc` | 最適化の条件付けを変えます。`cart` / `dlc`のどちらも一律に高速・堅牢ではないため、問題のseedで比較してください |
| 小さいモードのフラット化 | `--flatten` | 余分な虚振動数モードのフラット化ループを実行（`grad`: dimer ループ、`hess`: RS-P-RFO 後の処理）。`--no-flatten` は `flatten_max_iter = 0` を強制 |

モード変位と optimizer の停止理由を確認し、対応する精度・座標設定を選んで再実行します。`--flatten` は余分なモードにだけ使用します。例:

```bash
pdb2reaction tsopt -i ts_candidate.xyz -q -1 -m 1 \
    --precision fp64 --coord-type dlc --flatten -o result_tsopt
```

[よくあるエラーのレシピ → 収束・後処理で止まる](recipes-common-errors.md) も参照してください。

### 生成物側から開始したスキャンの障壁の読み方

この TS 候補を生成した `scan`（または経路）が**生成物**から始まった場合、報告される生のバリアは**逆方向**のバリア `E(TS) − E(product)` です。通常ほしい順方向のバリアは反応物から計算します。

| 実行内容 | 順方向バリア |
| --- | --- |
| 生成物始点のスキャン | `E(TS) − E(reactant)` — 生の生成物始点の値では**ない** |

これはフラグではなく読み取り時の解釈です。特に結晶構造の生成物複合体から開始した場合は、バリアを引用する前にスキャンがどちらの端点から始まったかを必ず確認してください。{ref}`scan: スキャン方向とバリアの符号 <ja-scan-direction-barrier-sign>` も参照してください。

### 条件を揃えた変異体 vs WT 比較

MEP の入力契約と系をまたぐ比較を混同しないでください。各 R→IM→P
経路の内部では、すべての構造が同じ原子を同じ順序で持つ必要があります。
一方、実際の WT→変異体置換では残基種や原子数が変わり得るため、WT と
変異体の全エネルギーを直接差し引いてはいけません。各系の内部で求めた
活性化エネルギーまたは自由エネルギーを比較します。

`ΔΔG‡ = (G_TS − G_R)_mutant − (G_TS − G_R)_WT`

- 意図した変異以外については、選択する残基**位置**とクラスター境界・cap
  の方針を化学的に妥当な範囲で揃えます。半径による独立抽出では境界残基が
  一方だけに入ることがあるため、選択結果を監査してください。
- protonation、電荷決定方法、backend/model、precision、拘束、熱化学条件を
  揃えます。変異が formal charge や protonation を変える場合、検証済み全電荷
  は正当に異なり得るので、見かけを対称にするため同じ `-q` を強制しません。
- **同じ組成**の機構どうしを比較する場合は、共通の原子集合と順序を使います。

```bash
# 各経路内では原子を一致させ、2つのclusterでは境界条件を揃える
pdb2reaction all -i wt_cluster.pdb     -l 'GPP:-3,SAM:1' --tsopt --thermo -o result_wt
pdb2reaction all -i mutant_cluster.pdb -l 'GPP:-3,SAM:1' --tsopt --thermo -o result_mutant
```

## YAML 設定

共通セクションについては [YAML リファレンス](yaml-reference.md) を参照してください。必要な値だけ変更してください。

### 共通設定（両モード共通）

`geom` と `calc` のキーは正規定義から変更ありません。詳細は YAML リファレンスの [`geom`](yaml-reference.md#geom) と [`calc`](yaml-reference.md#calc) を参照してください。

`opt` ブロックは [`opt`](yaml-reference.md#opt) と同じキーを使用し、以下が `tsopt` 固有のデフォルトです:

```yaml
opt:
 thresh: baker # tsopt のデフォルト（`opt` は `gau`）
 out_dir: ./result_tsopt/ # tsopt のデフォルト（`opt` は `./result_opt/`）
```

```{note}
**energy plateau stop（opt-in、デフォルト無効）。** Hessian-family TS optimizer（RS-P-RFO、
RS-I-RFO、TRIM）は共通の `energy_plateau` 設定を参照し、`--stop-plateau` で有効化します。
有効時、直近50 stepの energy rangeが `--stop-plateau-thresh`（default `1×10⁻⁴ au`）を下回ると、
`stalled` として停止し、終端 PHVA を実行します。未収束のまま `max_cycles` に到達した場合は
PHVA を実行しません。backend/model/system依存のforce floorが選択閾値への
到達を妨げる場合に無駄なcycleを避けられます。デフォルトで無効なのは、平坦なenergyで停止した
TS探索が余分な虚振動を残したままになりやすいためです。
```

### Dimer モード（`--opt-mode grad`）

`--opt-mode grad`（Hessian Guided Dimer + L-BFGS）で使用します。

`hessian_dimer` ブロック全体（内部 `dimer:` とそのネストされた `lbfgs:` を含む）は [`hessian_dimer`](yaml-reference.md#hessian_dimer) に記載されています。内部 `lbfgs:` は [`lbfgs`](yaml-reference.md#lbfgs) セクションを継承し、以下が `tsopt` 固有の上書きです:

```yaml
hessian_dimer:
 dimer:
   lbfgs:
     out_dir: ./result_tsopt/ # tsopt の上書き（defaults.py の値は ./result_opt/）
```

### RS-P-RFO / RS-I-RFO モード（`--opt-mode hess`、デフォルト → RS-P-RFO）

`--opt-mode hess`（RS-P-RFO、デフォルト）で使用します（`rsirfo` は RS-I-RFO、`trim` は TRIM を選択。3 つともこのブロックを共用します）。

`rsirfo` ブロック全体は [`rsirfo`](yaml-reference.md#rsirfo) に記載されています（trust-region と Hessian-update の各キーは [`rfo`](yaml-reference.md#rfo) からも継承）。`tsopt` 固有の上書きは以下のとおりです:

```yaml
rsirfo:
 trust_max: 0.10 # 最大信頼半径 (bohr)
 out_dir: ./result_tsopt/ # tsopt の上書き（defaults.py の値は ./result_opt/）
 hessian_recalc: 500 # N マクロステップごとに exact Hessian を再計算
 saddle_recovery_check_interval: 50 # 自動回復をYAMLで有効化した場合のexact PHVA間隔
 saddle_recovery_max_cycles: 0 # n_imag=0 自動回復はデフォルト無効
```

```{tip}
最適化中に TS モードが別のルートに切り替わる場合（例: 複数の虚振動数が存在する場合）は `rsirfo.track_mode_by_overlap: true` を設定してください。
```

```{tip}
TS 収束が遅い場合や最適化中に TS モードが失われる場合は、`rsirfo` セクションの `hessian_recalc` を小さくしてみてください（例: 50--200）。正確なHessian再計算の頻度を上げることで、追加のHessian評価コストと引き換えに堅牢性が向上します。
```

## 注記

- 絶対値が設定した閾値（デフォルト 5 cm⁻¹）未満の虚振動は、モードファイル出力と平坦化では無視します。最終的な TS 判定ではすべての負の振動数を数えます。Hessian-family optimizer は一次鞍点のrootを1個だけ追跡します。YAMLでは1要素のlist（例: `rsirfo.roots: [0]`）で設定し、空listまたは複数rootは拒否されます。Dimer は別の単数 key `hessian_dimer.root`（default `0`）を使います。`tsopt` に `--root` CLI flag はありません（[`irc`](irc.md) とは異なります）。
- `--opt-mode` はワークフロー選択用です（デフォルト: `rsprfo`）。YAML のモードマッピングを手動で変更するのではなく、目的のアルゴリズムに合ったモードを選択してください。
- Dimer方向、回転force、flatten、最終exact PHVA検証は`freq`と同じ固定の constrained 処理を使用します。Dimerは中心imageが変わるたびにこの基底を再構築します。全凍結anchorと両立する真の剛体null方向でない限り、active fragmentの並進を差し引きません。Hessian RFO最適化自体は、この射影を行わずactive-DOF Cartesian Hessian を扱います。詳細は[凍結原子](freeze-atoms.md#凍結境界での剛体モード)を参照してください。
- 設定の優先順位は {ref}`CLI 規約: 設定の優先順位 <ja-configuration-precedence>` を参照してください。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) — 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — 詳細な切り分け
- [path-search](path-search.md) — TS 候補（HEI）を特定する MEP 探索
- [irc](irc.md) — 最適化された TS からの反応経路追跡
- [freq](freq.md) — 完全な振動解析と熱化学補正（虚振動数チェックは `tsopt` が内部で実行済み）
- [all](all.md) — 抽出 → MEP → tsopt → IRC（→ オプションで freq/DFT）を連鎖する一気通貫ワークフロー
- [YAML リファレンス](yaml-reference.md) — `hessian_dimer`（Hessian Guided Dimer）と `rsirfo` の完全な設定オプション
- [用語集](glossary.md) — TS、Dimer、RS-I-RFO、Hessian の定義
