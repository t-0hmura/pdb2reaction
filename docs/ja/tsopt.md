# `tsopt`（遷移状態の構造最適化）

`tsopt` サブコマンドは、遷移状態（TS）の候補構造を 1 次の鞍点へ最適化し、final geometry で Hessian を計算して虚振動数の本数（n_imag）を数えます。TS 最適化が成功すると、反応モードの虚振動が 1 つ出ます。

## 主な用途

* **TS 候補の仕上げ**: [`path-opt`](path-opt.md) / [`path-search`](path-search.md) の最高エネルギーのイメージ（HEI）や [`scan`](scan.md) の頂点を、最適化した TS に仕上げる
* **自作の構造の検証**: 手で作った候補が TS か（n_imag = 1）を確かめ、反応モードをアニメーションで確認する
* **`all` の TS 段のやり直し**: [`all`](all.md) で得た TS を、設定を変えて単独で最適化し直す

デフォルトの計算バックエンドは、Meta が公開した学習済みの[機械学習原子間ポテンシャル（MLIP）](backends.md)の **UMA** です。`-b/--backend` で **ORB**、**MACE**、**AIMNet2**、**DFT** も選べます。

候補がまだ無い場合は、先に次のコマンドで作ってください。

| 手元にあるもの | 候補を作るコマンド |
| --- | --- |
| 反応物**と**生成物 | [`path-opt`](path-opt.md)（2 構造。`hei.xyz`）または [`path-search`](path-search.md)（2 構造以上。結合が変わるセグメントごとに `hei_seg_NN.xyz`） |
| 反応物だけ、または動かしたい結合がある | [`scan`](scan.md) で反応する距離を少しずつ動かし、ほかの自由度を緩和する |

---

## 基本的な実行例

### 1. 標準の実行（RS-P-RFO）

電荷とスピン多重度を明示して実行します。

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --out-dir ./result_tsopt
```

### 2. Dimer 法

完全な Hessian を繰り返し計算するのが重い場合や、難しい候補で別の方法を試したい場合に使います。

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode dimer --out-dir ./result_tsopt_dimer
```

### 3. 余分な虚振動の除去

候補に虚振動が 2 つ以上ある場合は `--flatten` を付けます。

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --flatten --out-dir ./result_tsopt_flatten
```

### 4. 保存した Hessian から開始

同じ構造で `freq` などの `--dump-hess` で保存した Hessian を読み込み、計算し直さずに始めます。

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --read-hess ts_cand_hess.npy --out-dir ./result_tsopt
```

---

## 処理の仕組みと計算仕様

1. **構造の読み込みと境界の凍結**: {ref}`電荷 <ja-charge-specification>`は `-q` または `-l` から決まります。`--freeze-links`（デフォルト有効）では、切り出したクラスターの{ref}`キャップ水素 <ja-link-hydrogen-and-frozen-atoms>`の親原子を凍結し、可動原子だけで Hessian を扱います（PHVA: 部分 Hessian 振動解析）。PDB のモデルから得た `.xyz` の候補では、その PDB を `--ref-pdb` で渡すと、`--freeze-links` で{ref}`境界を凍結 <ja-freeze-atoms-and-restraints>`できます。
2. **最適化法の選択**（`--opt-mode`）: `hess`（デフォルト）は完全な Hessian を使う **RS-P-RFO**（制限ステップ分割有理関数最適化）を実行し、`rsirfo` と `trim` はそれぞれ RS-I-RFO（restricted-step image RFO）と TRIM（trust-region image minimization）を選びます。`dimer`（または `grad`）は **Hessian-guided Dimer** 法で、勾配を使って最低固有モードを追い、ときどき厳密な Hessian で方向を更新します。
3. **反応モードに沿った探索**: 反応モードの方向にはエネルギーを上り、それ以外の方向には下りながら、収束条件（`--thresh`）を満たすまで構造を動かします。RS-P-RFO は Bofill 式で Hessian を更新し、1 ステップを信頼半径 0.1 bohr 以内に収めます。デフォルトの `baker` は、力の最大値 3 × 10⁻⁴ 未満、力の RMS 2 × 10⁻⁴ 未満、ステップの最大値 3 × 10⁻⁴ 未満、ステップの RMS 2 × 10⁻⁴ 未満（原子単位）、エネルギー変化 10⁻⁶ hartree 未満の 5 つをすべて同時に求めます。どれも Gaussian の既定（`gau`）より厳しい条件です。
4. **最後の確認**: 収束すると、final geometry で Hessian を計算して n_imag を数え、各虚振動モードをアニメーションとして書き出します。ν < −5.00 cm⁻¹ のモードを虚振動として数え、−5.00 以上 0 cm⁻¹ 未満の値は数値誤差として扱います。この閾値は YAML の `freq.zero_cutoff_cm` で変えられます。凍結原子の扱いは [`freq`](freq.md#凍結境界での剛体モード) と同じです。
5. **余分な虚振動の除去（`--flatten`）**: 虚振動が 2 つ以上残る場合は、余分なモードに沿って構造をずらして最適化し直し、1 つになるか回数の上限に達するまで繰り返します。Dimer 法では、各回でダイマー方向も更新し、短い Dimer + L-BFGS の区間を実行します。

---

## TS の判定

結果は、実行がどう終わったかで決まります。

| 終わり方 | 端末の `[tsopt]` の判定の行 | 終了コード | `tsopt` が残すもの | `all` の次の動作 |
| --- | --- | --- | --- | --- |
| 収束 | `[tsopt] Converged (n_imag=1).`。n_imag ≥ 2 では `[tsopt] WARNING: Higher-order stationary point (n_imag=N, …)`、n_imag = 0 では `[tsopt] No imaginary mode detected. …` | 0 | final geometry、n_imag、虚振動モード | n_imag ≥ 1 なら IRC へ進み、n_imag = 0 なら IRC の前で止まる |
| エネルギーが変わらなくなって停止 | `[tsopt] ERROR: Not converged (plateau stop, n_imag=N).` | 1 | final geometry と n_imag | IRC の前で止まる |
| `--max-cycles` に達して未収束 | `[tsopt] ERROR: Not converged.` | 1 | final geometry（Hessian なし） | IRC の前で止まる |
| `--skip-final-freq` を付けて収束 | `[tsopt] Converged; terminal PHVA is unavailable.` | 0 | final geometry（Hessian なし） | 反応モードを確かめられないため、IRC の前で止まる |
| 最後の Hessian の計算に失敗 | `[tsopt] Converged; terminal PHVA is unavailable.` | 1 | final geometry と、`hessian_status: failed` とその理由 | IRC の前で止まる |

n_imag は次のように読みます。

| n_imag | 意味 |
| --- | --- |
| 1 | 1 次の鞍点。モードが狙った原子を動かしているかを確かめてから、[`irc`](irc.md) を実行してください |
| 0 | 虚振動なし。構造が極小点の側へ緩和しています |
| 2 以上 | 高次の鞍点。余分な虚振動が残っています。`all` はそれでも、MEP の方向にいちばん近い虚振動のモードに沿って IRC を流すので、IRC の端点でそのモードがどこへつながるかを確かめられます |

n_imag は端末の `[tsopt]` の判定の行か、`--out-json` を付けたときの `result.json` の `n_imaginary_modes` で読みます。

(ja-wrong-imaginary-mode-count)=
### 最適化後に虚振動数の本数が誤っている場合

n_imag が 1 でない場合や、モードが狙った反応の原子を動かしていない場合は、次を試してください。これらは組み合わせて使えます。

| 結果 | 試すこと |
| --- | --- |
| n_imag = 0 | 候補が鞍点から遠い状態です。経路探索や scan でよりよい候補を作ってください。`all` では `--refine-path` で再帰的な `path-search` を実行し、HEI を細かく求め直せます。増えた素過程のそれぞれに TS 最適化と IRC がかかるため、計算量は増えます |
| n_imag ≥ 2 | 各モードの動きを確かめてください。`--flatten` を付けて最適化し直すか、この候補で `--precision fp32` / `fp64` や `--coord-type cart` / `dlc` を比べてください |
| モードは 1 つだが動きが違う | どの原子が動くかを確かめ、狙った反応に近い候補から始めてください |

例として、fp64 と DLC 座標で、flatten を有効にしてやり直す場合は次のようにします。

```bash
pdb2reaction tsopt -i ts_candidate.pdb -q -1 -m 1 \
    --precision fp64 --coord-type dlc --flatten -o result_tsopt
```

ほかの手は {ref}`TS が取れないとき <ja-ts-search-fails>` に、そのほかの失敗は [トラブルシューティング](troubleshooting.md) にあります。

---

## 主な出力ファイル

実行が終わると、`--out-dir` に次のファイルができます。

```text
result_tsopt/
├─ final_geometry.xyz             # final geometry（常に出力）
├─ final_geometry.pdb             # 同じ構造の PDB（PDB/mmCIF 入力。Gaussian 入力では .gjf）
├─ vib/
│  ├─ imag_-385.20cm-1_trj.xyz    # 虚振動モードごとのアニメーション
│  └─ imag_-385.20cm-1.pdb        # 同じアニメーションの PDB（PDB/mmCIF 入力）
├─ optimization_trj.xyz           # 最適化の軌跡（--dump。Dimer 法では optimization_all_trj.xyz）
└─ result.json                    # 結果の要約（--out-json）
```

{ref}`mmCIF の入力 <ja-mmcif-input>`と、PDB の欄に入りきらない大きな PDB の入力では、元の識別子を保った `.cif` も書きます。

* **final geometry**: `final_geometry.*` を、[`irc`](irc.md) に渡す TS として使います。
* **反応モード**: `vib/imag_*_trj.xyz` を PyMOL や VMD で開き、生成・切断される結合に沿って原子が動いているかを確かめてください。
* **要約**: `--out-json` を付けると、[`result.json`](json-output.md) に終わり方 `optimization_status` と `hessian_status` が記録されます。

---

## 主な CLI オプション

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | 入力構造ファイル（`.pdb`, `.cif`, `.mmcif`, `.xyz`, `.gjf`）。軌跡は 1 フレームを `.xyz` に切り出してから指定（{ref}`軌跡から 1 フレームを取り出す <ja-trajectory-one-frame>` を参照） |
| `-q, --charge` | 整数 | `None` | 系全体の総電荷。`-l` を使う場合と `.gjf` 入力のほかは必須 |
| `-m, --multiplicity` | 整数 | `1` | スピン多重度（2S+1） |
| `-l, --ligand-charge` | 文字列 | `None` | リガンドの総電荷（例: `-1`）または残基名ごとの電荷（例: `'GPP:-3,SAM:1'`）。`-q` を省いたときに使用（PDB/mmCIF 入力または `--ref-pdb`） |
| `--ref-pdb` | パス | `None` | `.xyz`・`.gjf` 入力に対応づける PDB/mmCIF のトポロジー（座標は `-i` のものを使用） |
| `-o, --out-dir` | パス | `./result_tsopt/` | 出力先ディレクトリ |
| `-b, --backend` | 文字列 | `uma` | バックエンド（`uma`, `orb`, `mace`, `aimnet2`, `dft`） |
| `--opt-mode` | `hess` / `dimer` / `rsirfo` / `trim` | `hess` | 最適化法: RS-P-RFO / Dimer / RS-I-RFO / TRIM（`rsprfo` = `hess`、`grad` = `dimer`）。`opt` では `grad` は L-BFGS を指す（{ref}`コマンドごとの --opt-mode <ja-opt-mode-semantics>` を参照） |
| `--ref-mode` | パス | `None` | 反応モードの参照方向（`.npz`, `.npy`, テキスト）。`all` が MEP から渡すもので、通常は指定しない。Dimer 法では使わない |
| `--hessian-calc-mode` | `FiniteDifference` / `Analytical` | `FiniteDifference` | Hessian の計算法（有限差分 / 解析的） |
| `--flatten/--no-flatten` | フラグ | `False` | 余分な虚振動を除く |
| `--freeze-links/--no-freeze-links` | フラグ | `True` | キャップ水素の親原子を凍結（PDB/mmCIF 入力または `--ref-pdb`） |
| `--freeze-atoms` | 文字列 | `None` | 凍結する原子（1 始まり、カンマ区切り: 例 `'1,3,5'`） |
| `--thresh` | プリセット | `baker` | 収束条件（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） |
| `--max-cycles` | 整数 | `100000` | 最適化サイクルの上限 |
| `--stop-plateau/--no-stop-plateau` | フラグ | `False` | エネルギーが変わらなくなったら（直近 50 サイクルの幅が 1e-4 hartree 未満）止め、Hessian を計算 |
| `--skip-final-freq/--no-skip-final-freq` | フラグ | `False` | 収束後の最後の Hessian を省く |
| `--read-hess` | パス | `None` | Hessian を計算せず `.npy` ファイルから読んで開始（Cartesian、Hartree/bohr²、全原子または可動原子だけ） |
| `--dump-hess` | パス | `None` | final geometry の Hessian を `.npy` ファイルに保存（`freq`・`tsopt`・`irc` の `--read-hess` 用）。最後の Hessian を計算したときだけ書く |
| `--precision` | `fp32` / `fp64` | バックエンドごと（`uma`: `fp32`、`orb`・`mace`: `fp64`） | バックエンドの精度。`aimnet2` は `fp64` を受け付けない（{ref}`MLIP バックエンド: 精度 <ja-precision-by-gpu-class>` を参照） |
| `--coord-type` | `cart` / `redund` / `dlc` / `tric` | `cart` | 最適化に使う座標系：デカルト座標 / 冗長内部座標 / 非局在化内部座標（DLC）/ 並進・回転を含む内部座標（TRIC） |
| `--config` | パス | `None` | コマンドラインのオプションより前に適用する YAML ファイル |
| `--dump` | フラグ | `False` | 最適化の軌跡を書き出す |
| `--out-json/--no-out-json` | フラグ | `False` | 結果の要約を `result.json` に出力（[JSON 出力リファレンス](json-output.md)） |

全オプションは `pdb2reaction tsopt --help-advanced` または [自動生成 CLI リファレンス](../reference/commands/tsopt.md) を参照してください。

> **補足:** YAML では、Dimer 法は `hessian_dimer:` ブロックを読み、RS-P-RFO・RS-I-RFO・TRIM は `rsirfo:` ブロックを共用します。キーの一覧は YAML リファレンスの [`rsirfo`](yaml-reference.md#rsirfo) と [`hessian_dimer`](yaml-reference.md#hessian_dimer) にあります。

> **補足:** 最適化の途中で反応モードが別の Hessian 固有ベクトル（root）に入れ替わる場合は、`rsirfo.track_mode_by_overlap: true` を設定してください。

> **補足:** 収束が遅い場合は、`rsirfo.hessian_recalc`（デフォルト `500`）を 50〜200 に下げてください。厳密な Hessian を計算し直す間隔が短くなり、計算は増えますが収束しやすくなります。

---

## 使用上の注意点

(ja-flatten-precedence-caveat)=
### `--flatten` を使うとき

flatten の回数は、Dimer 法でも RS-P-RFO・RS-I-RFO・TRIM でも、YAML の 1 つのキー `hessian_dimer.flatten_max_iter` で決まります。

| コマンドライン | flatten の回数 |
| --- | --- |
| `--flatten` も `--no-flatten` も付けない | `0`（無効）。YAML で `hessian_dimer.flatten_max_iter` を指定した場合はその値 |
| `--flatten` | YAML の値（正の値の場合）、無ければ `50` |
| `--no-flatten` | YAML に値があっても `0` |

`--flatten` は欠けている反応モードを作れません。n_imag = 0 の場合は、よりよい候補を作ってください。

### そのほかの注意

* **上り方向のステップは常に許可**: 鞍点探索では反応モードに沿ってエネルギーを上る必要があるため、YAML で指定しても `tsopt` は `reject_uphill: false` を保ちます。`--reject-uphill/--no-reject-uphill` は、`opt` と `all` の端点の最適化で使うフラグです。
* **生成物側から scan した障壁**: この候補を作った scan が生成物から始まった場合、障壁の読み方は {ref}`scan: スキャン方向とバリアの符号 <ja-scan-direction-barrier-sign>` を参照してください。
* **追う固有ベクトル（root）は 1 つ**: 最適化は 1 つの固有ベクトルに沿って上ります（`0` が最小の固有値）。`rsirfo.roots: [0]` のように 1 要素のリストで指定します。Dimer 法では `hessian_dimer.root` を使います。`tsopt` に `--root` フラグはありません。
* **そのほかの RS-P-RFO の設定**: `trust_norm: max_atom` はステップ全体ではなく原子ごとの変位を制限し（Cartesian 座標だけ）、`hessian_update: ts_bfgs` は Bofill の代わりに TS-BFGS で Hessian を更新します。どちらも信頼半径は変えません。
* **追加の探索は指定したときだけ**: 収束後は、n_imag が 1 でなくても、自動では追加の探索をしません。`--flatten` を使うか、`rsirfo.saddle_recovery_max_cycles` を `0` より大きくしてください（デフォルト `0`）。後者では、厳密な Hessian に虚振動が無いとき、RS-P-RFO・RS-I-RFO・TRIM がエネルギーを上る向きにステップを進めます。
* **併用できない組み合わせ**: `--skip-final-freq` と `--dump-hess`、2 以上の `--uma-workers` と `--hessian-calc-mode Analytical`。
* **`--skip-final-freq` と `--flatten`**: RS-P-RFO・RS-I-RFO・TRIM では、`--skip-final-freq` を付けると、最後の Hessian を使う `--flatten` も省かれます。
* **`--read-hess` を RS-P-RFO・RS-I-RFO・TRIM で使う場合**: ファイルの Hessian が最初の厳密な Hessian の代わりになるので、`rsirfo.hessian_init` はデフォルトの `calc` のままにしてください。ほかの値ではエラーで止まります。
* **Dimer の方向**: Dimer 法は今の方向を出力先の `.dimer_mode.dat` に書きます。
* **`--ref-mode` と凍結原子**: `--ref-mode` は MEP から反応の方向を与えるだけで、凍結境界の扱いは変えません。

---

## 関連ドキュメント

* [irc](irc.md) — 最適化した TS からの反応経路の追跡
* [freq](freq.md) — 完全な振動解析と熱化学補正
* [path-opt](path-opt.md) / [path-search](path-search.md) / [scan](scan.md) — TS 候補の作成
* [all](all.md) — 抽出・MEP・TS 最適化・IRC・振動解析を一度に実行するワークフロー
* [反応機構を調べるコツ](mechanism-tips.md) — TS が取れないときに試すこと
* [トラブルシューティング](troubleshooting.md) — 実行が失敗したときの切り分け
* [YAML リファレンス](yaml-reference.md) — `rsirfo` と `hessian_dimer` のすべての設定
* [用語集](glossary.md) — TS、Dimer、Hessian などの用語
* {ref}`終了コード <ja-exit-codes>` — 終了ステータスの意味
