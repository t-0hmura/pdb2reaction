# トラブルシューティング

症状ベースで当たりを付けたい場合は、先に [典型エラー別レシピ](recipes-common-errors.md) を参照してください。

## 実行前チェックリスト

長い計算を回す前に、最低限次を確認してください。

- UMA のモデルがダウンロードできる（Hugging Face のログイン/トークンが利用可能）
- 酵素系ワークフローでは、入力 PDB に **水素** と **元素記号（element column）** が入っている
- 複数の PDB を与える場合、**同じ原子が同じ順序** で並んでいる（座標だけが異なる）

---

(ts-input-extraction)=
## 入力 / 抽出の問題

### 「Element symbols are missing … add-elem-info を実行してください」
典型的なメッセージ:

```text
Element symbols are missing in '...'.
Please run `pdb2reaction add-elem-info -i...` to populate element columns before running extract.
```

対処:
- 次を実行して element 列（元素記号列）を補完します。

  ```bash
  pdb2reaction add-elem-info -i input.pdb -o input_with_elem.pdb
  ```

- その後、`extract` / `all` を補完後の PDB で再実行します。

原因:
- PDB の element 列が空だったり不統一だったりすることが多く、`extract` は原子タイプ判定のために元素記号を必要とします。

---

### 「[multi] Atom count mismatch …」「[multi] Atom order mismatch …」
典型的なメッセージ:

```text
[multi] Atom count mismatch between input #1 and input #2:...
[multi] Atom order mismatch between input #1 and input #2.
```

対処:
- **すべて** の構造を同じ前処理ワークフロー（同じプロトン化ツール、同じ設定）で作り直します。
- 水素付加を行う場合、全フレームで同一順序になる手順を選びます。

ヒント:
- MD 由来なら、同一のトポロジー/軌跡からフレーム抽出する方が安全です（異なるツールで生成した PDB を混ぜると順序がズレやすい）。

---

### 「活性部位モデル（バインディングポケット）が空っぽ / 必要な残基が落ちる」
症状:
- 抽出された活性部位モデルが想定より小さい
- 触媒残基が含まれない

対処の例:
- `--radius` を増やしてください（例: 2.6 → 3.5 Å）
- `--selected-resn` で残基を強制包含してください（例: `--selected-resn 'A:123,B:456'`）。残基名・ID・chain付きselectorの指定形式は CLI 規約の {ref}`ja-selected-resn-takes-ids` を参照。
- 以前に `--exclude-backbone` を明示的に有効化していて、その削除が強すぎる場合は当該フラグを外す（または `--no-exclude-backbone` を渡す）ことで主鎖原子を保持できます

---

### エネルギー・障壁の計算値が不正確

症状:
- 計算されたエネルギーや反応障壁が不合理に見える
- モデルサイズを大きくすると結果が大幅に変わる

対処:
- 抽出された活性部位モデルが小さすぎると、エネルギーや障壁の計算値が不正確になることがあります。抽出半径を大きくする（例: `-r 4.0` 以上）ことで、タンパク質環境をより多く含めて精度を改善できます:

  ```bash
  pdb2reaction extract -i complex.pdb -c 'SUB' -o model.pdb -r 4.0
  ```

---

### 非標準残基が正しく切断されない

抽出された活性部位モデルに非標準の 3 文字コードを持つ修飾アミノ酸残基（リン酸化セリン、メチル化リシンなど）が含まれている場合、デフォルトでは主鎖切断やキャップ水素付加が適用されません。`--modified-residue` で登録してください:

```bash
pdb2reaction extract -i complex.pdb -c PRE --modified-residue "SEP,TPO,MLY" -o pocket.pdb
```

`--modified-residue` で対応できない場合（残基の主鎖トポロジーが特殊な場合など）は、活性部位モデルを手動で構築し、下流のコマンド（`opt`、`tsopt`、`path-opt` など）に直接渡してください。

---

(ts-charge-spin)=
## 電荷 / スピンの問題

### 「電荷が必須」系のエラー（非 GJF 入力）
`.gjf` でない入力では、複数ステージで総電荷が必要になります。`-q/--charge` を省略した場合、workflow は `--ligand-charge/-l`（PDB/mmCIF、または `--ref-pdb` topology 付き XYZ/GJF）、YAML `calc.charge`、または `.gjf` template から電荷を解決します。いずれの経路でも解決できないと上記のエラーで停止します。

対処:
- 電荷と多重度を明示する:

  ```bash
  pdb2reaction path-search -i R.pdb P.pdb -q 0 -m 1
  ```

- あるいは（抽出ありの場合）残基名ごとの電荷マッピングを与える:

  ```bash
  pdb2reaction extract -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
  ```

---

(ts-install-env)=
## インストール / 環境の問題

### UMA のダウンロード/認証エラー（Hugging Face）
症状:
- モデルをダウンロードできない、認証が必要、といったエラー。

対処:
- 環境/マシンごとに一度ログインします。

  ```bash
  hf auth login
  ```

- HPC では、計算ノードから HF キャッシュ（ホームディレクトリ等）に書き込み可能か確認してください。

---

### CUDA / PyTorch の不整合
症状:
- GPU があるのに `torch.cuda.is_available()` が False
- import 時に CUDA runtime error が出る

対処:
- クラスターの CUDA と整合する PyTorch を入れます。
- GPU が見えているか確認します:

  ```bash
  nvidia-smi
  python -c "import torch; print(torch.version.cuda, torch.cuda.is_available())"
  ```

---

### ORB backend の import に失敗する

対処:
- 現行 extra を `pip install "pdb2reaction[orb]"` で導入し、`python -m pip check` を実行します。
- 現行 `orb-models` は `torch_scatter` を依存に持ちません。実際の error が、別途導入した `torch_scatter` 自体を名指ししていない限り、PyG index や source-build 手順を追加しないでください。

---

### DMF モードが動かない（`cyipopt` がない、または `No module named 'dmf'`）
DMF（`--mep-mode dmf`）を使うときに IPOPT/`cyipopt` の import エラーが出る場合:

対処:
- `pdb2reaction` を入れる前に conda-forge から `cyipopt` を入れるのが簡単です。

  ```bash
  conda install -c conda-forge cyipopt
  ```
- `pydmf` は `pdb2reaction` の依存として同梱されています。`No module named 'dmf'` が出る場合は `pip install --force-reinstall pdb2reaction` で再インストールしてください。

---

### 図のエクスポートが失敗する（Chrome がない）
Plotly/Chrome 系のエラーで静的画像が出ない場合:

対処:
- headless Chrome を一度入れます。

  ```bash
  plotly_get_chrome -y
  ```

---

(ts-calc-conv)=
## 計算 / 収束の問題

### 最適化が `max_cycles` に達し、`max(force)` が閾値をわずかに超える

症状:
- optimizerが `max_cycles` 上限まで回り、最終summaryで `max(force)` や `rms(force)` が選択閾値をわずかに上回る。
- 一方で、energyは選択した `energy_plateau_thresh` の範囲内で平坦化している。

原因:
- MLIPのforce noise/flatnessはbackend、model、precision、system、hardware stackに依存し、選択したforce閾値への到達を妨げる場合があります。

対処:
- 無限反復は `--max-cycles` が常に防ぎます。エネルギーが平坦化した時点で早く止めたい場合は、**エネルギー平坦化停止**を opt-in してください。`--stop-plateau`（YAML `opt.energy_plateau: true`）を指定すると、設定した力・ステップの収束条件を満たさないまま、直近 `--stop-plateau-window`（デフォルト 50）ステップのエネルギーレンジが `--stop-plateau-thresh`（デフォルト `1×10⁻⁴ au ≈ 0.06 kcal/mol`）を下回った時点で、計算は `stalled`（**未収束**）として停止します。平坦化だけでは停留点を受理できません。最終構造、力、結果statusを確認してから再実行してください。
- 再実行では、必要に応じて力の閾値を `--thresh gau`（`opt` のデフォルト）または `--thresh gau_loose` に変更してください。
- `--stop-plateau-thresh` / `--stop-plateau-window`（YAML `opt.energy_plateau_thresh` / `opt.energy_plateau_window`）で判定を調整できます。既定どおり無効のままにすれば、`--max-cycles` だけが上限になります。
- 注意: この平坦地形停止は **chain-of-states オプティマイザ**（`path-opt`、`path-search` の string/GSM（Growing String Method）/DMF（Direct Max Flux）段階）では**スキップ**されます（単一のスカラーエネルギー履歴ではなく、イメージごとのエネルギー配列を保持しているため）。

---

### TS 最適化が収束しない

症状:
- TS 最適化が多くのサイクルを回しても収束しない
- 最適化後も Hessian 行列に複数の負の固有値が残る（虚振動数が 2 本以上）

対処の例（CLI フラグと YAML キーは補完的、必要に応じて併用してください）:
- オプティマイザモードを切り替えてください: `--opt-mode grad`（Dimer 法）または `--opt-mode hess`（RS-P-RFO 法、デフォルト）
- 余分な虚振動数モードのフラット化を有効にしてください: `--flatten`（単独の `tsopt`、`opt`、および `pdb2reaction all` で利用可能。デフォルトは無効）
- coarse MEP のHEIが悪い場合は、`all`を`--refine-path`付きで再実行してください。ただし悪いpathを不要な複数segmentへ分割してcostを増大させることがあるためデフォルトOFFです。まずcoarse MEP を確認してください
- 最大サイクル数を増やしてください: `--max-cycles 20000`（単独の `tsopt` の場合）、`--tsopt-max-cycles 20000`（`all` の場合）
- より厳しい収束閾値を使ってください: `--thresh baker` または `--thresh gau_tight`
- YAML でステップサイズ / 信頼半径を縮小してください — L-BFGS/Dimer: `lbfgs.max_step` / `hessian_dimer.lbfgs.max_step`、RFO/RS-I-RFO: `rfo.trust_radius` / `rfo.trust_min` / `rfo.trust_max`（および `rsirfo` セクション）。セクション構成は [YAML リファレンス](yaml-reference.md) を参照

---

### IRC が正常に終了しない

症状:
- IRC が明確な極小構造に到達する前に停止する
- エネルギーが振動したり勾配ノルムが大きいままになる

対処の例:
- 単独の `irc`: `--step-size 0.05`（デフォルト: 0.10 bohr）、必要なら `--max-cycles 200`。
- `all`: `--irc-step-size 0.05`。上限は YAML `irc.max_cycles`（`all --max-cycles` は MEP 用）。
- 物理的な停止条件を無視するには、単独で `--never-stop`、`all` で `--irc-never-stop` を指定し、軌跡と端点を確認してください。
- IRC 実行前に Cartesian PHVA の負の振動数が **ちょうど 1 本** であることを確認してください。5 cm⁻¹の閾値はモードファイル、平坦化、任意の鞍点回復に使いますが、この本数は変えません。

---

### MEP（最小エネルギー経路）探索（GSM/DMF）が失敗または予期しない結果

症状:
- 経路探索が有効な MEP なしで終了
- 結合変化が正しく検出されない

対処の例:
- `--max-nodes` をデフォルトの 20 より増やしてください（複雑な反応には 30 や 40 など）
- 端点の事前最適化を有効にしてください: `--preopt`
- 別の MEP 手法を試してください: `--mep-mode dmf`（GSM が失敗した場合）またはその逆
- YAML で結合検出パラメータを調整してください（`bond.bond_factor`、`bond.delta_fraction`）

## パフォーマンス / 安定性のヒント

- **VRAM 不足**: 必要な残基が残ることを確認してから活性部位モデルを小さくする、`--max-nodes` を減らす、full-Hessian optimizer が不要なら `--opt-mode grad` を使う
- **解析 Hessian が遅いまたは OOM**: 対象 backend/model と代表的な原子数で試験するまでは、移植性のあるデフォルト `FiniteDifference` を維持してください。解析 autograd は通常メモリピークが大きくなりますが、一般化できる原子数・VRAM の境界値はありません
- **workers > 1**: hardware と workload によっては UMA throughput が向上しますが、並列 predictor は解析 Hessian を持ちません。`Analytical` を明示すると `BackendError`（`RuntimeError` のサブクラス）で停止します。解析 Hessian には `--workers 1`、並列実行には `FiniteDifference` を指定してください
- **大規模系**: 化学的に妥当な小さい活性部位モデルを作り、半径・境界位置への感度を確認してください。multi-GPU の対応範囲は backend/workflow ごとに異なり、GPU 数を増やしても worker 当たりのメモリが減るとは限りません
- **HPC で DFT を回すとき**: 選択した計算が temporary disk を使う場合は、容量と性能を確認した filesystem に `PYSCF_TMPDIR` を設定してください。`/tmp`、`$PBS_O_WORKDIR`、shared filesystem のどれが適切かは site ごとに異なります

## バックエンド選択ガイド

backend 名だけから速度・メモリ量・化学的精度を推定しないでください。
backend、model、precision、package version、GPU、原子数、workflow option
を記録し、実際の計算と同条件で pilot を行います。UMA はデフォルトで、
multi-worker inference を持つ唯一の built-in backend です。ORB は明示 fp32
の `float32-high`/TF32 matmul が有限差分 Hessian を noisy にし得るため fp64
がデフォルトです。MACE は対応 package stack が `fairchem-core` と競合する
ため、現状は専用環境が必要です。AIMNet2 の charge/spin 対応範囲は選択した
model に依存します。どの backend でも、optimizer の収束、意図した有意な
虚振動 1 本、正しい IRC endpoint を確認して TS を採用し、科学的主張に
必要なら cross-backend または DFT で検証してください。

## GPU メモリ (VRAM) 目安

VRAM は原子数だけでは決まりません。model architecture、neighbor 数、
precision、Hessian mode、凍結自由度、software version が影響します。同じ
設定の代表的な pilot でピークを測り、一時的な増加の余裕を残してください。
`torch.cuda.OutOfMemoryError` では `FiniteDifference` を維持または選択し、
対応している小さい model を使うか、化学と境界位置を確認して cluster を
縮小します。

## 不具合報告のときに添えると助かる情報

- 実行したコマンド（コピペ可能な形）
- `summary.log`（またはコンソール出力）
- 再現する最小入力（可能なら）
- OS / Python / CUDA / PyTorch バージョン
