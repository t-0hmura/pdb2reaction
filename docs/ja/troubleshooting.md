# トラブルシューティング

症状を早見表で探し、示された節で対処を読んでください。

(ja-troubleshooting-quick-table)=
## 早見表

| 症状 | 最初にやること | 詳細（節） |
| --- | --- | --- |
| **入力 / 抽出** | | |
| 元素列が空で `extract` が止まる（`Element symbols are missing in '...'`）。`all` は空の元素列を自分で埋め、割り当てられない原子が残ると止まる | 元の PDB に `add-elem-info` を適用してください | {ref}`入力 / 抽出の問題 <ja-input-extraction-problems>` |
| `[multi] Atom count mismatch` / `[multi] Atom order mismatch` | 同じ前処理ツール・同じ設定で全 PDB を作り直してください。原子の順序を決めた後は並べ替えません | {ref}`入力 / 抽出の問題 <ja-input-extraction-problems>` |
| ユーティリティのコマンド（`add-elem-info`、`fix-altloc`、`bond-summary`、`energy-diagram`、`trj2fig`）がエラーで止まる、または警告を出す | そのコマンドのページでメッセージを探してください。エラーは「使用上の注意点」に、`add-elem-info` の `[WARN] Could not confidently assign` は「主な出力ファイル」にあります | [add-elem-info](add-elem-info.md#主な出力ファイル)、[fix-altloc](fix-altloc.md#使用上の注意点)、[bond-summary](bond-summary.md#使用上の注意点)、[energy-diagram](energy-diagram.md#使用上の注意点)、[trj2fig](trj2fig.md#使用上の注意点) |
| `bond-summary` がある組で `ERROR: Atom types and ordering must be identical.` を出す | 同じ原子が同じ順に並ぶ構造どうしを比べてください | [bond-summary](bond-summary.md#使用上の注意点)、{ref}`入力 / 抽出の問題 <ja-input-extraction-problems>` |
| `fix-altloc` が `Output exists: <path> (use --overwrite to overwrite)` で止まる | `--overwrite` を付けるか、`-o` で別の出力先を指定してください | [fix-altloc](fix-altloc.md#主な-cli-オプション) |
| **電荷 / スピン** | | |
| `-q/--charge is required` / `Total charge could not be resolved` | `-q/--charge` または `-l/--ligand-charge` を明示してください | {ref}`電荷 / スピンの問題 <ja-charge-spin-problems>` |
| `Cluster electron count inconsistent`（`all --dry-run` では `--dry-run parity check failed`） | 電荷と多重度が電子数と合いません。`-q`・`-l` を直すか、`-m` を設定してください（電子数が奇数なら `-m 2` など） | {ref}`電荷 / スピンの問題 <ja-charge-spin-problems>` |
| 計算は通るが状態やエネルギーが不自然 | 渡した電荷と多重度を見直してください | {ref}`電荷 / スピンの問題 <ja-charge-spin-problems>` |
| **計算 / 収束** | | |
| `all` が `success` 以外の `Scientific status:` で終わり、`RESULT WARNING:` の行を出す | 各行に、失敗したセグメントと段（TS 最適化・MEP・IRC・端点の最適化）が出ます。その段を下の行で探してください。`summary.json` の `scientific_status_reasons` にも同じ理由が入ります（例：`all:segment_1:tsopt:ts_optimization_not_converged`） | [実行と要求段階の完了状況](json-output.md#実行と要求段階の完了状況) |
| UMA で `--uma-workers` を 2 以上にし、`--hessian-calc-mode Analytical` と併用すると `BackendError` | 解析 Hessian には `--uma-workers 1`、並列実行には `FiniteDifference` を指定してください | {ref}`パフォーマンス <ja-troubleshooting-performance>`、{ref}`workers と解析 Hessian <ja-workers-analytical-error>` |
| 実行時に CUDA のメモリ不足（`torch.cuda.OutOfMemoryError`） | 既定の `FiniteDifference` Hessian のままにする、`--max-nodes` を減らすか小さい MLIP モデル（`--backend-model`）を使う、VRAM の大きい GPU に移る、の順に試してください。`--radius` を小さくして切り出し直すのは最後の手です | {ref}`GPU メモリ <ja-troubleshooting-gpu-memory>` |
| TS 最適化は収束したが n_imag が 1 でない | n_imag ≥ 2（`TS imaginary-mode validation found n_imag=…`）：`--flatten` を付けてください（`tsopt`・`opt`・`all` で使えます）。n_imag = 0（`[tsopt] No imaginary mode detected.`）：別の候補から始めるか、`all` で `--refine-path` を使ってください | {ref}`TS 最適化 <ja-troubleshooting-ts>`、{ref}`TS が取れないとき <ja-ts-search-fails>` |
| TS 最適化が収束しない（`TS optimization did not converge`） | まず TS 候補を確かめ、次にオプティマイザを切り替え（`tsopt --opt-mode` / `all --opt-mode-post`）、それから YAML でステップサイズか信頼半径を小さくしてください | {ref}`TS 最適化 <ja-troubleshooting-ts>`、{ref}`TS が取れないとき <ja-ts-search-fails>` |
| IRC が正常に終了しない | 単独の `irc`：`--step-size` を小さく、`--max-cycles` を大きく。`all`：`--irc-step-size` / `--irc-max-cycles`。先に端点を確かめてください | {ref}`IRC <ja-troubleshooting-irc>` |
| エネルギーが変わらなくなり、opt / TS 最適化が `stalled` で止まる（`TS optimization status is stalled`） | 未収束として扱い、final geometry と力・ステップの条件を確かめてから、別の閾値やオプティマイザの設定で再実行してください | {ref}`max_cycles とプラトー停止 <ja-troubleshooting-max-cycles>` |
| MEP（最小エネルギー経路）探索（GSM / DMF）が失敗する（`MEP optimization did not converge`） | `--max-nodes` を既定の 20 より増やす、`--preopt` を有効のままにする（`all`・`path-search`・`path-opt` では既定で有効、`scan`・`scan2d`・`scan3d` では既定で無効）、別の `--mep-mode` を試す | {ref}`MEP 探索 <ja-troubleshooting-mep>` |
| `freq` がエラーで止まる | 動ける原子を 1 つ以上残してください | {ref}`freq のエラー <ja-troubleshooting-freq>` |
| DFT の SCF が収束しない | `--no-dft-low-memory`（密度フィッティング）か、YAML のレベルシフト（`mf: {level_shift: 0.2}`）を試してください。書く場所は、`dft` と `--dft` では `dft.pyscf`、`-b dft` では `calc.dft.pyscf` です | [dft の使用上の注意点](dft.md#使用上の注意点)、[DFT バックエンドの使用上の注意点](dft-backend.md#使用上の注意点) |
| DFT で GPU メモリが足りない | 小さい基底関数を使う、モデルを削る、メモリの大きい GPU に移る。`--dft` で足りないときは、`pdb2reaction dft` を別に実行してください | [DFT バックエンドの使用上の注意点](dft-backend.md#使用上の注意点) |
| **インストール / 環境** | | |
| DMF モードのインポートエラー（`cyipopt`）、または `No module named 'dmf'` | `conda install -c conda-forge cyipopt`（`pydmf` は `pdb2reaction` と一緒に入ります）。既定の GPU 版でそれでもインポートに失敗するときは、`pip install 'pydmf[torch]'` | {ref}`インストール / 環境の問題 <ja-installation-environment-problems>` |
| UMA モデルで 401 / 403 / アクセス制限付きリポジトリのエラー（`huggingface_hub.errors.GatedRepoError`） | `hf auth login` でログインし、UMA モデルのライセンスに同意してください | {ref}`インストール / 環境の問題 <ja-installation-environment-problems>` |
| `e3nn` / `fairchem-core` のインポートの競合（UMA の環境に MACE を入れた） | MACE 専用の環境を使ってください | {ref}`インストール / 環境の問題 <ja-installation-environment-problems>` |
| `ORB backend requires orb-models and torch`（AIMNet2 / MACE も同様） | バックエンドの追加パッケージを入れてください：`pip install "pdb2reaction[orb]"`。MACE は別の環境に入れます | {ref}`インストール / 環境の問題 <ja-installation-environment-problems>` |
| CUDA / GPU の実行時エラー | GPU、PyTorch のビルド、ドライバをまとめて確かめてください | {ref}`インストール / 環境の問題 <ja-installation-environment-problems>` |
| 図の出力に失敗する | `plotly_get_chrome -y` でヘッドレス Chrome を入れてください | {ref}`インストール / 環境の問題 <ja-installation-environment-problems>` |

## 実行前チェックリスト

長い計算を回す前に、次を確かめてください。

- このマシンで Hugging Face にログインできている（既定の UMA に必要）。
- 入力の PDB/mmCIF に **水素** と **元素記号** が入っている。
- 複数の PDB を与える場合、**同じ原子が同じ順序** で並んでいる。

---

(ja-input-extraction-problems)=
## 入力 / 抽出の問題

### `Element symbols are missing in '...'`

- **症状**：`extract` が `Element symbols are missing in '...'. For PDB input, run pdb2reaction add-elem-info -i ... before extract` で止まる。`all` は抽出の前に空の元素列を自分で埋め、割り当てられない原子が残ると同じメッセージで止まる。
- **原因**：PDB の元素列（77–78 桁）が空のことが多く、`extract` は原子の種類を決めるために元素記号を使います。mmCIF の入力では `_atom_site.type_symbol` が必要です。
- **対処**：`add-elem-info` で元素列を埋め、新しいファイルで再実行してください。`[WARN] Could not confidently assign N atoms` の下に並んだ原子は、元素記号を手で 77–78 桁に右詰めで書いてください。

  ```bash
  pdb2reaction add-elem-info -i input.pdb -o input_with_elem.pdb
  ```

### `[multi] Atom count mismatch` / `[multi] Atom order mismatch`

- **症状**：複数の入力を与えた実行が `[multi] Atom count mismatch between input #1 and input #2: ...` や `[multi] Atom order mismatch between input #1 and input #2.` で止まる。
- **原因**：構造ごとに別のツールや設定で前処理した、またはプロトン化やパラメータ化をやり直した後に原子の順序が変わった。
- **対処**：**すべて** の構造を、同じプロトン化ツール・同じ設定で作り直してください。MD のスナップショットなら、同じトポロジーと軌跡からフレームを取り出します。トポロジーを作った後は、PDB の原子を並べ替えません。

### 活性部位モデルが空になる・触媒残基が入らない

- **症状**：切り出したモデルが想定より小さい、または触媒残基が含まれない。
- **原因**：この部位には半径 `-r/--radius`（既定 2.6 Å）が小さすぎる、または `--exclude-backbone` で削りすぎた。
- **対処**：`--radius` を大きくするか（例：2.6 → 3.5 Å）、残基を足してください。`--selected-resn 'A:TYR:44'` は少なくともその側鎖を足し、`-c` に足すと、`-r` が 0 より大きければ残基が丸ごと残ります。詳しくは {ref}`モデルを広げる <ja-model-setup-larger>` と [残基セレクタ](cli-conventions.md#残基セレクタ) にあります。`--exclude-backbone` を付けていて削りすぎる場合は、`--no-exclude-backbone` を渡してください。

### エネルギーや障壁がモデルの大きさで変わる

- **症状**：エネルギーや障壁が不自然に見える、またはモデルを大きくすると大きく変わる。
- **原因**：切り出したモデルが小さすぎる。
- **対処**：半径を大きくし、結果がモデルの大きさと境界の位置でどう変わるかを確かめてください。

  ```bash
  pdb2reaction extract -i complex.pdb -c 'SUB' -o model.pdb -r 4.0
  ```

### 修飾残基が切断されない

- **症状**：登録されていない修飾アミノ酸の主鎖が切られず、キャップ水素も付かない。
- **原因**：主鎖の切断とキャップ水素の付加には、アミノ酸カタログへの登録が必要です。SEP、TPO、MLY は登録済みです。
- **対処**：未登録の残基だけを、`--modified-residue "XAA:0"` のように公称電荷を付けて登録してください。名前だけを書くと電荷 0 になるため、電荷を持つ登録済みの残基の再登録には使えません。主鎖のトポロジーが特殊な場合は、[活性部位モデルを手で組み](model-setup.md)、下流のコマンドに直接渡してください。

---

(ja-charge-spin-problems)=
## 電荷 / スピンの問題

まず、総電荷と多重度が対象の状態に合っているか、`-l/--ligand-charge` の各残基名が構造にあるかを確かめてください。重要な実行では `-q/--charge` か `-l/--ligand-charge` と `-m` を明示します。決まりは {ref}`電荷の指定 <ja-charge-specification>` にあります。

### `-q/--charge is required` / `Total charge could not be resolved`

- **症状**：`.gjf` 以外の入力で、`-q/--charge is required unless the input is a .gjf template with charge metadata.`、`all` では `[all] Total charge could not be resolved.` で止まる。
- **原因**：`-q/--charge` を省くと、電荷は `-l/--ligand-charge`（PDB/mmCIF、または `--ref-pdb` 付きの XYZ/GJF）、YAML の `calc.charge`、`.gjf` テンプレートから決まります。そのどれも使えなかった。
- **対処**：電荷と多重度を明示するか、抽出ありの場合は残基ごとの電荷を与えてください。

  ```bash
  pdb2reaction path-search -i R.pdb P.pdb -q 0 -m 1
  pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
  ```

---

(ja-installation-environment-problems)=
## インストール / 環境の問題

まず、使っている環境にオプションのパッケージが入っているか、PyTorch から GPU が見えるかを確かめてください。直した後は `pdb2reaction --version` と `python -c "import torch; print(torch.cuda.is_available())"` で確かめ、本番の前に一度 `--dry-run` を付けて実行し、オプションと入力を確かめます。

| 症状 | 原因 | 対処 |
| --- | --- | --- |
| UMA のダウンロードに失敗する（`huggingface_hub.errors.GatedRepoError`、`401`、`403`） | Hugging Face にログインしていない、または UMA モデルのライセンスに同意していない | 環境・マシンごとに一度 `hf auth login` を実行し、Hugging Face の UMA モデルのページでライセンスに同意してください。HPC では、計算ノードから Hugging Face のキャッシュディレクトリに書き込めるかを確かめます |
| `ORB backend requires orb-models and torch`、`AIMNet2 backend requires torch and aimnet`、`Could not import mace.calculators because mace-torch is not installed` | バックエンドのパッケージがこの環境に入っていない | ORB：`pip install "pdb2reaction[orb]"`。AIMNet2：`pip install "pdb2reaction[aimnet]"`。MACE：別の環境（次の行） |
| `e3nn` / `fairchem-core` のインポートの競合 | UMA の環境に MACE を入れた。`mace-torch` は `e3nn==0.4.4` に固定し、`fairchem-core` は `e3nn>=0.5` を必要とする | MACE 専用の conda 環境を使ってください：`pip uninstall -y fairchem-core && pip install 'mace-torch>=0.3.8'` |
| 追加パッケージを入れても ORB のインポートに失敗する | 環境の中のパッケージが `orb-models` と競合している | `python -m pip check` を実行してください。PyG / `torch_scatter` は、エラーがそれを名指ししたときだけ入れます。今の `orb-models` はこれを必要としません |
| `torch.cuda.is_available()` が `False`、またはインポート時に CUDA の実行時エラー | PyTorch のビルドが計算ノードの GPU・ドライバに合っていない | `nvidia-smi`、`python -m torch.utils.collect_env`、`python -m pip check` で、割り当てられた GPU、入っている wheel、ドライバを確かめてください。真偽値だけでは原因は分かりません。`nvidia-smi` が示す `CUDA Version` はドライバが扱える最も新しい CUDA です。それ以下の CUDA の wheel（`cu126`、`cu130`、`cu132`）を入れてください。ローカルの CUDA toolkit は要りません |
| `--mep-mode dmf` が `DMF mode (--mep-mode dmf) requires ase, cyipopt, and pydmf>=1.2` で止まる、または `No module named 'dmf'` | `cyipopt` が無い、または `pydmf` の GPU 用の部分が無い | `conda install -c conda-forge cyipopt` を、できれば `pdb2reaction` を入れる前に実行してください。`pydmf` は `pdb2reaction` と一緒に入ります。既定の GPU 版 DMF でそれでもインポートに失敗するときは、エラーの文が勧めるとおり `pip install 'pydmf[torch]'` を実行してください |
| 図の出力に失敗する（Plotly / Chrome） | ヘッドレス Chrome が無い | `plotly_get_chrome -y` を一度実行してください |

### DMF が IPOPT 内で極端に遅い

IPOPT/MUMPS が並列版 BLIS を使う環境では、入れ子の並列化で長い待ちが生じることがあります。ジョブスクリプトなどで、Python や CLI の起動前に `BLIS_NUM_THREADS=1` を設定してください。`OMP_NUM_THREADS` などほかのスレッド設定は変えずにおきます。`BLIS_JC_NT`、`BLIS_PC_NT`、`BLIS_IC_NT`、`BLIS_JR_NT`、`BLIS_IR_NT` の手動設定はこの制限より優先されるため、そのジョブの設定から外してください。起動済みの Notebook は、設定変更後にカーネルを再起動します。詳しくは [BLIS のスレッド設定](https://github.com/flame/blis/blob/2.0/docs/Multithreading.md)を参照してください。

---

(ja-calculation-convergence-problems)=
## 計算 / 収束の問題

まず TS 候補を確かめてください。TS 最適化が成功すると、反応モードの虚振動が 1 つ出ます。ν < −5.00 cm⁻¹ のモードを虚振動として数えます。閾値は YAML の `freq.zero_cutoff_cm` で変えられます。TS 最適化が収束しないときは {ref}`TS 最適化 <ja-troubleshooting-ts>` を、収束したが n_imag が 1 でないときは {ref}`TS が取れないとき <ja-ts-search-fails>` を読んでください。

(ja-troubleshooting-max-cycles)=
### 最適化が `max_cycles` に達し、`max(force)` が閾値をわずかに超える

- **症状**：オプティマイザが `max_cycles` まで回り、最後の要約で `max(force)` や `rms(force)` が選んだ閾値をわずかに上回る。一方でエネルギーはもう変わっていない。
- **原因**：MLIP の力のノイズや平坦さで、力の閾値に届かないことがあります。その程度はバックエンド、モデル、精度、系、ハードウェアで変わります。
- **対処**：final geometry と力を確かめてから、次の手で再実行してください。
  - 収束の基準を `--thresh gau_loose` に緩めます。既定は `opt` では `gau`、`tsopt` では [`baker`](tsopt.md#処理の仕組みと計算仕様) です。
  - `--max-cycles` まで回さずに早く止めたいときは、`--stop-plateau` を付けます。直近 `--stop-plateau-window`（既定 50）ステップのエネルギーの幅が `--stop-plateau-thresh`（既定 `1×10⁻⁴ au`）を下回ると、`stalled`（**未収束**）として止まります。

(ja-troubleshooting-ts)=
### TS 最適化が収束しない・虚振動が複数残る

- **症状**：TS 最適化が多くのサイクルを回しても収束しない、または最適化の後に n_imag が 2 以上になる。
- **報告される内容**：RS-P-RFO・RS-I-RFO・TRIM・Dimer の TS 最適化が max cycles に達して未収束のときは Hessian を計算しないので、n_imag は出ません。エネルギーが変わらなくなって止まったときは、必ず Hessian を計算して n_imag を出します。
- **最適化が収束しないときの対処**：次を順に試してください。
  1. オプティマイザを RS-P-RFO（既定）と Dimer 法の間で切り替える：単独では `tsopt --opt-mode hess` / `dimer`、`all` では `--opt-mode-post hess` / `grad`（Dimer）。
  2. YAML でステップサイズを小さくする：RS-P-RFO・RS-I-RFO・TRIM では `rsirfo.trust_radius` / `trust_min` / `trust_max`、Dimer では `hessian_dimer.lbfgs.max_step`。[YAML リファレンス](yaml-reference.md#ts-最適化セクション) を参照してください。
  3. 経路のよりよい HEI（最高エネルギーのイメージ）など、別の候補から始める。
- **n_imag ≥ 2 が残るときの対処**：`--flatten` を付けて最適化し直すか、単独では `tsopt --thresh`、`all` では `--thresh-post` で、収束の基準を既定の `baker` から `gau_tight` か `gau_vtight` に締めてください。`--refine-path`、反応を段に分ける、初期構造を作り直すなどのほかの手は、{ref}`TS が取れないとき <ja-ts-search-fails>` にあります。

(ja-troubleshooting-irc)=
### IRC が正常に終了しない

IRC が収束せずに止まっても、端点の最適化で狙った R と P に着けば使えます。まず最適化後の端点を確かめてください。

- **症状**：IRC が明確な極小構造に着く前に止まる、またはエネルギーが振動し勾配ノルムが大きいままになる。
- **原因**：この曲面にはステップが大きすぎる、サイクル上限が低すぎる、または開始構造に虚振動が 2 つ以上ある。
- **対処**：
  - 単独の `irc`：`--step-size 0.05`（既定 0.10 bohr）、必要なら `--max-cycles 200`（既定 125）。
  - `all`：`--irc-step-size 0.05`、必要なら `--irc-max-cycles 200`。
  - 開始構造が n_imag = 1 であることを確かめてください。
  - 物理的な停止条件を無視してサイクル上限まで追うには、単独で `irc --never-stop`、`all` で `--irc-never-stop` を指定し、軌跡と端点を確かめてください。

(ja-troubleshooting-mep)=
### MEP 探索（GSM / DMF）が失敗する・結合変化を取りこぼす

- **症状**：経路探索が使える MEP を作らずに終わる、または予想した結合変化が出ない。
- **対処**：
  - 複雑な反応では `--max-nodes`（既定 20）を 30 や 40 に増やしてください。
  - 端点の事前最適化は既定で有効です。`--no-preopt` を付けていたら外してください。
  - 別の手法を試してください：`--mep-mode dmf` ↔ `gsm`。
  - YAML の `bond.bond_factor` と `bond.delta_fraction` で結合変化の検出を調整してください。

(ja-troubleshooting-freq)=
### `freq` がエラーで止まる

- **すべての原子を固定した**：動ける原子が無いと解析する振動が無いので、`freq` はエラーで止まります。`--freeze-atoms`、YAML の `geom.freeze_atoms`、`--freeze-links` で固定されるキャップ水素の親原子を確かめてください。{ref}`原子の固定と距離の拘束 <ja-freeze-atoms-and-restraints>` を参照してください。
- **`--uma-workers` を 2 以上にして `--hessian-calc-mode Analytical` と併用した**：{ref}`パフォーマンス <ja-troubleshooting-performance>` を参照してください。

---

(ja-troubleshooting-performance)=
## パフォーマンス / 安定性のヒント

- **workers > 1**: ハードウェアと計算の内容によっては UMA のスループットが上がりますが、並列の推論では解析 Hessian を計算できません。`Analytical` を明示すると、`RuntimeError` のサブクラスの `BackendError` が `Analytical Hessian cannot be combined with UMA workers>1` を出して止まります。解析 Hessian には `--uma-workers 1`、並列実行には `FiniteDifference` を指定してください。複数のノードでワーカーを動かすときは、[HPC 実行例](hpc-example.md) のジョブスクリプトを使ってください。
- **大規模系**: 化学的に妥当な小さい活性部位モデルを作り、半径と境界の位置で結果がどう変わるかを確かめてください。複数 GPU への対応はバックエンドとワークフローごとに異なり、GPU を増やしてもワーカーあたりのメモリが減るとは限りません。
- **HPC で DFT を回すとき**: PySCF/GPU4PySCF が選んだ計算で一時ディスクを使う場合は、容量と速度を確かめたファイルシステムを `PYSCF_TMPDIR` に指定してください。ノードの `/tmp`、`$PBS_O_WORKDIR`、ほかの共有ディレクトリのどれが向くかは、計算機ごとに異なります。

(ja-troubleshooting-gpu-memory)=
## GPU メモリ (VRAM) 目安

VRAM は原子数だけでは決まらず、MLIP モデル、数値精度、Hessian の計算モード、固定した原子で変わるので、同じ設定の代表的な試行で測ってください。`torch.cuda.OutOfMemoryError` が出たら、次を順に試します。

1. 既定の `--hessian-calc-mode FiniteDifference` を維持するか、それに切り替える。解析 Hessian は、通常メモリのピークが大きくなります。
2. `--max-nodes` を減らすか、`--backend-model` で小さい MLIP モデルを使う。`opt` と `scan` では、`hess` ではなく、Hessian を使わない {ref}`--opt-mode grad <ja-opt-mode-semantics>`（L-BFGS）のままにする。
3. メモリの大きい GPU に移る。
4. クラスターモデルを縮小する。必要な残基が残ることと境界の位置を確かめてからにします。

## 不具合報告のときに添えると助かる情報

実行したとおりのコマンド、`summary.log`（またはコンソール出力）、再現する最小の入力、OS・Python・CUDA・PyTorch の版を添えてください。

## 関連ドキュメント

- [反応機構を調べるコツ](mechanism-tips.md) — TS が取れないときに試すこと
- [インストール](installation.md) — 環境の構築とオプションのバックエンド
- [MLIP バックエンド](backends.md) — バックエンドの選び方
- [クラスターモデルの組み方](model-setup.md) — 活性部位モデルを確かめる・削る・広げる
