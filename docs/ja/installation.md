# インストール

`pdb2reaction` は Linux 環境（ワークステーションや HPC）向けで、本番計算では通常 CUDA 対応 GPU を使用します。ビルド済みの **PyTorch** wheel は CUDA ランタイムを同梱するため、通常は互換性のある NVIDIA ドライバーだけが必要です。

## 最短のインストール

`nvidia-smi` の右上に出る `CUDA Version` は、ドライバーが扱える最も新しい CUDA です。PyTorch の wheel は、これ以下の `cu126`・`cu130`・`cu132` から選びます。以下のコマンドは推奨の `cu130` を使います。

### 必須

```bash
# 1) Python 3.12 の環境を作成してアクティブ化（例は conda。venv などでもよい）
# 2) CUDA 対応の PyTorch ビルドをインストール
# 3) pdb2reaction をインストール
# 4) Plotly 静的画像 (PNG) エクスポート用のヘッドレス Chrome をインストール
#    Chromium のバイナリをダウンロード（インターネット接続が必要）

conda create -n p2r -c conda-forge --override-channels python=3.12 -y
conda activate p2r

pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130  # 推奨。cu126 / cu132 なら末尾を置き換え
pip install pdb2reaction

plotly_get_chrome -y
```

UMA にはコンパイラが、PNG の出力には Chrome のライブラリが要ります。足りないときは先に入れてください: {ref}`システムの前提パッケージ <ja-system-prerequisites>`。

最後に、UMA モデルをダウンロードできるように **Hugging Face Hub** にログインします。無料の HF アカウントと読み取り専用トークンが要り、UMA モデルのページでライセンスへの同意が要ることもあります:

```bash
hf auth login
# またはスクリプト内でトークン指定する場合:
hf auth login --token '<YOUR_ACCESS_TOKEN>' --add-to-git-credential
```

ログインはマシン・環境ごとに 1 回行えば十分です。

`pdb2reaction --version` でバージョンが表示されれば、インストールは成功です。

### 任意

DMF（`--mep-mode dmf`）か `--solvent`（xTB）を使う場合は、代わりに cyipopt・xTB・PyTorch を conda-forge から入れて環境を作り、pdb2reaction を導入します:

```bash
conda create -n p2r -c conda-forge --override-channels python=3.12 cyipopt xtb "pytorch=2.13.0=cuda130*" cuda-nvcc -y
conda activate p2r
pip install pdb2reaction
plotly_get_chrome -y
```

`nvidia-smi` の `CUDA Version` が 13.0 より小さいときは、`cuda130` を `cuda129` に置き換えてください。ORB・AIMNet2・MACE・DFT の導入は {ref}`詳細なインストール手順 <ja-step-by-step-installation>` の手順 7 を参照してください。

(ja-system-prerequisites)=
## システムの前提パッケージ

pip と conda は、コンパイラと Chrome が使うライブラリを入れません。ワークステーションには通常ありますが、最小構成のコンテナや素のクラスターノードにはないことが多いです。

- **UMA 用のコンパイラ:** UMA は初回の実行時にコードをビルドします。GPU では C コンパイラ（`gcc`）、CPU では C++ コンパイラ（`g++`）が要ります。上の「任意」の conda 環境では、`cuda-nvcc` が C コンパイラを入れます。
- **Chrome 用のライブラリ:** PNG の図（`energy-diagram`、`trj2fig`、`all` の図）はヘッドレス Chrome で描きます。Chrome は `plotly_get_chrome -y` でダウンロードし、下の OS のライブラリを読み込みます。

Ubuntu/Debian（Ubuntu 22.04 と Debian 12 では `libasound2t64` の代わりに `libasound2` と書きます。root なら `sudo` は不要）:

```bash
sudo apt-get update
sudo apt-get install -y gcc g++
sudo apt-get install -y ca-certificates fonts-liberation libasound2t64 libatk-bridge2.0-0 \
  libatk1.0-0 libatspi2.0-0 libcairo2 libcups2 libcurl4 libdbus-1-3 libexpat1 libgbm1 \
  libglib2.0-0 libgtk-3-0 libnspr4 libnss3 libpango-1.0-0 libudev1 libvulkan1 libx11-6 libxcb1 \
  libxcomposite1 libxdamage1 libxext6 libxfixes3 libxkbcommon0 libxrandr2 wget xdg-utils
```

Rocky/RHEL:

```bash
sudo dnf install -y gcc gcc-c++
sudo dnf install -y alsa-lib at-spi2-atk at-spi2-core atk ca-certificates cairo cups-libs \
  dbus-libs expat glib2 gtk3 libX11 libXcomposite libXdamage libXext libXfixes libXrandr \
  libgcc liberation-fonts libxcb libxkbcommon mesa-libgbm nspr nss nss-util pango systemd-libs \
  vulkan-loader wget xdg-utils
```

sudo が使えない HPC では、`module load gcc` か、環境へのコンパイラの追加（`conda install -c conda-forge c-compiler cxx-compiler`）を使います。Chrome のライブラリは管理者に依頼してください。入っていなくても `all` と `path-search` は最後まで実行され、PNG だけが作られません。

どれかが足りないときは、コマンドが止まり、`Missing system prerequisite:` に続けて同じコマンドを表示します。

(ja-step-by-step-installation)=
## 詳細なインストール手順

環境を段階的に構築する場合:

1. **クラスターやビルドが必要とする場合のみ CUDA toolkit をロード**

    ビルド済みの PyTorch wheel に `nvcc` は不要です。依存パッケージをソースからビルドする場合は `module avail cuda` を確認し、クラスターが指定するコンパイラと CUDA toolkit の組み合わせをロードしてください:

    ```bash
    module load cuda/<your-version>   # 例: cuda/12.6 または cuda/12.9
    ```

2. **conda 環境を作成してアクティブ化**

    ```bash
    conda create -n <your-env> -c conda-forge --override-channels python=3.12 -y
    conda activate <your-env>
    ```

3. **cyipopt・xTB・PyTorch を conda-forge からインストール**
    `--mep-mode dmf`（DMF 法）で MEP を探索するときと、xTB の溶媒補正（`--solvent`）を使うときに必要です。この手順で PyTorch も入るので、手順 4 は飛ばします。GSM のみを使い `--solvent` も使わない場合はスキップできます。それでも `--mep-mode dmf` がインポートのエラーで止まるときは、{ref}`インストール / 環境の問題 <ja-installation-environment-problems>` を参照してください。

    ```bash
    conda install -c conda-forge --override-channels cyipopt xtb "pytorch=2.13.0=cuda130*" cuda-nvcc -y   # nvidia-smi の CUDA が 13.0 より小さいなら cuda129
    ```

4. **適切な CUDA ビルドの PyTorch をインストール**

    手順 3 で PyTorch を入れた場合は、この手順を飛ばします。推奨の例（`cu130`）:

    ```bash
    pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130
    ```

    PyTorch 2.13.0 の公式 wheel には `cu126`、`cu132`、`cpu` もあります。上の最短のインストールにある `nvidia-smi` の見方で wheel を選び、手順 8 で GPU が使えるかを確かめてください。[PyTorch の版の対応表](https://pytorch.org/get-started/previous-versions/) を参照してください。

5. **`pdb2reaction` 本体と可視化用 Chrome をインストール**

    ```bash
    pip install pdb2reaction
    plotly_get_chrome -y
    ```

    UMA が使うコンパイラと Chrome が使うライブラリはシステムのパッケージです: {ref}`システムの前提パッケージ <ja-system-prerequisites>` を参照してください。

    aarch64 では、UMA が `fairchem-core` を通して読み込む `warp-lang` 1.18.0 の wheel が glibc 2.35 を必要とし、それより古いシステム（`ldd --version` で確認できます）では import が `GLIBC_2.35' not found` で止まります。その場合は `warp-lang==1.17.0` をインストールしてください（`nvalchemi-toolkit-ops` には 1.13.0 以上が必要です）。

6. **Hugging Face Hub (UMA モデル) にログイン**

    ```bash
    hf auth login
    ```

    利用許諾と非対話型ログインは、上の「必須」を参照してください。

    詳細は上流プロジェクトを参照してください:

    - fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
    - Hugging Face トークンとセキュリティ: <https://huggingface.co/docs/hub/security-tokens>

7. **（任意）追加の MLIP バックエンドをインストール**

    pdb2reaction はデフォルトで UMA を使用します。他のバックエンドは対応する追加パッケージを導入し、`-b/--backend`（例: `-b orb`）で選択します:

    **ORB**: `orb-models` は `nvalchemi-toolkit-ops<0.4` に固定しているので、同じ環境に入れると UMA の `fairchem-core` が古い版に下がります。下の MACE と同じく別の conda 環境に入れ、UMA は元の環境で使います。ORB は Python 3.12 でだけ使えます（`orb-models` が固定する `dm-tree==0.1.8` に Python 3.13 用の wheel がないためです）。環境の名前 `orb-env` は変えても構いません。PyTorch の wheel は手順 4 と同じように選びます。

    ```bash
    conda create -n orb-env -c conda-forge --override-channels python=3.12 -y
    conda activate orb-env
    pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130
    pip install "pdb2reaction[orb]"
    ```

    **AIMNet2**:

    ```bash
    pip install "pdb2reaction[aimnet]"
    ```

    **MACE**: `mace-torch` が要求する `e3nn==0.4.4` が UMA の `fairchem-core` と衝突するので、別の conda 環境に入れます。環境の名前 `mace-env` は変えても構いません。PyTorch の wheel は手順 4 と同じように選びます。

    ```bash
    conda create -n mace-env -c conda-forge --override-channels python=3.12 -y
    conda activate mace-env
    pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130
    pip install pdb2reaction
    pip uninstall -y fairchem-core
    pip install 'mace-torch>=0.3.16'
    ```

    **DFT**（`-b dft`、`--dft`、`pdb2reaction dft`）: `[dft]` は Linux x86_64 で CUDA 13 用の GPU4PySCF を導入します（手順 4 の cu130 / cu132 と、手順 3 の `cuda130` の PyTorch 向け）。cu126 の wheel か `cuda129` のビルドを使うときは、代わりに `[dft-cuda12]` を導入します。aarch64 では [GPU4PySCF](https://github.com/pyscf/gpu4pyscf) をソースからビルドしてください。GPU のない環境でも `[dft]` で PySCF が入り、`--dft-engine cpu` で DFT を実行できます。

    ```bash
    pip install "pdb2reaction[dft]"
    ```

    DFT を使う場面と、MLIP の TS を DFT で確かめる方法は [MLIP の TS を DFT で確かめる](dft-backend.md) を参照してください。

8. **インストールの確認**

    ```bash
    pdb2reaction --version
    hf auth whoami
    ```

    1 行目でインストールされたバージョンが、2 行目で Hugging Face のユーザー名（UMA のダウンロード用）が表示されます。GPU アクセスの確認:

    ```bash
    python -c "import torch; print('CUDA:', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"
    ```

    `CUDA: False` の場合、バージョンを変える前に、インストールした wheel、スケジューラがジョブに GPU を見せているか、ドライバー、環境のライブラリを確認してください:

    ```bash
    python -m torch.utils.collect_env
    python -m pip check
    ```

## システム要件

**OS:** Linux、または Windows の WSL。

**GPU / CUDA:** 選んだ wheel に対応するドライバーの NVIDIA GPU（最短のインストールを参照）。新しい GPU アーキテクチャでは新しい wheel が必要なことがあります。CPU のみでも実行できますが、通常は大幅に遅くなります。

**VRAM・RAM・ディスク:** メモリはモデル、原子数、Hessian の計算方式とともに増え、ディスクには環境、モデルの重み、軌跡と Hessian が入ります。デフォルトの UMA で Hessian を計算するときの VRAM のおおよその目安は、デフォルトの有限差分（`FiniteDifference`）で 8 GB に約 900 原子、16 GB に約 2,000 原子、24 GB に約 3,000 原子、96 GB に約 1 万原子、解析 Hessian（`Analytical`）で 8 GB に約 200 原子、16 GB に約 400 原子、24 GB に約 600 原子、96 GB に約 1,500 原子です。VRAM に収まっても、1,000 原子を超えると Hessian 1 回に 1 時間以上かかり、反応機構解析の全体の時間が長くなります。巨大な系を扱う場合は、代表的な計算を 1 つ対象の計算ノードで実行し、最大使用量を見てください。

## 次のステップ

- [はじめに](getting-started.md) — 最短の実行と、次に読むページ
- [クイックスタート: `pdb2reaction all`](quickstart-all.md) — R と P から MEP を作る
- [クイックスタート: `pdb2reaction all --scan-lists`](quickstart-scan.md) — 1 つの構造から経路を作る
- [クイックスタート: TS-only モード](quickstart-tsopt.md) — TS 候補を最適化して確かめる
- [MLIP の TS を DFT で確かめる](dft-backend.md) — TS を DFT で精密化して確かめる
- [共通オプションと残基・原子の指定](cli-conventions.md) — 共通のオプションと、残基・原子の指定の書き方
- [HPC 実行例](hpc-example.md) — UMA のワーカーを複数の GPU ノードに広げるジョブスクリプト
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと試すこと
