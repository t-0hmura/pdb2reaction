# インストール

`pdb2reaction` は Linux 環境（ワークステーションや HPC）向けで、本番計算では通常 CUDA 対応 GPU を使用します。ビルド済みの **PyTorch** wheel は CUDA ランタイムを同梱するため、通常は互換性のある NVIDIA ドライバーだけが必要です。

## クイックスタート

`nvidia-smi` の右上に出る `CUDA Version` は、ドライバーが扱える最も新しい CUDA です。PyTorch の wheel は、これ以下の `cu126`・`cu130`・`cu132` から選びます。以下のコマンドは推奨の `cu130` を使います。

### 必須

```bash
# 1) conda 環境を作成してアクティブ化
# 2) CUDA 対応の PyTorchビルドをインストール
# 3) pdb2reactionをインストール
# 4) Plotly 静的画像 (PNG) エクスポート用のヘッドレス Chrome をインストール
#    Chromium のバイナリをダウンロード（インターネット接続が必要）

conda create -n p2r python=3.12 -y
conda activate p2r
TORCH_INDEX=cu130  # 推奨。cu126 / cu132 も選べる
pip install 'torch==2.13.0' --index-url "https://download.pytorch.org/whl/${TORCH_INDEX}"
pip install pdb2reaction
plotly_get_chrome -y
```

最後に、UMA モデルをダウンロードできるように **Hugging Face Hub** にログインします。無料の HF アカウントと読み取り専用トークンが要り、UMA モデルのページでライセンスへの同意が要ることもあります:

```bash
hf auth login
# またはスクリプト内でトークン指定する場合:
hf auth login --token '<YOUR_ACCESS_TOKEN>' --add-to-git-credential
```

これはマシン/環境ごとに 1 回だけ行う必要があります。

`pdb2reaction --version` でバージョンが表示されれば、インストールは成功です。

### 任意

DMF を使う場合は、環境をアクティブ化した直後に、{ref}`詳細なインストール手順 <ja-step-by-step-installation>` の手順 3 で cyipopt も導入してください。ORB・AIMNet2・MACE・DFT の導入は同じ節の手順 7 を参照してください。

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
    conda create -n <your-env> python=3.12 -y
    conda activate <your-env>
    ```

3. **cyipopt をインストール**
    `--mep-mode dmf`（DMF 法）で MEP を探索するときに必要です。GSM のみを使用する場合はスキップできます。それでも `--mep-mode dmf` がインポートのエラーで止まるときは、{ref}`インストール / 環境の問題 <ja-installation-environment-problems>` を参照してください。

    ```bash
    conda install -c conda-forge cyipopt -y
    ```

4. **適切な CUDA ビルドの PyTorch をインストール**

    推奨の例（`cu130`）:

    ```bash
    pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130
    ```

    PyTorch 2.13.0 の公式 wheel には `cu126`、`cu132`、`cpu` もあります。上のクイックスタートにある `nvidia-smi` の見方で wheel を選び、手順 8 で GPU が使えるかを確かめてください。[PyTorch の版の対応表](https://pytorch.org/get-started/previous-versions/) を参照してください。

5. **`pdb2reaction` 本体と可視化用 Chrome をインストール**

    ```bash
    pip install pdb2reaction
    plotly_get_chrome -y
    ```

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

    **ORB**（Python 3.11／3.12 が必要、3.12 推奨）:

    ```bash
    pip install "pdb2reaction[orb]"
    ```

    **AIMNet2**:

    ```bash
    pip install "pdb2reaction[aimnet]"
    ```

    **MACE**: `mace-torch` が要求する `e3nn==0.4.4` が UMA の `fairchem-core` と衝突するので、別の conda 環境に入れます。環境の名前 `mace-env` は変えても構いません。PyTorch の wheel は手順 4 と同じように選びます。

    ```bash
    conda create -n mace-env python=3.11 -y
    conda activate mace-env
    pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130
    pip install pdb2reaction
    pip uninstall -y fairchem-core
    pip install 'mace-torch>=0.3.8'
    ```

    **DFT**（`-b dft`、`--dft`、`pdb2reaction dft`）: `[dft]` は Linux x86_64 で CUDA 13 用の GPU4PySCF を導入します（手順 4 の cu130 / cu132 の PyTorch 向け）。cu126 の wheel を使うときは、代わりに `[dft-cuda12]` を導入します。aarch64 では [GPU4PySCF](https://github.com/pyscf/gpu4pyscf) をソースからビルドしてください。GPU のない環境でも `[dft]` で PySCF が入り、`--dft-engine cpu` で DFT を実行できます。

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

**GPU / CUDA:** 選んだ wheel に対応するドライバーの NVIDIA GPU（クイックスタートを参照）。新しい GPU アーキテクチャでは新しい wheel が必要なことがあります。CPU のみでも実行できますが、通常は大幅に遅くなります。

**VRAM・RAM・ディスク:** メモリはモデル、原子数、Hessian の計算方式とともに増え、ディスクには環境、モデルの重み、軌跡と Hessian が入ります。既定の UMA で Hessian を計算するときの VRAM のおおよその目安は、既定の有限差分（`FiniteDifference`）で 8 GB に約 900 原子、16 GB に約 2,000 原子、24 GB に約 3,000 原子、96 GB に約 1 万原子、解析 Hessian（`Analytical`）で 8 GB に約 200 原子、16 GB に約 400 原子、24 GB に約 600 原子、96 GB に約 1,500 原子です。代表的な計算を 1 つ対象の計算ノードで流し、最大使用量を見てください。

## 次のステップ

- [はじめに](getting-started.md) — 最短の実行と、次に読むページ
- [クイックスタート: `pdb2reaction all`](quickstart-all.md) — R と P から MEP を作る
- [クイックスタート: `pdb2reaction all --scan-lists`](quickstart-scan.md) — 1 つの構造から経路を作る
- [クイックスタート: TS-only モード](quickstart-tsopt.md) — TS 候補を最適化して確かめる
- [MLIP の TS を DFT で確かめる](dft-backend.md) — TS を DFT で詰めて確かめる
- [共通オプションと残基・原子の指定](cli-conventions.md) — 共通のオプションと、残基・原子の指定の書き方
- [HPC 実行例](hpc-example.md) — UMA のワーカーを複数の GPU ノードに広げるジョブスクリプト
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと試すこと
