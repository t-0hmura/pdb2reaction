# インストール

`pdb2reaction` は、CUDA 対応 GPU を備えた Linux 環境（ローカルワークステーションまたは HPC クラスター）が必要です。特に **PyTorch**、**fairchem-core (UMA)**、**gpu4pyscf-cuda12x** などの依存関係は、動作する CUDA インストールを前提としています。

詳細は上流プロジェクトを参照してください:

- fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
- Hugging Faceトークンとセキュリティ: <https://huggingface.co/docs/hub/security-tokens>

## クイックスタート

以下は多くの CUDA 12.9 クラスターで動作する最小限のセットアップ例です。モジュール名やバージョンはお使いの環境に合わせて調整してください。この例はデフォルトの GSM による MEP 探索（`--mep-mode gsm`）を前提としています。DMF（`--mep-mode dmf`）を使用する場合は、先に conda で cyipopt をインストールしてください。

### 必須

```bash
# 1) CUDA 対応の PyTorchビルドをインストール
# 2) pdb2reactionをインストール
# 3) Plotly 静的画像 (PNG) エクスポート用のヘッドレス Chrome をインストール
#    ~150 MB の Chromium バイナリをダウンロード（インターネット接続必要）

pip install torch --index-url https://download.pytorch.org/whl/cu129
pip install pdb2reaction
plotly_get_chrome -y
```

最後に、UMA モデルをダウンロードできるように **Hugging Face Hub** にログインします（無料の HF アカウントと読み取り専用トークンが必要。<https://huggingface.co/facebook/UMA> でモデルライセンスの承認が必要な場合あり）:

```bash
huggingface-cli login
# または
hf auth login --token '<YOUR_ACCESS_TOKEN>' --add-to-git-credential
```

これはマシン/環境ごとに1回だけ行う必要があります。

### 任意

- MEP 探索で Direct Max Flux（DMF）法を使用する場合は、conda 環境を作成し、pdb2reaction のインストール前に cyipopt をインストールしてください。

  ```bash
  # 専用のconda環境を作成してアクティブ化
  conda create -n pdb2reaction python=3.11 -y
  conda activate pdb2reaction

  # cyipoptをインストール（MEP 探索のDMF法に必要）
  conda install -c conda-forge cyipopt -y
  ```

- *環境モジュール*を使用する HPC クラスターでは、PyTorch をインストールする**前に** CUDA をロードしてください。`module avail cuda` で利用可能なバージョンを確認し、ターゲットの PyTorch wheel に合うバージョン（例: `cu126` ↔ CUDA 12.6、`cu129` ↔ CUDA 12.9）をロードしてください:

  ```bash
  module load cuda/<your-version>   # 例: cuda/12.6 または cuda/12.9
  ```

> **ヒント:** UMA がデフォルトの MLIP バックエンドです。ORB や AIMNet2 を使用するには、対応する extra をインストール（例: `pip install "pdb2reaction[orb]"`）し、コマンドに `-b/--backend orb` を渡してください。[詳細なインストール手順](#ja-step-by-step-installation)の手順 7 を参照してください。

```{warning}
**MACE:** MACE は `e3nn==0.4.4` を必要としますが、`fairchem-core`（UMA）と競合します。正準の MACE 導入手順は `pip uninstall -y fairchem-core && pip install mace-torch` です。UMA と MACE は同一環境で共存できないため、両方必要な場合は別々の conda 環境を使ってください。（古いメモにある `--no-deps mace-torch` 方式は torch-scatter / e3nn が pin されないため推奨しません。）
```


(ja-step-by-step-installation)=
## 詳細なインストール手順

環境を段階的に構築する場合:

1. **CUDAをロード（HPCで環境モジュールを使用する場合）**

    `module avail cuda` で利用可能なバージョンを確認し、ターゲットの
    PyTorch wheel に合うバージョン（例: `cu126` は CUDA 12.6、`cu129`
    は CUDA 12.9）をロードしてください:

    ```bash
    module load cuda/<your-version>
    ```

2. **conda環境を作成してアクティブ化**

    ```bash
    conda create -n pdb2reaction python=3.11 -y
    conda activate pdb2reaction
    ```

3. **cyipopt をインストール**
    MEP 探索で DMF 法（`--mep-mode dmf`）を使用する場合に必要です。GSM のみを使用する場合はスキップできます。

    ```bash
    conda install -c conda-forge cyipopt -y
    ```

4. **適切なCUDAビルドのPyTorchをインストール**

    CUDA 12.9の場合:

    ```bash
    pip install torch --index-url https://download.pytorch.org/whl/cu129
    ```

    PyTorch は CUDA ドライバーバージョンに合わせたビルドが必要です。互換性は [PyTorch Get Started](https://pytorch.org/get-started/locally/) で確認してください。CPU のみの実行もサポートされますが、大幅に遅くなります（10-100 倍）。

5. **`pdb2reaction` 本体と可視化用Chromeをインストール**

    ```bash
    pip install pdb2reaction
    plotly_get_chrome -y
    ```

6. **Hugging Face Hub (UMAモデル) にログイン**

    ```bash
    hf auth login
    ```

    参照:

    - <https://github.com/facebookresearch/fairchem>
    - <https://huggingface.co/facebook/UMA>
    - <https://huggingface.co/docs/hub/security-tokens>

7. **（任意）追加の MLIP バックエンドをインストール**

    pdb2reaction はデフォルトで UMA を使用します。他のバックエンドを使用する場合は、対応するオプション依存関係をインストールしてください:

    ```bash
    # ORB バックエンド
    pip install "pdb2reaction[orb]"

    # AIMNet2 バックエンド
    pip install "pdb2reaction[aimnet]"

    # MACE バックエンド (UMA の fairchem-core と競合 — 必ず別 conda 環境で実施。
    # 既定の UMA 環境で実行すると UMA が永久破壊される)
    # cd <separate_mace_env>
    # pip uninstall -y fairchem-core && pip install mace-torch

    # DFT 一点計算の後処理（`--dft` / `pdb2reaction dft`）
    # gpu4pyscf-cuda12x、PySCF、および関連依存をインストールします。
    # 注: gpu4pyscf-cuda12x は x86_64 wheel のみ。aarch64 では
    # `pdb2reaction dft --engine gpu` が ClickException を投げるため
    # `--engine cpu` (PySCF) を使うか、[dft] extras をスキップ。
    pip install "pdb2reaction[dft]"
    ```

    暗黙溶媒補正を使用するには、[xTB](https://github.com/grimme-lab/xtb) をインストールし、`xtb` コマンドが `PATH` 上で利用可能であることを確認してください。

    #### xTB のインストール

    **ALPB 溶媒和モデルの場合**（推奨の出発点）:

    ```bash
    conda install -c conda-forge xtb
    ```

    **CPCM-X 溶媒和モデルの場合**（ソースからのビルドが必要）:

    ```bash
    git clone --depth 1 https://github.com/grimme-lab/xtb.git
    cd xtb
    cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DWITH_CPCMX=ON
    make -C build -j8
    ```

    GCC >= 10 が必要です。実行時に `CPXHOME` を `build/_deps/cpcmx-src/` に設定してください。

    カスタム xTB バイナリを使用するには、YAML 設定で `xtb_cmd` キーを設定するか、Python で `calc.xtb_cmd` を使用してください。

8. **インストールの確認**

    ```bash
    pdb2reaction --version
    ```

    インストールされたバージョンが表示されます。GPU アクセスの確認:

    ```bash
    python -c "import torch; print('CUDA:', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"
    ```

    `CUDA: False` の場合、CUDA モジュールのロードと PyTorch ビルドの CUDA バージョンを確認してください。

## 次の導線

- [はじめに](getting-started.md) — プロジェクト概要、パイプラインの各ステージ、ワークフローモード
- [クイックスタート: `pdb2reaction all`](quickstart-all.md) — 2 つの PDB から end-to-end 実行
- [クイックスタート: 単一構造スキャン](quickstart-scan.md) — `--scan-lists` で 1 つの PDB から MEP
- [クイックスタート: TS 最適化](quickstart-tsopt-freq.md) — TS 候補の最適化と検証
- [CLI 規約](cli-conventions.md) — フラグの優先順位、原子/残基セレクタ、共通オプション
- [トラブルシューティング](troubleshooting.md) と [典型エラー別レシピ](recipes-common-errors.md)

