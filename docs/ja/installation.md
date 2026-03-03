# インストール

`pdb2reaction` は、CUDA 対応 GPU を備えた Linux 環境（ローカルワークステーションまたは HPC クラスター）向けに設計されています。特に **PyTorch**、**fairchem-core (UMA)**、**gpu4pyscf-cuda12x** などの依存関係は、動作する CUDA インストールを前提としています。

詳細は上流プロジェクトを参照してください:

- fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
- Hugging Faceトークンとセキュリティ: <https://huggingface.co/docs/hub/security-tokens>

## クイックスタート

以下は多くの CUDA 12.9 クラスターで動作する最小限のセットアップ例です。モジュール名やバージョンはお使いの環境に合わせて調整してください。この例はデフォルトの GSM による MEP 探索（`--mep-mode gsm`）を前提としています。DMF（`--mep-mode dmf`）を使用する場合は、先に conda で cyipopt をインストールしてください。

```bash
# 1) CUDA 対応の PyTorchビルドをインストール
# 2) GitHubからpdb2reactionをインストール
# 3) Plotly図表エクスポート用のヘッドレスChromeをインストール

pip install torch --index-url https://download.pytorch.org/whl/cu129
pip install git+https://github.com/t-0hmura/pdb2reaction.git
plotly_get_chrome -y
```

最後に、UMAモデルをダウンロードできるように **Hugging Face Hub** にログインします:

```bash
# Hugging Face CLI
hf auth login --token '<YOUR_ACCESS_TOKEN>' --add-to-git-credential
```

または

```bash
# クラシックCLI
huggingface-cli login
```

これはマシン/環境ごとに1回だけ行う必要があります。

- MEP 探索で Direct Max Flux（DMF）法を使用する場合は、conda 環境を作成し、pdb2reaction のインストール前に cyipopt をインストールしてください。
 ```bash
 # 専用のconda環境を作成してアクティブ化
 conda create -n pdb2reaction python=3.11 -y
 conda activate pdb2reaction

 # cyipoptをインストール（MEP 探索のDMF法に必要）
 conda install -c conda-forge cyipopt -y
 ```

- *環境モジュール*を使用する HPC クラスターでは、PyTorch をインストールする**前に** CUDA をロードしてください。
 ```bash
 module load cuda/12.9
 ```


## 詳細なインストール手順

環境を段階的に構築する場合:

1. **CUDAをロード（HPCで環境モジュールを使用する場合）**

 ```bash
 module load cuda/12.9
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

 （クラスターが推奨する場合は別の互換バージョンを使用できます。）

5. **`pdb2reaction` 本体と可視化用Chromeをインストール**

 ```bash
 pip install git+https://github.com/t-0hmura/pdb2reaction.git
 plotly_get_chrome -y
 ```

6. **Hugging Face Hub (UMAモデル) にログイン**

 ```bash
 huggingface-cli login
 ```

 参照:

 - <https://github.com/facebookresearch/fairchem>
 - <https://huggingface.co/facebook/UMA>
 - <https://huggingface.co/docs/hub/security-tokens>

7. **インストールの確認**

 ```bash
 pdb2reaction --version
 ```

 インストールされたバージョンが表示されます（例: `{{ version }}`）。
