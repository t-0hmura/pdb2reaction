# 再現性と決定性

GPU 上の MLIP 推論は、デフォルトではビット単位で再現しないことがあります。並列リダクション（atomic add、scatter 演算）の累積順序が hardware scheduling に依存するためです。project の UMA smoke benchmark では小さい drift が観測されていますが、その大きさは別 backend/model/GPU/software stack や長い optimizer trajectory に対する保証ではありません。科学的再現性は file identity だけでなく、化学的に意味のある observable で評価してください。

同一 stack での再実行性（例: golden-file 回帰 test）が必要な場合は、backend/model、package version、GPU、input、command を固定して `--deterministic` を使用してください。

## `--deterministic`

`--deterministic` はすべての計算サブコマンド（`opt`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, `all`, `sp`）で受け付けられます。`torch.use_deterministic_algorithms` と `index_reduce_` shim を有効化して PyTorch に strict determinism を要求しますが、任意の custom ASE calculator やすべての third-party custom kernel を制御することはできません。

```bash
pdb2reaction opt -i input.pdb -q 0 --deterministic
pdb2reaction all -i r.pdb p.pdb -q -1 --tsopt --deterministic
```

- これは**プロセス全体に作用します**。`all` で設定すると内部のすべてのステージに伝播するため、ステージごとに渡す必要はありません。
- これは**より低速です**。決定的な scatter/reduce カーネルはデフォルトのものよりスループットが低くなります。厳密な同一 stack 再実行性が必要な場合にのみ使用してください。
- PyTorch が非決定的と認識する演算に決定的 kernel がなければ**明示的に失敗します**。PyTorch の管理外にある外部/custom kernel まで検出する保証ではありません。
- 環境変数 `PDB2REACTION_STRICT_DETERMINISTIC=1` は、CI や直接の Python API（`create_calculator`）に対する同等のエントリポイントです。

### backend ごとの適用範囲

| バックエンド | `--deterministic` |
|---|---|
| `uma` | 同一 stack の end-to-end smoke gate で2実行の出力座標を比較 |
| `orb` / `mace` | PyTorch strict mode は有効になるが、installed third-party backend version ごとに smoke test が必要 |
| `aimnet2` | **非対応 — 拒否されます**（下記参照） |
| `custom` | user 提供 ASE calculator が determinism を所有し、この flag だけでは保証されない |

## 精度と再現性

`--precision` の変更は数値誤差や optimizer trajectory を変え得ますが、fp64 だけで GPU 実行がビット単位で同一になるわけではありません。リダクション順序に起因する非決定性は精度とは独立しています。strict mode が対象とするのは同一 stack の再実行性であり、version や hardware を跨ぐ identity ではありません。

モデル精度と Hessian dtype は独立した設定項目です。Hessian は fp64 が既定ですが、`calc.hessian_double: false` でモデルの native dtype を明示的に選べます。`--precision fp64` を渡すと Hessian も fp64 に強制されるため、optimizer の線形代数がモデルより低い精度で警告なく実行されることはありません。

(ja-precision-by-gpu-class)=
### backend と用途による精度の選択

`--precision` は MLIP 推論の浮動小数点精度（`fp32` | `fp64`、大文字小文字を区別しない）を選択します。バックエンド非依存であり、CLI は値を各バックエンド固有のキー（UMA `precision`、ORB `precision`、MACE `default_dtype`）に振り分けます。`aimnet2` では `fp32` は no-op で、`fp64` はモデル入力が上流で float32 にキャストされるため拒否されます。指定しない場合、UMA は fp32、ORB と MACE は fp64 が既定です（{doc}`バックエンド <backends>` の「精度」節を参照）。GPU class は cost に影響しますが、TS の妥当性基準は変えません。

| 用途 | 推奨 | 理由 |
| --- | --- | --- |
| 通常実行 | 未指定 (`auto`) | UMA/AIMNet2 fp32、ORB/MACE fp64 という tested default を保つ。 |
| 速度優先screening | 必要な場合だけ明示的 `--precision fp32` | ORB/MACE の既定を下げるため、その finite-difference Hessian を最終結果として信頼しない。 |
| 最終 TS/Hessian | ORB/MACE は fp64 を維持し、noise が問題なら UMA fp64 を検討 | precision にかかわらず独立 freq と IRC が必要。 |

OMol で学習された UMA バックエンドでは fp64 が TS 最適化と Hessian に影響し得ます。hardware cost を測定し、本番設定を記録してください。`--deterministic` は別の同一-stack再実行性 control であり、低精度 PES の精度を改善するものではありません。

AIMNet2 はこれらの機能には対応していません:

- **`--precision fp64`** — AIMNet2 のモデル入力は上流で float32 にキャストされるため、「fp64」実行は実際には fp64 になりません。
- **`--deterministic`** — AIMNet2 はカスタム CUDA カーネルを通じて原子間力を計算しますが、これは `torch.use_deterministic_algorithms` の制御外にあるため、原子間力はビット単位で再現可能ではありません（エネルギーは再現可能です）。PyTorch の決定的モードはこのカスタム演算を検出も制御もしないため、この制限は明示的に報告されます。

厳密な同一 stack 再実行性が必要なら、対応 backend に
`--deterministic` を付け、installed stack で2実行の smoke test を通して
ください。UMA には project の end-to-end gate がありますが、ORB/MACE は
installed third-party version ごとの gate が必要です。
