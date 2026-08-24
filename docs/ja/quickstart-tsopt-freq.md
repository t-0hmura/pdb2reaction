# クイックスタート: `pdb2reaction all --tsopt`（TS-only モード）

## 目的

既に手元にある TS 候補構造に対して、`pdb2reaction all --tsopt` で `tsopt → irc` を実行します。`--thermo` を加えると `freq`、`--dft` を加えると DFT 一点計算も実行します。上流の `extract` / `path-opt` はスキップされます。

## 事前に必要なもの

- pdb2reaction がインストール済み（[インストール](installation.md)を参照）
- TS 候補構造 1 つ: PDB/mmCIF（残基／電荷情報を持つため推奨）、XYZ、または GJF
- 電荷指定（次のいずれか **ひとつ**を必ず指定）:
  - `-q INT` — 全電荷を整数で直接指定
  - `-l 'RES:Q,...'` — 残基ごとのリガンド電荷（タンパク質–リガンド複合体 PDB の場合）
  - `.gjf` ヘッダから自動取得（Gaussian 入力を渡した場合のみ）
- `-m/--multiplicity` — デフォルトは `1`（一重項）。ラジカル種では明示が必要です
- `.xyz` 入力では電荷を `-q`、または `--ref-pdb cluster.pdb` と `-l` の組み合わせで解決します。多重度は一重項の `1` がデフォルトで、開殻系では `-m` を明示してください
- TS のみモードに入る条件は **3 つすべて成立**: (1) `-i` 入力がちょうど 1 つ、(2) `--scan-lists` が無い、(3) `--tsopt` が指定されている。そうでない場合 CLI は入力ゲートで `BadParameter` を送出します（`Provide at least two structures with -i/--input in reaction order, or use a single structure with --scan-lists, or a single structure with --tsopt.`）。

## 最小コマンド

```bash
pdb2reaction all -i ts_candidate.pdb -l 'SAM:1,GPP:-3' \
    --tsopt --thermo -o ./result_ts_only
```

XYZ 一重項では `-q` を明示し、多重度は省略できます:

```bash
pdb2reaction all -i ts_candidate.xyz -q -1 -b uma \
    --tsopt --thermo -o ./result_ts_only
```

`--tsopt` で検証チェーンを起動し、`--thermo` で freq から ZPE / Gibbs 補正を追加します。両ステージとも同一バックエンド（デフォルト UMA）で実行されます。

### （任意）DFT 一点計算を追加

```bash
pdb2reaction all -i ts_candidate.pdb -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft --dft-func-basis 'wb97m-v/def2-tzvpd' \
    -o ./result_ts_only
```

> **VRAM 注意:** `--dft` は GPU4PySCF で一点計算を実行します。必要メモリは
> 構造・基底・汎関数・精度・software stackに依存するため、対象nodeで代表構造を
> pilot実行し、peak memoryを測定してください。OOM時は `--dft` を外し、より小さい
> 基底または縮小clusterで `pdb2reaction dft` を単独実行するか、より大きいnodeへ
> 移してください。`[dft]` extraのinstallも必要です
> （[インストール](installation.md) を参照）。

## 期待される出力

成功時の出力ツリー:

```text
result_ts_only/
├── summary.log                                # 人間可読サマリー
├── summary.json                               # status: success | partial | failed
└── segments/
    └── seg_01/                                # TS のみモードの成果物
        ├── reactant.pdb                       # 正準 R/TS/P は seg_01/ 直下（TS のみモード）
        ├── ts.pdb
        ├── product.pdb
        ├── ts/
        │   ├── final_geometry.{xyz,pdb}
        │   └── vib/imag_*_trj.xyz             # 虚振動モード軌跡
        ├── irc/
        │   └── {forward,backward,finished}_irc_trj.xyz
        ├── freq/{R,TS,P}/
        │   ├── frequencies_cm-1.txt
        │   └── thermoanalysis.yaml
        └── dft/{R,TS,P}/                      # --dft 時のみ
            └── result.yaml                    # 常に出力（--dft 時）
```

**確認ポイント（実行順）:**

1. `summary.json` の `scientific_status` が `"success"` であること。`"partial"` または `"failed"` なら `scientific_status_reasons` を確認します。互換用の `status` は科学的な可否判定には使いません。
2. `post_segments[0].ts_imag.n_imag == 1` — 一次鞍点の必要条件です。振動数の大きさだけで反応性を判定せず、modeを可視化して IRC 接続を確認します。系ごとのnoise評価で柔らかい非反応modeを特定した場合だけ、YAMLの opt-in filter `irc.imag_below` をデフォルトの `0.0` より負側へ設定します。IRC が受理するのは `ν <= imag_below` のmodeです。
3. `segments/seg_01/irc/{forward,backward}_irc_trj.xyz` を PyMOL で開き、R 端・P 端まで到達していることを確認。
4. `segments[0].bond_changes` が空でなく、想定どおりの結合切断・形成が記録されていること。
5. `segments/seg_01/freq/{R,TS,P}/frequencies_cm-1.txt` を診断用に確認します。一次鞍点の条件として TS は虚振動がちょうど 1 つ必要です。R/P の虚振動は熱化学解析を妨げません。

**トラブルシュート:**

| 症状 | 原因 | 対処 |
|---|---|---|
| `post_segments[0].ts_imag.n_imag == 0` | TS 候補が極小に落ちてしまう | 経路情報のない通常の TS-only mode は目的の隣接鞍点を特定できず、自動 saddle recovery の default budget も 0 です。endpoint がある場合は `path-search` で TS 候補を取り直します |
| `n_imag >= 2` | 高次鞍点候補または TS 未収束 | TS の認定には虚振動がちょうど 1 つ必要です。周波数が小さいことだけを理由に余分なモードを除外せず、各虚振動モードの変位を確認してください。必要に応じて `--thresh-post` を厳しくし、`--flatten` と再最適化で余分なモードを除去してから認定します。 |
| `segments[0].bond_changes` が空（`""` または `(no covalent changes detected)`）、または IRC が想定と違う終点に到達 | 虚振動が反応座標方向と一致していない、または TS が同じ井戸同士を結んでいる（反応物側と生成物側が同一極小） | `segments/seg_01/ts/vib/imag_*_trj.xyz` を PyMOL で可視化し、虚振動が想定の反応方向か確認。違う場合は TS 候補を取り直す |

## 補足

- `tsopt` 単独でのパラメータ調整（`--opt-mode`、`--max-cycles`、Hessian オプション）は [tsopt](tsopt.md) を参照。
- デフォルトの `FiniteDifference` を基準にし、選択した backend/model と対象系で速度・メモリ・結果を検証できた場合にだけ `--hessian-calc-mode Analytical` を明示してください。
- 全オプションを確認するには `pdb2reaction all --help-advanced`。

## 次のステップ

- 複数構造からの MEP 経路: [クイックスタート: `pdb2reaction all`](quickstart-all.md)
- 単一構造からのスキャン駆動: [クイックスタート: `pdb2reaction all --scan-lists`](quickstart-scan.md)
- 全オプションリファレンス: [all](all.md) / [tsopt](tsopt.md) / [irc](irc.md) / [freq](freq.md) / [dft](dft.md)
