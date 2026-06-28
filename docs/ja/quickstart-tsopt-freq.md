# クイックスタート: `pdb2reaction all --tsopt`（TS-only モード）

## 目的

既に手元にある TS 候補構造に対して、`tsopt → irc → freq →（任意で dft）` の検証・熱力学チェーンだけを `pdb2reaction all` で一括実行します。上流の `extract` / `path-opt` はスキップされます。

## 事前に必要なもの

- pdb2reaction がインストール済み（[インストール](installation.md)を参照）
- TS 候補構造 1 つ: `.pdb`（残基／電荷情報を持つため推奨）または `.xyz`
- 電荷指定（次のいずれか **ひとつ**を必ず指定）:
  - `-q INT` — 全電荷を整数で直接指定
  - `-l 'RES:Q,...'` — 残基ごとのリガンド電荷（タンパク質–リガンド複合体 PDB の場合）
  - `.gjf` ヘッダから自動取得（Gaussian 入力を渡した場合のみ）
- `-m/--multiplicity` — デフォルトは `1`（一重項）。ラジカル種では明示が必要です
- `.xyz` 入力では `-q` と `-m` が必須です。`-l` 形式で電荷を渡したい場合は `--ref-pdb cluster.pdb` を併用してください
- TS のみモードに入る条件は **3 つすべて成立**: (1) `-i` 入力がちょうど 1 つ、(2) `--scan-lists` が無い、(3) `--tsopt`（または `--tsopt True`）が指定されている。そうでない場合 CLI は入力ゲートで `BadParameter` を送出します（`Provide at least two structures with -i/--input in reaction order, or use a single structure with --scan-lists, or a single structure with --tsopt True.`）。

## 最小コマンド

```bash
pdb2reaction all -i ts_candidate.pdb -l 'SAM:1,GPP:-3' \
    --tsopt --thermo -o ./result_ts_only
```

XYZ 入力では `-q` と `-m` を明示します:

```bash
pdb2reaction all -i ts_candidate.xyz -q -1 -m 1 -b uma \
    --tsopt --thermo -o ./result_ts_only
```

`--tsopt` で検証チェーンを起動し、`--thermo` で freq から ZPE / Gibbs 補正を追加します。両ステージとも同一 backend（デフォルト UMA）で実行されます。

### （任意）DFT 一点計算を追加

```bash
pdb2reaction all -i ts_candidate.pdb -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft --dft-func-basis 'wb97m-v/def2-tzvpd' \
    -o ./result_ts_only
```

> **VRAM 注意:** `--dft` は GPU4PySCF で一点計算を実行します。クラスター > 200 原子 / VRAM < 24 GB では `CUDA out of memory` の可能性あり。OOM 時は `--dft` を外して `pdb2reaction dft` を小さい基底または縮小クラスターで単独実行するか、VRAM の大きい node に移してください。`[dft]` extra のインストールも必要です（[インストール](installation.md) を参照）。

## 期待される出力

成功時の出力ツリー:

```text
result_ts_only/
├── summary.log                                # 人間可読サマリー
├── summary.json                               # status: success | partial | failed
└── segments/
    └── seg_01/                                # TS のみモードの成果物
        ├── structures/                        # 正準 R/TS/P（TS のみモード）
        │   ├── reactant.pdb
        │   ├── ts.pdb
        │   └── product.pdb
        ├── ts/
        │   ├── final_geometry.{xyz,pdb}
        │   └── vib/imag_*_trj.xyz             # 虚振動モードアニメーション
        ├── irc/
        │   └── {forward,backward,finished}_irc_trj.xyz
        ├── freq/{R,TS,P}/
        │   ├── frequencies_cm-1.txt
        │   └── thermoanalysis.yaml
        └── dft/{R,TS,P}/                      # --dft 時のみ
            ├── result.yaml                    # 常に出力（--dft 時）
            └── result.json                    # --out-json で opt-in
```

**確認ポイント（実行順）:**

1. `summary.json` の `status` が `"success"` であること。`segments[0].barrier_kcal` と `segments[0].delta_kcal` も確認。
2. `post_segments[0].ts_imag.n_imag == 1` — 一次鞍点の必要条件。
3. `segments/seg_01/irc/{forward,backward}_irc_trj.xyz` を PyMOL で開き、R 端・P 端まで到達していることを確認。
4. `segments[0].bond_changes` が空でなく、想定どおりの結合切断・形成が記録されていること。
5. `segments/seg_01/freq/{R,TS,P}/frequencies_cm-1.txt` で R / P に虚振動が無く、TS で 1 つだけあること。

**トラブルシュート:**

| 症状 | 原因 | 対処 |
|---|---|---|
| `post_segments[0].ts_imag.n_imag == 0` | TS 候補が極小に落ちてしまう | TS 候補を取り直す（例: `path-search`）。TS のみモードでは失った鞍点を回復できません |
| `n_imag >= 2` | 縮退した負モードあり | `--flatten` で余剰モードを除去。`hessian_dimer.flatten_max_iter` は [tsopt](tsopt.md) を参照 |
| `segments[0].bond_changes` が空（`""` または `(no covalent changes detected)`）、または IRC が想定と違う終点に到達 | 虚振動が反応座標方向と一致していない、または TS が同一井戸を結んでいる | `segments/seg_01/ts/vib/imag_*_trj.xyz` を PyMOL で可視化し、虚振動が想定の反応方向か確認。違う場合は TS 候補を取り直す |
| `freq/{R,P}/frequencies_cm-1.txt` に虚振動が残る | IRC の終点が真の極小に達していない | 収束を厳しくする（`--thresh-post baker`）か YAML で IRC max cycles を伸ばす。[freq](freq.md) を参照 |

## 補足

- `tsopt` 単独でのパラメータ調整（`--opt-mode`、`--max-cycles`、Hessian オプション）は [tsopt](tsopt.md) を参照。
- VRAM に余裕がある場合は `--hessian-calc-mode Analytical` を推奨します（デフォルトは `FiniteDifference`）。
- 全オプションを確認するには `pdb2reaction all --help-advanced`。

## 次のステップ

- 複数構造からの MEP 経路: [クイックスタート: `pdb2reaction all`](quickstart-all.md)
- 単一構造からのスキャン駆動: [クイックスタート: `pdb2reaction all --scan-lists`](quickstart-scan.md)
- 全オプションリファレンス: [all](all.md) / [tsopt](tsopt.md) / [irc](irc.md) / [freq](freq.md) / [dft](dft.md)
