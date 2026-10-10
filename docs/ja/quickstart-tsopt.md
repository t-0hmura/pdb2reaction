# クイックスタート: `pdb2reaction all --tsopt`（TS-only モード）

TS-only モードは、手元の遷移状態（TS）の候補 1 つを、最小エネルギー経路（MEP）の探索を省いて確かめます。`pdb2reaction all --tsopt` は TS を最適化し、固有反応座標（IRC）を両方向に追跡して、両端の構造（反応物 R と生成物 P）を最適化します。`--thermo` で振動解析と熱化学補正を、`--dft` で R・TS・P の DFT 一点計算を追加できます。TS の候補があるなら、このモードに渡して、直に TS の探索を始められます。

---

## 主な用途

* **スキャンや MEP の候補を TS に詰める**: スキャンの最高点や、MEP の最高エネルギーのイメージ（HEI）を TS まで最適化
* **別の方法で作った候補を確かめる**: 他のプログラムで得た構造や手で組んだ構造が、狙った R と P を結ぶ TS（n_imag = 1）かを確認
* **DFT に渡す前に MLIP の TS を確かめる**: 機械学習原子間ポテンシャル（MLIP）で得た TS を、[DFT](dft-backend.md) で詰める前に確認

## 最小コマンド

TS 候補を 1 つ渡し、`--tsopt` を付けます。同梱例には TS 候補が無いので、次のコマンドは [`all` のクイックスタート](quickstart-all.md) の実行で得た HEI を使います。自分の反応では、自分で用意した候補を渡してください。

```bash
pdb2reaction all -i result_all/_work/path_opt/hei_seg_01.pdb -l 'SAM:1,GPP:-3' \
    --tsopt --thermo -o ./result_ts_only
```

端末の最後のほうの `====== Pipeline summary ======` の下に `Scientific status: success` と出れば成功で、`summary.json` の `scientific_status` にも同じ値が入ります。

他のプログラムで得た候補などの XYZ では、モデルの総電荷を `-q` で明示します。`all` のクイックスタートのモデルは 0 です。

```bash
pdb2reaction all -i ts_candidate.xyz -q 0 \
    --tsopt --thermo -o ./result_ts_only
```

### （任意）DFT 一点計算を追加

`--dft` で R・TS・P の DFT 一点計算を追加し、`--func-basis` で汎関数と基底を指定します。デフォルトは `wb97m-v/def2-svp` です。

```bash
pdb2reaction all -i result_all/_work/path_opt/hei_seg_01.pdb -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft --func-basis 'wb97m-v/def2-tzvpd' \
    -o ./result_ts_only
```

## 実行の前に

* **入力**: PDB/mmCIF・XYZ・GJF の TS 候補 1 つ。クラスターモデルの切り出しは `-c` を指定したときだけ行い、省略すると構造をそのまま使います。
* **電荷と多重度**: {ref}`電荷 <ja-charge-specification>`は `-q`、`-l`（PDB/mmCIF のみ）、YAML の `calc.charge`、GJF ヘッダーのいずれかで与えます。多重度は `-m` → YAML の `calc.spin` → GJF ヘッダー → 1 の順で決まります。
* **TS-only モードになる条件**: 入力が 1 つで、`--tsopt` を付け、`--scan-lists` を付けないとき。入力が 2 つ以上なら [MEP 探索](quickstart-all.md)、入力 1 つに `--scan-lists` を付けると[スキャン](quickstart-scan.md)になります。

## 期待される出力

成功すると、次のファイルが書き出されます。

```text
result_ts_only/
├── summary.log                     # 実行の要約
├── summary.json                    # 結果（scientific_status を含む）
└── segments/
    └── seg_01/
        ├── reactant.pdb            # R/TS/P の構造（XYZ 入力では .xyz、GJF 入力では .gjf）
        ├── ts.pdb
        ├── product.pdb
        ├── energy_diagram_MLIP.png # R–TS–P のエネルギー図（--thermo で energy_diagram_G_MLIP.png も）
        ├── ts/
        │   ├── final_geometry.{xyz,pdb}
        │   └── vib/imag_*_trj.xyz  # 虚振動モードごとのアニメーション
        ├── irc/
        │   └── {forward,backward,finished}_irc_trj.xyz
        ├── freq/{R,TS,P}/          # --thermo のとき
        │   ├── frequencies_cm-1.txt
        │   └── thermoanalysis.yaml
        └── dft/{R,TS,P}/           # --dft のとき
            └── result.yaml
```

## 結果の確認

1. **完了状況**: `scientific_status` には、指定した段がすべて収束すると `success`、そうでなければ `partial` か `failed` が入り、[理由](json-output.md#実行の完了と指定した段の完了)は `scientific_status_reasons` に出ます。虚振動のモードが、できる結合と切れる結合を動かしているか、端点が狙った R と P か、の 2 点は自分で確かめてください。
2. **TS のモード**: TS 最適化が成功すると、反応モードの虚振動が 1 つ出ます。このとき端末に `[tsopt] Converged (n_imag=1).` と出て、`summary.json` の `post_segments[0].tsopt.n_imaginary_modes` に本数が記録されます。`segments/seg_01/ts/vib/imag_*_trj.xyz` をビューアで開き、できる結合と切れる結合に沿って原子が動くかを確認してください。
3. **端点**: `segments/seg_01/irc/finished_irc_trj.xyz` と R/TS/P の構造 `reactant.pdb`・`ts.pdb`・`product.pdb` を開き、`segments[0].bond_changes` を読みます。端点は狙った R と P のはずです。IRC が収束しなくても、端点の最適化で狙った R と P に着けば、その結果は使えます。
4. **端点の振動数**: `--thermo` のとき、`segments/seg_01/freq/{R,TS,P}/frequencies_cm-1.txt` に符号つきの全振動数が出ます。R と P には −5.00 cm⁻¹ より小さい値（虚振動）が無いはずです。

| 結果 | 次に試すこと |
|---|---|
| n_imag = 0 | MEP の HEI やスキャンの最高点など、よりよい候補から始めます。TS-only モードには手がかりにする経路がありません。 |
| n_imag ≥ 2 | 各虚振動のモードを見ます。`--flatten` を付けて最適化し直すか、`all --thresh-post gau_tight`（デフォルトの [`baker`](tsopt.md#処理の仕組みと計算仕様) より厳しい）または `tsopt --thresh gau_tight` で収束を厳しくします。 |
| `bond_changes` が空、または端点が狙いと違う | TS のモードと IRC を確認します。経路が別の極小点どうしを結んでいる可能性があります。 |
| R や P に虚振動が残る | 端点の構造とモードを確認します。`--thresh-post gau_tight` で端点の最適化を厳しくするか、`--irc-max-cycles`（デフォルト 125）で IRC を延ばします。 |

それでも TS が取れないときは、{ref}`TS を確かめる <ja-mechanism-check-ts>` と {ref}`TS が取れないとき <ja-ts-search-fails>` を参照してください。

## 使用上の注意点

* **XYZ・GJF のキャップ水素**: 同じ原子の PDB を `--ref-pdb` で渡したときだけ、キャップ水素の親原子が固定されます。{ref}`原子の固定と距離の拘束 <ja-freeze-atoms-and-restraints>` を参照してください。
* **Hessian の計算法**: デフォルトの `--hessian-calc-mode FiniteDifference` のまま使ってください。`--hessian-calc-mode Analytical` は、使うバックエンド・モデルと対象の系で速度・メモリ・結果を確かめてから指定します。
* **エネルギー**: `post_segments[0].mlip.barrier_kcal` が ΔE‡（TS − R）、`.delta_kcal` が ΔE（P − R）で、単位は kcal/mol です。どちらも最適化した TS と端点から求めます。MEP が無いので、`segments[0].barrier_kcal` と `.delta_kcal` にも同じ値が入ります。`--thermo` のときは `post_segments[0].gibbs_mlip.barrier_kcal` と `.delta_kcal` が ΔG‡ と ΔG、`--dft` のときは `post_segments[0].dft.barrier_kcal` と `.delta_kcal` が DFT の値です。
* **R と P の名付け**: MEP が無いと反応の向きが分からないため、TS-only モードはエネルギーが高いほうの IRC 端点を R、低いほうを P と呼び、この決まりを `summary.json` の `endpoint_assignment` に記録します。この名前は化学的な反応の向きではありません。P からの障壁は `barrier_kcal − delta_kcal` です。
* **R・P の虚振動**: R や P に虚振動が残っても熱化学は計算され、そのモードは ZPE と G から除かれます。極小点として扱う前にモードを確認してください。
* **IRC に進む条件**: `all` が IRC に進むのは、TS 最適化が収束し、最後の Hessian の計算が終わり、n_imag が 1 以上のときだけです。n_imag が 2 以上のときの IRC は、虚振動の 1 つに沿った診断用の計算で、構造が一次の鞍点になるわけではありません。
* **`tsopt` の細かい設定**: `--opt-mode`、`--max-cycles`、Hessian のオプションを変えるときは、[`tsopt`](tsopt.md) を単独で実行してください。`all` の全オプションは `pdb2reaction all --help-advanced` で確認できます。
* **DFT**: TS そのものを DFT で最適化する方法と、DFT 用の追加パッケージのインストールや GPU メモリについては [MLIP の TS を DFT で確かめる](dft-backend.md) を参照してください。

## 次のステップ

- [MLIP の TS を DFT で確かめる](dft-backend.md): TS を DFT で詰めて確かめる
- [反応機構を調べるコツ](mechanism-tips.md): TS の確かめ方と、TS が取れないときに試すこと
- [`tsopt`](tsopt.md)・[`irc`](irc.md)・[`freq`](freq.md): 各段を単独で実行する
- [クイックスタート: `pdb2reaction all`](quickstart-all.md): R と P から MEP を作る
- [クイックスタート: scan](quickstart-scan.md): 1 つの構造から経路を作る
- [`all`](all.md)・[`dft`](dft.md): オプションの説明
- [トラブルシューティング](troubleshooting.md): エラーメッセージや症状から対処を探す
