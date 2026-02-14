# ドキュメント校正・再構成 実行計画（v7 / 2026-02-14）

この計画は、`pdb2reaction` のドキュメント全体（EN/JA）を対象に、前後関係の矛盾をなくし、導線を単純化し、配布時に再現できる形へ整えるための実行手順です。

対象は `README.md`, `docs/*.md`, `docs/ja/*.md`, `test/test.md`, `docs/_static/**`, `docs/_templates/**` とし、`docs/_build/**` は生成物として精読・検索対象から除外します。

---

## 1. 目標（Definition of Done）

1. README から各ガイドへ一貫して辿れる。
2. TS 記述が「候補」と「検証済み」を混同しない。
3. EN/JA の章構成と注意点が対応している。
4. チュートリアルが 1 ページ構成で再現可能。
5. ドキュメント構造が簡素で、不要な階層がない。
6. リンク切れ・旧名称・内部ID漏れが 0 件。

---

## 2. ディレクトリ簡素化方針（固定ルール）

最終的な構造は次を維持する。

```text
docs/
├─ index.md
├─ getting-started.md
├─ tutorial.md
├─ concepts.md
├─ troubleshooting.md
├─ (subcommand docs: all/extract/opt/tsopt/path_search/... )
├─ _static/
│  └─ bezA_case_study/
├─ _templates/
└─ ja/
   ├─ index.md
   ├─ getting-started.md
   ├─ tutorial.md
   ├─ concepts.md
   ├─ troubleshooting.md
   └─ (subcommand docs)
```

禁止:
- `docs/tutorials/`, `docs/ja/tutorials/`
- `docs/_static/tutorials/`
- `tutorial-*` の分割ページ再導入

---

## 3. 用語・表記ルール（全ページ強制）

- `first-pass`: 1コマンド結果は初期案
- `TS candidate`: 未検証のTS候補
- `validated TS`: `freq` + `irc` で妥当性確認済み

表記統一:
- CLI 名: `path-search`, `path-opt`
- 出力ディレクトリ: `path_search/`, `path_opt/`
- 内部ID（例: `test3`）はユーザー向け本文で禁止

禁止表現:
- 「1コマンドでTS確定」の含意
- 「収束すれば必ず虚数振動数1本」の断定

---

## 4. ページ役割の再定義

- `README.md`: 価値提案、最小実行例、期待値コントロール、主要導線
- `docs/index.md`, `docs/ja/index.md`: 目標別ナビ + サブコマンド索引
- `docs/getting-started.md`, `docs/ja/getting-started.md`: 導入、環境、実行モード入口
- `docs/tutorial.md`, `docs/ja/tutorial.md`: ケーススタディ + smoke test マトリクス（1ファイル）
- `docs/concepts.md`, `docs/ja/concepts.md`: 用語と段階の全体像
- `docs/troubleshooting.md`, `docs/ja/troubleshooting.md`: エラー対処の実用導線

---

## 5. 実行フェーズ（順次実施）

### Phase 0: ベースライン固定
- [ ] 検索対象を `docs/_build/**` 除外で統一
- [ ] 現在のファイル構造を記録（`docs`, `docs/ja`, `docs/_static`）

### Phase 1: 入口導線の整合
- [ ] `README.md` を短く保ち、詳細を `getting-started` へ委譲
- [ ] `index/getting-started/tutorial/troubleshooting` の EN/JA で導線と注意文を一致させる
- [ ] tutorial 参照を `tutorial.md` へ一本化する

### Phase 2: TS 文脈の統一
- [ ] `all`, `path_search`, `path_opt`, `tsopt` に TS candidate 定義を明記
- [ ] `freq`, `irc` に validated TS の検証条件を明記
- [ ] TS に関する See Also を双方向で整える

### Phase 3: サブコマンドページ精読
- [ ] EN/JA ペア単位で文脈矛盾を解消
- [ ] CLI 例の再現性（入力条件、出力先）を明示
- [ ] 長文段落を分割し、判断基準は箇条書き化

### Phase 4: チュートリアル配布性の保証
- [ ] `docs/_static/bezA_case_study/r_input.pdb` を再現入力として明示
- [ ] 同梱ログのパスを相対化（環境依存情報を除去）
- [ ] `docs/_static/...`（配布データ）と実行出力 `./bezA_case_study/` を明確に区別

### Phase 5: 機械検証ゲート
- [ ] Markdown リンク切れ 0
- [ ] 旧階層/旧ファイル参照 0
- [ ] 内部ID（`test3` 等）の本文残存 0
- [ ] Sphinx ビルド成功（可能環境のみ）

### Phase 6: 最終レビュー
- [ ] EN/JA ペアの差分をレビュー
- [ ] 変更サマリを `improve_plan.md` 実行ログへ追記

---

## 6. 精読対象（全ページ）

入口:
- `README.md`
- `docs/index.md`, `docs/ja/index.md`
- `docs/getting-started.md`, `docs/ja/getting-started.md`
- `docs/tutorial.md`, `docs/ja/tutorial.md`
- `docs/concepts.md`, `docs/ja/concepts.md`
- `docs/troubleshooting.md`, `docs/ja/troubleshooting.md`

メイン/TS/解析:
- `docs/all.md`, `docs/ja/all.md`
- `docs/path_search.md`, `docs/ja/path_search.md`
- `docs/path_opt.md`, `docs/ja/path_opt.md`
- `docs/tsopt.md`, `docs/ja/tsopt.md`
- `docs/freq.md`, `docs/ja/freq.md`
- `docs/irc.md`, `docs/ja/irc.md`

前処理/スキャン/補助:
- `docs/extract.md`, `docs/ja/extract.md`
- `docs/opt.md`, `docs/ja/opt.md`
- `docs/scan.md`, `docs/ja/scan.md`
- `docs/scan2d.md`, `docs/ja/scan2d.md`
- `docs/scan3d.md`, `docs/ja/scan3d.md`
- `docs/dft.md`, `docs/ja/dft.md`
- `docs/trj2fig.md`, `docs/ja/trj2fig.md`
- `docs/energy-diagram.md`, `docs/ja/energy-diagram.md`
- `docs/yaml-reference.md`, `docs/ja/yaml-reference.md`
- `docs/uma_pysis.md`, `docs/ja/uma_pysis.md`
- `docs/cli-conventions.md`, `docs/ja/cli-conventions.md`
- `docs/glossary.md`, `docs/ja/glossary.md`
- `docs/add_elem_info.md`, `docs/ja/add_elem_info.md`
- `docs/fix_altloc.md`, `docs/ja/fix_altloc.md`
- `test/test.md`

---

## 7. 検証コマンド（実行用）

```bash
# 1) 旧構造参照（生成物除外）
rg -n --glob '!docs/_build/**' \
  'docs/tutorials/|docs/ja/tutorials/|docs/_static/tutorials/|tutorial-case-study|tutorial-smoke-tests' \
  README.md docs test

# 2) 内部ID/旧命名の残存
rg -n --glob '!docs/_build/**' \
  'test3|bezA-test3|tutorial-bezA' \
  README.md docs test/test.md

# 3) Markdownリンク検証（docs + README）
perl <<'PERL'
use strict; use warnings;
use File::Basename qw(dirname); use File::Spec; use File::Find;
my $fail = 0;
sub check_file {
  my ($file, $content, $base) = @_;
  while ($content =~ /\]\(([^)]+)\)/g) {
    my $link = $1; $link =~ s/\s+$//;
    my ($target) = split(/\s+/, $link, 2);
    next if !defined($target) || $target eq '';
    next if $target =~ /^#/;
    next if $target =~ m{^(?:https?://|mailto:)};
    $target =~ s/[#?].*$//;
    next if $target eq '';
    my $path = File::Spec->catfile($base, $target);
    if (!-e $path) { print "BROKEN: $file -> $target\n"; $fail = 1; }
  }
}
find({ no_chdir => 1, wanted => sub {
  return if -d $_;
  return if $_ !~ /\.md\z/;
  return if $File::Find::name =~ m{\Adocs/_build/};
  open my $fh, '<', $_ or do { print "BROKEN: $_ -> (unreadable)\n"; $fail = 1; return; };
  local $/; my $txt = <$fh>; close $fh;
  check_file($_, $txt, dirname($_));
}}, 'docs');
if (-e 'README.md') {
  open my $rfh, '<', 'README.md' or die "Failed to open README.md: $!\n";
  local $/; my $r = <$rfh>; close $rfh;
  check_file('README.md', $r, '.');
}
exit($fail);
PERL

# 4) HTMLビルド（環境が整っている場合）
make -C docs html
```

---

## 8. 実行ログ（更新用）

- 2026-02-14: `improve_plan.md` を v7 として全面改訂（1ファイルtutorial前提、簡素化方針、全ページ精読フロー、機械検証ゲートを再定義）。
