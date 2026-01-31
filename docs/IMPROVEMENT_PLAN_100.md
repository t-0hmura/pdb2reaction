# pdb2reaction Documentation - 100点達成プラン

## 現状評価（改善後）

| カテゴリ | 現在スコア | 目標 | ギャップ |
|---------|-----------|------|---------|
| 正確性 (Accuracy) | 94 | 100 | -6 |
| わかりやすさ (Clarity) | 94 | 100 | -6 |
| 構造 (Structure) | 93 | 100 | -7 |
| レイアウト (Layout) | 92 | 100 | -8 |
| **総合** | **93** | **100** | **-7** |

---

## フェーズ1: 正確性の改善 (94→100)

### 1.1 README.md のリンク修正（重大）
**問題**: README.md のドキュメントリンクが `docs/` を指しているが、実際は `docs/`
```markdown
# 現状（誤）
[docs/getting-started.md](docs/getting-started.md)

# 修正後
[docs/getting-started.md](docs/getting-started.md)
```

**対象ファイル**: `README.md` (全リンク)

### 1.2 オプション名の統一
**問題**: `--mult` vs `--multiplicity` の表記が混在

| ファイル | 現状 | 修正 |
|---------|------|------|
| `all.md` | `-m, --mult` | そのまま（all固有） |
| 他サブコマンド | `-m, --multiplicity` | 明確な注記追加 |

**追加修正**: `getting-started.md` と `index.md` に以下の注記追加：
```markdown
> **Note:** The `all` command uses `--mult` while other subcommands use `--multiplicity`.
```

### 1.3 YAML の重複キー修正
**問題**: `calc` セクションに `freeze_atoms` が2回出現

**対象ファイル**: `yaml-reference.md`, `path_search.md`, `tsopt.md`, `freq.md`, `irc.md`, `scan.md`, `opt.md`, `path_opt.md`

```yaml
# 現状（重複）
calc:
  freeze_atoms: null         # calculator-level frozen atoms  ← 削除
  ...
  hessian_calc_mode: ...

# 修正後（コメントで説明）
calc:
  ...
  # Note: freeze_atoms is set via geom.freeze_atoms and propagated to calculator
```

### 1.4 デフォルト値の整合性チェック
**確認項目**:
- `--thresh` のデフォルト値（`gau` vs `baker`）
- `--opt-mode` のデフォルト値（`light` vs `heavy`）
- `--hessian-calc-mode` のデフォルト値（`FiniteDifference`）

**対象**: 全サブコマンドのCLIオプション表

### 1.5 バージョン・更新日の追加
各ファイルのフッターに:
```markdown
---
*Last updated: January 2026*
```

---

## フェーズ2: わかりやすさの改善 (94→100)

### 2.1 クイックサマリーボックスの追加
各サブコマンドドキュメントの先頭に「TL;DR」セクション追加

**例** (`tsopt.md`):
```markdown
## Overview

> **TL;DR:** `tsopt` finds and optimizes transition states using either Dimer (`--opt-mode light`) or RS-I-RFO (`--opt-mode heavy`) methods. Use `--hessian-calc-mode Analytical` for best performance when VRAM permits.
```

**対象**: `tsopt.md`, `path_search.md`, `freq.md`, `irc.md`, `scan.md`, `dft.md`, `opt.md`, `path_opt.md`

### 2.2 ワークフロー図の追加
`concepts.md` と `all.md` にMermaid形式のフローチャート追加

```markdown
```{mermaid}
graph TD
    A[Input PDBs] --> B{Has -c/--center?}
    B -->|Yes| C[extract: Pocket Extraction]
    B -->|No| D[Use Full Structure]
    C --> E{Multiple Inputs?}
    D --> E
    E -->|≥2| F[path-search: MEP Search]
    E -->|1 + scan-lists| G[scan: Staged Scan]
    E -->|1 + tsopt| H[tsopt: TS Optimization]
    G --> F
    F --> I{--tsopt True?}
    I -->|Yes| J[tsopt + IRC]
    I -->|No| K[Done]
    J --> L{--thermo True?}
    L -->|Yes| M[freq: Vibrational Analysis]
    M --> N{--dft True?}
    N -->|Yes| O[dft: Single-Point DFT]
    O --> K
    H --> L
```
```

### 2.3 用語集への参照追加
技術用語の初出時にリンク追加

**例**:
```markdown
# 現状
The command runs MEP search using GSM.

# 修正後
The command runs [MEP](glossary.md#mep) search using [GSM](glossary.md#gsm).
```

**対象**: `getting-started.md`, `concepts.md`, `all.md`

### 2.4 YAMLコメントの充実
最小限のYAML例にも目的を説明するコメント追加

```yaml
# Minimal configuration for enzyme reaction pathway
calc:
  model: uma-s-1p1               # Default UMA model (S variant)
  hessian_calc_mode: Analytical  # Faster but requires more VRAM
gs:
  max_nodes: 12                  # More nodes = smoother path
  climb: true                    # Enable climbing image for TS location
```

### 2.5 サブコマンド間の関係明確化
`index.md` にサブコマンドの依存関係テーブル追加

```markdown
### Subcommand Dependencies

| Subcommand | Typical Input | Typical Output | Often Followed By |
|------------|---------------|----------------|-------------------|
| `extract` | Full PDB | Pocket PDB | `path-search`, `scan` |
| `path-search` | Pocket PDBs | MEP trajectory | `tsopt` |
| `tsopt` | TS candidate | Optimized TS | `irc`, `freq` |
| `irc` | Optimized TS | IRC trajectory | `freq` |
| `freq` | Any structure | Frequencies, thermo | `dft` |
```

---

## フェーズ3: 構造の改善 (93→100)

### 3.1 セクション順序の標準化
全サブコマンドドキュメントで以下の順序に統一:

1. `## Overview` (TL;DR含む)
2. `## Usage`
3. `## Examples`
4. `## Workflow`
5. `## CLI Options`
6. `## Outputs`
7. `## Notes`
8. `## YAML Configuration`
9. `## See Also`

**対象**: 全サブコマンドドキュメント（現状でセクション順が異なるもの）

### 3.2 "See Also" セクションの追加
全サブコマンドドキュメントの末尾に関連リンク追加

**例** (`tsopt.md`):
```markdown
## See Also

- [path-search](path_search.md) — MEP search that identifies TS candidates
- [irc](irc.md) — IRC from optimized TS
- [freq](freq.md) — Vibrational analysis to confirm single imaginary frequency
- [YAML Reference](yaml-reference.md) — Complete YAML options for `hessian_dimer` and `rsirfo`
- [Glossary](glossary.md) — Definitions of TS, Dimer, RS-I-RFO
```

### 3.3 CLIオプション表の論理グループ化
長いオプション表をカテゴリ別に分割

**例** (`all.md`):
```markdown
### Input/Output Options
| Option | Description | Default |
|--------|-------------|---------|
| `-i, --input` | ... | ... |
| `--out-dir` | ... | ... |

### Extraction Options
| Option | Description | Default |
|--------|-------------|---------|
| `-c, --center` | ... | ... |
| `--radius` | ... | ... |

### MEP Search Options
| Option | Description | Default |
|--------|-------------|---------|
| `--max-nodes` | ... | ... |
| `--mep-mode` | ... | ... |

### Post-Processing Options
| Option | Description | Default |
|--------|-------------|---------|
| `--tsopt` | ... | ... |
| `--thermo` | ... | ... |
| `--dft` | ... | ... |
```

### 3.4 見出し階層の修正
一部ファイルで見出しレベルがスキップされている問題を修正

```markdown
# 問題（h2からh4へスキップ）
## Section
#### Subsection

# 修正後
## Section
### Subsection
```

---

## フェーズ4: レイアウトの改善 (92→100)

### 4.1 Admonition の統一
全ファイルで以下の形式に統一:

```markdown
# Sphinx/MyST形式（推奨）
```{important}
Text here.
```

```{tip}
Text here.
```

```{warning}
Text here.
```

```{note}
Text here.
```
```

**対象**: 現在 `> **Note:**` 形式を使用しているファイル

### 4.2 大きなテーブルの分割
55行を超えるCLIオプション表を論理グループに分割（3.3と連携）

### 4.3 水平線による視覚的区切りの統一
主要セクション間に `---` を挿入

```markdown
## Overview
...

---

## Usage
...

---

## Examples
```

### 4.4 コードブロックのコンテキスト追加
コードブロックの前に1行の説明を追加

```markdown
# 現状（説明なし）
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP'
```

# 修正後
Run a basic MEP search between reactant and product:
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP'
```
```

### 4.5 出力構造図のフォーマット改善
`tree` 形式の出力に説明コメント追加

```text
out_dir/
├─ summary.log          # Quick inspection summary
├─ summary.yaml         # Machine-readable summary
├─ pockets/             # Extracted cluster models (one per input)
│   ├─ pocket_R.pdb
│   └─ pocket_P.pdb
└─ path_search/         # MEP results
    ├─ mep.trj          # Complete MEP trajectory
    ├─ mep.pdb          # PDB version (when template available)
    └─ seg_000_*/       # Per-segment details
```

---

## 日本語版の同期

上記の全修正を `docs/ja/` 内の対応ファイルにも適用する。

| 英語ファイル | 日本語ファイル |
|-------------|---------------|
| `all.md` | `ja/all.md` |
| `getting-started.md` | `ja/getting-started.md` |
| `concepts.md` | `ja/concepts.md` |
| `troubleshooting.md` | `ja/troubleshooting.md` |
| `glossary.md` | `ja/glossary.md` |
| `yaml-reference.md` | `ja/yaml-reference.md` |
| ... | ... |

---

## 実装優先順位

### 高優先度（即座に改善効果あり）
1. README.md のリンク修正 (Phase 1.1)
2. TL;DR サマリーの追加 (Phase 2.1)
3. See Also セクションの追加 (Phase 3.2)
4. Admonition の統一 (Phase 4.1)

### 中優先度
5. オプション表の論理グループ化 (Phase 3.3)
6. YAML 重複キーの修正 (Phase 1.3)
7. 用語集参照の追加 (Phase 2.3)
8. セクション順序の標準化 (Phase 3.1)

### 低優先度（細部の磨き上げ）
9. ワークフロー図の追加 (Phase 2.2)
10. コードブロックのコンテキスト追加 (Phase 4.4)
11. バージョン・更新日の追加 (Phase 1.5)

---

## 予想スコア改善

| カテゴリ | 現在 | Phase 1 後 | Phase 2 後 | Phase 3 後 | Phase 4 後 |
|---------|-----|-----------|-----------|-----------|-----------|
| 正確性 | 94 | 98 | 99 | 99 | 100 |
| わかりやすさ | 94 | 95 | 99 | 99 | 100 |
| 構造 | 93 | 94 | 96 | 100 | 100 |
| レイアウト | 92 | 93 | 95 | 97 | 100 |
| **総合** | **93** | **95** | **97** | **99** | **100** |

---

## 作業量見積もり

| フェーズ | 対象ファイル数 | 主な変更内容 |
|---------|---------------|-------------|
| Phase 1 | ~15 | 正確性修正、リンク修正 |
| Phase 2 | ~12 | TL;DR、図、参照追加 |
| Phase 3 | ~14 | セクション再構成 |
| Phase 4 | ~14 | フォーマット統一 |
| 日本語同期 | ~14 | 上記の翻訳・適用 |

---

*計画作成日: 2026年1月31日*
