# Documentation Improvement Plan
## Target: Score 90+ in Accuracy, Clarity, Structure, and Layout

---

## Current Assessment (docs)

| Category | Current Score | Target | Gap |
|----------|--------------|--------|-----|
| Accuracy (正確性) | 90/100 | 92+ | +2 |
| Clarity (わかりやすさ) | 86/100 | 92+ | +6 |
| Structure (構造) | 85/100 | 92+ | +7 |
| Layout (レイアウト) | 88/100 | 92+ | +4 |
| **Overall** | **88/100** | **92+** | **+4** |

---

## Phase 1: Critical Improvements (Highest Impact)

### 1.1 Remove YAML Duplication from all.md
**Problem:** all.md contains ~140 lines of redundant YAML configuration (lines 197-333) that duplicates yaml-reference.md.

**Solution:**
- Replace with concise summary (5-10 lines) + link to yaml-reference.md
- Keep only a minimal illustrative example

**Impact:** Structure +3, Clarity +2

### 1.2 Centralize Boolean Option Guidance
**Problem:** Boolean option warnings (`True`/`False` not `true`/`false`) are scattered across multiple files.

**Solution:**
- Add a "CLI Conventions" section to getting-started.md or concepts.md
- Reference this section from other docs instead of repeating

**Impact:** Clarity +2, Structure +1

### 1.3 Simplify `--ligand-charge` Explanation
**Problem:** The charge precedence logic in all.md (lines 75-76) is dense and hard to parse.

**Solution:**
- Use a decision table or flowchart-style explanation
- Add concrete examples for each scenario

**Impact:** Clarity +3

---

## Phase 2: Structural Enhancements

### 2.1 Add Quick Reference Card to getting-started.md
**Solution:**
- Add a "Quick Reference" section with common command patterns
- Include a concise cheat sheet for first-time users

**Impact:** Clarity +2

### 2.2 Improve concepts.md Navigation
**Solution:**
- Add ASCII workflow diagram (already exists, but enhance)
- Add cross-links to each subcommand's documentation

**Impact:** Structure +2

### 2.3 Standardize Document Structure
**Solution:** Ensure all subcommand docs follow this template:
1. Overview (1-2 sentences)
2. Usage (basic syntax)
3. Examples (2-3 realistic examples)
4. Workflow (numbered steps)
5. CLI Options (table)
6. Outputs (directory structure)
7. Notes (tips and caveats)

**Impact:** Structure +2, Clarity +1

---

## Phase 3: Layout and Polish

### 3.1 Consistent Admonition Usage
**Solution:**
- Use `> **Note:**` for tips
- Use `> **Warning:**` for critical caveats
- Use `> **Important:**` for must-know information

### 3.2 Table Formatting Consistency
**Solution:**
- Ensure all CLI option tables have consistent column widths
- Use `_None_` or `—` consistently for empty defaults

### 3.3 Code Block Language Tags
**Solution:**
- Ensure all code blocks have appropriate language tags (`bash`, `yaml`, etc.)

---

## Implementation Order

1. **all.md** - Remove YAML duplication (highest impact)
2. **getting-started.md** - Add CLI conventions section
3. **all.md** - Simplify charge explanation
4. **concepts.md** - Improve structure
5. **getting-started.md** - Add Quick Reference
6. **Apply same changes to ja/** - Mirror improvements

---

## Expected Final Scores

| Category | Before | After |
|----------|--------|-------|
| Accuracy | 90 | 93 |
| Clarity | 86 | 93 |
| Structure | 85 | 93 |
| Layout | 88 | 92 |
| **Overall** | **88** | **93** |

---

## Files to Modify

| File | Changes |
|------|---------|
| `all.md` | Remove YAML duplication, simplify charge explanation |
| `getting-started.md` | Add CLI conventions, Quick Reference section |
| `concepts.md` | Improve navigation, enhance workflow diagram |
| `ja/all.md` | Mirror English changes |
| `ja/getting-started.md` | Mirror English changes |
| `ja/concepts.md` | Mirror English changes |
