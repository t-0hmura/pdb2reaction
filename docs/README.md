# pdb2reaction Documentation

This directory contains the documentation source files for pdb2reaction.

## Building the Documentation

### Prerequisites

Install the documentation dependencies (from the repository root):

```bash
pip install -r docs/requirements-docs.txt
```

Or, if you're already inside the `docs/` directory:

```bash
pip install -r requirements-docs.txt
```

Or install the optional documentation extras (if available in pyproject.toml):

```bash
pip install pdb2reaction[docs]
```

### Build HTML Documentation

```bash
cd docs
make html
```

The built documentation will be available at `docs/_build/html/index.html`.

### Live Preview (Development)

For live preview during documentation development:

```bash
cd docs
make livehtml
```

This starts a local server with auto-reload at `http://127.0.0.1:8000`.

### Build Other Formats

```bash
# PDF (requires LaTeX)
make latexpdf

# EPUB
make epub

# Single HTML file
make singlehtml
```

## Documentation Structure

```
docs/
├── conf.py              # Sphinx configuration
├── index.md             # Main documentation index
├── ja/                  # Japanese documentation sources
├── Makefile             # Build script (Linux/macOS)
├── make.bat             # Build script (Windows)
├── requirements-docs.txt # Documentation dependencies
├── _static/             # Static assets (CSS, images)
├── _templates/          # Custom templates
├── yaml-reference.md    # YAML configuration reference
└── *.md                 # Subcommand documentation
```

## Writing Documentation

### Markdown Support

Documentation is written in Markdown with MyST extensions. Supported features include:

- Standard Markdown syntax
- GitHub Flavored Markdown tables
- Code blocks with syntax highlighting
- Admonitions (notes, warnings, tips)
- Cross-references between documents

### Admonitions

```markdown
:::{note}
This is a note.
:::

:::{warning}
This is a warning.
:::

:::{tip}
This is a tip.
:::
```

### Cross-References

Link to other documents:
```markdown
See [YAML Reference](yaml-reference.md) for details.
```

Link to sections:
```markdown
See [Geometry Settings](yaml-reference.md#geom) for details.
```

## Deployment

### GitHub Pages (Recommended)

The workflow file is already configured at `.github/workflows/docs.yml`.

**Setup Steps:**

1. **Enable GitHub Pages in your repository:**
   - Go to **Settings** → **Pages**
   - Under "Build and deployment", select **GitHub Actions** as the source

2. **Push to main branch:**
   - The documentation will automatically build and deploy when you push changes to the `pdb2reaction/docs/` directory

3. **Access your documentation:**
   - URL: `https://<username>.github.io/pdb2reaction/`
   - Example: `https://t-0hmura.github.io/pdb2reaction/`

**Manual trigger:**
You can also manually trigger the build from the Actions tab → Documentation → Run workflow.

**Workflow details:**
- Builds on push to `main`/`master` branches
- Only triggers when files in `pdb2reaction/docs/` are changed
- Pull requests will build but not deploy (for preview)

### Read the Docs (Alternative)

Add a `.readthedocs.yaml` file to the repository root:

```yaml
version: 2

build:
  os: ubuntu-22.04
  tools:
    python: "3.11"

sphinx:
  configuration: pdb2reaction/docs/conf.py

python:
  install:
    - requirements: pdb2reaction/docs/requirements-docs.txt
```

Then connect your repository at [readthedocs.org](https://readthedocs.org/).
