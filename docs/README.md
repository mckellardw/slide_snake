# Documentation

This directory contains the documentation for slide_snake, built using [Sphinx](https://www.sphinx-doc.org/) and hosted on [Read the Docs](https://readthedocs.org/).

## Building the Documentation Locally

To build the documentation locally, you'll need to install the required dependencies:

```bash
pip install -r requirements.txt
```

Then, from the `docs/` directory, run:

```bash
# For HTML output
sphinx-build -b html . _build/html

# For PDF output (requires LaTeX)
sphinx-build -b latex . _build/latex
```

The built documentation will be available in the `_build/` directory.

### Test Read the Docs Compatibility
```bash
# Run the test script (Unix/Linux/macOS)
../test_docs.sh

# Or on Windows
../test_docs.bat

# Or manually test with same flags as Read the Docs
python -m sphinx -T -W --keep-going -b html -d _build/doctrees -D language=en . _build/html
```

## Read the Docs Integration

The documentation is automatically built and deployed by Read the Docs when changes are pushed to the main branch. The configuration is defined in the `.readthedocs.yaml` file in the root of the repository.

## Documentation Structure

- `index.rst` - Main documentation index
- `*.md` - Markdown documentation files
- `conf.py` - Sphinx configuration
- `requirements.txt` - Python dependencies for building docs
- `_static/` - Static files (images, CSS, etc.)
- `_templates/` - Custom Sphinx templates (if needed)

## Supported Formats

The documentation supports:
- **Markdown** files (`.md`) via MyST parser
- **reStructuredText** files (`.rst`) 
- **Jupyter Notebooks** (`.ipynb`) if needed
- **Mermaid diagrams** via sphinxcontrib-mermaid
