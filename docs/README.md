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
