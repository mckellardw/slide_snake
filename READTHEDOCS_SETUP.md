# Read the Docs Setup Summary

This file summarizes the changes made to prepare the slide_snake repository for Read the Docs compatibility.

## Files Created/Modified

### 1. `.readthedocs.yaml` (modified)
- Enabled Python requirements installation
- Set `fail_on_warning: true` for better error detection
- Added build jobs for dependency management
- Enabled PDF and EPUB format generation

### 2. `docs/requirements.txt` (modified)
- Added Sphinx core dependencies with version constraints
- Added `sphinx-rtd-theme` for consistent styling
- Added `myst-parser` for Markdown file support
- Added `sphinx-copybutton` for code block copy functionality
- Added optional `sphinxcontrib-mermaid` for diagrams

### 3. `docs/conf.py` (modified)
- Updated to use `sphinx_rtd_theme`
- Added support for Markdown files via MyST parser
- Configured proper extensions for documentation generation
- Added intersphinx mapping for external documentation links
- Set up proper HTML theme options
- Added logo path and copyright information
- Excluded unnecessary files from documentation build

### 4. `docs/index.rst` (modified)
- Improved structure with organized sections
- Added logo display
- Better organization of documentation sections

### 5. New Files Created
- `docs/_static/` - Directory for static assets
- `docs/.gitignore` - Ignore build artifacts
- `docs/README.md` - Documentation for the docs setup
- `docs/Makefile` - Unix build automation
- `docs/make.bat` - Windows build automation

## How to Use

### For Read the Docs (Automatic)
1. Connect your GitHub repository to Read the Docs
2. The `.readthedocs.yaml` configuration will automatically handle the build
3. Documentation will be built and deployed automatically on commits

### For Local Development
```bash
# Install dependencies
cd docs/
pip install -r requirements.txt

# Build documentation
make html
# or on Windows: make.bat html

# Clean build artifacts
make clean
```

### For Live Development
```bash
# Install sphinx-autobuild for live reloading
pip install sphinx-autobuild

# Start live server
make livehtml
```

## Features Enabled

- ✅ **Markdown Support**: All `.md` files in docs/ are processed
- ✅ **Copy Code Buttons**: Easy copying of code blocks
- ✅ **Cross-references**: Links to Python and Snakemake documentation
- ✅ **Modern Theme**: Clean, responsive Read the Docs theme
- ✅ **Logo Integration**: Project logo displayed in documentation
- ✅ **Multi-format Output**: HTML, PDF, and EPUB generation
- ✅ **Error Detection**: Build fails on warnings to catch issues early
- ✅ **GitHub Pages**: Automatic `.nojekyll` file generation

## Troubleshooting

### Common Issues Fixed

1. **CSV Syntax Highlighting Warnings**: 
   - **Issue**: Pygments doesn't recognize 'csv' as a valid lexer
   - **Solution**: Added `suppress_warnings = ['misc.highlighting_failure']` to `conf.py`

2. **Missing _static Directory**:
   - **Issue**: Sphinx expects `_static` directory to exist
   - **Solution**: Created `_static/.gitkeep` to ensure directory exists in version control

3. **Build Failures on Warnings**:
   - **Issue**: `fail_on_warning: true` causes build to fail on minor warnings
   - **Solution**: Properly configured warning suppression for known issues

## Next Steps

1. **Connect to Read the Docs**: 
   - Sign up at readthedocs.org
   - Connect your GitHub repository
   - Configure the project settings

2. **Set up webhooks** (optional but recommended):
   - Enable automatic builds on Git pushes
   - Configure branch-specific documentation

3. **Customize further** (optional):
   - Add custom CSS in `docs/_static/`
   - Configure additional Sphinx extensions
   - Set up automated testing of documentation builds
