#!/bin/bash
# Test script to verify Read the Docs configuration

echo "Testing Sphinx documentation build..."

# Change to docs directory
cd "$(dirname "$0")/docs" || exit 1

# Clean previous builds
echo "Cleaning previous builds..."
rm -rf _build/

# Test build with warnings as errors (same as Read the Docs)
echo "Building documentation with strict warning checking..."
python -m sphinx -T -W --keep-going -b html -d _build/doctrees -D language=en . _build/html

if [ $? -eq 0 ]; then
    echo "✅ Documentation build successful!"
    echo "📄 Generated files:"
    ls -la _build/html/*.html | head -5
    echo "🌐 Open _build/html/index.html in your browser to view the documentation"
else
    echo "❌ Documentation build failed!"
    exit 1
fi
