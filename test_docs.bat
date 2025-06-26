@echo off
REM Test script to verify Read the Docs configuration

echo Testing Sphinx documentation build...

REM Change to docs directory
cd /d "%~dp0docs"

REM Clean previous builds
echo Cleaning previous builds...
if exist "_build" rmdir /s /q "_build"

REM Test build with warnings as errors (same as Read the Docs)
echo Building documentation with strict warning checking...
python -m sphinx -T -W --keep-going -b html -d _build/doctrees -D language=en . _build/html

if %errorlevel% equ 0 (
    echo ✅ Documentation build successful!
    echo 📄 Generated files:
    dir _build\html\*.html | findstr /V "Directory"
    echo 🌐 Open _build\html\index.html in your browser to view the documentation
) else (
    echo ❌ Documentation build failed!
    exit /b 1
)
