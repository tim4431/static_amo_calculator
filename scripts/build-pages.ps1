$ErrorActionPreference = "Stop"

$panel = Join-Path $PSScriptRoot "..\\.venv\\Scripts\\panel.exe"
$root = Resolve-Path (Join-Path $PSScriptRoot "..")

Push-Location $root
try {
    & $panel convert app.py --to pyodide-worker --out docs --title "Physics Calculators" --requirements pyodide-requirements.txt
    Copy-Item docs\app.html docs\index.html -Force
}
finally {
    Pop-Location
}

