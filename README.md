# Physics Calculators with Panel + Pyodide

This repository is a starting point for a browser-side physics calculator site:

- `Panel` provides the UI, tabs, widgets, plots, and layout.
- `Pyodide` lets the same Python app run statically on GitHub Pages.
- Each calculator is split into `core.py` and `ui.py`.
- `physics_calculators/registry.py` is the shared wrapper that mounts calculators into tabs and shows their source code.

## Current calculators

- `Cavity Mode`: ported from the original `test.py`
- `Trap Frequencies`: a lightweight Gaussian dipole-trap example

## Local development

Create a virtual environment and install dependencies:

```powershell
uv venv .venv --python 3.12
uv pip install --python .venv\Scripts\python.exe -r requirements.txt
```

Run the live-reload development server:

```powershell
.venv\Scripts\panel serve app.py --dev --show
```

## Static build for GitHub Pages

Build a pure static site that runs Python in the browser:

```powershell
.\scripts\build-pages.ps1
```

Then serve the generated `docs/` locally if you want to inspect the static output.

## Architecture pattern for future calculators

When you add a new calculator:

1. Add `physics_calculators/calculators/<name>/core.py`
2. Add `physics_calculators/calculators/<name>/ui.py`
3. Register it in `physics_calculators/registry.py`

That is enough for it to appear as a new tab, together with its Python source panes.

## ARC / Pyodide caveat

Your idea of using ARC for rubidium polarizability is reasonable, but ARC and some of its dependencies may not work inside Pyodide. If that happens, the clean path is:

- keep the same Panel UI structure
- keep the same `core.py`/`ui.py` split
- run that specific calculator on a Python server instead of in-browser Pyodide

So this repo is a good foundation either way.
