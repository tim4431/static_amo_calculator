from __future__ import annotations

import panel as pn

from physics_calculators.registry import build_tabs

pn.extension("mathjax", sizing_mode="stretch_width")


def build_app() -> pn.viewable.Viewable:
    intro = pn.pane.Markdown(
        """
# Physics Calculators

This prototype uses a shared Python wrapper to register multiple calculators, expose them as tabs,
and show both the live UI and the underlying source code. The same codebase can be served locally
with `panel serve` or exported to a static Pyodide site for GitHub Pages.

### Development model

- `core.py`: physics-only computation logic
- `ui.py`: widgets, plots, layout, and interaction rules
- `registry.py`: one place to wire calculators into the site

### Notes

- Edit any Python file and run `panel serve app.py --dev` for live-reload during local development.
- The deployed GitHub Pages site runs Python in the browser through Pyodide.
- Packages such as ARC may not be Pyodide-compatible; in that case keep the same UI wrapper but switch that calculator to a server-backed mode.
"""
    )

    return pn.Column(
        intro,
        build_tabs(),
        sizing_mode="stretch_width",
        max_width=1400,
    )


app = build_app()
app.servable(title="Physics Calculators")

