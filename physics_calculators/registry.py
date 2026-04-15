from __future__ import annotations

import inspect
from dataclasses import dataclass
from importlib import import_module

import panel as pn


@dataclass(frozen=True)
class CalculatorSpec:
    slug: str
    title: str
    description: str
    ui_module: str
    core_module: str


CALCULATORS = [
    CalculatorSpec(
        slug="cavity",
        title="Cavity Mode",
        description="Two-mirror cavity eigenmode visualizer ported from the original marimo file.",
        ui_module="physics_calculators.calculators.cavity.ui",
        core_module="physics_calculators.calculators.cavity.core",
    ),
    CalculatorSpec(
        slug="trap",
        title="Trap Frequencies",
        description="Gaussian beam optical trap frequencies in the harmonic approximation.",
        ui_module="physics_calculators.calculators.trap.ui",
        core_module="physics_calculators.calculators.trap.core",
    ),
]


def _source_view(module_name: str, label: str) -> pn.viewable.Viewable:
    module = import_module(module_name)
    source = inspect.getsource(module)
    module_path = getattr(module, "__file__", module_name)
    editor = pn.widgets.TextAreaInput(
        value=source,
        sizing_mode="stretch_both",
        min_height=620,
    )
    return pn.Column(
        pn.pane.Markdown(
            f"### {label}\n\n`{module_path}`\n\n"
            "This editor is editable in the page. Edits change the in-browser copy only and do not write back to the file on disk."
        ),
        editor,
        sizing_mode="stretch_both",
    )


def _calculator_view(spec: CalculatorSpec) -> pn.viewable.Viewable:
    ui_module = import_module(spec.ui_module)
    calculator = ui_module.build()
    return pn.Column(
        pn.pane.Markdown(spec.description),
        pn.Tabs(
            ("Calculator", calculator),
            ("Core Python", _source_view(spec.core_module, "Core source")),
            ("UI Python", _source_view(spec.ui_module, "UI source")),
            dynamic=True,
            sizing_mode="stretch_width",
        ),
        sizing_mode="stretch_width",
    )


def build_tabs() -> pn.Tabs:
    return pn.Tabs(
        *((spec.title, _calculator_view(spec)) for spec in CALCULATORS),
        dynamic=True,
        sizing_mode="stretch_width",
    )
