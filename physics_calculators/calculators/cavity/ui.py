from __future__ import annotations

import numpy as np
import panel as pn
from matplotlib import pyplot as plt

from physics_calculators.calculators.cavity.core import (
    cavity_mode_from_geometry,
    hg_intensity_2d,
)


def _summary_markdown(mode, probe_z_mm: float, probe_radius_mm: float) -> str:
    return f"""
### Mode summary

- Stability: `g1 = {mode.g1:.5f}`, `g2 = {mode.g2:.5f}`, `g1*g2 = {mode.g_product:.5f}`
- Waist position: `{mode.waist_position_mm:.4f} mm`
- Rayleigh range: `{mode.rayleigh_range_mm:.4f} mm`
- Waist radius: `{mode.waist_radius_mm * 1e3:.3f} um`
- Beam radius on mirror 1: `{mode.mirror1_beam_radius_mm * 1e3:.3f} um`
- Beam radius on mirror 2: `{mode.mirror2_beam_radius_mm * 1e3:.3f} um`
- Probe plane: `z = {probe_z_mm:.4f} mm`, `w(z) = {probe_radius_mm * 1e3:.3f} um`
- FSR: `{mode.fsr_mhz:.3f} MHz`
"""


def _build_figure(
    mode,
    probe_z_mm: float,
    probe_radius_mm: float,
    transverse_halfwidth_um: float,
    m_index: int,
    n_index: int,
):
    z_mm = np.linspace(0.0, mode.cavity_length_mm, 800)
    w_um = 1e3 * mode.w_of_z(z_mm)

    halfwidth_mm = transverse_halfwidth_um * 1e-3
    x_mm = np.linspace(-halfwidth_mm, halfwidth_mm, 250)
    y_mm = np.linspace(-halfwidth_mm, halfwidth_mm, 250)
    intensity = hg_intensity_2d(
        x_mm=x_mm,
        y_mm=y_mm,
        beam_radius_mm=probe_radius_mm,
        m_index=m_index,
        n_index=n_index,
    )

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    cavity_ax = axes[0]
    cavity_ax.plot(z_mm, +w_um, lw=2, label="+w(z)")
    cavity_ax.plot(z_mm, -w_um, lw=2, label="-w(z)")
    cavity_ax.fill_between(z_mm, +w_um, -w_um, alpha=0.18)

    mirror_height = max(1.15 * np.max(w_um), 1.1 * transverse_halfwidth_um)
    cavity_ax.plot([0, 0], [-mirror_height, mirror_height], lw=4, color="black")
    cavity_ax.plot(
        [mode.cavity_length_mm, mode.cavity_length_mm],
        [-mirror_height, mirror_height],
        lw=4,
        color="black",
    )
    cavity_ax.axvline(mode.waist_position_mm, ls="--", lw=1.5, label="waist")
    cavity_ax.axvline(probe_z_mm, ls=":", lw=1.5, label="probe plane")
    cavity_ax.set_xlabel("z (mm)")
    cavity_ax.set_ylabel("transverse coordinate (um)")
    cavity_ax.set_title("Cavity eigenmode envelope")
    cavity_ax.set_xlim(-0.02 * mode.cavity_length_mm, 1.02 * mode.cavity_length_mm)
    cavity_ax.set_ylim(-mirror_height, mirror_height)
    cavity_ax.grid(alpha=0.25)
    cavity_ax.legend(loc="upper right")

    transverse_ax = axes[1]
    image = transverse_ax.imshow(
        intensity,
        origin="lower",
        extent=[
            -transverse_halfwidth_um,
            transverse_halfwidth_um,
            -transverse_halfwidth_um,
            transverse_halfwidth_um,
        ],
        aspect="equal",
    )
    transverse_ax.set_xlabel("x (um)")
    transverse_ax.set_ylabel("y (um)")
    transverse_ax.set_title(
        rf"$\mathrm{{TEM}}_{{{m_index}{n_index}}}$ at $z={probe_z_mm:.3f}$ mm"
    )
    fig.colorbar(image, ax=transverse_ax, fraction=0.046, pad=0.04, label="normalized intensity")
    fig.tight_layout()
    return fig


def build() -> pn.viewable.Viewable:
    wavelength_nm = pn.widgets.FloatInput(name="Wavelength (nm)", value=780.0, step=1, start=300, end=2000)
    cavity_length_mm = pn.widgets.FloatInput(name="Cavity length L (mm)", value=20.0, step=0.1, start=1)
    mirror1_radius_mm = pn.widgets.FloatInput(name="Mirror 1 radius R1 (mm)", value=25.0, step=1, start=1)
    mirror2_radius_mm = pn.widgets.FloatInput(name="Mirror 2 radius R2 (mm)", value=25.0, step=1, start=1)
    m_index = pn.widgets.IntInput(name="Transverse mode index m", value=0, step=1, start=0, end=6)
    n_index = pn.widgets.IntInput(name="Transverse mode index n", value=0, step=1, start=0, end=6)
    probe_percent = pn.widgets.FloatSlider(
        name="Probe plane (% of cavity length)",
        value=50.0,
        start=0.0,
        end=100.0,
        step=1.0,
    )
    transverse_halfwidth_um = pn.widgets.FloatInput(
        name="Transverse view half-width (um)",
        value=250.0,
        step=10.0,
        start=20.0,
    )

    controls = pn.Column(
        pn.pane.Markdown(
            """
### Inputs

Use a very large mirror radius such as `1e9 mm` to approximate a planar mirror.
""",
            margin=(0, 0, 10, 0),
        ),
        wavelength_nm,
        cavity_length_mm,
        mirror1_radius_mm,
        mirror2_radius_mm,
        m_index,
        n_index,
        probe_percent,
        transverse_halfwidth_um,
        sizing_mode="stretch_width",
    )

    @pn.depends(
        wavelength_nm,
        cavity_length_mm,
        mirror1_radius_mm,
        mirror2_radius_mm,
        m_index,
        n_index,
        probe_percent,
        transverse_halfwidth_um,
    )
    def render(
        wavelength_nm_value: float,
        cavity_length_mm_value: float,
        mirror1_radius_mm_value: float,
        mirror2_radius_mm_value: float,
        m_index_value: int,
        n_index_value: int,
        probe_percent_value: float,
        transverse_halfwidth_um_value: float,
    ):
        try:
            mode = cavity_mode_from_geometry(
                wavelength_nm=wavelength_nm_value,
                cavity_length_mm=cavity_length_mm_value,
                mirror1_radius_mm=mirror1_radius_mm_value,
                mirror2_radius_mm=mirror2_radius_mm_value,
            )
        except Exception as exc:
            return pn.pane.Alert(str(exc), alert_type="danger", sizing_mode="stretch_width")

        probe_z_mm = mode.cavity_length_mm * probe_percent_value / 100.0
        probe_radius_mm = float(mode.w_of_z(probe_z_mm))

        figure = _build_figure(
            mode=mode,
            probe_z_mm=probe_z_mm,
            probe_radius_mm=probe_radius_mm,
            transverse_halfwidth_um=transverse_halfwidth_um_value,
            m_index=int(m_index_value),
            n_index=int(n_index_value),
        )

        return pn.Column(
            pn.pane.Markdown(_summary_markdown(mode, probe_z_mm, probe_radius_mm)),
            pn.pane.Matplotlib(figure, tight=True, dpi=144, sizing_mode="stretch_width"),
            sizing_mode="stretch_both",
        )

    return pn.Column(
        pn.pane.Markdown(
            """
## Cavity Mode Visualizer

This tab reproduces the original `test.py` cavity app in Panel. The left side is pure Python control logic,
the right side is a live view rendered by Panel and suitable for Pyodide export.
"""
        ),
        pn.FlexBox(
            pn.Card(
                controls,
                title="Controls",
                sizing_mode="stretch_height",
                styles={"flex": "1 1 320px", "min-width": "280px"},
            ),
            pn.Card(
                render,
                title="Visualization",
                sizing_mode="stretch_both",
                styles={"flex": "2 1 720px", "min-width": "320px"},
            ),
            flex_wrap="wrap",
            gap="1rem",
            sizing_mode="stretch_width",
        ),
        sizing_mode="stretch_width",
    )
