from __future__ import annotations

import panel as pn
from matplotlib import pyplot as plt

from physics_calculators.calculators.trap.core import (
    ATOMIC_MASSES_AMU,
    gaussian_trap_frequencies,
    normalized_trap_profiles,
)


def _summary_markdown(result) -> str:
    return f"""
### Trap summary

- Atom: `{result.species}` (`{result.atomic_mass_amu:.6f} amu`)
- Beam waist: `{result.waist_um:.3f} um`
- Wavelength: `{result.wavelength_nm:.1f} nm`
- Trap depth: `{result.trap_depth_mk:.4f} mK`
- Rayleigh range: `{result.rayleigh_range_um:.3f} um`
- Radial trap frequency: `{result.radial_frequency_hz:.2f} Hz`
- Axial trap frequency: `{result.axial_frequency_hz:.2f} Hz`
"""


def _build_figure(waist_um: float, wavelength_nm: float):
    radial_um, radial_profile, axial_um, axial_profile = normalized_trap_profiles(
        waist_um=waist_um,
        wavelength_nm=wavelength_nm,
        radial_extent_um=max(4.0 * waist_um, 20.0),
        axial_extent_um=max(3.0 * waist_um**2 / wavelength_nm * 1e3, 100.0),
    )

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].plot(radial_um, radial_profile, lw=2)
    axes[0].axvline(-waist_um, ls="--", lw=1, color="gray")
    axes[0].axvline(+waist_um, ls="--", lw=1, color="gray")
    axes[0].set_xlabel("radial position (um)")
    axes[0].set_ylabel("normalized potential")
    axes[0].set_title("Radial Gaussian trap profile")
    axes[0].grid(alpha=0.25)

    axes[1].plot(axial_um, axial_profile, lw=2, color="tab:orange")
    axes[1].set_xlabel("axial position (um)")
    axes[1].set_ylabel("normalized potential")
    axes[1].set_title("Axial trap profile")
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    return fig


def build() -> pn.viewable.Viewable:
    species = pn.widgets.Select(name="Atom species", options=list(ATOMIC_MASSES_AMU.keys()), value="Rb-87")
    wavelength_nm = pn.widgets.FloatInput(name="Trap wavelength (nm)", value=1064.0, step=1.0, start=400.0)
    waist_um = pn.widgets.FloatInput(name="Beam waist w0 (um)", value=40.0, step=1.0, start=1.0)
    trap_depth_mk = pn.widgets.FloatInput(name="Trap depth U0 / kB (mK)", value=1.0, step=0.05, start=0.001)

    controls = pn.Column(
        pn.pane.Markdown(
            """
### Inputs

This example uses the harmonic approximation near the beam focus.
It is a lightweight placeholder for a future ARC-backed polarizability calculator.
"""
        ),
        species,
        wavelength_nm,
        waist_um,
        trap_depth_mk,
        sizing_mode="stretch_width",
    )

    @pn.depends(species, wavelength_nm, waist_um, trap_depth_mk)
    def render(species_value: str, wavelength_nm_value: float, waist_um_value: float, trap_depth_mk_value: float):
        try:
            result = gaussian_trap_frequencies(
                species=species_value,
                wavelength_nm=wavelength_nm_value,
                waist_um=waist_um_value,
                trap_depth_mk=trap_depth_mk_value,
            )
        except Exception as exc:
            return pn.pane.Alert(str(exc), alert_type="danger", sizing_mode="stretch_width")

        figure = _build_figure(waist_um=waist_um_value, wavelength_nm=wavelength_nm_value)
        return pn.Column(
            pn.pane.Markdown(_summary_markdown(result)),
            pn.pane.Matplotlib(figure, tight=True, dpi=144, sizing_mode="stretch_width"),
            sizing_mode="stretch_both",
        )

    return pn.Column(
        pn.pane.Markdown(
            """
## Gaussian Beam Trap Frequencies

This tab demonstrates the second calculator path in the shared wrapper. It is intentionally separated into
`core.py` and `ui.py`, so later you can swap in ARC-based polarizability or rubidium-specific physics without
rewriting the page scaffolding.
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
