import marimo

__generated_with = "0.23.1"
app = marimo.App()


@app.cell(hide_code=True)
def _():
    import math
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    from numpy.polynomial.hermite import hermval

    def cavity_mode_from_geometry(
        wavelength_nm: float,
        L_mm: float,
        R1_mm: float,
        R2_mm: float,
    ) -> dict:
        """
        Two-mirror linear cavity, positive radii = concave mirrors.
        Use a very large radius (e.g. 1e9 mm) to approximate a plane mirror.

        Returns:
            z0_mm : waist position measured from mirror 1
            zR_mm : Rayleigh range
            w0_mm : waist radius
            w1_mm : beam radius on mirror 1
            w2_mm : beam radius on mirror 2
            g1, g2, gprod
            FSR_MHz
            plus helper functions w(z), R(z)
        """
        if wavelength_nm <= 0:
            raise ValueError("Wavelength must be positive.")
        if L_mm <= 0:
            raise ValueError("Cavity length must be positive.")
        if R1_mm <= 0 or R2_mm <= 0:
            raise ValueError("This version assumes positive concave radii.")

        lam_mm = wavelength_nm * 1e-6

        # Treat very large radius as approximately planar
        R1 = math.inf if R1_mm >= 1e8 else float(R1_mm)
        R2 = math.inf if R2_mm >= 1e8 else float(R2_mm)

        def g(R):
            return 1.0 if math.isinf(R) else 1.0 - L_mm / R

        g1 = g(R1)
        g2 = g(R2)
        gprod = g1 * g2

        # Standard stability region, with confocal allowed.
        if not (0.0 <= gprod <= 1.0):
            raise ValueError(f"Unstable cavity: g1*g2 = {gprod:.6f} is outside [0, 1].")
        if math.isinf(R1) and math.isinf(R2):
            raise ValueError("Plane-plane cavity is not stable.")

        # Waist position z0 and Rayleigh range zR
        # z0 is measured from mirror 1 toward mirror 2.
        if math.isinf(R1) and not math.isinf(R2):
            if R2 <= L_mm:
                raise ValueError("For a plano-concave cavity, require R2 > L.")
            z0 = 0.0
            zR = math.sqrt(L_mm * (R2 - L_mm))

        elif math.isinf(R2) and not math.isinf(R1):
            if R1 <= L_mm:
                raise ValueError("For a concave-plano cavity, require R1 > L.")
            z0 = L_mm
            zR = math.sqrt(L_mm * (R1 - L_mm))

        else:
            # Finite-finite cavity
            # Special handling for confocal: R1 = R2 = L
            if abs(R1 - L_mm) < 1e-12 and abs(R2 - L_mm) < 1e-12:
                z0 = L_mm / 2.0
                zR = L_mm / 2.0
            else:
                denom = R1 + R2 - 2.0 * L_mm
                if abs(denom) < 1e-14:
                    raise ValueError(
                        "Geometry is too close to a singular marginal case."
                    )

                z0 = L_mm * (R2 - L_mm) / denom
                zR2 = z0 * (R1 - z0)

                if zR2 <= 0:
                    raise ValueError(
                        "No real Gaussian eigenmode found for this geometry."
                    )
                zR = math.sqrt(zR2)

        w0 = math.sqrt(lam_mm * zR / math.pi)

        def w_of_z(z_mm):
            z = np.asarray(z_mm, dtype=float)
            return w0 * np.sqrt(1.0 + ((z - z0) / zR) ** 2)

        def R_of_z(z_mm):
            z = np.asarray(z_mm, dtype=float)
            dz = z - z0
            out = np.empty_like(dz, dtype=float)
            near_waist = np.isclose(dz, 0.0)
            out[near_waist] = np.inf
            out[~near_waist] = dz[~near_waist] * (1.0 + (zR / dz[~near_waist]) ** 2)
            return out

        c = 299_792_458.0
        FSR_MHz = c / (2.0 * L_mm * 1e-3) / 1e6

        return {
            "wavelength_nm": wavelength_nm,
            "L_mm": L_mm,
            "R1_mm": R1,
            "R2_mm": R2,
            "g1": g1,
            "g2": g2,
            "gprod": gprod,
            "z0_mm": z0,
            "zR_mm": zR,
            "w0_mm": w0,
            "w1_mm": float(w_of_z(0.0)),
            "w2_mm": float(w_of_z(L_mm)),
            "FSR_MHz": FSR_MHz,
            "w_of_z": w_of_z,
            "R_of_z": R_of_z,
        }

    def hg_intensity_2d(
        x_mm: np.ndarray,
        y_mm: np.ndarray,
        w_mm: float,
        m: int,
        n: int,
    ) -> np.ndarray:
        X, Y = np.meshgrid(x_mm, y_mm)

        xi = np.sqrt(2.0) * X / w_mm
        yi = np.sqrt(2.0) * Y / w_mm

        coeff_m = np.zeros(m + 1)
        coeff_n = np.zeros(n + 1)
        coeff_m[m] = 1.0
        coeff_n[n] = 1.0

        Hm = hermval(xi, coeff_m)
        Hn = hermval(yi, coeff_n)

        intensity = (Hm**2) * (Hn**2) * np.exp(-2.0 * (X**2 + Y**2) / w_mm**2)

        max_val = np.max(intensity)
        if max_val > 0:
            intensity = intensity / max_val
        return intensity

    return cavity_mode_from_geometry, hg_intensity_2d, mo, np, plt


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Cavity mode visualizer

    This app visualizes the Gaussian eigenmode of a two-mirror cavity.

    **Conventions**
    - `L`: cavity length
    - `R1`, `R2`: positive concave mirror radii
    - To approximate a plane mirror, set its radius to a very large value such as `1e9 mm`

    The right panel shows the transverse intensity of the Hermite-Gaussian mode
    $\mathrm{TEM}_{mn}$ at a user-selected longitudinal position.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    wavelength_nm = mo.ui.number(
        start=300,
        stop=2000,
        step=1,
        value=780,
        label="Wavelength (nm)",
    )
    L_mm = mo.ui.number(
        start=1,
        stop=500,
        step=0.1,
        value=20.0,
        label="Cavity length L (mm)",
    )
    R1_mm = mo.ui.number(
        start=1,
        stop=1e9,
        step=1,
        value=25.0,
        label="Mirror 1 radius R1 (mm)",
    )
    R2_mm = mo.ui.number(
        start=1,
        stop=1e9,
        step=1,
        value=25.0,
        label="Mirror 2 radius R2 (mm)",
    )

    m = mo.ui.number(
        start=0,
        stop=6,
        step=1,
        value=0,
        label="Transverse mode index m",
    )
    n = mo.ui.number(
        start=0,
        stop=6,
        step=1,
        value=0,
        label="Transverse mode index n",
    )

    probe_percent = mo.ui.slider(
        start=0,
        stop=100,
        step=1,
        value=50,
        show_value=True,
        include_input=True,
        label="Probe plane (% of cavity length)",
    )

    transverse_halfwidth_um = mo.ui.number(
        start=20,
        stop=5000,
        step=10,
        value=250,
        label="Transverse view half-width (µm)",
    )

    controls = mo.vstack(
        [
            wavelength_nm,
            L_mm,
            R1_mm,
            R2_mm,
            m,
            n,
            probe_percent,
            transverse_halfwidth_um,
        ],
        gap=1.0,
    )
    return (
        L_mm,
        R1_mm,
        R2_mm,
        controls,
        m,
        n,
        probe_percent,
        transverse_halfwidth_um,
        wavelength_nm,
    )


@app.cell(hide_code=True)
def _(
    L_mm,
    R1_mm,
    R2_mm,
    cavity_mode_from_geometry,
    probe_percent,
    transverse_halfwidth_um,
    wavelength_nm,
):
    try:
        mode = cavity_mode_from_geometry(
            wavelength_nm=wavelength_nm.value,
            L_mm=L_mm.value,
            R1_mm=R1_mm.value,
            R2_mm=R2_mm.value,
        )
        error_message = None
    except Exception as e:
        mode = None
        error_message = str(e)

    probe_z_mm = None
    w_probe_mm = None
    halfwidth_mm = transverse_halfwidth_um.value * 1e-3

    if mode is not None:
        probe_z_mm = mode["L_mm"] * probe_percent.value / 100.0
        w_probe_mm = float(mode["w_of_z"](probe_z_mm))
    return halfwidth_mm, mode, probe_z_mm, w_probe_mm


@app.cell(hide_code=True)
def _(
    controls,
    halfwidth_mm,
    hg_intensity_2d,
    m,
    mo,
    mode,
    n,
    np,
    plt,
    probe_z_mm,
    w_probe_mm,
):


    # Longitudinal beam envelope
    z = np.linspace(0.0, mode["L_mm"], 800)
    w_um = 1e3 * mode["w_of_z"](z)
    z0 = mode["z0_mm"]
    L = mode["L_mm"]

    # Transverse mode image at chosen z
    x = np.linspace(-halfwidth_mm, halfwidth_mm, 250)
    y = np.linspace(-halfwidth_mm, halfwidth_mm, 250)
    intensity = hg_intensity_2d(
        x_mm=x,
        y_mm=y,
        w_mm=w_probe_mm,
        m=int(m.value),
        n=int(n.value),
    )
    halfwidth_um = halfwidth_mm * 1e3

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    # --- Left panel: cavity geometry + Gaussian envelope ---
    ax = axes[0]
    ax.plot(z, +w_um, lw=2, label="+w(z)")
    ax.plot(z, -w_um, lw=2, label="-w(z)")
    ax.fill_between(z, +w_um, -w_um, alpha=0.18)

    mirror_height = max(1.15 * np.max(w_um), 1.1 * halfwidth_um)
    ax.plot([0, 0], [-mirror_height, mirror_height], lw=4, color="black")
    ax.plot([L, L], [-mirror_height, mirror_height], lw=4, color="black")

    ax.axvline(z0, ls="--", lw=1.5, label="waist")
    ax.axvline(probe_z_mm, ls=":", lw=1.5, label="probe plane")

    ax.set_xlabel("z (mm)")
    ax.set_ylabel("transverse coordinate (µm)")
    ax.set_title("Cavity eigenmode envelope")
    ax.set_xlim(-0.02 * L, 1.02 * L)
    ax.set_ylim(-mirror_height, mirror_height)
    ax.legend(loc="upper right")
    ax.grid(alpha=0.25)

    # --- Right panel: transverse intensity ---
    ax = axes[1]
    im = ax.imshow(
        intensity,
        origin="lower",
        extent=[
            -halfwidth_um,
            halfwidth_um,
            -halfwidth_um,
            halfwidth_um,
        ],
        aspect="equal",
    )
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.set_title(
        rf"$\mathrm{{TEM}}_{{{int(m.value)}{int(n.value)}}}$ at $z={probe_z_mm:.3f}$ mm"
    )
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="normalized intensity")

    plt.tight_layout()
    mo.hstack([controls, fig], widths=[1.0, 2.5], align="start")
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
