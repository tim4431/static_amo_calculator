from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.hermite import hermval


@dataclass(frozen=True)
class CavityModeResult:
    wavelength_nm: float
    cavity_length_mm: float
    mirror1_radius_mm: float
    mirror2_radius_mm: float
    g1: float
    g2: float
    g_product: float
    waist_position_mm: float
    rayleigh_range_mm: float
    waist_radius_mm: float
    mirror1_beam_radius_mm: float
    mirror2_beam_radius_mm: float
    fsr_mhz: float

    def w_of_z(self, z_mm: float | np.ndarray) -> np.ndarray:
        z = np.asarray(z_mm, dtype=float)
        return self.waist_radius_mm * np.sqrt(
            1.0 + ((z - self.waist_position_mm) / self.rayleigh_range_mm) ** 2
        )


def cavity_mode_from_geometry(
    wavelength_nm: float,
    cavity_length_mm: float,
    mirror1_radius_mm: float,
    mirror2_radius_mm: float,
) -> CavityModeResult:
    """Solve the Gaussian eigenmode of a two-mirror linear cavity."""
    if wavelength_nm <= 0:
        raise ValueError("Wavelength must be positive.")
    if cavity_length_mm <= 0:
        raise ValueError("Cavity length must be positive.")
    if mirror1_radius_mm <= 0 or mirror2_radius_mm <= 0:
        raise ValueError("This calculator assumes positive concave radii.")

    wavelength_mm = wavelength_nm * 1e-6
    mirror1 = math.inf if mirror1_radius_mm >= 1e8 else float(mirror1_radius_mm)
    mirror2 = math.inf if mirror2_radius_mm >= 1e8 else float(mirror2_radius_mm)

    def g_parameter(radius_mm: float) -> float:
        return 1.0 if math.isinf(radius_mm) else 1.0 - cavity_length_mm / radius_mm

    g1 = g_parameter(mirror1)
    g2 = g_parameter(mirror2)
    g_product = g1 * g2

    if not (0.0 <= g_product <= 1.0):
        raise ValueError(f"Unstable cavity: g1*g2 = {g_product:.6f} is outside [0, 1].")
    if math.isinf(mirror1) and math.isinf(mirror2):
        raise ValueError("Plane-plane cavity is not stable.")

    if math.isinf(mirror1) and not math.isinf(mirror2):
        if mirror2 <= cavity_length_mm:
            raise ValueError("For a plano-concave cavity, require R2 > L.")
        waist_position_mm = 0.0
        rayleigh_range_mm = math.sqrt(cavity_length_mm * (mirror2 - cavity_length_mm))
    elif math.isinf(mirror2) and not math.isinf(mirror1):
        if mirror1 <= cavity_length_mm:
            raise ValueError("For a concave-plano cavity, require R1 > L.")
        waist_position_mm = cavity_length_mm
        rayleigh_range_mm = math.sqrt(cavity_length_mm * (mirror1 - cavity_length_mm))
    else:
        if (
            abs(mirror1 - cavity_length_mm) < 1e-12
            and abs(mirror2 - cavity_length_mm) < 1e-12
        ):
            waist_position_mm = cavity_length_mm / 2.0
            rayleigh_range_mm = cavity_length_mm / 2.0
        else:
            denominator = mirror1 + mirror2 - 2.0 * cavity_length_mm
            if abs(denominator) < 1e-14:
                raise ValueError("Geometry is too close to a singular marginal case.")

            waist_position_mm = cavity_length_mm * (mirror2 - cavity_length_mm) / denominator
            rayleigh_range_squared = waist_position_mm * (mirror1 - waist_position_mm)
            if rayleigh_range_squared <= 0:
                raise ValueError("No real Gaussian eigenmode exists for this geometry.")
            rayleigh_range_mm = math.sqrt(rayleigh_range_squared)

    waist_radius_mm = math.sqrt(wavelength_mm * rayleigh_range_mm / math.pi)
    c = 299_792_458.0
    fsr_mhz = c / (2.0 * cavity_length_mm * 1e-3) / 1e6

    result = CavityModeResult(
        wavelength_nm=wavelength_nm,
        cavity_length_mm=cavity_length_mm,
        mirror1_radius_mm=mirror1,
        mirror2_radius_mm=mirror2,
        g1=g1,
        g2=g2,
        g_product=g_product,
        waist_position_mm=waist_position_mm,
        rayleigh_range_mm=rayleigh_range_mm,
        waist_radius_mm=waist_radius_mm,
        mirror1_beam_radius_mm=0.0,
        mirror2_beam_radius_mm=0.0,
        fsr_mhz=fsr_mhz,
    )
    return CavityModeResult(
        **{
            **result.__dict__,
            "mirror1_beam_radius_mm": float(result.w_of_z(0.0)),
            "mirror2_beam_radius_mm": float(result.w_of_z(cavity_length_mm)),
        }
    )


def hg_intensity_2d(
    x_mm: np.ndarray,
    y_mm: np.ndarray,
    beam_radius_mm: float,
    m_index: int,
    n_index: int,
) -> np.ndarray:
    """Return a normalized HG_mn intensity image."""
    x_grid, y_grid = np.meshgrid(x_mm, y_mm)
    xi = np.sqrt(2.0) * x_grid / beam_radius_mm
    yi = np.sqrt(2.0) * y_grid / beam_radius_mm

    coefficients_m = np.zeros(m_index + 1)
    coefficients_n = np.zeros(n_index + 1)
    coefficients_m[m_index] = 1.0
    coefficients_n[n_index] = 1.0

    hermite_m = hermval(xi, coefficients_m)
    hermite_n = hermval(yi, coefficients_n)
    intensity = (hermite_m**2) * (hermite_n**2) * np.exp(
        -2.0 * (x_grid**2 + y_grid**2) / beam_radius_mm**2
    )

    max_intensity = np.max(intensity)
    if max_intensity > 0:
        intensity = intensity / max_intensity
    return intensity

