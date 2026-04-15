from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


ATOMIC_MASSES_AMU = {
    "Rb-87": 86.9091805,
    "Rb-85": 84.9117897,
    "Cs-133": 132.90545196,
    "Na-23": 22.98976928,
}

AMU_TO_KG = 1.66053906660e-27
BOLTZMANN_CONSTANT = 1.380649e-23


@dataclass(frozen=True)
class TrapFrequencies:
    species: str
    atomic_mass_amu: float
    wavelength_nm: float
    waist_um: float
    trap_depth_mk: float
    rayleigh_range_um: float
    radial_frequency_hz: float
    axial_frequency_hz: float


def gaussian_trap_frequencies(
    species: str,
    wavelength_nm: float,
    waist_um: float,
    trap_depth_mk: float,
) -> TrapFrequencies:
    """
    Harmonic approximation for a red-detuned Gaussian dipole trap.

    Inputs use a supplied trap depth instead of deriving it from power/polarizability,
    which keeps the example Pyodide-compatible and independent of ARC.
    """
    if species not in ATOMIC_MASSES_AMU:
        raise ValueError(f"Unsupported atom species: {species}")
    if wavelength_nm <= 0:
        raise ValueError("Wavelength must be positive.")
    if waist_um <= 0:
        raise ValueError("Beam waist must be positive.")
    if trap_depth_mk <= 0:
        raise ValueError("Trap depth must be positive.")

    mass_kg = ATOMIC_MASSES_AMU[species] * AMU_TO_KG
    wavelength_m = wavelength_nm * 1e-9
    waist_m = waist_um * 1e-6
    trap_depth_joule = trap_depth_mk * 1e-3 * BOLTZMANN_CONSTANT
    rayleigh_range_m = math.pi * waist_m**2 / wavelength_m

    radial_angular_frequency = math.sqrt(4.0 * trap_depth_joule / (mass_kg * waist_m**2))
    axial_angular_frequency = math.sqrt(2.0 * trap_depth_joule / (mass_kg * rayleigh_range_m**2))

    return TrapFrequencies(
        species=species,
        atomic_mass_amu=ATOMIC_MASSES_AMU[species],
        wavelength_nm=wavelength_nm,
        waist_um=waist_um,
        trap_depth_mk=trap_depth_mk,
        rayleigh_range_um=rayleigh_range_m * 1e6,
        radial_frequency_hz=radial_angular_frequency / (2.0 * math.pi),
        axial_frequency_hz=axial_angular_frequency / (2.0 * math.pi),
    )


def normalized_trap_profiles(
    waist_um: float,
    wavelength_nm: float,
    radial_extent_um: float,
    axial_extent_um: float,
    point_count: int = 400,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return normalized radial and axial 1D Gaussian trap profiles."""
    waist_m = waist_um * 1e-6
    wavelength_m = wavelength_nm * 1e-9
    rayleigh_range_m = math.pi * waist_m**2 / wavelength_m

    radial_um = np.linspace(-radial_extent_um, radial_extent_um, point_count)
    axial_um = np.linspace(-axial_extent_um, axial_extent_um, point_count)

    radial_profile = -np.exp(-2.0 * (radial_um * 1e-6) ** 2 / waist_m**2)
    axial_profile = -1.0 / (1.0 + (axial_um * 1e-6 / rayleigh_range_m) ** 2)
    return radial_um, radial_profile, axial_um, axial_profile

