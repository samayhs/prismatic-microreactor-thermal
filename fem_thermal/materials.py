"""
Material property models for the prismatic gas-cooled reactor block.

All temperatures are in Kelvin. Every conductivity model returns BOTH the value
k(T) and its derivative dk/dT, because the FEM solver builds an analytic Newton
tangent -- the temperature dependence of k is what makes the conduction problem
nonlinear, and differentiating it by hand (rather than numerically) is a
deliberate demonstration of FEM depth.

Correlations are representative of nuclear-grade materials (IG-110-class graphite,
a TRISO-in-graphite fuel compact) rather than a specific vendor datasheet; the
functional shapes -- graphite conductivity falling with temperature, a lower and
flatter compact conductivity -- are the physically important features.
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
STEFAN_BOLTZMANN = 5.670374419e-8  # W / (m^2 K^4)

# TRISO fuel design limit. TRISO particles retain fission products reliably up to
# ~1600 C; peak fuel temperature staying below this is the central thermal margin
# argument for a gas-cooled microreactor.
TRISO_LIMIT_K = 1600.0 + 273.15  # ~1873 K


# ---------------------------------------------------------------------------
# Thermal conductivity models  ->  return (k, dk/dT)
# ---------------------------------------------------------------------------
def k_graphite(T: np.ndarray | float):
    """Nuclear-grade graphite matrix conductivity [W/m-K].

    Decreases with temperature (phonon-transport dominated). Linear fit from
    ~130 W/m-K near room temperature to ~50 W/m-K near 1500 K, floored so the
    solver stays well-posed at high temperature.
    """
    T = np.asarray(T, dtype=float)
    k0 = 130.0      # W/m-K near 300 K
    alpha = 0.065   # W/m-K per K
    k_lin = k0 - alpha * (T - 300.0)
    k_min = 30.0
    k = np.maximum(k_lin, k_min)
    # derivative is -alpha where above the floor, 0 where clamped
    dk = np.where(k_lin > k_min, -alpha, 0.0)
    return k, dk


def k_fuel_compact(T: np.ndarray | float):
    """Effective conductivity of a TRISO-in-graphite fuel compact [W/m-K].

    Lower and flatter than the surrounding graphite because the TRISO particles
    and packing degrade bulk conduction.
    """
    T = np.asarray(T, dtype=float)
    k0 = 25.0
    alpha = 0.005
    k_lin = k0 - alpha * (T - 300.0)
    k_min = 12.0
    k = np.maximum(k_lin, k_min)
    dk = np.where(k_lin > k_min, -alpha, 0.0)
    return k, dk


# Registry indexed by the integer material tag the mesh assigns to each element.
MAT_GRAPHITE = 0
MAT_FUEL = 1
CONDUCTIVITY = {
    MAT_GRAPHITE: k_graphite,
    MAT_FUEL: k_fuel_compact,
}


# ---------------------------------------------------------------------------
# Surface / boundary properties
# ---------------------------------------------------------------------------
EMISSIVITY_GRAPHITE = 0.80  # oxidized graphite surface, dimensionless


class Scenario:
    """Bundles the load case: heat generation, coolant convection, and the
    outer-boundary heat rejection path (convection + radiation).

    Two scenarios tell the reactor-design story:

      * NORMAL   -- full fission power, forced helium convection in the channels.
                    Question: peak fuel temperature and its margin to 1600 C.

      * PASSIVE  -- forced cooling lost. Only decay heat is produced, channel
                    convection collapses to a weak natural-convection value, and
                    the block sheds heat mainly by RADIATION from its outer
                    surface to the vessel/ambient. Question: does passive
                    radiative rejection keep the fuel below the TRISO limit?
    """

    def __init__(
        self,
        name: str,
        q_fuel: float,        # volumetric heat generation in fuel [W/m^3]
        h_coolant: float,     # convective coefficient on channel walls [W/m^2-K]
        T_coolant: float,     # bulk helium temperature [K]
        h_outer: float,       # convective coefficient on outer surface [W/m^2-K]
        T_inf: float,         # vessel/ambient sink temperature [K]
        radiation: bool,      # enable sigma-eps-(T^4 - Tinf^4) on outer surface
        emissivity: float = EMISSIVITY_GRAPHITE,
    ):
        self.name = name
        self.q_fuel = q_fuel
        self.h_coolant = h_coolant
        self.T_coolant = T_coolant
        self.h_outer = h_outer
        self.T_inf = T_inf
        self.radiation = radiation
        self.emissivity = emissivity


# Nominal volumetric power density in the fuel compacts [W/m^3]. Representative of
# a gas-cooled microreactor fuel compact (~tens of MW/m^3 in the compact itself).
Q_NOMINAL = 30.0e6

# Decay heat a short time after shutdown is a few percent of nominal.
DECAY_FRACTION = 0.03

NORMAL = Scenario(
    name="Normal operation (forced helium cooling)",
    q_fuel=Q_NOMINAL,
    h_coolant=2500.0,          # forced convection, helium (Dittus-Boelter order)
    T_coolant=600.0,           # ~327 C bulk helium
    h_outer=8.0,               # near-insulated outer boundary during operation
    T_inf=400.0,
    radiation=False,
)

PASSIVE = Scenario(
    name="Loss of forced cooling (decay heat, passive radiation)",
    q_fuel=Q_NOMINAL * DECAY_FRACTION,
    h_coolant=25.0,            # weak natural convection in stagnant helium
    T_coolant=600.0,
    h_outer=6.0,               # weak natural convection to vessel gap
    T_inf=400.0,
    radiation=True,            # radiative rejection now carries the load
)
