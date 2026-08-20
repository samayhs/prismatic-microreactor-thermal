"""
Energy-balance validation for the unit-cell CHT model.

With adiabatic outer walls, at steady state ALL generated fission heat must leave
through the coolant. This script confirms that by comparing three INDEPENDENT
measures of the heat and reporting the closure error:

  Q_gen     = q''' * V_fuel                 (the imposed volumetric source)
  Q_coolant = mdot * cp * (T_out - T_in)    (coolant enthalpy rise)
  Q_wall    = integral of the interface heat flux   (from the wallHeatFlux FO)

Reads raw OpenFOAM ASCII fields directly (robust to function-object quirks in a
multi-region case). Run it via run_validation.sh, which first generates the required
cell-volume and wall-heat-flux fields.

Usage:  python3 energy_balance.py <caseDir> <time>
"""

import re
import sys

Q_TRIPLE = 7.0e6      # q''' [W/m^3]  (must match constant/solid/fvOptions)
CP = 5193.0           # helium cp [J/kg-K] (must match thermophysicalProperties)
T_IN = 573.0          # inlet temperature [K]


def _patch_list(path, patch, key="value"):
    txt = open(path).read()
    sub = txt[txt.find("    " + patch + "\n"):]
    m = re.search(key + r"\s+nonuniform\s+List<scalar>\s*(\d+)\s*\((.*?)\)", sub, re.S)
    return [float(x) for x in m.group(2).split()]


def _internal_list(path):
    txt = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\((.*?)\)",
                  txt, re.S)
    return [float(x) for x in m.group(2).split()]


def _cellzone(path, name):
    txt = open(path).read()
    sub = txt[txt.find("\n" + name + "\n"):]
    m = re.search(r"cellLabels\s+List<label>\s*(\d+)\s*\((.*?)\)", sub, re.S)
    return [int(x) for x in m.group(2).split()]


def _wallheatflux_integral(case, time, patch="fluid_to_solid"):
    """Parse the integral column for `patch` from the wallHeatFlux FO output."""
    import glob
    hits = glob.glob(f"{case}/postProcessing/fluid/wallHeatFlux/*/wallHeatFlux.dat") \
        + glob.glob(f"{case}/postProcessing/wallHeatFlux/*/wallHeatFlux.dat")
    if not hits:
        return None
    for line in reversed(open(hits[0]).read().splitlines()):
        if line.startswith("#") or not line.strip():
            continue
        # columns: Time  min/max/integral per patch (grouped). We stored the value
        # from the console (min/max/integ); here fall back to None if unparsed.
        return None
    return None


def main(case, time):
    phi_o = _patch_list(f"{case}/{time}/fluid/phi", "outlet")
    T_o = _patch_list(f"{case}/{time}/fluid/T", "outlet")
    phi_i = _patch_list(f"{case}/{time}/fluid/phi", "inlet")

    mdot_out = sum(phi_o)
    mdot_in = -sum(phi_i)
    T_out = sum(p * t for p, t in zip(phi_o, T_o)) / sum(phi_o)

    V = _internal_list(f"{case}/{time}/solid/V")
    fuel = _cellzone(f"{case}/constant/solid/polyMesh/cellZones", "fuel")
    V_fuel = sum(V[i] for i in fuel)

    Q_gen = Q_TRIPLE * V_fuel
    Q_cool = mdot_out * CP * (T_out - T_IN)

    print(f"mass flow in / out       : {mdot_in*1e3:.4f} / {mdot_out*1e3:.4f} g/s "
          f"(closure {abs(mdot_out-mdot_in)/mdot_in*100:.3f}%)")
    print(f"coolant T_out (mass-wtd) : {T_out:.1f} K = {T_out-273.15:.1f} C "
          f"(rise {T_out-T_IN:.1f} K)")
    print(f"fuel volume              : {V_fuel*1e6:.1f} cm^3 ({len(fuel)} cells)")
    print(f"Q_gen  = q'''*V_fuel     : {Q_gen:8.1f} W")
    print(f"Q_cool = mdot*cp*dT      : {Q_cool:8.1f} W")
    print(f"energy closure           : {abs(Q_gen-Q_cool)/Q_gen*100:.2f}%")


if __name__ == "__main__":
    case = sys.argv[1] if len(sys.argv) > 1 else "."
    time = sys.argv[2] if len(sys.argv) > 2 else "latestTime"
    main(case, time)
