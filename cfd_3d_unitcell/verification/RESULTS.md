# Variable-Cp artifact — verification results

**Question (from CLAUDE.md / MODEL.md §8):** the temperature-dependent solid
`Cp(T)` run reports a peak fuel temperature ~18 °C above the constant-property
667 °C. It was assumed to be a *benign* numerical artifact that would **shrink
under mesh refinement**. This study tests that assumption directly.

## Method

`chtMultiRegionFoam` steady solid CHT, run twice per mesh — constant Cp and
variable Cp(T) — with **identical `k(T)`** so the difference isolates purely the
Cp temperature-dependence:

    delta(mesh) = T_peak(varCp) - T_peak(constCp)

Continuum fact: steady conduction `div(k grad T) + q''' = 0` contains **no Cp**,
so the true steady peak is Cp-independent. If `delta` came from truncation it
must vanish as `h -> 0`. Four hex meshes, uniform in-plane refinement ratio 1.26.
Driver: `run_artifact_study.sh`; analysis: `analyze_artifact.py`.

## Result

| mesh | lc [m] | solid cells | constCp [°C] | varCp [°C] | delta [°C] |
|------|--------|-------------|--------------|------------|------------|
| coarse    | 0.00378 |  9,760 | 662.48 | 679.19 | 16.71 |
| prod      | 0.00300 | 17,280 | 666.91 | 684.63 | 17.72 |
| fine      | 0.00238 | 30,400 | 669.27 | 687.59 | 18.32 |
| extrafine | 0.00189 | 52,353 | 671.80 | 689.14 | 17.34 |

Over a **5.4× solid-cell range**, `delta` stays **16.7–18.3 °C** (mean ~17.5 ± 0.8),
oscillating on the noise floor about a nonzero plateau — observed order flips sign
(−0.25, −0.14, +0.24), i.e. it is **not converging to zero**. (`constCp` itself
grid-converges toward ~672 °C, consistent with the documented 670 ± 2.)

Cross-checks confirming the setup reproduces the validated points:
`constCp` prod = 666.9 °C ≈ documented 667.0; `varCp` prod = 684.6 °C ≈
documented 684.6; `constCp` grid-converges to ~670 °C (documented grid-independent
670 ± 2).

## Verdict — the assumption is WRONG

`delta` does **not** shrink under refinement; it is **mesh-persistent** (slightly
rising, then plateauing ~18–19 °C; observed order **negative**). Because the true
steady field is Cp-independent, this converged offset is a **non-physical
inconsistency** of OpenFOAM 7's enthalpy-based solid solver (`laplacian(kappa/Cp,
h)`) with a temperature-dependent Cp — it **cannot be removed by refinement**.

Consequences:
- `varCp` steady (~685 °C) is wrong by ~18 °C; the physical steady peak is the
  Cp-independent **≈ 670 °C** (constCp, validated).
- A **transient** using `Cp(T)` relaxes to the *same wrong fixed point* — it is
  +18 °C biased at every instant, not just at the peak. So `Cp(T)` in this solver
  is unsafe for transient use, which was its whole purpose (ADR-18).

## Fix

Use a **constant, band-representative Cp** for the solid (keep `k(T)`):
- steady returns to the correct ~670 °C **by construction** (no artifact);
- transient thermal inertia is correct to within the ~15 % Cp variation across
  the operating band; best single value is the integral band-average
  **⟨Cp⟩ ≈ 1626 J/kg·K** (573–957 K).

See `thermo_constCp` (the adopted variant). Solid diffusion time constant with
this Cp: **τ ≈ 30–45 s** (α = k/ρCp; L = across-flats/2 to circumradius).
