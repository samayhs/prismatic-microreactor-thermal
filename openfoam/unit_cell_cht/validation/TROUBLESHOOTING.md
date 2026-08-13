# Troubleshooting Log — Post-Processing Function Objects in a Multi-Region Case

**Context:** extracting energy-balance quantities (outlet bulk temperature, mass flow,
fuel-zone volume) from the converged `chtMultiRegionFoam` solution. Several attempts
failed before the cause was found. Recorded here because the debugging is instructive
and the fix is non-obvious.

## Symptom

Custom `surfaceFieldValue` / `volFieldValue` function objects (FOs) either errored on
the object registry or ran **silently** (no console output, no `postProcessing/` files),
while the built-in `wallHeatFlux` preset worked fine.

## The failure chain (and how each was resolved)

| # | Error / symptom | Cause | Resolution |
|---|-----------------|-------|-----------|
| 1 | inline `-func "volFieldValue(...)"` silent | inline `-func` expects a *preset template* whose arg slots match; a hand-written spec configures a no-op | move to a `-dict` file |
| 2 | `problem while reading header for object …` | a `-dict` file must be a valid OpenFOAM dictionary | add a `FoamFile{}` header |
| 3 | `request for objectRegistry region0 … available: fluid` | ran `postProcess -region fluid`, but the FO defaulted to the `region0` registry | add `region fluid;` to the FO |
| 4 | `request for objectRegistry fluid … available: region0` | ran standalone `postProcess` (no `-region`) → only `region0` loaded, but FO now asked for `fluid` | **root cause below** |
| 5 | preset `patchAverage(patch=outlet,...)` printed `<patchName>` | preset arg not substituted under the expected key | superseded by the fix |
| 6 | `Cannot find functionObject file patchFlowRate` | wrong preset name (OF7 uses `flowRatePatch`) | not needed after the fix |
| 7 | `onEnd is not in enumeration` | `writeControl onEnd` is an OpenFOAM.com/ESI keyword, absent in OF7 | use `writeControl timeStep` (or `writeTime`) |

## Root cause

**The standalone `postProcess` utility loads only ONE mesh region into the object
registry.** In a split conjugate case there is no single mesh — there is `fluid` *and*
`solid`. So a region-tagged FO can never find its region (only one is loaded), and an
un-tagged FO defaults to `region0`, which also isn't what a named region is called.
That is the chicken-and-egg in rows 3–4: fixing one half just flips the error.

## The fix

Define the FOs in `system/controlDict`'s `functions{}` block, tag each with its
`region` (`fluid` / `solid`), use a valid `writeControl` (`timeStep`), and run them
through the **solver's** post-process mode:

```bash
chtMultiRegionFoam -postProcess -latestTime      # NOT: postProcess -region ...
```

The solver loads **all** regions into the registry, so every `region` tag resolves.
See `energyBalance.functions` for the working block.

## Verification

With the fix, the FOs returned `weightedAverage(outlet) of T = 723.83 K` and
`sum(outlet) of phi = 0.0050472 kg/s` — **identical** to the values obtained by parsing
the raw fields in `energy_balance.py` (723.8 K, 5.047 g/s). The two independent methods
agreeing is itself a cross-check on the energy balance.

## Takeaway

None of these were physics or solution errors — they were all post-processing
addressing issues specific to multi-region cases. The raw-field Python parser
(`energy_balance.py`) is kept as an independent, tool-agnostic cross-check, and is
arguably more transparent for a report than a chain of function objects.
