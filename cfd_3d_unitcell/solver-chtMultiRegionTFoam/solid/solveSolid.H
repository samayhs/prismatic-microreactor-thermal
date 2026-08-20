{
    // Temperature-based solid energy equation (replaces the stock enthalpy form).
    //
    //     rho*Cp(T) dT/dt = div(kappa(T) grad T) + q'''
    //
    // Solving directly in T removes the variable-Cp bias of the enthalpy solver:
    // the stock form div(alpha grad h), alpha=kappa/Cp, inserts Cp as a nodal value
    // in alpha but as an interval-mean in grad h, and the two cancel only for
    // constant Cp -- leaving a mesh-independent offset (~17 C here). In the T form,
    // steady conduction div(kappa grad T)+q'''=0 is Cp-INDEPENDENT (correct by
    // construction), and Cp(T) enters ONLY the transient inertia term, where it
    // physically belongs. For a steadyState ddt scheme the ddt term is identically
    // zero, so the steady result carries no Cp dependence at all.
    //
    // We solve on thermo.T() in place so the coupled interface BC (which exchanges
    // temperature + kappa) updates directly; the enthalpy field is resynced to the
    // new T afterwards so the thermo state (he, rho, alpha) stays consistent.
    volScalarField& T = const_cast<volScalarField&>(thermo.T());

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix TEqn
        (
            fvm::ddt(betav*rho*cp, T)
          - fvm::laplacian(betav*kappa, T, "laplacian(alpha,h)")
          ==
            betav*Q            // volumetric fission heat [W/m^3] (see setRegionSolidFields.H)
        );

        TEqn.relax();

        fvOptions.constrain(TEqn);

        TEqn.solve();

        fvOptions.correct(T);
    }

    // Resync enthalpy to the solved temperature, then update rho/alpha via the
    // thermo (its calculate() inverts he->T, which returns the same T).
    h = thermo.he(p, T);
}

thermo.correct();

Info<< "Min/max T:" << min(thermo.T()).value() << ' '
    << max(thermo.T()).value() << endl;
