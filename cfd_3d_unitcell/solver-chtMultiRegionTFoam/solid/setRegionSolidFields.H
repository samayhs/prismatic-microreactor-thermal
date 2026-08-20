solidThermo& thermo = thermos[i];

tmp<volScalarField> trho = thermo.rho();
const volScalarField& rho = trho();

tmp<volScalarField> tcp = thermo.Cp();
const volScalarField& cp = tcp();

// Isotropic conductivity kappa(T) and pressure, for the temperature-based solid
// energy equation (see solveSolid.H). Solving in T instead of enthalpy removes the
// variable-Cp discretisation bias of the stock enthalpy solver: steady conduction
// div(kappa grad T)+q'''=0 is Cp-independent, and Cp(T) enters only the transient
// inertia term rho*Cp*dT/dt.
tmp<volScalarField> tkappa = thermo.kappa();
const volScalarField& kappa = tkappa();

const volScalarField& p = thermo.p();

// Volumetric fission heat source Q [W/m^3] on the "fuel" cellZone, for the
// temperature energy equation. Read from constant/solid/heatSource (qVol, zone).
// This replaces the rho-weighted fvOptions source, which is dimensioned for the
// ENTHALPY equation (rho*field/time) and is dimensionally inconsistent with a T
// equation (which needs W/m^3). The stock enthalpy solver keeps using fvOptions;
// this T solver uses Q instead.
IOdictionary heatSourceDict
(
    IOobject
    (
        "heatSource",
        solidRegions[i].time().constant(),
        solidRegions[i],
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);
const scalar qVol = heatSourceDict.lookupOrDefault<scalar>("qVol", 0.0);
const word qZone = heatSourceDict.lookupOrDefault<word>("zone", "fuel");

volScalarField Q
(
    IOobject
    (
        "Q",
        solidRegions[i].time().timeName(),
        solidRegions[i],
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    solidRegions[i],
    dimensionedScalar("Q", dimPower/dimVolume, 0.0)
);
{
    const label zid = solidRegions[i].cellZones().findZoneID(qZone);
    if (zid != -1)
    {
        const labelList& zc = solidRegions[i].cellZones()[zid];
        forAll(zc, ci)
        {
            Q[zc[ci]] = qVol;
        }
    }
    else
    {
        WarningInFunction << "cellZone '" << qZone
            << "' not found in solid region " << solidRegions[i].name()
            << " -- no volumetric heat source applied." << endl;
    }
}

tmp<volSymmTensorField> taniAlpha;
if (!thermo.isotropic())
{
    volSymmTensorField& aniAlpha = aniAlphas[i];
    tmp<volVectorField> tkappaByCp = thermo.Kappa()/cp;
    const coordinateSystem& coodSys = coordinates[i];

    aniAlpha.primitiveFieldRef() =
        coodSys.R().transformVector(tkappaByCp());
    aniAlpha.correctBoundaryConditions();

    taniAlpha = tmp<volSymmTensorField>
    (
        new volSymmTensorField(aniAlpha)
    );
}


volScalarField& h = thermo.he();

const volScalarField& betav = betavSolid[i];

fv::options& fvOptions = solidHeatSources[i];

solidNoLoopControl& pimple = pimples.solid(i);
