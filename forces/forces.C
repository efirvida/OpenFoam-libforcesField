/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "forces.H"
#include "fvcGrad.H"
#include "porosityModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace functionObjects
    {
        defineTypeNameAndDebug(forcesFields, 0);
        addToRunTimeSelectionTable(functionObject, forcesFields, dictionary);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
void Foam::functionObjects::forcesFields::initialise()
{
    if (initialised_)
    {
        return;
    }

    if (directForceDensity_)
    {
        if (!obr_.foundObject<volVectorField>(fDName_))
        {
            FatalErrorInFunction
                << "Could not find " << fDName_ << " in database."
                << exit(FatalError);
        }
    }
    else
    {
        if (!obr_.foundObject<volVectorField>(UName_) || !obr_.foundObject<volScalarField>(pName_))
        {
            FatalErrorInFunction
                << "Could not find " << UName_ << ", " << pName_
                << exit(FatalError);
        }

        if (rhoName_ != "rhoInf" && !obr_.foundObject<volScalarField>(rhoName_))
        {
            FatalErrorInFunction
                << "Could not find " << rhoName_
                << exit(FatalError);
        }
    }

    initialised_ = true;
}


void Foam::functionObjects::forcesFields::resetFields()
{
    force_[0] = Zero;
    force_[1] = Zero;
    force_[2] = Zero;

    volVectorField &force =
        const_cast<volVectorField &>(
            lookupObject<volVectorField>("Force"));

    force == dimensionedVector("Force", force.dimensions(), Zero);

}


Foam::tmp<Foam::volSymmTensorField>
Foam::functionObjects::forcesFields::devRhoReff() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (obr_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const cmpTurbModel& turb =
            obr_.lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        return turb.devRhoReff();
    }
    else if (obr_.foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const incompressible::turbulenceModel& turb =
            obr_.lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        return rho()*turb.devReff();
    }
    else if (obr_.foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const fluidThermo& thermo =
            obr_.lookupObject<fluidThermo>(fluidThermo::dictName);

        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        return -thermo.mu()*dev(twoSymm(fvc::grad(U)));
    }
    else if
    (
        obr_.foundObject<transportModel>("transportProperties")
    )
    {
        const transportModel& laminarT =
            obr_.lookupObject<transportModel>("transportProperties");

        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        return -rho()*laminarT.nu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (obr_.foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
             obr_.lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu
        (
            "nu",
            dimViscosity,
            transportProperties.lookup("nu")
        );

        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        return -rho()*nu*dev(twoSymm(fvc::grad(U)));
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        return volSymmTensorField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::forcesFields::mu() const
{
    if (obr_.foundObject<fluidThermo>(basicThermo::dictName))
    {
        const fluidThermo& thermo =
             obr_.lookupObject<fluidThermo>(basicThermo::dictName);

        return thermo.mu();
    }
    else if
    (
        obr_.foundObject<transportModel>("transportProperties")
    )
    {
        const transportModel& laminarT =
            obr_.lookupObject<transportModel>("transportProperties");

        return rho()*laminarT.nu();
    }
    else if (obr_.foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
             obr_.lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu
        (
            "nu",
            dimViscosity,
            transportProperties.lookup("nu")
        );

        return rho()*nu;
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for dynamic viscosity calculation"
            << exit(FatalError);

        return volScalarField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::forcesFields::rho() const
{
    if (rhoName_ == "rhoInf")
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("rho", dimDensity, rhoRef_)
            )
        );
    }
    else
    {
        return(obr_.lookupObject<volScalarField>(rhoName_));
    }
}


Foam::scalar Foam::functionObjects::forcesFields::rho(const volScalarField& p) const
{
    if (p.dimensions() == dimPressure)
    {
        return 1.0;
    }
    else
    {
        if (rhoName_ != "rhoInf")
        {
            FatalErrorInFunction
                << "Dynamic pressure is expected but kinematic is provided."
                << exit(FatalError);
        }

        return rhoRef_;
    }
}


void Foam::functionObjects::forcesFields::addToFields(
    const label patchi,
    const vectorField &fN,
    const vectorField &fT,
    const vectorField &fP)
{
    volVectorField &force =
        const_cast<volVectorField &>(
            lookupObject<volVectorField>("Force"));

    vectorField &pf = force.boundaryFieldRef()[patchi];
    pf += fN + fT + fP;
}


void Foam::functionObjects::forcesFields::addToFields(
    const labelList &cellIDs,
    const vectorField &fN,
    const vectorField &fT,
    const vectorField &fP)
{
    volVectorField &force =
        const_cast<volVectorField &>(
            lookupObject<volVectorField>("Force"));

    forAll(cellIDs, i)
    {
        label celli = cellIDs[i];
        force[celli] += fN[i] + fT[i] + fP[i];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::functionObjects::forcesFields::forcesFields(
    const word &name,
    const Time &runTime,
    const dictionary &dict,
    bool readFields)
    : fvMeshFunctionObject(name, runTime, dict),
      force_(3),
      patchSet_(),
      pName_(word::null),
      UName_(word::null),
      rhoName_(word::null),
      directForceDensity_(false),
      fDName_(""),
      rhoRef_(vGreat),
      pRef_(0),
      porosity_(false),
      initialised_(false)
{
    if (readFields)
    {
        read(dict);
    }
}

Foam::functionObjects::forcesFields::forcesFields(
    const word &name,
    const objectRegistry &obr,
    const dictionary &dict,
    bool readFields)
    : fvMeshFunctionObject(name, obr, dict),
      force_(3),
      patchSet_(),
      pName_(word::null),
      UName_(word::null),
      rhoName_(word::null),
      directForceDensity_(false),
      fDName_(""),
      rhoRef_(vGreat),
      pRef_(0),
      porosity_(false),
      initialised_(false)
{
    if (readFields)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::functionObjects::forcesFields::~forcesFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::functionObjects::forcesFields::calcForcesMoment()
{
    initialise();
    resetFields();

    force_[0] = Zero;
    force_[1] = Zero;
    force_[2] = Zero;

    if (directForceDensity_)
    {
        const volVectorField& fD = obr_.lookupObject<volVectorField>(fDName_);

        const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();

            scalarField sA(mag(Sfb[patchi]));

            // Normal force = surfaceUnitNormal*(surfaceNormal & forceDensity)
            vectorField fN(Sfb[patchi]/sA*(Sfb[patchi] & fD.boundaryField()[patchi]));

            // Tangential force (total force minus normal fN)
            vectorField fT(sA*fD.boundaryField()[patchi] - fN);

            //- Porous force
            vectorField fP(fT.size(), Zero);

            addToFields(patchi, fN, fT, fP);

        }
    }
    else
    {
        const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

        const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();

        tmp<volSymmTensorField> tdevRhoReff = devRhoReff();
        const volSymmTensorField::Boundary& devRhoReffb = tdevRhoReff().boundaryField();

        // Scale pRef by density for incompressible simulations
        scalar pRef = pRef_/rho(p);

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();

            vectorField fN(rho(p)*Sfb[patchi]*(p.boundaryField()[patchi] - pRef));

            vectorField fT(Sfb[patchi] & devRhoReffb[patchi]);

            vectorField fP(fT.size(), Zero);

            addToFields(patchi, fN, fT, fP);

        }
    }

    if (porosity_)
    {
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        const volScalarField rho(this->rho());
        const volScalarField mu(this->mu());

        const HashTable<const porosityModel*> models = obr_.lookupClass<porosityModel>();

        if (models.empty())
        {
            WarningInFunction
                << "Porosity effects requested, but no porosity models found "
                << "in the database"
                << endl;
        }

        forAllConstIter(HashTable<const porosityModel*>, models, iter)
        {
            // non-const access required if mesh is changing
            porosityModel& pm = const_cast<porosityModel&>(*iter());

            vectorField fPTot(pm.force(U, rho, mu));

            const labelList& cellZoneIDs = pm.cellZoneIDs();

            forAll(cellZoneIDs, i)
            {
                label zoneI = cellZoneIDs[i];
                const cellZone& cZone = mesh_.cellZones()[zoneI];

                const vectorField d(mesh_.C(), cZone);
                const vectorField fP(fPTot, cZone);

                const vectorField fDummy(fP.size(), Zero);

                addToFields(cZone, fDummy, fDummy, fP);

            }
        }
    }

}


bool Foam::functionObjects::forcesFields::read(const dictionary &dict)
{
    fvMeshFunctionObject::read(dict);

    initialised_ = false;

    Info << "\nForces: " << nl;

    directForceDensity_ = dict.lookupOrDefault("directForceDensity", false);

    const polyBoundaryMesh &pbm = mesh_.boundaryMesh();

    patchSet_ = pbm.patchSet(wordReList(dict.lookup("patches")));

    if (directForceDensity_)
    {
        // Optional entry for fDName
        fDName_ = dict.lookupOrDefault<word>("fD", "fD");
    }
    else
    {
        // Optional entries U and p
        pName_ = dict.lookupOrDefault<word>("p", "p");
        UName_ = dict.lookupOrDefault<word>("U", "U");
        rhoName_ = dict.lookupOrDefault<word>("rho", "rho");

        // Reference density needed for incompressible
        if (obr_.foundObject<transportModel>("transportProperties"))
        {
            rhoName_ = "rhoInf";
            const dictionary &transportProperties =
                obr_.lookupObject<dictionary>("transportProperties");

            dimensionedScalar rhoRef_(
                "rho",
                dimDensity,
                transportProperties.lookup("rho"));
        }

        // Reference pressure, 0 by default
        pRef_ = dict.lookupOrDefault<scalar>("pRef", 0.0);
    }

    dict.readIfPresent("porosity", porosity_);
    if (porosity_)
    {
        Info << "    Including porosity effects" << endl;
    }
    else
    {
        Info << "    Not including porosity effects" << endl;
    }

    Info << "    Force field will be written\n" << endl;

    volVectorField *forcePtr(
        new volVectorField(
            IOobject(
                "Force",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh_,
            dimensionedVector("Force", dimForce, Zero)));

    mesh_.objectRegistry::store(forcePtr);

    return true;
}


bool Foam::functionObjects::forcesFields::write()
{
    calcForcesMoment();

    lookupObject<volVectorField>("Force").write();

    return true;
}


bool Foam::functionObjects::forcesFields::execute()
{
    return true;
}
// ************************************************************************* //
