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

Class
    Foam::functionObjects::forces

Description
    Calculates the forces by integrating the pressure and
    skin-friction forces over a given list of patches.

    Member function forces::write() store and write volume field 
    representations of forces.

    Example of function object specification:
    \verbatim
    forceField_1
    {
        type          forcesField;
        libs          ("libforcesField.so");
        ...
        patches       (<patch list>);
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: forces       | yes         |
        patches      | Patches included in the forces calculation | yes |
        p            | Pressure field name     | no          | p
        U            | Velocity field name     | no          | U
        rho          | Density field name (see below) | no   | rho
        CofR         | Centre of rotation (see below) | no   |
        directForceDensity | Force density supplied directly (see below)|no|no
        fD           | Name of force density field (see below) | no | fD
    \endtable

Note
  - For incompressible cases, set \c rho in the \c "constant/transportProperties"
    dictionary.
  - If the force density is supplied directly, set the \c directForceDensity
    flag to 'yes', and supply the force density field using the \c
    fDName entry


SourceFiles
    forces.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_forces_H
#define functionObjects_forces_H

#include "fvMeshFunctionObject.H"
#include "coordinateSystem.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace functionObjects
    {

        /*---------------------------------------------------------------------------*\
                                Class forces Declaration
        \*---------------------------------------------------------------------------*/

        class forcesFields: public fvMeshFunctionObject
        {
            protected:

                //- Pressure, viscous and porous force per bin
                List<Field<vector>> force_;

                //- Patches to integrate forces over
                labelHashSet patchSet_;

                //- Name of pressure field
                word pName_;

                //- Name of velocity field
                word UName_;

                //- Name of density field (optional)
                word rhoName_;

                //- Is the force density being supplied directly?
                Switch directForceDensity_;

                //- The name of the force density (fD) field
                word fDName_;

                //- Reference density needed for incompressible calculations
                scalar rhoRef_;

                //- Reference pressure
                scalar pRef_;

                //- Coordinate system used when evaluating forces/moments
                coordinateSystem coordSys_;

                //- Flag to indicate whether we are using a local co-ordinate sys
                bool localSystem_;

                //- Flag to include porosity effects
                bool porosity_;

                //- Initialised flag
                bool initialised_;

                //- Initialise the fields
                void initialise();

                //- Reset the fields prior to accumulation of force/moments
                void resetFields();

                //- Return the effective viscous stress (laminar + turbulent).
                tmp<volSymmTensorField> devRhoReff() const;

                //- Dynamic viscosity field
                tmp<volScalarField> mu() const;

                //- Return rho if specified otherwise rhoRef
                tmp<volScalarField> rho() const;

                //- Return rhoRef if the pressure field is dynamic, i.e. p/rho
                //  otherwise return 1
                scalar rho(const volScalarField &p) const;

                //- Add patch contributions to force and moment fields
                void addToFields(
                    const label patchi,
                    const vectorField &fN,
                    const vectorField &fT,
                    const vectorField &fP);

                //- Add cell contributions to force and moment fields
                void addToFields(
                    const labelList &cellIDs,
                    const vectorField &fN,
                    const vectorField &fT,
                    const vectorField &fP);

            public:
                //- Runtime type information
                TypeName("forcesField");

                // Constructors

                //- Construct from Time and dictionary
                forcesFields(
                    const word &name,
                    const Time &runTime,
                    const dictionary &dict,
                    const bool readFields = true);

                //- Construct from objectRegistry and dictionary
                forcesFields(
                    const word &name,
                    const objectRegistry &obr,
                    const dictionary &dict,
                    const bool readFields = true);

                //- Destructor
                virtual ~forcesFields();

                //- Read the forces data
                virtual bool read(const dictionary &);

                //- Calculate the forces and moments
                virtual void calcForcesMoment();

                //- Execute, currently does nothing
                virtual bool execute();

                //- Write the forces
                virtual bool write();
        };


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
