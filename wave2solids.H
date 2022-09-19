/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    wave2solids

Description
    interFoam fluid model, allowing for a dynamic mesh.
    This should be equivalent to interDyMFoam

Author
    Philip Cardiff, UCD. All rights reserved.
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    wave2solids.C

\*---------------------------------------------------------------------------*/

#ifndef wave2solids_H
#define wave2solids_H

#include "fluidModel.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"
#endif
#include "pimpleControl.H"
#include "fvm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

/*---------------------------------------------------------------------------*\
                        Class wave2solids Declaration
\*---------------------------------------------------------------------------*/

class wave2solids
:
    public fluidModel
{
    // Private data
    IOdictionary waveDictionary;
	dimensionedVector referencePoint;        
	scalar sL;

	//- pd
	volScalarField pd;

	//- Volume fraction of phase 1
	volScalarField alpha1;

	//- Transport model
	twoPhaseMixture twoPhaseProperties;

        //- Interface
        interfaceProperties interface;

        //- Density of phase 1
        const dimensionedScalar& rho1;

        //- Density of phase 2
        const dimensionedScalar& rho2;

        //- Average density field
        volScalarField rho;

        //- Mass flux
        surfaceScalarField rhoPhi;

        //- Gravity by cell-centre field
        volScalarField gh;

        //- Gravity by face-centre field
        surfaceScalarField ghf;

        //- pd reference cell index
        label pdRefCell;

        //- pd reference value
        scalar pdRefValue;

        //- p reference value
        scalar pRefValue;

        //- Turbulence model
        autoPtr<incompressible::turbulenceModel> turbulence;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        wave2solids(const wave2solids&);

        //- Disallow default bitwise assignment
        void operator=(const wave2solids&);

public:

    //- Runtime type information
    TypeName("wave2solids");

    // Constructors

        //- Construct from components
        wave2solids
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~wave2solids()
        {}


    // Member Functions

        // Access

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> patchViscousForce
            (
                const label patchID
            ) const;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> patchPressureForce
            (
                const label patchID
            ) const;

        // Edit

            //- Evolve the fluid model
            virtual bool evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
