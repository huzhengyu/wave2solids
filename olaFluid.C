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

\*---------------------------------------------------------------------------*/

#include "olaFluid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(olaFluid, 0);
addToRunTimeSelectionTable(physicsModel, olaFluid, fluid);
addToRunTimeSelectionTable(fluidModel, olaFluid, dictionary);


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

olaFluid::olaFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),

    pd
    (
        IOobject
        (
            "pd",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),

    alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),

    twoPhaseProperties(U(), phi(), "alpha1"),

    rho1(twoPhaseProperties.rho1()),
    rho2(twoPhaseProperties.rho2()),

    // Need to store rho for ddt(rho, U)
    rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + (scalar(1) - alpha1)*rho2,
        alpha1.boundaryField().types()
    ),

    rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phi()
    ),
    
    gh("gh", g() & mesh().C()),
    ghf("ghf", g() & mesh().Cf()),
    pdRefCell(0),
    pdRefValue(0.0),
    pRefValue(0.0),

    interface(alpha1, U(), twoPhaseProperties),
    
    turbulence
    (
        incompressible::turbulenceModel::New(U(), phi(), twoPhaseProperties)
    )

{
    UisRequired();
    // Reset p dimensions
    Info<< "Resetting the dimensions of p" << endl;
    p().dimensions().reset(dimPressure);
    p() = pd + rho*gh;

    rho.oldTime();
    setRefCell(pd, pimple().dict(), pdRefCell, pdRefValue);
    //setRefCell(p(), fluidProperties(), pdRefCell_, pdRefValue_);

    mesh().schemesDict().setFluxRequired(pd.name());

    if (pd.needReference())
    {
        pRefValue = readScalar(pimple().dict().lookup("pRefValue"));

        p() += dimensionedScalar
        (
            "p",
            p().dimensions(),
            pRefValue - getRefCellValue(p(), pdRefCell)
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> olaFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF() =
        (
            mesh().boundary()[patchID].nf()
          & (-turbulence->devReff()().boundaryField()[patchID])
        );


    return tvF;
}


tmp<scalarField> olaFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );
    tpF() = p().boundaryField()[patchID];

    return tpF;
}


bool olaFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    fvc::makeAbsolute(phi(), U());
    
    dynamicFvMesh& mesh = fluidModel::mesh();

    // Take a reference to the pimple control
    pimpleControl& pimple = fluidModel::pimple();

    const Time& runTime = fluidModel::runTime();

    // For now, we check for FSI mesh update here; a better way will be to
    // create a FSI dynamicFvMesh: to-do
    bool meshChanged;

    if (fluidModel::fsiMeshUpdate())
    {
        // The FSI interface is in charge of calling mesh.update()
        meshChanged =  fluidModel::fsiMeshUpdateChanged();
    }
    else
    {
        meshChanged =  mesh.update();
        reduce(meshChanged, orOp<bool>());
    }

    // Update gh fields as the mesh may have moved
    gh = g() & mesh.C();
    ghf = g() & mesh.Cf();
    // gh = g() & (mesh.C() - referencePoint);
    // ghf = g() & (mesh.Cf()- referencePoint);

    //if (correctPhi && meshChanged)
    if (meshChanged)
    {
        #include "correctPhi.H"
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi(), U());

    // CourantNo
    //fluidModel::CourantNo();

    // Pressure-velocity corrector
    while (pimple.loop())
    {
        twoPhaseProperties.correct();        

        #include "alphaEqnSubCycle.H"

        #include "UEqn.H"

        // --- PISO loop
        while (pimple.correct())
        {
            #include "pEqn.H"
        }

        p() = pd + rho*gh;

        if (pd.needReference())
        {
            p() +=
            dimensionedScalar
            (
                "p",
                p().dimensions(),
                pRefValue - getRefCellValue(p(), pdRefCell)
            );
        }

        gradp() = fvc::grad(p());

        gradU() = fvc::grad(U());
        
        turbulence->correct();
        //turbulence().correct();
    }

    // Make the fluxes absolute for when runTime++ is called
    //fvc::makeAbsolute(phi(), U());

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
