/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "mixturekEpsilon.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mixturekEpsilon, 0);
addToRunTimeSelectionTable(RASModel, mixturekEpsilon, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mixturekEpsilon::mixturekEpsilon
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const fluidThermo& thermophysicalModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, rho, U, phi, thermophysicalModel, turbulenceModelName),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.92
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            coeffDict_,
            -0.33
        )
    ),
    C4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C4",
            coeffDict_,
            1.00
        )
    ),

    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),
    Prt_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prt",
            coeffDict_,
            1.0
        )
    ),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateEpsilon("epsilon", mesh_)
    ),
    mut_
    (
        IOobject
        (
            "mut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateMut("mut", mesh_)
    ),
    alphat_
    (
        IOobject
        (
            "alphat",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateAlphat("alphat", mesh_)
    )
{
    bound(k_, kMin_);
    dimensionedScalar kMin2_("kMin2_", k_.dimensions(), 1.0E-10);
    bound(k_, kMin2_);
    bound(epsilon_, epsilonMin_);
    dimensionedScalar epsilonMin2_("epsilonMin2_", epsilon_.dimensions(), 1.0E-10);
    bound(epsilon_, epsilonMin2_);

    mut_ = Cmu_*rho_*sqr(k_)/epsilon_;
    mut_.correctBoundaryConditions();

    alphat_ = mut_/Prt_;
    alphat_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> mixturekEpsilon::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - (mut_/rho_)*dev(twoSymm(fvc::grad(U_))),
            k_.boundaryField().types()
	//overline(u'v'),same as the incompressiblekEpsilonmodel
        )
    );
}


tmp<volSymmTensorField> mixturekEpsilon::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -muEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> mixturekEpsilon::divDevRhoReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(muEff(), U)
      - fvc::div(muEff()*dev2(T(fvc::grad(U))))
    );
}


bool mixturekEpsilon::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        C3_.readIfPresent(coeffDict());
        C4_.readIfPresent(coeffDict());
        sigmak_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        Prt_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void mixturekEpsilon::correct()
{
    if (!turbulence_)
    {
        // Re-calculate viscosity
        mut_ = rho_*Cmu_*sqr(k_)/epsilon_;
        mut_.correctBoundaryConditions();

        // Re-calculate thermal diffusivity
        alphat_ = mut_/Prt_;
        alphat_.correctBoundaryConditions();

        return;
    }

    RASModel::correct();

    volScalarField divU(fvc::div(phi_/fvc::interpolate(rho_)));

    if (mesh_.moving())
    {
        divU += fvc::div(mesh_.phi());
    }

    tmp<volTensorField> tgradU = fvc::grad(U_);
    volScalarField G(GName(), mut_*(tgradU() && dev(twoSymm(tgradU()))));
    tgradU.clear();
    const volScalarField& mix_ract_ = mesh_.lookupObject<volScalarField>("uppgradp");
    // Update epsilon and G at the wall
    epsilon_.boundaryField().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(rho_, epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
	fvm::SuSp(C1_*G/k_,epsilon_)
      - fvm::SuSp(((2.0/3.0)*C1_ + C3_)*rho_*divU, epsilon_)
      - fvm::Sp(C2_*rho_*epsilon_/k_, epsilon_)
      - fvm::SuSp(C4_*mix_ract_/k_,epsilon_)
    );

    epsEqn().relax();

    epsEqn().boundaryManipulate(epsilon_.boundaryField());


    solve(epsEqn);
    bound(epsilon_, epsilonMin_);
    dimensionedScalar epsilonMin2_("epsilonMin2_", epsilon_.dimensions(), 1.0E-10);
    bound(epsilon_, epsilonMin2_);

    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(rho_, k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        fvm::SuSp(G/k_,k_)
      - fvm::SuSp((2.0/3.0)*rho_*divU, k_)
      - fvm::Sp(rho_*epsilon_/k_, k_)
      + fvm::SuSp(mix_ract_/k_,k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);
    dimensionedScalar kMin2_("kMin2_", k_.dimensions(), 1.0E-10);
    bound(k_, kMin2_);


    // Re-calculate viscosity
    mut_ = rho_*Cmu_*sqr(k_)/epsilon_;
    mut_.correctBoundaryConditions();

    // Re-calculate thermal diffusivity
    alphat_ = mut_/Prt_;
    alphat_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
