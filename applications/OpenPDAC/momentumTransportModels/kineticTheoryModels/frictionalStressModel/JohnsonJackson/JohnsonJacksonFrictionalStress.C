/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "JohnsonJacksonFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{
    defineTypeNameAndDebug(JohnsonJackson, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        JohnsonJackson,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
readCoeffs
(
    const dictionary& coeffDict
)
{
    Fr_.read(coeffDict);
    eta_.read(coeffDict);
    p_.read(coeffDict);

    phi_.read(coeffDict);
    phi_ *= constant::mathematical::pi/180.0;

    alphaDeltaMin_.read(coeffDict);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
JohnsonJackson
(
    const dictionary& coeffDict
)
:
    frictionalStressModel(coeffDict),
    Fr_("Fr", dimensionSet(1, -1, -2, 0, 0), coeffDict),
    eta_("eta", dimless, coeffDict),
    p_("p", dimless, coeffDict),
    phi_("phi", dimless, coeffDict),
    alphaDeltaMin_("alphaDeltaMin", dimless, coeffDict)
{
    phi_ *= constant::mathematical::pi/180.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
~JohnsonJackson()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
frictionalPressure
(
    const phaseModel& phase,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax
) const
{
    const volScalarField alphas = 1.0 - continuousPhase;
        
    return
        Fr_*pow(max(alphas - alphaMinFriction, scalar(0)), eta_)
       /pow(max(alphasMax - alphas, alphaDeltaMin_), p_);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
frictionalPressurePrime
(
    const phaseModel& phase,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax
) const
{
    const volScalarField alphas = 1.0 - continuousPhase;

    // TODO: CHECK IF THIS MAKE THE CONVERGENCE WORST OR BETTER
    // ADDED LINES 117 and 119, COMMENTED 116           
    return Fr_*
    (
        eta_*pow(max(alphas - alphaMinFriction, scalar(0)), eta_ - 1)
       //*(alphasMax - alphas)
       *max(alphasMax - alphas, alphaDeltaMin_)
      + p_*pow(max(alphas - alphaMinFriction, scalar(0)), eta_)
      *pos(alphasMax - alphas - alphaDeltaMin_)
    )/pow(max(alphasMax - alphas, alphaDeltaMin_), p_ + 1);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::nu
(
    const phaseModel& phase,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax,
    const volScalarField& pf,
    const volScalarField& rho,
    const volSymmTensorField& D
) const
{
    return volScalarField::New
    (
        IOobject::groupName
        (
            Foam::typedName<frictionalStressModel>("nu"),
            phase.group()
        ),
        dimensionedScalar(dimTime, 0.5)*pf/rho*sin(phi_)
    );
}


// ************************************************************************* //
