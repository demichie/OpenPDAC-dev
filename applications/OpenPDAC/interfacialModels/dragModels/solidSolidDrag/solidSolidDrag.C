/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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

#include "solidSolidDrag.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(solidSolidDrag, 0);
    addToRunTimeSelectionTable(dragModel, solidSolidDrag, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



Foam::tmp<Foam::volScalarField>
Foam::dragModels::solidSolidDrag::KSolidSolid
(
    const phaseModel& gas,
    const phaseModel& solid1,
    const phaseModel& solid2
) const
{

    const phaseSystem& fluid = gas.fluid();
    const volScalarField& alphag = gas;
    const volScalarField& alphas1 = solid1;
    const volScalarField& alphas2 = solid2;
    const volScalarField& d1 = solid1.d();
    const volScalarField& d2 = solid2.d();
    const volScalarField& rho1 = solid1.rho();
    const volScalarField& rho2 = solid2.rho();
    const scalar Pi = constant::mathematical::pi;
        
    volScalarField g0 = 1.0 / alphag; 
        
    forAll(fluid.phases(), phasei)
    {
        const phaseModel& phase = fluid.phases()[phasei];
        if (phase.incompressible())
        {
            const volScalarField& alphas = phase;
    	    g0 += ( 3.0 * ( alphas / phase.d() ) * d1 * d2 ) / 
    	    ( sqr(alphag) * ( d1 + d2 ) );
        }

    }  

    return
        ( 3.0 * ( 1.0 + E_ ) * ( Pi / 2.0 + Cf_ * sqr(Pi) / 8.0 ) 
        * alphas1 * rho1 * alphas2 * rho2 * sqr( d1 + d2 )
        * g0 ) / ( 2.0 * Pi * ( rho1 * pow(alphas1, 3.0) 
        + rho2 * pow(alphas2, 3.0) ) );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::solidSolidDrag::solidSolidDrag
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dragModel(dict, interface, registerObject),
    interface_(interface),
    gasName_(dict.lookup("gas")),
    solid1Name_(dict.lookup("solid1")),
    solid2Name_(dict.lookup("solid2")),
    E_("E", dimless, dict),
    Cf_("Cf", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::solidSolidDrag::~solidSolidDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::solidSolidDrag::K() const
{
    const phaseModel& gas = interface_.fluid().phases()[gasName_];
    const phaseModel& solid1 = interface_.fluid().phases()[solid1Name_];
    const phaseModel& solid2 = interface_.fluid().phases()[solid2Name_];

    if (interface_.contains(solid1) && interface_.contains(solid2))
    {
        return KSolidSolid(gas, solid1, solid2);
    }

    FatalErrorInFunction
        << "The interface " << interface_.name() << " does not contain two "
        << "out of the gas, liquid and solid phase models."
        << exit(FatalError);

    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::dragModels::solidSolidDrag::Kf() const
{
    return fvc::interpolate(K());
}


// ************************************************************************* //
