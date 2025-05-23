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

#include "LunPressure.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace granularPressureModels
{
    defineTypeNameAndDebug(Lun, 0);

    addToRunTimeSelectionTable
    (
        granularPressureModel,
        Lun,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::Lun::Lun
(
    const dictionary& coeffDict
)
:
    granularPressureModel(coeffDict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::Lun::~Lun()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Lun::granularPressureCoeff
(
    const volScalarField& alpha1,
    const phaseModel& continuousPhase,
    const volScalarField& g0,
    const volScalarField& rho1,
    const dimensionedScalar& e
) const
{

    return rho1*alpha1*(1 + 2*(1 + e)*alpha1*g0);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Lun::
granularPressureCoeffPrime
(
    const volScalarField& alpha1,
    const phaseModel& continuousPhase,
    const volScalarField& g0,
    const volScalarField& g0prime,
    const volScalarField& rho1,
    const dimensionedScalar& e
) const
{
    return rho1*(1 + alpha1*(1 + e)*(4*g0 + 2*g0prime*alpha1));
}

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Lun::
granularPressureCoeff
(
    const phaseModel& phase1,
    const phaseModel& continuousPhase,
    const PtrList<volScalarField>& g0_im,
    const volScalarField& rho1,
    const dimensionedScalar& e
) const
{
    // pressure coeff from Eq. B3 MFIX Equations 2012

    volScalarField alpha1 = phase1;
    const phaseSystem& fluid = phase1.fluid();
        
    volScalarField pCoeff = alpha1*rho1;
    dimensionedScalar eta = 0.5*( 1.0 + e );

    forAll(fluid.phases(), phasei)
    {
        const phaseModel& phase = fluid.phases()[phasei];
        
        if (&phase != &continuousPhase)
        {
            const volScalarField& alpha = phase;
	
            pCoeff += alpha1*rho1*(4.0 * eta * alpha * g0_im[phasei]);  

        }
    }

    return 1.0*pCoeff;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Lun::
granularPressureCoeffPrime
(
    const phaseModel& phase1,
    const phaseModel& continuousPhase,
    const PtrList<volScalarField>& g0_im,
    const PtrList<volScalarField>& g0prime_im,
    const volScalarField& rho1,
    const dimensionedScalar& e
) const
{
    // pressure coeff from Eq. B3 MFIX Equations 2012

    volScalarField alpha1 = phase1;
    const phaseSystem& fluid = phase1.fluid();

    // Eq. B12 MFIX Equations 2012
    dimensionedScalar eta = 0.5*( 1.0 + e );

    volScalarField pCoeffprime = rho1;
            
    forAll(fluid.phases(), phasei)
    {
        const phaseModel& phase = fluid.phases()[phasei];
        
        if (&phase != &continuousPhase)
        {
            const volScalarField& alpha = phase;
	
            pCoeffprime += rho1*(4.0 * eta * alpha * g0_im[phasei]);  

        }
    }    

    pCoeffprime += alpha1*rho1*(4.0*eta*g0prime_im[phase1.index()]);
    
    return 1.0*pCoeffprime;
}


// ************************************************************************* //
