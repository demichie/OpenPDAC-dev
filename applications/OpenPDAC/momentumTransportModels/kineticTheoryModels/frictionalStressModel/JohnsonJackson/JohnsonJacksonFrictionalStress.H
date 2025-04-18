/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson

Description

SourceFiles
    JohnsonJacksonFrictionalStress.C

\*---------------------------------------------------------------------------*/

#ifndef JohnsonJacksonFrictionalStress_H
#define JohnsonJacksonFrictionalStress_H

#include "frictionalStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{

/*---------------------------------------------------------------------------*\
                       Class JohnsonJackson Declaration
\*---------------------------------------------------------------------------*/

class JohnsonJackson
:
    public frictionalStressModel
{
    // Private Data

        //- Material constant for frictional normal stress
        dimensionedScalar Fr_;

        //- Material constant for frictional normal stress
        dimensionedScalar eta_;

        //- Material constant for frictional normal stress
        dimensionedScalar p_;

        //- Angle of internal friction
        dimensionedScalar phi_;

        //- Lower limit for (alphaMax - alpha1)
        dimensionedScalar alphaDeltaMin_;


    // Private member functions

        //- Read the coefficients from coeffDict
        virtual bool readCoeffs(const dictionary& coeffDict);


public:

    //- Runtime type information
    TypeName("JohnsonJackson");


    // Constructors

        //- Construct from the coefficients dictionary
        JohnsonJackson(const dictionary& coeffDict);


    //- Destructor
    virtual ~JohnsonJackson();


    // Member Functions

        virtual tmp<volScalarField> frictionalPressure
        (
            const phaseModel& phase,
            const phaseModel& continuousPhase,
            const dimensionedScalar& alphaMinFriction,
            const volScalarField& alphasMax
        ) const;

        virtual tmp<volScalarField> frictionalPressurePrime
        (
            const phaseModel& phase,
            const phaseModel& continuousPhase,
            const dimensionedScalar& alphaMinFriction,
            const volScalarField& alphasMax
        ) const;

        virtual tmp<volScalarField> nu
        (
            const phaseModel& phase,
            const phaseModel& continuousPhase,
            const dimensionedScalar& alphaMinFriction,
            const volScalarField& alphasMax,
            const volScalarField& pf,
            const volScalarField& rho,
            const volSymmTensorField& D
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frictionalStressModels
} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
