/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2025 OpenFOAM Foundation
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
    Foam::OneResistanceHeatTransferPhaseSystem

Description
    Class which models interfacial heat transfer between a number of phases. A
    single heat transfer model is used for each interface.

See also
    TwoResistanceHeatTransferPhaseSystem

SourceFiles
    OneResistanceHeatTransferPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef OneResistanceHeatTransferPhaseSystem_H
#define OneResistanceHeatTransferPhaseSystem_H

#include "HeatTransferPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blendedHeatTransferModel;

/*---------------------------------------------------------------------------*\
            Class OneResistanceHeatTransferPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class OneResistanceHeatTransferPhaseSystem
:
    public HeatTransferPhaseSystem<BasePhaseSystem>
{
protected:

    // Protected typedefs

        typedef HashPtrTable
        <
            blendedHeatTransferModel,
            phaseInterfaceKey,
            phaseInterfaceKey::hash
        > heatTransferModelTable;


    // Protected data

        // Sub Models

            //- Heat transfer models
            const heatTransferModelTable heatTransferModels_;


public:

    // Constructors

        //- Construct from fvMesh
        OneResistanceHeatTransferPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~OneResistanceHeatTransferPhaseSystem();


    // Member Functions

        //- Return the heat transfer matrices
        virtual autoPtr<phaseSystem::heatTransferTable> heatTransfer() const;

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "OneResistanceHeatTransferPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
