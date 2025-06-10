/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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

#include "OpenPDAC.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSup.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::solvers::OpenPDAC::compositionPredictor()
{
    autoPtr<HashPtrTable<fvScalarMatrix>> popBalSpecieTransferPtr =
        populationBalanceSystem_.specieTransfer();
    HashPtrTable<fvScalarMatrix>& popBalSpecieTransfer =
        popBalSpecieTransferPtr();

    fluid_.correctReactions();

    forAll(fluid.multicomponentPhases(), multicomponentPhasei)
    {
        phaseModel& phase = fluid_.multicomponentPhases()[multicomponentPhasei];

        UPtrList<volScalarField>& Y = phase.YRef();
        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();

        forAll(Y, i)
        {
            if (phase.solveSpecie(i))
            {
                fvScalarMatrix YiEqn
                (
                    phase.YiEqn(Y[i])
                 ==
                    *popBalSpecieTransfer[Y[i].name()]
                  + fvModels().source(alpha, rho, Y[i])
                );

                YiEqn.relax();
                
                fvConstraints().constrain(YiEqn);
                
                YiEqn.solve("Yi");
                
                Y[i] = max(0.0,min(Y[i],1.0));
                
                fvConstraints().constrain(Y[i]);
            }
            else
            {
                Y[i].correctBoundaryConditions();
            }
        }
    }

    fluid_.correctSpecies();
}


void Foam::solvers::OpenPDAC::energyPredictor()
{
    autoPtr<HashPtrTable<fvScalarMatrix>> heatTransferPtr =
        heatTransferSystem_.heatTransfer();
    HashPtrTable<fvScalarMatrix>& heatTransfer =
        heatTransferPtr();

    autoPtr<HashPtrTable<fvScalarMatrix>> popBalHeatTransferPtr =
        populationBalanceSystem_.heatTransfer();
    HashPtrTable<fvScalarMatrix>& popBalHeatTransfer =
        popBalHeatTransferPtr();

    forAll(fluid.thermalPhases(), thermalPhasei)
    {
        phaseModel& phase = fluid_.thermalPhases()[thermalPhasei];

        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();

        fvScalarMatrix EEqn
        (
            phase.heEqn()
         ==
            *heatTransfer[phase.name()]
          + *popBalHeatTransfer[phase.name()]
          + fvModels().source(alpha, rho, phase.thermo().he())
        );

        if (pimple.dict().lookupOrDefault<Switch>("totalEnergy", false) &&
            pimple.dict().lookupOrDefault<Switch>("dragEnergyCorrection", false)
           )
        {
            PtrList<volScalarField> dragEnergyTransfers(movingPhases.size());
            momentumTransferSystem_.dragEnergy(dragEnergyTransfers);
            EEqn -= dragEnergyTransfers[thermalPhasei];
        }   

        EEqn.relax();
        fvConstraints().constrain(EEqn);
        EEqn.solve();
        fvConstraints().constrain(phase.thermo().he());
    }

    fluid_.correctThermo();
    fluid_.correctContinuityError(populationBalanceSystem_.dmdts());
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::OpenPDAC::thermophysicalPredictor()
{
    for (int Ecorr=0; Ecorr<nEnergyCorrectors; Ecorr++)
    {

        const word&continuousPhaseName_ = fluid.continuousPhaseName();
        const phaseModel& continuousPhase = fluid.phases()[continuousPhaseName_];


        forAll(fluid.thermalPhases(), thermalPhasei)
        {
            const phaseModel& phase = fluid.thermalPhases()[thermalPhasei];

            if ( (&phase != &continuousPhase) && correctTdispersed )
            {
                volScalarField& heNew = const_cast<volScalarField&>(phase.thermo().he());
                volScalarField heTcont = phase.thermo().he(p_,continuousPhase.thermo().T());
                        
                heNew = pos(heNew)*heNew 
                        + neg0(heNew)*heTcont;

                // volScalarField he1 = phase.thermo().he(p_,phase.thermo().T());
                        
                heNew = pos0(phase-phase.residualAlpha())*heNew 
                        + neg(phase-phase.residualAlpha())*heTcont;

                fluid_.correctThermo();
            }
        }

        fluid_.correctThermo();

        fluid_.predictThermophysicalTransport();
        compositionPredictor();
        energyPredictor();

        forAll(fluid.thermalPhases(), thermalPhasei)
        {
            const phaseModel& phase = fluid.thermalPhases()[thermalPhasei];

            Info<< phase.name() << " min/max T "
                << min(phase.thermo().T()).value()
                << " - "
                << max(phase.thermo().T()).value()
                << endl;
        }
            
        bool checkResidual(true);
        bool doCheck(false);

        forAll(fluid.thermalPhases(), thermalPhasei)
        {
            const phaseModel& phase = fluid.thermalPhases()[thermalPhasei];

            word name(phase.thermo().he().name());
            const DynamicList<SolverPerformance<scalar>>& sp
            (
                Residuals<scalar>::field(mesh, name)
            );
            label n = sp.size();
                
            scalar r0 = cmptMax(sp[n-1].initialResidual());
            Info << name << " initial residual " << r0 << endl;
            if (energyControlDict.found(name))
            {
                doCheck = true;
                scalar residual(energyControlDict.lookup<scalar>(name));
                checkResidual = checkResidual && ( r0 <= residual);;          
            } 
        }
        Info << "Iteration " 
             << Ecorr+1 
             << " Check for intial Energy Residual " 
             << checkResidual 
             << endl;
        if (doCheck) 
        {
            if ( checkResidual )
            {
                convergenceFlag = true;
                break;
            }
            else
            {                   
                convergenceFlag = false;
            }
            Info << "convergenceFlag = " << convergenceFlag << endl;
        }            
    }
}


// ************************************************************************* //
