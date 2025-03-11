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

\*---------------------------------------------------------------------------*/

#include "ThermalPhaseModel.H"
#include "phaseSystem.H"
#include "fvcMeshPhi.H"
#include "fvcDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::ThermalPhaseModel<BasePhaseModel>::filterPressureWork
(
    const tmp<volScalarField>& pressureWork
) const
{
    const volScalarField& alpha = *this;

    scalar pressureWorkAlphaLimit =
        this->thermo_->properties()
       .lookupOrDefault("pressureWorkAlphaLimit", 0.0);

    if (pressureWorkAlphaLimit > 0)
    {
        return
        (
            max(alpha - pressureWorkAlphaLimit, scalar(0))
           /max(alpha - pressureWorkAlphaLimit, pressureWorkAlphaLimit)
        )*pressureWork;
    }
    else
    {
        return pressureWork;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::ThermalPhaseModel<BasePhaseModel>::ThermalPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    ThermophysicalTransportPhaseModel<BasePhaseModel>
    (
        fluid,
        phaseName,
        referencePhase,
        index
    ),
    g_(fluid.mesh().lookupObject<uniformDimensionedVectorField>("g")),
    totalEnergy(fluid.mesh().solution().dict().subDict("PIMPLE").subDict("energyControl").lookupOrDefault<Switch>("totalEnergy", false))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::ThermalPhaseModel<BasePhaseModel>::~ThermalPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::ThermalPhaseModel<BasePhaseModel>::correctThermo()
{
    BasePhaseModel::correctThermo();

    this->thermo_->correct();
}


template<class BasePhaseModel>
bool Foam::ThermalPhaseModel<BasePhaseModel>::isothermal() const
{
    return false;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::ThermalPhaseModel<BasePhaseModel>::heEqn()
{
    const volScalarField& alpha = *this;
    const volScalarField& rho = this->rho();

    const tmp<volVectorField> tU(this->U());
    const volVectorField& U(tU());

    const tmp<surfaceScalarField> talphaRhoPhi(this->alphaRhoPhi());
    const surfaceScalarField& alphaRhoPhi(talphaRhoPhi());

    const tmp<surfaceScalarField> talphaPhi(this->alphaPhi());
    const surfaceScalarField& alphaPhi(talphaPhi());


    const tmp<surfaceScalarField> tPhi(this->phi());
    const surfaceScalarField& phi(tPhi());


    const tmp<volScalarField> tcontErr(this->continuityError());
    const volScalarField& contErr(tcontErr());

    tmp<volScalarField> tK(this->K());
    const volScalarField& K(tK());

    volScalarField& he = this->thermo_->he();

    tmp<fvScalarMatrix> tEEqn
    (
        fvm::ddt(alpha, rho, he)
      + fvm::div(alphaRhoPhi, he)
      - fvm::Sp(contErr, he)
      + this->divq(he)
     ==
        alpha*this->Qdot()
    );

    Info << "totalEnergy: " << totalEnergy << endl;        
    
    if (totalEnergy)
    {
        tEEqn.ref() += ( fvc::ddt(alpha, rho, K) + fvc::div(alphaRhoPhi, K)
                 - contErr*K - alpha*rho*(U&g_));

        // Add the appropriate pressure-work term
        if (he.name() == this->thermo_->phasePropertyName("e"))
        {
            tEEqn.ref() += filterPressureWork
            (
                fvc::div
                (
                    fvc::absolute(alphaRhoPhi, alpha, rho, U),
                    this->fluidThermo().p()/rho
                )
              + (fvc::ddt(alpha) - contErr/rho)*this->fluidThermo().p()
            );
        }
        else if (this->thermo_->dpdt())
        {
            tEEqn.ref() -= filterPressureWork(alpha*this->fluid().dpdt());
        }        
    }
    else
    {
        // Add the appropriate pressure-work term
        if (he.name() == this->thermo_->phasePropertyName("e"))
        {
            tEEqn.ref() -= filterPressureWork
            (        
                this->fluidThermo().p() * 
                ( fvc::ddt(alpha)
                + fvc::div
                    (
                        fvc::absolute(phi,U),
                        alpha
                    )
                )
            );
        }
        else if (this->thermo_->dpdt())
        {
            volScalarField DpDt(alpha*this->fluid().dpdt() + 
            fvc::div
            (
                fvc::absolute(alphaPhi, alpha, U),
                this->fluidThermo().p()
            ) 
            - fvc::div
            (
                fvc::absolute(phi,U),
                alpha
            )*this->fluidThermo().p()
            );            
            tEEqn.ref() -= filterPressureWork(DpDt);
        }        
    }

    return tEEqn;
}


// ************************************************************************* //
