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

#include "TwoResistanceHeatTransferPhaseSystem.H"
#include "heatTransferModel.H"
#include "fvmSup.H"
#include "rhoFluidMulticomponentThermo.H"

// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::addDmdtHefs
(
    const phaseSystem::dmdtfTable& dmdtfs,
    const phaseSystem::dmdtfTable& Tfs,
    const latentHeatScheme scheme,
    const latentHeatTransfer transfer,
    phaseSystem::heatTransferTable& eqns
) const
{
    HeatTransferPhaseSystem<BasePhaseSystem>::addDmdtHefsWithoutL
    (
        dmdtfs,
        Tfs,
        scheme,
        eqns
    );

    // Loop the pairs
    forAllConstIter(phaseSystem::dmdtfTable, dmdtfs, dmdtfIter)
    {
        const phaseInterface interface(*this, dmdtfIter.key());

        const volScalarField& dmdtf = *dmdtfIter();

        const volScalarField& Tf = *Tfs[interface];

        const phaseModel& phase1 = interface.phase1();
        const phaseModel& phase2 = interface.phase2();
        const rhoThermo& thermo1 = phase1.thermo();
        const rhoThermo& thermo2 = phase2.thermo();

        // Transfer coefficients
        const sidedBlendedHeatTransferModel& heatTransferModel =
            *heatTransferModels_[interface];
        const volScalarField H1(heatTransferModel.modelInThe(phase1).K());
        const volScalarField H2(heatTransferModel.modelInThe(phase2).K());
        const volScalarField H1Fac(H1/(H1 + H2));
        const volScalarField HEff(H1Fac*H2);

        // Latent heat contribution
        switch (transfer)
        {
            case latentHeatTransfer::heat:
            {
                *eqns[phase1.name()] +=
                  - HEff*(thermo2.T() - thermo1.T()) + H1*(Tf - thermo1.T());

                *eqns[phase2.name()] +=
                  - HEff*(thermo1.T() - thermo2.T()) + H2*(Tf - thermo2.T());

                break;
            }
            case latentHeatTransfer::mass:
            {
                const volScalarField L(this->L(interface, dmdtf, Tf, scheme));

                *eqns[phase1.name()] += H1Fac*dmdtf*L;
                *eqns[phase2.name()] += (1 - H1Fac)*dmdtf*L;

                break;
            }
        }
    }
}


template<class BasePhaseSystem>
void Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::addDmidtHefs
(
    const phaseSystem::dmidtfTable& dmidtfs,
    const phaseSystem::dmdtfTable& Tfs,
    const latentHeatScheme scheme,
    const latentHeatTransfer transfer,
    phaseSystem::heatTransferTable& eqns
) const
{
    HeatTransferPhaseSystem<BasePhaseSystem>::addDmidtHefsWithoutL
    (
        dmidtfs,
        Tfs,
        scheme,
        eqns
    );

    // Loop the pairs
    forAllConstIter(phaseSystem::dmidtfTable, dmidtfs, dmidtfIter)
    {
        const phaseInterface interface(*this, dmidtfIter.key());

        const volScalarField& Tf = *Tfs[interface];

        const phaseModel& phase1 = interface.phase1();
        const phaseModel& phase2 = interface.phase2();
        const rhoThermo& thermo1 = phase1.thermo();
        const rhoThermo& thermo2 = phase2.thermo();

        // Transfer coefficients
        const sidedBlendedHeatTransferModel& heatTransferModel =
            *heatTransferModels_[interface];
        const volScalarField H1(heatTransferModel.modelInThe(phase1).K());
        const volScalarField H2(heatTransferModel.modelInThe(phase2).K());
        const volScalarField H1Fac(H1/(H1 + H2));
        const volScalarField HEff(H1Fac*H2);

        // Loop the species
        forAllConstIter(HashPtrTable<volScalarField>, *dmidtfIter(), dmidtfJter)
        {
            const word& specie = dmidtfJter.key();

            const volScalarField& dmidtf = *dmidtfJter();

            // Latent heat contribution
            switch (transfer)
            {
                case latentHeatTransfer::heat:
                {
                    // Do nothing. This term is handled outside the specie loop.

                    break;
                }
                case latentHeatTransfer::mass:
                {
                    const volScalarField Li
                    (
                        this->Li(interface, specie, dmidtf, Tf, scheme)
                    );

                    *eqns[phase1.name()] += H1Fac*dmidtf*Li;
                    *eqns[phase2.name()] += (1 - H1Fac)*dmidtf*Li;

                    break;
                }
            }
        }

        // Latent heat contribution
        switch (transfer)
        {
            case latentHeatTransfer::heat:
            {
                *eqns[phase1.name()] +=
                  - HEff*(thermo2.T() - thermo1.T()) + H1*(Tf - thermo1.T());

                *eqns[phase2.name()] +=
                  - HEff*(thermo1.T() - thermo2.T()) + H2*(Tf - thermo2.T());

                break;
            }
            case latentHeatTransfer::mass:
            {
                // Do nothing. This term is handled inside the specie loop.

                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
TwoResistanceHeatTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    HeatTransferPhaseSystem<BasePhaseSystem>(mesh),
    heatTransferModels_
    (
        this->template
        generateInterfacialModels<sidedBlendedHeatTransferModel>()
    )
{
    // Check that models have been specified on both sides of the interfaces
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const phaseInterface& interface = heatTransferModelIter()->interface();

        forAllConstIter(phaseInterface, interface, iter)
        {
            if (!heatTransferModelIter()->haveModelInThe(iter()))
            {
                FatalErrorInFunction
                    << "A heat transfer model for the " << iter().name()
                    << " side of the " << interface.name()
                    << " interface is not specified"
                    << exit(FatalError);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
~TwoResistanceHeatTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr
    (
        new phaseSystem::heatTransferTable()
    );

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        eqns.insert
        (
            phase.name(),
            new fvScalarMatrix(phase.thermo().he(), dimEnergy/dimTime)
        );
    }

    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const sidedBlendedHeatTransferModel& model = *heatTransferModelIter();

        const phaseModel& phase1 = model.interface().phase1();
        const phaseModel& phase2 = model.interface().phase2();
        const volScalarField& he1 = phase1.thermo().he();
        const volScalarField& he2 = phase2.thermo().he();
        const volScalarField Cpv1(phase1.thermo().Cpv());
        const volScalarField Cpv2(phase2.thermo().Cpv());

        const volScalarField H1(model.modelInThe(phase1).K());
        const volScalarField H2(model.modelInThe(phase2).K());
        const volScalarField HEff(H1*H2/(H1 + H2));

        *eqns[phase1.name()] +=
            HEff*(phase2.thermo().T() - phase1.thermo().T())
          + H1/Cpv1*he1 - fvm::Sp(H1/Cpv1, he1);
        *eqns[phase2.name()] +=
            HEff*(phase1.thermo().T() - phase2.thermo().T())
          + H2/Cpv2*he2 - fvm::Sp(H2/Cpv2, he2);
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
predictThermophysicalTransport()
{
    BasePhaseSystem::predictThermophysicalTransport();
    correctInterfaceThermo();
}


template<class BasePhaseSystem>
void Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
correctThermophysicalTransport()
{
    BasePhaseSystem::correctThermophysicalTransport();
}


template<class BasePhaseSystem>
bool Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
