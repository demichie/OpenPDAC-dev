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
    Foam::phaseSystem

Description
    Class to represent a system of phases and model interfacial transfers
    between them.

SourceFiles
    phaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef phaseSystem_H
#define phaseSystem_H

#include "IOdictionary.H"

#include "phaseModel.H"
#include "phaseInterface.H"
#include "phaseInterfaceKey.H"
#include "HashPtrTable.H"
#include "PtrListDictionary.H"
#include "hashedWordList.H"

#include "pimpleNoLoopControl.H"

#include "IOMRFZoneList.H"
#include "fvModels.H"
#include "fvConstraints.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class interfaceSurfaceTensionModel;
class pressureReference;
class nonOrthogonalSolutionControl;

/*---------------------------------------------------------------------------*\
                         Class phaseSystem Declaration
\*---------------------------------------------------------------------------*/

class phaseSystem
:
    public IOdictionary
{
public:

    // Public Types

        struct alphaControl
        {
            //- Function to calculate the number of explicit MULES sub-cycles
            //  from the alpha Courant number
            autoPtr<Function1<scalar>> nAlphaSubCyclesPtr;

            label nAlphaSubCycles;

            //- Number of alpha correctors
            //  Usually used to improve the accuracy at high Courant numbers
            //  with semi-implicit MULES, MULESCorr = true
            label nAlphaCorr;

            scalar vDotResidualAlpha;

            void correct(const scalar CoNum)
            {
                nAlphaSubCycles = ceil(nAlphaSubCyclesPtr->value(CoNum));
            }

            void read(const dictionary& dict)
            {
                nAlphaSubCyclesPtr = Function1<scalar>::New
                (
                    "nAlphaSubCycles",
                    dimless,
                    dimless,
                    dict
                );

                nAlphaCorr = dict.lookupOrDefault<label>("nAlphaCorr", 1);

                vDotResidualAlpha =
                    dict.lookupOrDefault("vDotResidualAlpha", 1e-4);
            }
        };


    // Public Typedefs

        typedef HashPtrTable<fvVectorMatrix> momentumTransferTable;

        typedef HashPtrTable<fvScalarMatrix> heatTransferTable;

        typedef HashPtrTable<fvScalarMatrix> specieTransferTable;

        typedef PtrListDictionary<phaseModel> phaseModelList;

        typedef UPtrList<phaseModel> phaseModelPartialList;

        typedef
            HashPtrTable
            <
                volScalarField,
                phaseInterfaceKey,
                phaseInterfaceKey::hash
            >
            dmdtfTable;

        typedef
            HashPtrTable
            <
                HashPtrTable<volScalarField>,
                phaseInterfaceKey,
                phaseInterfaceKey::hash
            >
            dmidtfTable;



protected:

    // Protected typedefs

        typedef
            HashPtrTable
            <
                interfaceSurfaceTensionModel,
                phaseInterfaceKey,
                phaseInterfaceKey::hash
            >
            interfaceSurfaceTensionModelTable;

        typedef
            HashTable
            <
                scalar,
                phaseInterfaceKey,
                phaseInterfaceKey::hash
            >
            cAlphaTable;


    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to pimpleNoLoopControl
        const pimpleNoLoopControl& pimple_;

        //- Optional MRF zones
        IOMRFZoneList MRF_;

        //- Name of continuous phase
        word continuousPhaseName_;

        //- Name of optional reference phase which is not solved for
        // but obtained from the sum of the other phases
        word referencePhaseName_;

        //- Phase models
        phaseModelList phaseModels_;

        //- Moving phase models
        phaseModelPartialList movingPhaseModels_;

        //- Stationary phase models
        phaseModelPartialList stationaryPhaseModels_;

        //- Anisothermal phase models
        phaseModelPartialList anisothermalPhaseModels_;

        //- Multi-component phase models
        phaseModelPartialList multicomponentPhaseModels_;

        //- Total volumetric flux
        surfaceScalarField phi_;

        //- Rate of change of pressure
        volScalarField dpdt_;

        //- Interface compression coefficients
        const cAlphaTable cAlphas_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;


        // Sub Models

            //- Surface tension models
            const interfaceSurfaceTensionModelTable
                interfaceSurfaceTensionModels_;


    // Protected member functions

        //- Return the sum of the phase fractions of the moving phases
        tmp<volScalarField> sumAlphaMoving() const;

        //- Re-normalise the velocity of the phases
        //  around the specified mixture mean
        void setMixtureU(const volVectorField& Um);


        // Functions required for interface compression

            //- Normal to interface between two phases
            //  Used for interface compression
            tmp<surfaceVectorField> nHatfv
            (
                const volScalarField& alpha1,
                const volScalarField& alpha2
            ) const;

            //- Normal to interface between two phases dotted with face areas
            //  Used for interface compression
            tmp<surfaceScalarField> nHatf
            (
                const volScalarField& alpha1,
                const volScalarField& alpha2
            ) const;

            //- Curvature of interface between two phases
            //  Used for interface compression
            tmp<volScalarField> K
            (
                const phaseModel& alpha1,
                const phaseModel& alpha2
            ) const;


        // Sub-model construction

            //- Return the dictionary containing interfacial model or value
            //  settings for the given name. Performs backwards compatibility
            //  conversions.
            template<class Type>
            dictionary interfacialDict
            (
                const word& name,
                const bool defaultIsEmpty
            ) const;


public:

    //- Runtime type information
    TypeName("phaseSystem");

    //- Default name of the phase properties dictionary
    static const word propertiesName;


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseSystem,
            dictionary,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from fvMesh
        phaseSystem(const fvMesh& mesh);


    // Selectors

        static autoPtr<phaseSystem> New
        (
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~phaseSystem();


    // Member Functions

        // Access

            //- Return the mesh
            inline const fvMesh& mesh() const;

            //- Return pimpleNoLoopControl
            inline const pimpleNoLoopControl& pimple() const;

            //- Return the phase models
            inline const phaseModelList& phases() const;

            //- Access the phase models
            inline phaseModelList& phases();

            //- Return the models for phases that are moving
            inline const phaseModelPartialList& movingPhases() const;

            //- Access the models for phases that are moving
            inline phaseModelPartialList& movingPhases();

            //- Return the models for phases that are stationary
            inline const phaseModelPartialList& stationaryPhases() const;

            //- Access the models for phases that are stationary
            inline phaseModelPartialList& stationaryPhases();

            //- Return the models for phases that have variable temperature
            inline const phaseModelPartialList& anisothermalPhases() const;

            //- Access the models for phases that have variable temperature
            inline phaseModelPartialList& anisothermalPhases();

            //- Return the models for phases that have multiple species
            inline const phaseModelPartialList& multicomponentPhases() const;

            //- Access the models for phases that have multiple species
            inline phaseModelPartialList& multicomponentPhases();

            //- Return the phase not given as an argument in a two-phase system
            //  An error is generated if the system is not two-phase
            inline const phaseModel& otherPhase(const phaseModel& phase) const;

            //- Return the mixture flux
            inline const surfaceScalarField& phi() const;

            //- Access the mixture flux
            inline surfaceScalarField& phi();

            //- Return the rate of change of the pressure
            inline const volScalarField& dpdt() const;

            //- Access the rate of change of the pressure
            inline volScalarField& dpdt();

            //- Return MRF zones
            inline const IOMRFZoneList& MRF() const;

            //- Access the fvModels
            inline Foam::fvModels& fvModels(fvMesh& mesh);

            //- Access the fvModels
            inline const Foam::fvModels& fvModels() const;

            //- Access the fvConstraints
            inline Foam::fvConstraints& fvConstraints(fvMesh& mesh);

            //- Access the fvConstraints
            inline const Foam::fvConstraints& fvConstraints() const;


        // Sub-model construction

            //- Return the model name. This is the same as the model's typename
            //  but without "Model" on the end.
            template<class ModelType>
            word modelName() const;

            //- Generate interfacial-model lists
            template<class ModelType, class ... InterfaceTypes>
            void generateInterfacialModels
            (
                const dictionary& dict,
                const phaseInterface& interface,
                PtrList<phaseInterface>& interfaces,
                PtrList<ModelType>& models
            ) const;

            //- Generate interfacial-model tables
            template<class ModelType>
            HashPtrTable
            <
                ModelType,
                phaseInterfaceKey,
                phaseInterfaceKey::hash
            > generateInterfacialModels() const;

            //- Generate interfacial-value tables
            template<class ValueType>
            HashTable
            <
                ValueType,
                phaseInterfaceKey,
                phaseInterfaceKey::hash
            > generateInterfacialValues(const word& valueName) const;

            //- Return the dictionary from which to construct a low-level
            //  sub-model. Checks that there is just one sub-dictionary then
            //  returns it.
            template<class ModelType>
            static const dictionary& modelSubDict(const dictionary& dict);

            //- Check that mass transfer is supported across the given interface
            template<class ModelType>
            void validateMassTransfer(const phaseInterface& interface) const;


        // Sub-model lookup

            //- Check availability of a sub model for a given interface
            template<class ModelType>
            bool foundInterfacialModel(const phaseInterface& interface) const;

            //- Return a sub model for an interface
            template<class ModelType>
            const ModelType& lookupInterfacialModel
            (
                const phaseInterface& interface
            ) const;


        // Properties

            //- Return the mixture density
            tmp<volScalarField> rho() const;

            //- Return the maximum packing
            tmp<volScalarField> alfasMax() const;

            word continuousPhaseName() const;

            //- Return the mixture velocity
            tmp<volVectorField> U() const;

            //- Return the surface tension coefficient for an interface
            tmp<volScalarField> sigma(const phaseInterfaceKey& key) const;

            //- Return the surface tension coefficient for an interface on a
            //  patch
            tmp<scalarField> sigma
            (
                const phaseInterfaceKey& key,
                const label patchi
            ) const;

            //- Indicator of the proximity of the interface
            //  Field values are 1 near and 0 away for the interface.
            tmp<volScalarField> nearInterface() const;

            //- Stabilisation for normalisation of the interface normal
            inline const dimensionedScalar& deltaN() const;

            //- Return the mass transfer rate for an interface
            virtual tmp<volScalarField> dmdtf
            (
                const phaseInterfaceKey& key
            ) const;

            //- Return the mass transfer rates for each phase
            virtual PtrList<volScalarField> dmdts() const;

            //- Return the mass transfer pressure implicit coefficients
            //  for each phase
            virtual PtrList<volScalarField> d2mdtdps() const;

            //- Return incompressibility
            bool incompressible() const;


        // Transfers

            //- Return the momentum transfer matrices for the cell-based
            //  algorithm
            virtual autoPtr<momentumTransferTable> momentumTransfer() = 0;

            //- Return the momentum transfer matrices for the face-based
            //  algorithm
            virtual autoPtr<momentumTransferTable> momentumTransferf() = 0;

            //- Return the force fluxes for the cell-based algorithm
            virtual PtrList<surfaceScalarField> Fs() const = 0;

            //- Return the force fluxes for the face-based algorithm
            virtual PtrList<surfaceScalarField> Ffs() const = 0;

            //- Return the inverse of the central + drag + virtual mass
            //  coefficient matrix
            virtual void invADVs
            (
                const PtrList<volScalarField>& As,
                PtrList<volVectorField>& HVms,
                PtrList<PtrList<volScalarField>>& invADVs,
                PtrList<PtrList<surfaceScalarField>>& invADVfs
            ) const = 0;

            //- Return the inverse of the face central + drag + virtual mass
            //  coefficient matrix
            virtual PtrList<PtrList<surfaceScalarField>> invADVfs
            (
                const PtrList<surfaceScalarField>& Afs,
                PtrList<surfaceScalarField>& HVmfs
            ) const = 0;

            //- Returns true if the phase pressure is treated implicitly
            //  in the phase fraction equation
            virtual bool implicitPhasePressure(const phaseModel& phase) const;

            //- Returns true if the phase pressure is treated implicitly
            //  in the phase fraction equation for any phase
            virtual bool implicitPhasePressure() const;

            //- Return the phase diffusivity
            //  divided by the momentum central coefficient
            virtual tmp<surfaceScalarField> alphaDByAf
            (
                const PtrList<volScalarField>& rAs
            ) const = 0;

            //- Return the flux corrections for the cell-based algorithm
            virtual PtrList<surfaceScalarField> ddtCorrs() const = 0;

            //- Set the cell and faces drag correction fields
            virtual void dragCorrs
            (
                PtrList<volVectorField>& dragCorrs,
                PtrList<surfaceScalarField>& dragCorrf
            ) const = 0;

            //- Set the cell and faces drag correction fields
            virtual void dragEnergy
            (
                PtrList<volScalarField>& dragEnergyTransfer
            ) const = 0;


            //- Re-normalise the flux of the phases
            //  around the specified mixture mean
            void setMixturePhi
            (
                const PtrList<surfaceScalarField>& alphafs,
                const surfaceScalarField& phim
            );

            //- Return the heat transfer matrices
            virtual autoPtr<heatTransferTable> heatTransfer() const = 0;

            //- Return the specie transfer matrices
            virtual autoPtr<specieTransferTable> specieTransfer() const = 0;

            //- Return the surface tension force
            tmp<surfaceScalarField> surfaceTension
            (
                const phaseModel& phase
            ) const;


        // Evolution

            //- Solve for the phase fractions
            virtual void solve
            (
                const alphaControl& alphaControls,
                const PtrList<volScalarField>& rAs
            );

            //- Correct the fluid properties other than those listed below
            virtual void correct();

            //- Correct the continuity errors
            virtual void correctContinuityError();

            //- Correct the kinematics
            virtual void correctKinematics();

            //- Correct the thermodynamics
            virtual void correctThermo();

            //- Correct the reactions
            virtual void correctReactions();

            //- Correct the species mass fractions
            virtual void correctSpecies();

            //- Predict the momentumTransport
            virtual void predictMomentumTransport();

            //- Predict the energy transport e.g. alphat
            virtual void predictThermophysicalTransport();

            //- Correct the momentumTransport
            virtual void correctMomentumTransport();

            //- Correct the energy transport e.g. alphat
            virtual void correctThermophysicalTransport();

            //- Update the fluid properties for mesh changes
            virtual void meshUpdate();

            //- Correct fixed-flux BCs to be consistent with the velocity BCs
            void correctBoundaryFlux();

            void correctPhi
            (
                const volScalarField& p_rgh,
                const autoPtr<volScalarField>& divU,
                const pressureReference& pressureReference,
                nonOrthogonalSolutionControl& pimple
            );


        // IO

            //- Read base phaseProperties dictionary
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volScalarField> byDt(const volScalarField& vf);
tmp<surfaceScalarField> byDt(const surfaceScalarField& sf);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseSystemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "phaseSystemTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
