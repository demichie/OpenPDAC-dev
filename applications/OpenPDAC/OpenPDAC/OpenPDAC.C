/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
#include "localEulerDdtScheme.H"
#include "surfaceFields.H"
#include "fvcDiv.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcMeshPhi.H"
#include "addToRunTimeSelectionTable.H"
#include "myHydrostaticInitialisation.H"

#include "IOobjectList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(OpenPDAC, 0);
    addToRunTimeSelectionTable(solver, OpenPDAC, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::OpenPDAC::readControls()
{
    fluidSolver::readControls();

    faceMomentum =
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false);

    dragCorrection =
        pimple.dict().lookupOrDefault<Switch>("dragCorrection", false);

    partialElimination =
        pimple.dict().lookupOrDefault<Switch>("partialElimination", false);

    nEnergyCorrectors =
        pimple.dict().lookupOrDefault<int>("nEnergyCorrectors", 1);
}


void Foam::solvers::OpenPDAC::correctCoNum()
{
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
    );

    forAll(phases, phasei)
    {
        sumPhi = max
        (
            sumPhi,
            fvc::surfaceSum(mag(phases[phasei].phi()))().primitiveField()
        );
    }

    CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    const scalar meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();

    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::OpenPDAC::OpenPDAC(fvMesh& mesh)
:
    fluidSolver(mesh),

    faceMomentum
    (
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false)
    ),

    dragCorrection
    (
        pimple.dict().lookupOrDefault<Switch>("dragCorrection", false)
    ),

    partialElimination
    (
        pimple.dict().lookupOrDefault<Switch>("partialElimination", false)
    ),

    nEnergyCorrectors
    (
        pimple.dict().lookupOrDefault<int>("nEnergyCorrectors", 1)
    ),

    trDeltaT
    (
        LTS
      ? new volScalarField
        (
            IOobject
            (
                fv::localEulerDdt::rDeltaTName,
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, 1),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
      : nullptr
    ),

    trDeltaTf
    (
        LTS && faceMomentum
      ? new surfaceScalarField
        (
            IOobject
            (
                fv::localEulerDdt::rDeltaTfName,
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, 1)
        )
      : nullptr
    ),

    buoyancy(mesh),

    fluidPtr(phaseSystem::New(mesh)),

    fluid(fluidPtr()),

    phases(fluid.phases()),

    phi(fluid.phi()),

    p(phases[0].thermoRef().p()),

    p_rgh(buoyancy.p_rgh),

    rho
    (
        IOobject
        (
            "rho",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.rho()
    ),
    
    carrierIdx(0),

    muC(phases[carrierIdx].thermo().mu()),

    muMix
    (
        IOobject
        (
            "muMix",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
        

    U
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    // Initialize cloud
    clouds(rho, U, muMix, buoyancy.g),
    
    pressureReference
    (
        p,
        p_rgh,
        pimple.dict(),
        fluid.incompressible()
    ),

    MRF(fluid.MRF())
{
    // Read the controls
    readControls();

    mesh.schemes().setFluxRequired(p_rgh.name());

    // create ph_rgh (p_rgh for hydrostatic pressure)
    volScalarField& ph_rgh = regIOobject::store
    (
        new volScalarField
        (
            IOobject
            (
                "ph_rgh",
                "0",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    );

    
    // Initialization of hydrostatic pressure profile
    hydrostaticInitialisation
    (
        p_rgh,
        ph_rgh,
        p,
        buoyancy.g,
        buoyancy.hRef,
        buoyancy.gh,
        buoyancy.ghf,
        fluid,
        pimple.dict()
    );
    

    // Correct mixture thermodynamics with new pressure    
    fluid.correctThermo();
    rho = fluid.rho();

    Info << "hRef " << buoyancy.hRef.value() << endl;

    Info<< "min p " << min(p).value() <<
  	               " max p " << max(p).value() << endl;
    Info<< "min p_rgh " << min(p_rgh).value() <<
   	               " max p_rgh " << max(p_rgh).value() << endl;
    Info<< "min rho " << min(rho).value() <<
   	               " max rho " << max(rho).value() << endl;

    // Search for carrier phase
    carrierIdx = 0;

    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];
        if (!phase.incompressible())
        {
    	    Info << phasei << " compressible" << endl;
    	    carrierIdx = phasei;
        }
	    else
        {
    	    Info << phasei << " incompressible" << endl;
        }

    }

    // Carrier phase viscosity
    muC = phases[carrierIdx].thermo().mu();

    Info<< "min muC " << min(muC).value() << " max muC " << max(muC).value() << endl;
    
    // Mixture viscosity
    muMix = muC * pow( 1.0 - ( 1.0 - phases[carrierIdx] ) / 0.62 , -1.55);

    // Compute mass-weighted mixture velocity
    U = 0.0* phases[0].U();
    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];
        U += phase * phase.rho() * phase.U() / rho;

    }

    Info<< "min muMix " << min(muMix).value() << " max muMix " << max(muMix).value() << endl;

    clouds.info();
    
    if (pimple.dict().lookupOrDefault<bool>("hydrostaticInitialisation", false))
    {
 
        const Time& runTime = mesh().time();
        scalar startTime_ = runTime.startTime().value();
        scalar deltaT = runTime.deltaT().value();

        // set small value for deltaT to evolve particles    
        const_cast<Time&>(runTime).setDeltaT(1.e-5*deltaT);

        // increase time iterator            
        const_cast<Time&>(runTime)++;
    
        // evolve particle cloud
        clouds.evolve();

        // restore startTime
        const_cast<Time&>(runTime).setTime(startTime_,startTime_);
                
        // restore deltaT            
        const_cast<Time&>(runTime).setDeltaT(deltaT);

        // write everything (including lagrangian)
        const_cast<Time&>(runTime).writeNow();
    }  
    
    if (transient())
    {
        correctCoNum();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::OpenPDAC::~OpenPDAC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::OpenPDAC::preSolve()
{
    // Read the controls
    readControls();

    if (transient())
    {
        correctCoNum();
    }
    else if (LTS)
    {
        setRDeltaT();
    }

    // Store divU from the previous mesh so that it can be
    // mapped and used in correctPhi to ensure the corrected phi
    // has the same divergence
    if (correctPhi || mesh.topoChanging())
    {
        // Construct and register divU for mapping
        divU = new volScalarField
        (
            "divU0",
            fvc::div
            (
                fvc::absolute(phi, fluid.movingPhases()[0].U())
            )
        );
    }

    fvModels().preUpdateMesh();

    // Update the mesh for topology change, mesh to mesh mapping
    mesh_.update();
}


void Foam::solvers::OpenPDAC::prePredictor()
{
    if (pimple.models())
    {
        fvModels().correct();
    }

    if (pimple.thermophysics() || pimple.flow())
    {
        fluid.solve(rAUs, rAUfs);
        fluid.correct();
        fluid.correctContinuityError();
    }

    if (pimple.flow() && pimple.predictTransport())
    {
        fluid.predictMomentumTransport();
    }
}


void Foam::solvers::OpenPDAC::postCorrector()
{
    if (pimple.flow() && pimple.correctTransport())
    {
        fluid.correctMomentumTransport();
        fluid.correctThermophysicalTransport();
    }
}


void Foam::solvers::OpenPDAC::postSolve()
{
    divU.clear();
    
    
    muMix = muC * pow( max(0.0, 1.0 - ( 1.0 - phases[carrierIdx] )) / 0.62 , -1.55);
    rho = fluid.rho();
    
    U *= 0.0;
    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];
        U += phase * phase.rho() * phase.U() / rho;

    }

    Info<< "min mu " << min(muMix).value() << " max mu " << max(muMix).value() << endl;

    clouds.evolve();
        
}


// ************************************************************************* //