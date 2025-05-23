/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2025 OpenFOAM Foundation
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

#include "phaseSystem.H"

#include "MULES.H"
#include "subCycle.H"

#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcMeshPhi.H"
#include "fvcSup.H"

#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseSystem::solve
(
    const alphaControl& alphaControls,
    const PtrList<volScalarField>& rAs
)
{
    const label nAlphaSubCycles = alphaControls.nAlphaSubCycles;
    const bool LTS = fv::localEulerDdt::enabled(mesh_);

    // Optional reference phase which is not solved for
    // but obtained from the sum of the other phases
    phaseModel* referencePhasePtr = nullptr;

    // The phases which are solved
    // i.e. the moving phases less the optional reference phase
    phaseModelPartialList solvePhases;

    if (referencePhaseName_ != word::null)
    {
        referencePhasePtr = &phases()[referencePhaseName_];

        solvePhases.setSize(movingPhases().size() - 1);
        label solvePhasesi = 0;
        forAll(movingPhases(), movingPhasei)
        {
            if (&movingPhases()[movingPhasei] != referencePhasePtr)
            {
                solvePhases.set(solvePhasesi++, &movingPhases()[movingPhasei]);
            }
        }
    }
    else
    {
        solvePhases = movingPhases();
    }

    forAll(phases(), phasei)
    {
        phases()[phasei].correctBoundaryConditions();
    }

    // Calculate the void fraction
    volScalarField alphaVoid
    (
        IOobject
        (
            "alphaVoid",
            mesh_.time().name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 1)
    );
    forAll(stationaryPhases(), stationaryPhasei)
    {
        alphaVoid -= stationaryPhases()[stationaryPhasei];
    }

    // Calculate the effective flux of the moving phases
    tmp<surfaceScalarField> tphiMoving(phi_);
    if (stationaryPhases().size())
    {
        tphiMoving = phi_/upwind<scalar>(mesh_, phi_).interpolate(alphaVoid);
    }
    const surfaceScalarField& phiMoving = tphiMoving();

    bool dilatation = false;
    forAll(movingPhases(), movingPhasei)
    {
        if (movingPhases()[movingPhasei].divU().valid())
        {
            dilatation = true;
            break;
        }
    }

for (int acorr=0; acorr<alphaControls.nAlphaCorr; acorr++)
        {
        PtrList<volScalarField::Internal> Sps(phases().size());
        PtrList<volScalarField::Internal> Sus(phases().size());

        forAll(movingPhases(), movingPhasei)
        {
            const phaseModel& phase = movingPhases()[movingPhasei];
            const volScalarField& alpha = phase;
            const label phasei = phase.index();

            Sps.set
            (
                phasei,
                new volScalarField::Internal
                (
                    IOobject
                    (
                        IOobject::groupName("Sp", phase.name()),
                        mesh_.time().name(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar(dimless/dimTime, 0)
                )
            );

            Sus.set
            (
                phasei,
                new volScalarField::Internal
                (
                    IOobject::groupName("Su", phase.name()),
                    min(alpha.v(), scalar(1))
                   *fvc::div(fvc::absolute(phi_, phase.U()))->v()
                )
            );

            if (dilatation)
            {
                // Construct the dilatation rate source term
                volScalarField::Internal vDot
                (
                    volScalarField::Internal::New
                    (
                        "vDot",
                        mesh_,
                        dimensionedScalar(dimless/dimTime, 0)
                    )
                );

                forAll(phases(), phasej)
                {
                    const phaseModel& phase2 = phases()[phasej];
                    const volScalarField& alpha2 = phase2;

                    if (&phase2 != &phase)
                    {
                        if (!phase.stationary() && phase.divU().valid())
                        {
                            vDot += alpha2()*phase.divU()()();
                        }

                        if (!phase2.stationary() && phase2.divU().valid())
                        {
                            vDot -= alpha()*phase2.divU()()();
                        }
                    }
                }

                volScalarField::Internal& Sp = Sps[phasei];
                volScalarField::Internal& Su = Sus[phasei];

                forAll(vDot, celli)
                {
                    if (vDot[celli] > 0)
                    {
                        Sp[celli] -=
                            vDot[celli]
                           /max
                            (
                                1 - alpha[celli],
                                alphaControls.vDotResidualAlpha
                            );
                        Su[celli] +=
                            vDot[celli]
                           /max
                            (
                                1 - alpha[celli],
                                alphaControls.vDotResidualAlpha
                            );
                    }
                    else if (vDot[celli] < 0)
                    {
                        Sp[celli] +=
                            vDot[celli]
                           /max
                            (
                                alpha[celli],
                                alphaControls.vDotResidualAlpha
                            );
                    }
                }
            }
        }

        tmp<volScalarField> trSubDeltaT;

        if (LTS && nAlphaSubCycles > 1)
        {
            trSubDeltaT =
                fv::localEulerDdt::localRSubDeltaT(mesh_, nAlphaSubCycles);
        }

        List<volScalarField*> alphaPtrs(phases().size());
        forAll(phases(), phasei)
        {
            alphaPtrs[phasei] = &phases()[phasei];
        }

        for
        (
            subCycle<volScalarField, subCycleFields> alphaSubCycle
            (
                alphaPtrs,
                nAlphaSubCycles
            );
            !(++alphaSubCycle).end();
        )
        {
            // Create correction fluxes
            PtrList<surfaceScalarField> alphaPhis(phases().size());

            tmp<surfaceScalarField> alphaDByAf;
            if (implicitPhasePressure() && (rAs.size()))
            {
                alphaDByAf = this->alphaDByAf(rAs);
            }

            forAll(movingPhases(), movingPhasei)
            {
                const phaseModel& phase = movingPhases()[movingPhasei];
                const volScalarField& alpha = phase;

                alphaPhis.set
                (
                    phase.index(),
                    new surfaceScalarField
                    (
                        IOobject::groupName("alphaPhiCorr", phase.name()),
                        fvc::flux
                        (
                            phase.phi()(),
                            alpha,
                            "div(phi," + alpha.name() + ')'
                        )
                    )
                );

                surfaceScalarField& alphaPhi = alphaPhis[phase.index()];

                if (!cAlphas_.empty())
                {
                    forAll(phases(), phasei)
                    {
                        const phaseModel& phase2 = phases()[phasei];
                        const volScalarField& alpha2 = phase2;

                        if (&phase2 == &phase) continue;

                        cAlphaTable::const_iterator cAlpha
                        (
                            cAlphas_.find(phaseInterface(phase, phase2))
                        );

                        if (cAlpha != cAlphas_.end())
                        {
                            const surfaceScalarField phir
                            (
                                phase.phi() - phase2.phi()
                            );

                            const surfaceScalarField phic
                            (
                                (mag(phi_) + mag(phir))/mesh_.magSf()
                            );

                            const surfaceScalarField phirc
                            (
                                min(cAlpha()*phic, max(phic))
                               *nHatf(alpha, alpha2)
                            );

                            const word phirScheme
                            (
                                "div(phir,"
                              + alpha2.name() + ',' + alpha.name()
                              + ')'
                            );

                            alphaPhi += fvc::flux
                            (
                                -fvc::flux(-phirc, alpha2, phirScheme),
                                alpha,
                                phirScheme
                            );
                        }
                    }
                }

                if (alphaDByAf.valid())
                {
                    alphaPhi +=
                        alphaDByAf()
                       *fvc::snGrad(alpha, "bounded")*mesh_.magSf();
                }

                phase.correctInflowOutflow(alphaPhi);

                MULES::limit
                (
                    geometricOneField(),
                    alpha,
                    phase.phi()(),
                    alphaPhi,
                    Sps[phase.index()],
                    Sus[phase.index()],
                    min(alphaVoid.primitiveField(), phase.alphaMax())(),
                    zeroField(),
                    false
                );
            }

            // Limit the flux corrections to ensure the phase fractions sum to 1
            {
                // Generate alphas for the moving phases
                UPtrList<const volScalarField> alphasMoving
                (
                    movingPhases().size()
                );

                UPtrList<surfaceScalarField> alphaPhisMoving
                (
                    movingPhases().size()
                );

                forAll(movingPhases(), movingPhasei)
                {
                    const phaseModel& phase = movingPhases()[movingPhasei];

                    alphasMoving.set(movingPhasei, &phase);

                    alphaPhisMoving.set
                    (
                        movingPhasei,
                        &alphaPhis[phase.index()]
                    );
                }

                MULES::limitSum(alphasMoving, alphaPhisMoving, phiMoving);
            }

            forAll(solvePhases, solvePhasei)
            {
                phaseModel& phase = solvePhases[solvePhasei];
                volScalarField& alpha = phase;

                surfaceScalarField& alphaPhi = alphaPhis[phase.index()];
                phase.correctInflowOutflow(alphaPhi);

                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha,
                    alphaPhi,
                    Sps[phase.index()],
                    Sus[phase.index()]
                );

                if (alphaSubCycle.index() == 1)
                {
                    phase.alphaPhiRef() = alphaPhi;
                }
                else
                {
                    phase.alphaPhiRef() += alphaPhi;
                }
            }

            if (alphaDByAf.valid())
            {
                // Update alphaDByAf due to changes in alpha
                alphaDByAf = this->alphaDByAf(rAs);

                forAll(solvePhases, solvePhasei)
                {
                    phaseModel& phase = solvePhases[solvePhasei];
                    volScalarField& alpha = phase;

                    fvScalarMatrix alphaEqn
                    (
                        fvm::ddt(alpha) - fvc::ddt(alpha)
                      - fvm::laplacian(alphaDByAf(), alpha, "bounded")
                    );

                    alphaEqn.solve();

                    phase.alphaPhiRef() += alphaEqn.flux();
                }
            }

            // TODO: update alphasmax and rescale sum(alfas) if > alfasmax
            // and then set alfac = 1 - sum(alfas)

            // Report the phase fractions and the phase fraction sum
            forAll(solvePhases, solvePhasei)
            {
                phaseModel& phase = solvePhases[solvePhasei];
                fvConstraints().constrain(phase);

                Info<< phase.name() << " fraction mean, min, max = "
                    << phase.weightedAverage(mesh_.V()).value()
                    << ' ' << min(phase).value()
                    << ' ' << max(phase).value()
                    << endl;
            }

            if (referencePhasePtr)
            {
                volScalarField& referenceAlpha = *referencePhasePtr;
                referenceAlpha = alphaVoid;

                forAll(solvePhases, solvePhasei)
                {
                    referenceAlpha -= solvePhases[solvePhasei];
                }
            }
            else
            {
                volScalarField sumAlphaMoving
                (
                    IOobject
                    (
                        "sumAlphaMoving",
                        mesh_.time().name(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar(dimless, 0)
                );
                forAll(movingPhases(), movingPhasei)
                {
                    sumAlphaMoving += movingPhases()[movingPhasei];
                }

                Info<< "Phase-sum volume fraction, min, max = "
                    << (sumAlphaMoving + 1 - alphaVoid)()
                      .weightedAverage(mesh_.V()).value()
                    << ' ' << min(sumAlphaMoving + 1 - alphaVoid).value()
                    << ' ' << max(sumAlphaMoving + 1 - alphaVoid).value()
                    << endl;

                // Correct the sum of the phase fractions to avoid drift
                forAll(movingPhases(), movingPhasei)
                {
                    movingPhases()[movingPhasei] *= alphaVoid/sumAlphaMoving;
                }
            }
        }

        if (nAlphaSubCycles > 1)
        {
            forAll(solvePhases, solvePhasei)
            {
                phaseModel& phase = solvePhases[solvePhasei];

                phase.alphaPhiRef() /= nAlphaSubCycles;
            }
        }

        forAll(solvePhases, solvePhasei)
        {
            phaseModel& phase = solvePhases[solvePhasei];

            phase.alphaRhoPhiRef() =
                fvc::interpolate(phase.rho())*phase.alphaPhi();

            phase.maxMin(0, 1);
        }

        if (referencePhasePtr)
        {
            phaseModel& referencePhase = *referencePhasePtr;

            referencePhase.alphaPhiRef() = phi_;

            forAll(solvePhases, solvePhasei)
            {
                phaseModel& phase = solvePhases[solvePhasei];
                referencePhase.alphaPhiRef() -= phase.alphaPhi();
            }

            referencePhase.alphaRhoPhiRef() =
                fvc::interpolate(referencePhase.rho())
               *referencePhase.alphaPhi();

            volScalarField& referenceAlpha = referencePhase;
            referenceAlpha = alphaVoid;

            forAll(solvePhases, solvePhasei)
            {
                referenceAlpha -= solvePhases[solvePhasei];
            }
        }
    }
}


// ************************************************************************* //
