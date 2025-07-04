#!/bin/sh

# =============================================================================
# Allrun Script for the Column Collapse Tutorial
# -----------------------------------------------------------------------------
# This script orchestrates the entire simulation workflow, based on the
# numerical experiments described in Valentine (2020).
# It models the collapse of a gas-particle mixture to study the initiation
# of pyroclastic currents.
#
# The workflow includes:
# 1. Cleaning the case and setting up the mesh.
# 2. Running a first simulation with an ".init" configuration.
# 3. Running a second simulation with a ".run" configuration.
#
# Usage: ./Allrun
# =============================================================================

# Change to the script's directory for robust execution
cd "${0%/*}" || exit 1

# Source the OpenFOAM functions for running applications
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"


# --- PHASE 0: CLEANING & MESHING ---

echo "--> Cleaning the case from previous runs..."
./Allclean

echo "--> Creating the background mesh with blockMesh..."
runApplication blockMesh

echo "--> Performing mesh quality check..."
runApplication checkMesh -allTopology -allGeometry

echo "--> Setting 2D empty boundary conditions..."
runApplication changeDictionary


# --- PHASE 1: INITIALIZATION RUN ---
# This stage runs the simulation with the first set of parameters, defined
# by the ".init" configuration files.

echo "--> Preparing for the initialization run..."
# Set up the dictionaries for this phase
cp ./system/controlDict.init ./system/controlDict
cp ./system/fvSolution.init ./system/fvSolution

# Copy base fields from org.0 and rename them for the ".init" run
rm -rf 0
cp -r org.0 0
echo "--> Setting up fields for the '.init' run..."
for field in alpha.air T.air U.air; do
    mv "0/${field}.init" "0/${field}"
done
for particle in particles1 particles2; do
    for field in alpha T U; do
        mv "0/${field}.${particle}.init" "0/${field}.${particle}"
    done
done

echo "--> Starting the initialization run (foamRun0)..."
runApplication $(getApplication)


# --- PHASE 2: MAIN SIMULATION RUN ---
# This stage runs the simulation with the second set of parameters, defined
# by the ".run" configuration files. This likely represents a different
# physical scenario (e.g., different particle concentration or properties).

echo "--> Preparing for the main simulation run..."
# Set up the dictionaries for this phase
cp ./system/controlDict.run ./system/controlDict
cp ./system/fvSolution.run ./system/fvSolution

# Rename the fields in the '0' directory for the ".run" configuration.
# NOTE: This setup overwrites the previous initial conditions in the 0 folder.
echo "--> Setting up fields for the '.run' run..."
for field in alpha.air T.air U.air; do
    mv "0/${field}.run" "0/${field}"
done
for particle in particles1 particles2; do
    for field in alpha T U; do
        mv "0/${field}.${particle}.run" "0/${field}.${particle}"
    done
done

echo "--> Starting the main simulation (foamRun1)..."
runApplication $(getApplication)

# -----------------------------------------------------------------------------
echo
echo "Allrun script finished successfully."
# =============================================================================
