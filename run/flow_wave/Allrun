#!/bin/sh

# =============================================================================
# Allrun Script for the Tutorial Case
# -----------------------------------------------------------------------------
# This script orchestrates the entire simulation workflow:
# 1. Cleans the case directory.
# 2. Generates the mesh via a multi-step process:
#    - Python script to create STL geometry.
#    - blockMesh for a background mesh.
#    - snappyHexMesh to conform to the STL.
#    - extrudeMesh to create a 2D domain.
# 3. Runs a two-stage simulation:
#    - An initialization run to set up initial fields.
#    - The main simulation run.
#
# Usage: ./Allrun
# =============================================================================

# Change to the script's directory for robust execution
cd "${0%/*}" || exit 1

# Source the OpenFOAM functions for running applications
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"


# --- PHASE 0: CLEANING ---

echo "--> Cleaning the case from previous runs..."
./Allclean.sh


# --- PHASE 1: MESHING ---
# In this phase, we generate the 2D computational mesh.

echo "--> Running Python script to generate STL geometry..."
(
    cd preprocessing || exit 1
    python3 createSTL.py
)

cp ./system/controlDict.init ./system/controlDict

echo "--> Creating the background mesh with blockMesh..."
runApplication blockMesh

echo "--> Performing initial mesh quality check..."
runApplication checkMesh -allTopology -allGeometry

mv log.checkMesh log.checkMesh0

echo "--> Conforming mesh to STL geometry with snappyHexMesh..."
runApplication snappyHexMesh -overwrite

echo "--> Extruding the mesh to create a 2.5D domain..."
runApplication extrudeMesh

echo "--> Performing final mesh quality check..."
runApplication checkMesh -allTopology -allGeometry

echo "--> Modifying dictionary entries with changeDictionary..."
runApplication changeDictionary


# --- PHASE 2: INITIALIZATION RUN ---
# Here, we run a special setup to initialize the fields.

echo "--> Preparing for the initialization run..."
# Set up the dictionaries for the initialization phase
cp ./system/controlDict.init ./system/controlDict
cp ./system/fvSolution.init ./system/fvSolution

# Copy base fields from org.0 and rename them for initialization
rm -rf 0
cp -r org.0 0
echo "--> Setting up fields for initialization..."
for field in alpha.air alpha.particles T.air T.particles U.air U.particles; do
    mv "0/${field}.init" "0/${field}"
done

echo "--> Starting the initialization run..."
# getApplication reads the solver name from the system/controlDict
runApplication $(getApplication)

mv log.foamRun log.foamRun0

# --- PHASE 3: MAIN SIMULATION RUN ---
# Now, we run the primary simulation with the final settings.

echo "--> Preparing for the main simulation run..."
# Set up the dictionaries for the main simulation phase
cp ./system/controlDict.run ./system/controlDict
cp ./system/fvSolution.run ./system/fvSolution

# Rename the fields in the '0' directory for the main run.
# This overwrites the initial state with a new set of conditions.
echo "--> Setting up fields for the main run..."
for field in alpha.air alpha.particles T.air T.particles U.air U.particles; do
    mv "0/${field}.run" "0/${field}"
done

echo "--> Starting the main simulation..."
# The solver will likely start from time 0 again, using the new fields.
runApplication $(getApplication)

# -----------------------------------------------------------------------------
echo
echo "Allrun script finished successfully."
# =============================================================================
