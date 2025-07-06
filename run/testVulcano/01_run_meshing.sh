#!/bin/sh

# =============================================================================
# Script 1: Mesh Generation
# -----------------------------------------------------------------------------
# This script performs all meshing operations:
#   - Cleans the case
#   - Creates a base mesh with blockMesh
#   - Decomposes the domain for parallel processing
#   - Modifies and deforms the mesh using topoSet and topoGrid
#
# After running, you can inspect the mesh quality and view it in ParaView.
#
# Usage: ./01_run_meshing.sh
# =============================================================================

# Change to the script's directory for robust execution
cd "${0%/*}" || exit 1

# Source the OpenFOAM functions for running applications
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

# --- MESHING ---

echo "--> Cleaning the case from previous runs..."
# Use the Allclean script to ensure a fresh start
./Allclean

echo "--> Running Python script for geometry preparation..."
python3 smoothCraterArea.py

cp ./system/controlDict.init ./system/controlDict

echo "--> Creating the base mesh with blockMesh..."
runApplication blockMesh

echo "--> Setting up initial conditions from 'org.0' directory..."
rm -rf 0
cp -r org.0 0

echo "--> Decomposing the domain for parallel execution..."
runApplication decomposePar

echo "--> Performing initial mesh quality check in parallel..."
runParallel checkMesh -allGeometry -allTopology -writeSets

mv log.checkMesh log.checkMesh0

echo "--> Applying topological changes to the mesh with topoSet..."
runParallel topoSet

echo "--> Deforming the mesh with topoGrid..."
runParallel topoGrid

echo "--> Performing final mesh quality check on the deformed mesh..."
runParallel checkMesh -allGeometry -allTopology -writeSets

# -----------------------------------------------------------------------------
echo
echo "MESHING SCRIPT COMPLETE."
echo "The mesh has been generated in the processor* directories."
echo "To view the decomposed mesh, open ParaView and open the case by selecting"
echo "the file with the '.foam' extension (e.g., case.foam)."
echo "Next, run ./02_run_fieldInitialization.sh"
# =============================================================================
