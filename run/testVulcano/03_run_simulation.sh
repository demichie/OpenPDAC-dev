#!/bin/sh

# =============================================================================
# Script 3: Main Simulation Run
# -----------------------------------------------------------------------------
# This script executes the main simulation.
# It assumes that the mesh and fields have been prepared by the previous scripts.
#
# It performs these steps:
#   - Swaps in the '.run' dictionaries for the main simulation.
#   - Runs the solver in parallel.
#   - Reconstructs the results for final visualization.
#
# Usage: ./03_run_simulation.sh
# =============================================================================

# Change to the script's directory for robust execution
cd "${0%/*}" || exit 1

# Source the OpenFOAM functions for running applications
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

# --- MAIN SIMULATION ---

echo "--> Preparing case for the main simulation run..."
# Swap in the final dictionaries for the simulation.
cp ./system/controlDict.run ./system/controlDict
cp ./system/fvSolution.run ./system/fvSolution
cp ./constant/cloudProperties.run ./constant/cloudProperties

echo "--> Starting the main simulation in parallel with $(getApplication)..."
# This command will block until the simulation is finished.
runParallel $(getApplication)

# --- POST-PROCESSING ---

echo "--> Reconstructing the case results..."
# Merges the processor* directories into the main time directories.
runApplication reconstructPar

# -----------------------------------------------------------------------------
echo
echo "SIMULATION SCRIPT COMPLETE."
echo "Results have been reconstructed and are ready for post-processing."
echo "You can now view the final results with 'paraFoam' or ParaView."
# =============================================================================
