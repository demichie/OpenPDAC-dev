#!/bin/sh
# This script runs all steps of the tutorial in sequence.
set -e # Exit immediately if a command exits with a non-zero status.

echo "--- Starting Step 1: Meshing ---"
./01_run_meshing.sh

echo "\n--- Starting Step 2: Field Initialization ---"
./02_run_fieldInitialization.sh

echo "\n--- Starting Step 3: Main Simulation ---"
./03_run_simulation.sh

echo "\n--- All steps completed successfully! ---"
