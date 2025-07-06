#!/bin/sh
# This script cleans the case, restoring it to its original state.

cd "${0%/*}" || exit 1

# Use the standard OpenFOAM utility to remove logs, time steps, etc.
foamCleanCase

# Explicitly remove the '0' directory
rm -rf 0
rm output_file.asc
rm modified_points.png

# Remove dictionaries that were copied by the run scripts
rm -f system/controlDict
rm -f system/fvSolution
rm -f constant/cloudProperties

echo "Clean-up complete."
