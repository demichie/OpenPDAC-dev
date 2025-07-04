#!/bin/sh

# =============================================================================
# Allclean Script for the Tutorial Case
# -----------------------------------------------------------------------------
# This script removes all generated files, restoring the case to its
# original, clean state.
#
# Usage: ./Allclean
# =============================================================================

cd "${0%/*}" || exit 1

# Use the standard OpenFOAM utility to remove logs, time steps, etc.
foamCleanCase

# Explicitly remove the '0' directory
rm -rf 0

# Remove dictionaries that were copied by the Allrun script
rm -f system/controlDict
rm -f system/fvSolution

# Clean up files generated in the preprocessing directory
rm -f preprocessing/*.stl
rm -f preprocessing/log.*
rm -rf preprocessing/__pycache__

echo "Clean-up complete."
