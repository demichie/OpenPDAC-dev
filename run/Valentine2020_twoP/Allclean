#!/bin/sh

# =============================================================================
# Allclean Script for the Column Collapse Tutorial
# -----------------------------------------------------------------------------
# This script removes all generated files, restoring the case to its
# original, clean state.
#
# Usage: ./Allclean
# =============================================================================

cd "${0%/*}" || exit 1

# Use the standard OpenFOAM utility to remove logs, time steps, etc.
foamCleanTutorials

# Explicitly remove the '0' directory
rm -rf 0

# Remove dictionaries that were copied by the Allrun script
rm -f system/controlDict
rm -f system/fvSolution

echo "Clean-up complete."
