#!/bin/sh

# =============================================================================
# Allclean Script for the Polydisperse Fluidized Bed Tutorial
# -----------------------------------------------------------------------------
# This script removes all generated files, restoring the case to its
# original, clean state.
#
# Usage: ./Allclean
# =============================================================================

cd "${0%/*}" || exit 1

# Use the standard OpenFOAM utility to remove logs, time steps,
# processor directories, etc.
foamCleanCase

# Explicitly remove the '0' directory to be certain of a clean start
rm -rf 0

echo "Clean-up complete."
