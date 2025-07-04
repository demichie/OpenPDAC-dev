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

echo "Clean-up complete."
