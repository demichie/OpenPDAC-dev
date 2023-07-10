# OpenPDAC-dev

OpenPDAC is an OpenFOAM module based on the module multiphaseEuler, distributed with OpenFOAM. 

With respect to the origianl module, in OpenPDAC the equations from the kinetic theory for granualar flows are modified to model multiple dispersed solid phases.

In addition, a lagrangian library is included in the model (one-way coupling with the gas-solid mixture). 

The module also implement an initialization of the hydrostatic pressure profile, which is needed for simulations on large domain. This allows you to use boundary conditions which are appropriate for inflow/outflow.

Four test cases are provided:

- a 3D explosion simulation with a synthetic topography;
- a 2D dilute flow over a wavy surface;
- a 2D fluidezed bed with two solid phases;
- a 3D impinging flow with two solid phases.

This version is based on *openfoam-dev_20230707_amd64.deb*

This code is not approved not endorsed by the OpenFOAM Foundation or by ESI Ltd, the owner of OpenFOAM.
