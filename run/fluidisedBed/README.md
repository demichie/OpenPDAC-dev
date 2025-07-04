# Polydisperse Fluidized Bed Simulation Tutorial

## 1. Case Summary

This tutorial demonstrates a **2D simulation of a fluidized bed with a polydisperse granular mixture**. Specifically, this case is a modification of the standard `fluidisedBed` tutorial available in the OpenFOAM repository, but it has been adapted to include **two distinct solid phases** (`particles1` and `particles2`) instead of one.

The primary goal of this tutorial is to **test and demonstrate the capabilities of the polydisperse kinetic theory models for granular flows**, as implemented in the **OpenPDAC** solver.

The simulation setup is as follows:
1.  A bed of two different types of particles is initialized at the bottom of a rectangular domain.
2.  A gas (air) is injected from an inlet at the bottom of the domain.
3.  As the gas flows upward through the bed, it exerts drag on the particles. When the drag force is sufficient to counteract the weight of the particles, the bed "fluidizes," behaving like a fluid.
4.  The simulation captures the mixing, segregation, and bubbling phenomena characteristic of polydisperse fluidized beds.

---

## 2. Case Features

-   **Phenomenon:** Fluidization of a polydisperse (two-phase) granular bed.
-   **Solver:** `OpenPDAC` (or a similar multiphase solver capable of handling polydisperse kinetic theory).
-   **Meshing:** A simple rectangular domain created with `blockMesh`.
-   **Initialization:** The `setFields` utility is used to define the initial location and volume fractions of both `particles1` and `particles2`.
-   **Key Physics:** This case allows for the study of phenomena like particle segregation (where different particle types separate due to differences in size or density) and bubbling, all governed by the polydisperse granular flow models.

---

## 3. How to Run the Case

This tutorial is designed to be run in parallel.

### Automated Execution

First, make the scripts executable:
```bash
chmod +x Allrun Allclean
```

To run the entire workflow from start to finish, simply execute:
```bash
./Allrun
```

The `Allrun` script will automatically:
1.  Create the mesh (`blockMesh`).
2.  Set the initial particle bed (`setFields`).
3.  Decompose the case for parallel running (`decomposePar`).
4.  Run the simulation (`runParallel`).
5.  Reconstruct the results for easy visualization (`reconstructPar`).

---

## 4. Cleaning the Case

To remove all generated data and reset the case to its original state, run the `Allclean` script:
```bash
./Allclean
```

---

## 5. Description of Key Files

-   **`Allrun` / `Allclean`**: Master scripts for running the case or cleaning it.
-   **`system/setFieldsDict`**: This is a critical file for this tutorial. It defines a `box` region at the bottom of the domain and specifies the initial volume fractions for `alpha.particles1` and `alpha.particles2` within that box.
-   **`system/decomposeParDict`**: Defines the number of processors and the method used for parallel decomposition.
-   **`constant/kineticTheoryProperties`**: Another crucial file. This is where the physical properties of the two solid phases are defined, such as their diameter, density, and coefficient of restitution. The differences between `particles1` and `particles2` in this file will drive the polydisperse behavior.
-   **`0/U`**: Defines the inlet velocity of the fluidizing gas at the bottom boundary.
