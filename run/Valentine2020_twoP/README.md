# Column Collapse and Pyroclastic Current Initiation Tutorial

## 1. Case Summary

This tutorial is based on the numerical experiments described in the paper **"Initiation of dilute and concentrated pyroclastic currents from collapsing mixtures and origin of their proximal deposits"** by G. A. Valentine (2020), *Bulletin of Volcanology*, 82:20.

The case models a **2D collapsing column of a gas-particle mixture** impacting a flat surface, similar to the setup shown in Figure 1b of the paper. The primary goal is to explore the dynamics of the "impact zone" and understand the conditions under which the collapsing mixture forms either a **concentrated, granular underflow** or a **dilute, turbulent current** (as illustrated in Figure 2 of the paper).

The simulation features:
-   **One gas phase** (air).
-   **Two distinct solid phases** (`particles1` and `particles2`), representing the "coarse" and "fine" particles discussed in the paper. This allows for the study of particle segregation and differential coupling with the gas phase.
-   A **two-stage simulation process**, managed by swapping configuration files.
    -   **Stage 1: Hydrostatic Initialization.** A preliminary run establishes a stable atmosphere. Using separate `.init` files is necessary because the inlet conditions would otherwise create an incorrect atmospheric profile.
    -   **Stage 2: Main Surge Simulation.** The main run uses different `.run` files to define the inlet conditions.

-    The `Allrun` script executes two separate runs (`.init` and `.run`) for the two stages.

This tutorial is an excellent example of how to use OpenFOAM to replicate and explore the physics described in a scientific publication.

---

## 2. Case Features

-   **Phenomenon:** 2D simulation of a collapsing eruptive column impacting a flat surface.
-   **Physics:** Multiphase flow with one gas and two solid phases (coarse and fine particles).
-   **Methodology:** Two-stages simulation to initialize a stable atmosphere and then to simulate the process. This is managed by swapping dictionary files (`controlDict`, `fvSolution`) and renaming field files in the `0` directory.
-   **Meshing:** A simple 2D domain created with `blockMesh`. `changeDictionary` is used to apply the `empty` boundary condition for 2D simulations.

---

## 3. How to Run the Case

This tutorial can be executed via a single, comprehensive `Allrun` script, which will run the two simulation scenarios sequentially.

### Automated Execution

First, make the scripts executable:
```bash
chmod +x Allrun Allclean
```

To run the entire workflow from start to finish, simply execute:
```bash
./Allrun
```

The script is divided into logical phases:
1.  **Meshing:** The simple 2D mesh is created.
2.  **Initialization Run:** The first simulation is executed using the `.init` configuration files.
3.  **Main Simulation Run:** The second simulation is executed using the `.run` configuration files. The initial conditions in the `0` directory are reset before this run.

---

## 4. Cleaning the Case

To remove all generated data and reset the case to its original state, run the `Allclean` script:
```bash
./Allclean
```

---

## 5. Description of Key Files

-   **`Allrun` / `Allclean`**: Master scripts for running the two simulation scenarios or cleaning the case.
-   **`system/controlDict.init` & `.run`**: Control dictionaries for the two different simulation runs. They likely differ in end time, write intervals, or other run-time parameters.
-   **`system/fvSolution.init` & `.run`**: Solution dictionaries, which might specify different numerical schemes or solver tolerances for each run.
-   **`org.0/`**: A "template" directory for the initial conditions. It contains two complete sets of boundary conditions for each field, distinguished by `.init` and `.run` suffixes. The `Allrun` script renames the appropriate files before each simulation stage. This is a common and effective method for managing multiple case setups within a single directory.
