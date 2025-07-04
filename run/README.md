# Advanced Multiphase Flow Tutorials for OpenPDAC

## 1. Overview

This repository contains a collection of advanced OpenFOAM tutorials designed to demonstrate and test the capabilities of the **OpenPDAC** solver for geophysical and multiphase flows. The cases cover a range of applications, from volcanic explosions to fluidized beds, and showcase sophisticated meshing and simulation workflows.

Each tutorial is self-contained in its own directory and includes all necessary scripts to run the simulation from start to finish.

---

## 2. Prerequisites

Before running these tutorials, please ensure you have the following installed and configured:
-   A working installation of **OpenFOAM**.
-   The **OpenPDAC** solver and its associated custom utilities (e.g., `topoGrid`) must be compiled and accessible in your environment.
-   **Python 3** is required for several preprocessing scripts.
-   **ParaView** (or `paraFoam`) for visualizing the results.

---

## 3. Tutorial Descriptions

Below is a brief description of each tutorial included in this collection.

### 1. `testVulcano` - 3D Phreatic Explosion with Real Topography
This is the most complex case, simulating a phreatic explosion using the real-world topography of Vulcano Island, Italy.
-   **Key Feature:** Uses the custom `topoGrid` utility (from OpenPDAC) to deform the mesh so it conforms to a real ESRI topography file.
-   **Physics:** Simulates an explosion originating from a sub-surface conduit beneath the topographic surface.
-   **Workflow:** Features a detailed, step-by-step execution (`01_...`, `02_...`, etc.) to guide the user through meshing, initialization, and the final run.

### 2. `synthTopo2D` - 2D Explosion with Conduit/Crater Geometry
This tutorial models a 2D explosion using a synthetic, flat topography but with complex sub-surface geometry.
-   **Key Feature:** Demonstrates local mesh refinement using `snappyHexMesh` to capture a distinct surface **crater** and a sub-surface **conduit**.
-   **Physics:** The explosion is triggered by applying overpressure **only within the conduit**, which then erupts through the crater.
-   **Workflow:** Uses a two-stage initialization process: first to establish a stable atmosphere, then to set the explosion conditions with `setFields`.

### 3. `flow_wave` - 2D Pyroclastic Surge over Topography
This case simulates a 2D pyroclastic surge flowing over a synthetic, non-flat topography.
-   **Key Feature:** Showcases an advanced technique to create a 2D mesh for complex geometries: a full 3D mesh is first generated with `snappyHexMesh` and then a 2D slice is extracted using `extrudeMesh`.
-   **Physics:** Models a high-energy flow over a landscape, with a robust two-stage initialization to correctly establish the atmospheric profile before introducing the surge.
-   **Workflow:** Ideal for learning how to handle 2D simulations with complex boundaries in OpenFOAM.

### 4. `Valentine2020_twoP` - 2D Column Collapse (Valentine, 2020)
This tutorial is a direct numerical replication of the experiments described in the scientific paper by G. A. Valentine (2020), *Bulletin of Volcanology*.
-   **Key Feature:** Models a collapsing column of a gas-particle mixture with **two solid phases** (coarse and fine particles).
-   **Physics:** Explores the conditions that lead to the formation of either dilute, turbulent currents or concentrated, granular underflows upon impact with the ground.
-   **Workflow:** Uses a two-stage run to test different simulation parameters by swapping configuration and field files.

### 5. `fluidisedBed` - Polydisperse Fluidized Bed
This is a classic fluid mechanics case, adapted to test polydisperse models.
-   **Key Feature:** A modification of the standard OpenFOAM `fluidisedBed` tutorial, but with **two distinct solid phases**.
-   **Physics:** Designed specifically to test and validate the **polydisperse kinetic theory models** implemented in OpenPDAC.
-   **Workflow:** A standard, single-run case that demonstrates phenomena like particle segregation and bubbling in a polydisperse system.

---

## 4. How to Run a Tutorial

All tutorials follow a similar execution pattern.

1.  Navigate into the directory of the tutorial you wish to run:
    ```bash
    cd name_of_the_tutorial
    [CODE_END]
2.  Make the scripts executable:
    ```bash
    chmod +x Allrun Allclean
    [CODE_END]
3.  Execute the main `Allrun` script to run the entire simulation from start to finish:
    ```bash
    ./Allrun
    [CODE_END]
4.  For more complex cases like `testVulcano`, you can also run the step-by-step scripts (`./01_run_...`, etc.) to inspect the results at each stage.
5.  To clean the case and remove all generated data, run:
    ```bash
    ./Allclean
    [CODE_END]

---

## 5. Author / Contact

[Your Name / Institution]
[your.email@example.com]
