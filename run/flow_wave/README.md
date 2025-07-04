# 2D Pyroclastic Surge Simulation Tutorial

## 1. Case Summary

This tutorial demonstrates how to set up and run a **2D simulation of a pyroclastic surge** flowing over a synthetic topography. The case showcases an advanced meshing technique and a robust two-stage initialization process, which are common in environmental and geophysical flows.

The workflow is structured as follows:
1.  **2D Mesh Generation:** The process uses a standard OpenFOAM technique to create a 2D computational domain. First, `snappyHexMesh` generates a fully 3D mesh, refining it around the topography. Then, `extrudeMesh` extracts a 2D slice from this 3D grid and extrudes it by a single cell. The result is a valid 3D mesh that is one cell thick, which is how OpenFOAM handles 2D simulations (by using `empty` boundary conditions on the front and back patches).
2.  **Two-Stage Initialization:**
    *   **Stage 1: Atmospheric Initialization.** A preliminary run is executed using specific `.init` boundary conditions. This is crucial for establishing a stable, hydrostatic pressure profile for the atmosphere *before* introducing the surge.
    *   **Stage 2: Main Simulation.** The boundary conditions are then swapped to the `.run` files, which define the high-energy inlet conditions for the pyroclastic surge. The main simulation is then executed. This two-step process prevents the high-energy inlet from disrupting the initial atmospheric equilibrium.

This case is an excellent example of advanced meshing for 2D cases and a robust initialization procedure for atmospheric flows.

---

## 2. Case Features

-   **Phenomenon:** 2D simulation of a pyroclastic surge over synthetic topography.
-   **Meshing:**
    -   `python`: Automates the creation of the input STL geometry for the topography.
    -   `snappyHexMesh`: Generates an initial 3D mesh refined around the topography.
    -   `extrudeMesh`: Extracts a 2D slice and extrudes it by one cell to create the final 2D computational domain.
-   **Simulation:**
    -   **Stage 1: Hydrostatic Initialization.** A preliminary run establishes a stable atmosphere. Using separate `.init` files is necessary because the final surge inlet conditions would otherwise create an incorrect atmospheric profile.
    -   **Stage 2: Main Surge Simulation.** The main run uses different `.run` files to define the energetic inlet conditions of the surge.

---

## 3. How to Run the Case

This tutorial can be executed via a single, comprehensive `Allrun` script.

### Automated Execution

First, make the scripts executable:
``bash
chmod +x Allrun Allclean
```

To run the entire workflow from start to finish, simply execute:
```bash
./Allrun
```

The script is divided into logical phases:
1.  **Meshing:** All steps related to 2D mesh generation are performed. After this phase, you can inspect the final mesh in ParaView to verify its structure.
2.  **Initialization Run:** The first simulation is executed to establish the hydrostatic atmospheric profile.
3.  **Main Simulation Run:** The final simulation of the pyroclastic surge is executed.

---

## 4. Cleaning the Case

To remove all generated data and reset the case to its original state, run the `Allclean` script:
```bash
./Allclean
```

---

## 5. Description of Key Files

-   **`Allrun` / `Allclean`**: Master scripts for running the entire case or cleaning it.
-   **`preprocessing/createSTL.py`**: A Python script that generates the `geometry.stl` file representing the synthetic topography.
-   **`system/snappyHexMeshDict`**: Controls the `snappyHexMesh` process, creating a 3D refined mesh.
-   **`system/extrudeMeshDict`**: Controls the `extrudeMesh` utility, which extracts a 2D slice from the 3D snappy-generated mesh and creates the final one-cell-thick domain for the 2D simulation.
-   **`system/controlDict.init` & `.run`**: Control dictionaries for the initialization and main simulation phases, respectively. They may differ in timings, write controls, etc.
-   **`system/fvSolution.init` & `.run`**: Solution dictionaries, potentially with different solvers or tolerances optimized for the quasi-static initialization versus the highly transient surge simulation.
-   **`org.0/`**: A "template" directory for the initial conditions. It contains two sets of boundary conditions for each field, distinguished by suffixes.
    -   The **`.init`** files are used to create a stable hydrostatic atmosphere.
    -   The **`.run`** files define the high-energy inlet conditions for the pyroclastic surge itself.
