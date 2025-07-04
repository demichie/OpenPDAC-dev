# 2D Synthetic Explosion Simulation with Conduit/Crater Geometry

## 1. Case Summary

This tutorial demonstrates a **2D simulation of a volcanic explosion** featuring a distinct sub-surface **conduit** and a surface **crater**. The case showcases a sophisticated workflow for meshing complex geometries and a precise, multi-stage initialization process.

The workflow is structured as follows:
1.  **2D Mesh Generation with Local Refinement:** A 2D mesh is generated using an advanced technique. First, `snappyHexMesh` creates a 3D mesh, refining it around the STL geometries of both the surface crater and the ground level. Then, `extrudeMesh` is used to extract a 2D slice, resulting in a one-cell-thick domain suitable for 2D simulation in OpenFOAM.
2.  **Two-Stage Initialization:** The simulation is initialized in two steps for physical accuracy.
    *   **Stage 1: Hydrostatic Initialization.** A preliminary run establishes a stable, quiescent atmosphere across the entire domain.
    *   **Stage 2: Explosion Trigger in the Conduit.** `topoSet` is used to select the cells corresponding to the sub-surface **conduit**. Then, `setFields` applies high pressure and temperature **only to this conduit region**, setting the trigger for the explosion.
3.  **Main Simulation:** The final simulation captures the dynamics of the high-pressure mixture erupting from the conduit, passing through the crater, and expanding into the atmosphere.

---

## 2. Case Features

-   **Phenomenon:** 2D simulation of a volcanic explosion propagating from a sub-surface conduit through a surface crater.
-   **Geometry:** Synthetic flat topography featuring a distinct surface crater and a deeper sub-surface conduit.
-   **Meshing:**
    -   `snappyHexMesh`: Creates a 3D mesh with local refinement around the crater and ground-level STL geometries.
    -   `extrudeMesh`: Extracts a 2D slice from the 3D mesh to create the final computational domain.
    -   `topoSet`: Crucially, it defines a cell set for the **conduit**, which is the source of the explosion.
-   **Initialization:**
    -   A robust two-stage process: first, a run to establish a hydrostatic atmosphere; second, `setFields` applies overpressure **only within the conduit**, not the crater.

---

## 3. How to Run the Case

This tutorial can be executed via a single, comprehensive `Allrun` script.

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
1.  **Meshing:** All steps related to 2D mesh generation are performed. After this phase, inspect the mesh in ParaView to verify the refinement around the crater and the correct geometry.
2.  **Initialization Run:** The first simulation establishes the hydrostatic atmosphere. Then, `topoSet` and `setFields` apply the explosion source conditions to the **conduit**.
3.  **Main Simulation Run:** The final simulation of the explosion is executed.

---

## 4. Cleaning the Case

To remove all generated data and reset the case to its original state, run the `Allclean` script:
```bash
./Allclean
```

---

## 5. Description of Key Files

-   **`Allrun` / `Allclean`**: Master scripts for running the entire case or cleaning it.
-   **`preprocessing/ASCtoSTL.py`**: A Python script that generates the necessary STL geometry files for the crater and conduit.
-   **`system/snappyHexMeshDict`**: Controls `snappyHexMesh`. This is key for defining the refinement levels around the crater geometry.
-   **`system/extrudeMeshDict`**: Controls `extrudeMesh`, which creates the final one-cell-thick domain for the 2D simulation.
-   **`system/topoSetDict-conduit`**: A dedicated dictionary used by `topoSet` to create a cell set named `conduit`. This set isolates the **sub-surface conduit cells** that will serve as the explosion source.
-   **`system/setFieldsDict`**: Defines the high-pressure and high-temperature values that are applied by `setFields` **only to the `conduit` cell set**.
-   **`system/controlDict.init` & `.run`**: Control dictionaries for the initialization and main simulation phases.
-   **`org.0/`**: A template directory for the initial conditions.
