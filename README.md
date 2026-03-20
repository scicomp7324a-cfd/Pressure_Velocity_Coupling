# Pressure_Velocity_Coupling
PressureVelocityCoupling
========================

Overview
--------
This project contains a C++ implementation of a modified SIMPLE-based
pressure–velocity coupling solver for the two-dimensional incompressible
lid-driven cavity problem.

The code is written using a finite-volume formulation on a face-addressed
Cartesian mesh. The implementation includes:
- mesh and boundary handling,
- sparse matrix and sparse linear-system support,
- Gauss–Seidel smoothing and solving,
- discretisation of convection, diffusion, transient, and source terms,
- pressure correction and face-flux correction,
- automated runs for selected mesh sizes and Reynolds numbers,
- output files for convergence monitoring and post-processing.

How to run
----------
Run the code from the project root using:

    ./run.sh

This is the intended way to build and execute the project.

What run.sh does
----------------
The script:
1. removes any existing build directory,
2. configures the project with CMake,
3. builds the executable,
4. enters the build directory,
5. runs the solver executable.

So the normal workflow is simply:

    ./run.sh

Requirements
------------
The project requires:
- a C++17-compatible compiler,
- CMake 3.16 or later,
- a Unix-like shell environment.

Project structure
-----------------
The main folders are:

- Mesh/
  Mesh data structures and mesh readers.

- SparseOperations/
  Sparse address, sparse matrix, sparse vector, and sparse linear-system classes.

- FVCore/
  Core finite-volume data structures such as fields, boundary conditions,
  gradients, flux fields, pressure Laplacian, and under-relaxation.

- Discretisation/
  Discretisation components for convection, diffusion, transient, and source terms.

- Solver/
  Gauss–Seidel smoother and solver.

- Files/
  Input files required by the code.
  This contains:
  - FieldFiles/
    Initial field data for Ux, Uy, and pressure.
  - MeshFiles/
    Mesh definitions for different mesh sizes.

Important run information
-------------------------
The main program runs the solver for the following mesh sizes:

- 20 x 20
- 40 x 40
- 60 x 60
- 80 x 80

For each mesh, it runs the following Reynolds numbers:

- 100
- 400
- 1000

Output
------
During the build, the Files/ directory is copied into the build directory.
An Output/ directory is also created automatically.

The solver writes case outputs into the Output/ folder, organised by
mesh size and Reynolds number.

Main source files
-----------------
- main.cpp
  Main driver that loops over all required mesh-size and Reynolds-number cases.

- run.sh
  Recommended script for building and running the whole project.

- CMakeLists.txt
  CMake build configuration.

Notes
-----
- The project should be run from the PressureVelocityCoupling project root.
- Any macOS metadata files such as .DS_Store can be ignored.

Summary
-------
This project is a modular finite-volume CFD code for studying pressure–velocity
coupling through a modified SIMPLE procedure for the lid-driven cavity problem.
The recommended way to execute it is:

    ./run.sh
