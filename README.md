# Conjugate Gradient Solver for 2D Laplace Equation (MPI Parallelized)

This repository provides an MPI-parallelized implementation of the conjugate gradient (CG) method to solve the linear system derived from discretizing the 2D Laplace equation using the finite difference method. The computation is performed in a Cartesian coordinate system with Dirichlet boundary conditions.

## Key Features
- **Finite Difference Discretization**: Discretizes the 2D Laplace equation using a uniform grid.
- **Conjugate Gradient Method**: Solves the generated sparse linear system efficiently.
- **MPI Parallelization**: Distributes the computation across multiple processes in a 2D Cartesian topology using MPI.

## Requirements
- An MPI library (e.g., MPICH, OpenMPI)
- A C/C++ compiler with MPI support (e.g., `mpic++`)

## Compilation
Use the provided `Makefile` or compile directly with:
```bash
mpic++ -std=c++17 main.cpp Grid.cpp -o cgsolver
```

## How to Run
The executable is named `cgsolver`. Use the following syntax to run the program:
```bash
mpirun -np N ./cgsolver nx ny iter np_x np_y [write_option]
```
### Parameters:
1. `N`: Number of MPI processes.
2. `nx`: Number of grid points in the x-direction.
3. `ny`: Number of grid points in the y-direction.
4. `iter`: Maximum number of iterations for the CG solver.
5. `np_x`: Number of MPI processes in the x-direction.
6. `np_y`: Number of MPI processes in the y-direction.
7. `[write_option]` (optional): Provide any value to enable saving the solution to a file.

### Example:
```bash
mpirun -np 4 ./cgsolver 100 100 1000 2 2 write
```
This example runs the solver on a 100x100 grid for up to 1000 iterations, using a 2x2 MPI process grid. The solution will be saved to `solution.txt` if the `write_option` is provided.

### Output:
- The results are output to the terminal. If `write_option` is enabled, the computed solution for the entire domain is written to `solution.txt`.

## Important Notes
- Ensure the product of `np_x` and `np_y` matches the total number of MPI processes specified with `N` in the `mpirun` command.
- The problem is solved in a two-dimensional domain [0, 1]x[0, 1], with initial conditions set to 0 at all points. The implementation assumes Dirichlet boundary conditions with predefined boundary values. The L2 norm is used to calculate the error, the residual, and to determine the criteria for early stopping. 