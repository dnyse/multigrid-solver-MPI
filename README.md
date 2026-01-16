# Multigrid Solver MPI

A 3D Navier-Stokes CFD solver with MPI parallelization, developed at NHR@FAU, University Erlangen-Nuremberg.

## Features

- 3D incompressible Navier-Stokes equations
- Pressure Poisson equation solvers:
  - Red-Black SOR (`rb`)
  - Multigrid (`mg`)
- MPI domain decomposition with 3D Cartesian topology
- VTK output for visualization (ParaView compatible)
- LIKWID marker support for performance analysis

## Project Structure

```
├── 3D-seq/          # Sequential version
│   ├── src/         # Source files
│   ├── *.par        # Parameter files (canal, dcavity)
│   └── Makefile
├── 3D-mpi/          # MPI parallel version
│   ├── src/         # Source files
│   ├── *.par        # Parameter files
│   └── Makefile
└── README.md
```

## Building

### Prerequisites

- C compiler (GCC, Clang, or Intel ICX)
- MPI implementation (OpenMPI, Intel MPI, etc.)
- Optional: LIKWID for performance monitoring

### Configuration

Edit `config.mk` to configure the build:

```makefile
TAG ?= ICX              # Compiler: GCC, CLANG, ICX
ENABLE_MPI ?= true      # Enable MPI parallelization
ENABLE_OPENMP ?= false  # Enable OpenMP (hybrid)
SOLVER ?= mg            # Solver: rb (Red-Black SOR), mg (Multigrid)
DEBUG ?= false          # Debug build
```

### Compile

```bash
cd 3D-mpi
make
```

This creates the executable `exe-<TAG>` (e.g., `exe-ICX`).

### Clean

```bash
make clean      # Remove build artifacts
make distclean  # Remove build artifacts and executable
```

## Running

### Sequential

```bash
cd 3D-seq
./exe-<TAG> canal.par
```

### MPI Parallel

```bash
cd 3D-mpi
mpirun -np <num_procs> ./exe-<TAG> canal.par
```

### With LIKWID Performance Monitoring

```bash
likwid-mpirun -np 4 -g FLOPS_DP ./exe-ICX canal.par
```

## Parameter Files

Parameter files (`.par`) configure the simulation. Key parameters:

| Parameter | Description |
|-----------|-------------|
| `name` | Simulation name |
| `xlength`, `ylength`, `zlength` | Domain dimensions |
| `imax`, `jmax`, `kmax` | Grid cells in each direction |
| `re` | Reynolds number |
| `te` | End time |
| `dt` | Time step size |
| `tau` | Safety factor for adaptive time stepping |
| `itermax` | Max pressure solver iterations |
| `eps` | Convergence tolerance |
| `omg` | SOR relaxation parameter |
| `levels` | Multigrid levels |
| `presmooth`, `postsmooth` | Multigrid smoothing iterations |

### Boundary Conditions

```
bcLeft, bcRight, bcBottom, bcTop, bcFront, bcBack
  1 = no-slip
  2 = free-slip
  3 = outflow
  4 = periodic
```

### Example: Laminar Canal Flow

```
name canal
bcLeft   3    bcRight  3
bcBottom 1    bcTop    1
bcFront  1    bcBack   1
re       100.0
xlength  30.0
ylength  4.0
zlength  4.0
imax     200
jmax     50
kmax     50
te       40.0
```

## Output

The solver outputs VTK files for visualization in ParaView:
- Pressure field
- Velocity vector field (u, v, w components)

## License

MIT License - see LICENSE file.

## Acknowledgments

Copyright (C) 2022-2024 NHR@FAU, University Erlangen-Nuremberg.
