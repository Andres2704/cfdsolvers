
# CFD Solvers in Julia

This repository is dedicated to implementing various **Computational Fluid Dynamics (CFD)** solvers in Julia. Each solver is designed to be as modular, readable, and educational as possible. The project begins with a classic test case: **2D lid-driven cavity flow**, solved using a finite volume scheme with **Rhie-Chow interpolation** for pressure-velocity coupling on a collocated grid and also **staggered grids**.

## Current Solver: 2D Lid-Driven Cavity Flow (Rhie-Chow)

### Features

- Collocated grid layout
- Rhie-Chow interpolation to avoid pressure-velocity decoupling
- Explicit time advancement
- Pressure Poisson equation solved using **SOR**
- Configurable physical and numerical parameters
- Output fields exported as `.txt` for easy visualization

---

## File Structure

```
.
├── main_simulation_rhiechow.jl         # Main script to configure and run the simulation
├── cavity_functions_rhiechow.jl        # Module with all solver functions
├── *.txt                               # Output result files
└── README.md
```

---

## How to Run

### 1. Requirements

Make sure you have [Julia](https://julialang.org/downloads/) installed (version 1.6+ recommended).

No external libraries are required, only the standard library `DelimitedFiles`.

### 2. Run the Simulation

```bash
julia main_simulation_rhiechow.jl
```

This will:
- Set up the cavity domain
- Run the time-stepping solver until the specified simulation time
- Export the velocity, pressure, and divergence fields to `.txt` files

---

## Configurable Parameters

You can change the simulation setup directly in `main_simulation_rhiechow.jl`:

```julia
# Mesh resolution
NX = 32
NY = 32

# Physical properties
RHO = 1.0         # Density [kg/m³]
MU  = 0.01        # Viscosity [Pa.s]

# Lid velocity (top wall)
ut = 1.0

# Solver settings
CFL = 0.1
END_TIME = 4.0
SOR_TOL = 1e-11
SOR_BETA = 1.872
MAXSOR_ITER = 300
```

---

## Output Files

After simulation, the following files will be generated:

| File Name                        | Description                     |
|----------------------------------|----------------------------------|
| `u_velocity_field_rhiechow.txt` | u-component of velocity field   |
| `v_velocity_field_rhiechow.txt` | v-component of velocity field   |
| `pressure_field_rhiechow.txt`   | Pressure field                  |
| `divergence_field_rhiechow.txt` | Divergence (∇·u) field          |
| `x_coordinate_rhiechow.txt`     | x-grid coordinates              |
| `y_coordinate_rhiechow.txt`     | y-grid coordinates              |
| `X_coordinate_field_rhiechow.txt` | 2D meshgrid for X              |
| `y_coordinate_field_rhiechow.txt` | 2D meshgrid for Y              |

These files can be post-processed with Python (e.g., using `matplotlib`), Julia (e.g., `Plots.jl`), MATLAB, or ParaView after format conversion.

---
