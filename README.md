# 3D Stable Fluid and More!

An implementation of a 3D fluid simulator based on Jos Stam's [Stable Fluids](https://pages.cs.wisc.edu/~chaol/data/cs777/stam-stable_fluids.pdf) method. The solver is suitable for visual effects involving smoke, gas, or any fluid-like motion.

![Fluid Simulation Demo](fluid/example-meshes/example.gif)

## Features

- **3D Stable Fluid Solver:** Real-time fluid simulation using an unconditionally stable, semi-Lagrangian approach
- **Vorticity Confinement:** Enhanced swirling patterns and detailed flow features
- **Solid Object Collisions:** Fluid interacts with solid obstacles (cubes, spheres) with proper flow dynamics
- **Multiple Visualization Modes:**
  - Volume ray marching for realistic fluid density rendering
  - Voxel-based visualization
  - Shell-only rendering for surface visualization
- **Interactive Controls:** Real-time addition of density and velocity through mouse interaction
- **Customizable Parameters:** Adjust viscosity, diffusion, iterations and more
- **Efficient Implementation:** OpenMP parallelization for improved performance

## Requirements

- C++17 or higher
- OpenGL 3.3+
- GLM for mathematics (included)
- Eigen for vector operations
- CMake 3.10 or higher

## Building the Project

```bash
mkdir build
cd build
cmake ..
make
```

## Running the Simulation 
```bash
./Fruits
```

## User Interaction
 - Left Mouse Drag: Add density and velocity in the direction of movement
 - Right Mouse Drag: Rotate camera
 - WASD: Move camera
 - P: Pause/start simulation
 - V: Increase vorticity
 - B: Decrease vorticity
 - 1/2/3/4: Change different colors

## Parameter Tuning
You can modify these parameters in the GUI or directly in the code:

 - Diffusion: Controls how quickly the density spreads (0.0001 - 0.01)
 - Viscosity: Controls the "thickness" of the fluid (0.0001 - 0.01)
 - Vorticity Strength: Controls the intensity of swirling motion (0.0 - 5.0)
 - Iterations: Controls accuracy of the simulation (5 - 40)

## Physics Implementation Details
This simulation solves the incompressible Navier-Stokes equations using:

 - Velocity Diffusion: Simulates viscous behavior
 - Projection: Enforces incompressibility through pressure solving
 - Advection: Semi-Lagrangian scheme for stable transport
 - Vorticity Confinement: Counteracts numerical dissipation
 - Solid Handling: Ensures fluid flows along solid boundaries
   
## Extensions
The simulation supports multiple interacting features:

Solid Objects: Add static obstacles like cubes and spheres
Multiple Fluid Types: Simulate different fluid types with varying properties
Density Fading: Realistic dissipation over time


## Credits
Based on Jos Stam's "Stable Fluids" paper
Implementation by [Unstable Fruits Team]
Visualization using OpenGL and custom GLSL shaders





