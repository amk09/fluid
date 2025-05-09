# 3D Stable Fluid and More!

An implementation of a 3D fluid simulator based on Jos Stam's [Stable Fluids](https://pages.cs.wisc.edu/~chaol/data/cs777/stam-stable_fluids.pdf) method. The solver is suitable for visual effects involving smoke, gas, or any fluid-like motion.

![Fluid Simulation Demo](./example-meshes/example.gif)

*Real-time fluid simulation with default color.*

![Fluid Simulation Demo](./example-meshes/olexample.gif)

*Offline precomputed fluid simulation with default color.*

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
- **Color Mixture:** Different fluid colors can appear at the same time
- **Obstacle Detection:** Obstacles can be placed freely in the domain, and the fluid automatically avoids diffusing into them.
- **Higher Resolution Computation:** Offline precomputing feature for higher resolution, makes it up to 256 cubes

## Requirements

- C++17 or higher
- OpenGL 3.3+
- GLM for mathematics (included)
- Eigen for vector operations
- CMake 3.10 or higher

### Eigen

After cloning the repository, run
```
git submodule update --init
```
to pull in the Eigen dependency

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
- W/A/S/D/R/F: Move camera
- C: Toggle camera orbit mode
- T: Toggle wireframe mode
- P: Pause/start simulation
- V: Increase vorticity strength
- B: Decrease vorticity strength
- 1~9, 0, L: Switch between different color maps
- M: Toggle between Volume and Shell rendering
- O: Add obstacle
- I: Clear obstacles
- ↑ / ↓ / ← / →: Move obstacle in Y/X plane
- = / -: Move obstacle forward/backward along Z axis
- Space/Delete/Return: Clear all fluids
- Esc: Quit application

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

- Solid Objects: Add static obstacles like cubes, or *fruits!
- Multiple Fluid Types: Simulate different fluid types with varying properties
- Density Fading: Better visualization over time

## Tried Extensions But Failed
There are some tried extensions but somehow cannot work out:

 - We tried assigning a dyeColor to each cell and uploading it as a second 3D texture alongside densityTex to achieve multicolor fluid rendering. While the effect looked great, the performance cost was too high—uploading and sampling two large 3D textures per frame significantly reduced the frame rate, even at 48×48×48 resolution. We ultimately dropped this approach and may explore lighter alternatives.
 
- We tested a bit-pattern encoded Lookup Table approach inspired by Unity's fluid techniques for complex mesh boundaries. While effective for simple geometries, it produced artifacts at corners and required excessive memory as complexity increased. We reverted to a direct boundary method that proved more reliable.
## Credits
Based on Jos Stam's "Stable Fluids" paper

Implementation by **Unstable Fruits Team**: 
[Jue Han](https://github.com/amk09),
[TianXing Ji](https://github.com/TianxingJi), 
[Xiaoxi Yang](https://github.com/yangxiaoxi65),
[Akash Singirikonda](https://github.com/AkashSingirikonda)

Visualization using OpenGL and custom GLSL shaders





