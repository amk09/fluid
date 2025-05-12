# 3D Stable Fluid and More!

An implementation of an interactive 3D fluid simulator based on Jos Stam's [Stable Fluids](https://pages.cs.wisc.edu/~chaol/data/cs777/stam-stable_fluids.pdf) method. The solver is suitable for visual effects involving smoke, gas, or any fluid-like motion. 

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
- **Higher Resolution Computation:** Offline precomputing feature for higher resolution, makes it up to 256 cubes. We have two methods of offline rendering, one adds all grid values for each timestep into a single binary and loads that into ram before rendering. The other stores these values to separate files after each update, after which each file is rendered from in real time.

## Examples
![Fluid Simulation Demo](./example/example.gif)

*Real-time fluid simulation with default color.*

![UnstableFruits](https://github.com/user-attachments/assets/fa42f7d6-cd5a-4cd1-b5ff-77babd87683f)

*Real-time fluid simulation with unstable fruits.*

![Fluid Simulation Demo](./example/olexample.gif)

*Offline precomputed fluid simulation with default color.*

https://github.com/user-attachments/assets/f8089f12-4efb-434c-b1ec-5d79d3cddaca

*Offline precomputed fluid simulation with changing color.*

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
- CUDA Implementation: We implemented the solver features in CUDA. Additional Information included below.

## Tried Extensions But Failed
There are some tried extensions but somehow cannot work out:

 - We tried assigning a dyeColor to each cell and uploading it as a second 3D texture alongside densityTex to achieve multicolor fluid rendering. While the effect looked great, the performance cost was too high—uploading and sampling two large 3D textures per frame significantly reduced the frame rate, even at 48×48×48 resolution. We ultimately dropped this approach and may explore lighter alternatives.
 
- We tested a bit-pattern encoded Lookup Table approach inspired by Unity's fluid techniques for complex mesh boundaries. While effective for simple geometries, it produced artifacts at corners and required excessive memory as complexity increased. We reverted to a direct boundary method that proved more reliable.

## Cuda Implementation
Here are the results that were generated using the binaries along with the frame-by-frame offline renderer in this repository.

https://github.com/user-attachments/assets/821f1772-0c38-4cd1-897e-1a8353022fb8

The requirements for this project are an Nvidia GPU and the CUDA library. This can be setup by following the instructions [here]([url](https://docs.nvidia.com/cuda/cuda-installation-guide-microsoft-windows/)). If you have access to the Brown Hydra cluster, after [ssh'ing]([url](https://cs.brown.edu/courses/cs017/content/docs/using-ssh.pdf)) into the CS department, you can run ``` .\interact.sh ```. This workflow is borrowed from CSCI 1390 and we would like to thank the course staff immensely for this.

To build, from build directory:
```
cmake .. && make -j
```

To run:
```
./fluid_sim
```
This creates binaries for density (density_XXXX.bin), velocity magnitude (velocity_XXXX.bin), and color (which is just the density values color_XXXX.bin), which are available in the ```visualizer/renderData``` folder.

To visualize using matplotlib:
```
python3 visualizer.py
```
A sequence of cross sections in X-Y, Y-Z, and Z-X planes will be generated in the ```visualizer/plots``` folder.

To visualize from binaries using our renderer, place the density, velocity magnitude, and color files in the build/offline_data folder, and set the ```offlineFileLoadedFF``` value at the top of ```fluidcube.cpp``` to true. This will generate a visualization similar to the one above.


The only methodological difference between this and amk09/fluid is that instead of using a Gauss-Seidel solver, we use a Jacobi Kernel. This is because the Gauss-Seidel solver typically writes from data it reads to, which a Jacobi kernel does not. The Jacobi kernel however converges slower and is run for 10 solver step rather than the 4 that's used in our CPU implementation. 

We see an up to **6.7x** speedup in the update step on CUDA. We compared the update/step methods of the cuda version to the OPEN_MP parallelized CPU version. The CUDA version runs on a GTX 1080 and the CPU version runs on a M2 Max MacBook Pro with 32GB of Unified Memory. The current cuda implementation makes use of global memory coalescing, but doesn't make use of tiling or shared memory coalescing. We hope to extend this with these optimizations in the future along with the interactive features.

## Performance

| Grid Size | CPU Time (ms) | CUDA Time (ms) | Speed-up (CPU / CUDA) |
| --------- | ------------- | -------------- | --------------------- |
| 16        | 0.9           | 1.9            | 0.47×                 |
| 32        | 3.0           | 2.0            | **1.50×**             |
| 64        | 23.1          | 4.8            | **4.81×**             |
| 128       | 174           | 30             | **5.80×**             |
| 256       | 1508          | 225            | **6.70×**             |


## Credits
Based on Jos Stam's "Stable Fluids" paper

Implementation by **Unstable Fruits Team**: 
[Xiaoxi Yang](https://github.com/yangxiaoxi65),
[Jue Han](https://github.com/amk09),
[TianXing Ji](https://github.com/TianxingJi), 
[Akash Singirikonda](https://github.com/AkashSingirikonda)

Presentation slides for CSCI2240 can be viewed [here](/example/Final_Unstable_Fruits.pptx)

Visualization using OpenGL and custom GLSL shaders







