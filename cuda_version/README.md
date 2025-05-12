# fluidCuda

This is an implementation of Jo Stam's stable fluids, based off of amk09/fluid using cuda. (This only implements the simulation mechanism and not any of the interactive features)

Here are the results that were generated using the binaries along with the frame-by-frame offline renderer in amk09/fluid.

https://github.com/user-attachments/assets/821f1772-0c38-4cd1-897e-1a8353022fb8

Here are the results from the cross sections visualized from MatPlotLib code in visualizer.py.

The requirements for this project are an Nvidia GPU and the CUDA library. This can be setup by following the instructions [here]([url](https://docs.nvidia.com/cuda/cuda-installation-guide-microsoft-windows/)). If you have access to the Brown Hydra cluster, after [ssh'ing]([url](https://cs.brown.edu/courses/cs017/content/docs/using-ssh.pdf)) into the CS department, you can run ``` .\interact.sh ```. This workflow is borrowed from CSCI 1390 and I would like to thank the course staff immensely for this.

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

To viualize from binaries using our renderer, place the density, velocity magnitude, and color files in the build/offline_data folder, and set the ```offlineFileLoadedFF``` value at the top of ```fluidcube.cpp``` to true. This will generate a visualization similar to the one above.


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
