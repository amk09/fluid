#include "fluid_simulation.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <vector>
#include <cmath> // For std::floor
#include <algorithm> // For std::swap
#include <numeric>   // For std::accumulate (potential use)

// CUDA error checking macro
#define CUDA_CHECK(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", \
                    __FILE__, __LINE__, cudaGetErrorString(error)); \
            exit(1); \
        } \
    } while(0)

// Forward declaration for CPU debugging function (using float*)
void CPUlinSolveJacobi(int b, float* x, const float* rhs, float a, float c, int iter, int N);

// Constructor
FluidSimulation::FluidSimulation(int width, int height, int depth) {
    params.width = width;
    params.height = height;
    params.depth = depth;
    params.dt = 0.1f;
    params.visc = 0.0000001f;
    params.diff = 0.00001f;
    
    allocateMemory();
}

// Destructor
FluidSimulation::~FluidSimulation() {
    freeMemory();
}

// Memory allocation
void FluidSimulation::allocateMemory() {
    size_t size = params.width * params.height * params.depth * sizeof(float);
    CUDA_CHECK(cudaMalloc(&d_density0, size));
    CUDA_CHECK(cudaMalloc(&d_density, size));

    CUDA_CHECK(cudaMalloc(&d_Vx, size));
    CUDA_CHECK(cudaMalloc(&d_Vy, size));
    CUDA_CHECK(cudaMalloc(&d_Vz, size));
    
    CUDA_CHECK(cudaMalloc(&d_Vx0, size));
    CUDA_CHECK(cudaMalloc(&d_Vy0, size));
    CUDA_CHECK(cudaMalloc(&d_Vz0, size));

    CUDA_CHECK(cudaMalloc(&d_tmp, size));
    
    resetFields();
}

// Memory cleanup
void FluidSimulation::freeMemory() {
    CUDA_CHECK(cudaFree(d_density0));
    CUDA_CHECK(cudaFree(d_density));

    CUDA_CHECK(cudaFree(d_Vx));
    CUDA_CHECK(cudaFree(d_Vy));
    CUDA_CHECK(cudaFree(d_Vz));
    
    CUDA_CHECK(cudaFree(d_Vx0));
    CUDA_CHECK(cudaFree(d_Vy0));
    CUDA_CHECK(cudaFree(d_Vz0));

    CUDA_CHECK(cudaFree(d_tmp));
}

// Reset all fields to zero
void FluidSimulation::resetFields() {
    size_t size = params.width * params.height * params.depth * sizeof(float);

    CUDA_CHECK(cudaMemset(d_density0, 0, size));
    CUDA_CHECK(cudaMemset(d_density, 0, size));

    CUDA_CHECK(cudaMemset(d_Vx, 0, size));
    CUDA_CHECK(cudaMemset(d_Vy, 0, size));
    CUDA_CHECK(cudaMemset(d_Vz, 0, size));
    
    CUDA_CHECK(cudaMemset(d_Vx0, 0, size));
    CUDA_CHECK(cudaMemset(d_Vy0, 0, size));
    CUDA_CHECK(cudaMemset(d_Vz0, 0, size));

    CUDA_CHECK(cudaMemset(d_tmp, 0, size));

}

// Initialize simulation
void FluidSimulation::initialize() {
    resetFields();
}

// Main simulation step
void FluidSimulation::step() {
    static int current_simulation_step = 0; // Keep track of the simulation step number

    int solver_iterations = 20; // Default iterations for diffuse/project
 
    // Sync threads after each step // This comment seems misplaced, cudaDeviceSynchronize is used specifically.
    
    // Velocity Step
    // 1. Diffuse Velocities (input: d_Vx, d_Vy, d_Vz from prev step; output: d_Vx0, d_Vy0, d_Vz0)
    diffuse(1, d_Vx0, d_Vx, params.visc, params.dt, solver_iterations, params.width);
    diffuse(2, d_Vy0, d_Vy, params.visc, params.dt, solver_iterations, params.width);
    diffuse(3, d_Vz0, d_Vz, params.visc, params.dt, solver_iterations, params.width);

    // Project diffused velocities (d_Vx0, d_Vy0, d_Vz0)
    // Using d_Vx as pressure buffer, d_Vy as divergence buffer.
    // Result of projection is in d_Vx0, d_Vy0, d_Vz0.
    project(d_Vx0, d_Vy0, d_Vz0, d_Vx, d_Vy, solver_iterations, params.width);
 
    // 2. Advect Velocities 
    // Input for advection: d_Vx0, d_Vy0, d_Vz0 (diffused and projected)
    // Advecting field: d_Vx0, d_Vy0, d_Vz0
    // Output: d_Vx, d_Vy, d_Vz
    advect(1, d_Vx, d_Vx0, d_Vx0, d_Vy0, d_Vz0, params.dt, params.width);
    advect(2, d_Vy, d_Vy0, d_Vx0, d_Vy0, d_Vz0, params.dt, params.width);
    advect(3, d_Vz, d_Vz0, d_Vx0, d_Vy0, d_Vz0, params.dt, params.width);
 
    // 3. Project advected velocities (d_Vx, d_Vy, d_Vz)
    // Using d_Vx0 as pressure buffer, d_Vy0 as divergence buffer.
    // Result of projection is in d_Vx, d_Vy, d_Vz.
    project(d_Vx, d_Vy, d_Vz, d_Vx0, d_Vy0, solver_iterations, params.width);

    // Density Step
    // 1. Diffuse Density (input: d_density from prev step; output: d_density0)
    //printf("[Step %d] Density sum before diffuse: %f (reading from d_density)\n", current_simulation_step, sumFieldHost(d_density));
    diffuse(0, d_density0, d_density, params.diff, params.dt, solver_iterations, params.width);
    //printf("[Step %d] Density sum after diffuse:  %f (result in d_density0)\n", current_simulation_step, sumFieldHost(d_density0));
    
    // 2. Advect Density 
    // Input for advection: d_density0 (diffused density)
    // Advecting field: d_Vx, d_Vy, d_Vz (divergence-free)
    // Output: d_density
    // printf("[Step %d] Density sum before advect: %f (reading from d_density0)\n", current_simulation_step, sumFieldHost(d_density0)); // Optional, should be same as after diffuse
    advect(0, d_density, d_density0, d_Vx, d_Vy, d_Vz, params.dt, params.width);
    CUDA_CHECK(cudaDeviceSynchronize());
    //printf("[Step %d] Density sum after advect:   %f (result in d_density)\n", current_simulation_step, sumFieldHost(d_density));

    current_simulation_step++;
}

// Device version of IX function for CUDA kernels
__device__ int IX_device(int x, int y, int z, int width, int height) {
    return x + width * (y + height * z);  // x-major order
}

// Host version of IX function
int FluidSimulation::IX(int x, int y, int z) {
    return x + params.width * (y + params.height * z);  // x-major order
}

// CUDA kernel for adding velocity
__global__ void addVelocityKernel(float* d_Vx, float* d_Vy, float* d_Vz, 
                                int x, int y, int z, float vx, float vy, float vz,
                                int width, int height, int depth) {
    int idx = (z * height + y) * width + x;
    if (x >= 0 && x < width && y >= 0 && y < height && z >= 0 && z < depth) {
        d_Vx[idx] += vx;
        d_Vy[idx] += vy;
        d_Vz[idx] += vz;
    }
}

// CUDA kernel for adding density
__global__ void addDensityKernel(float* d_density, int x, int y, int z, 
                               float amount, int width, int height, int depth) {
    int idx = (z * height + y) * width + x;
    if (x >= 0 && x < width && y >= 0 && y < height && z >= 0 && z < depth) {
        d_density[idx] += amount;
    }
}

__device__ float computeGaussianWeightForOffset(int dx, int dy, int dz, float sigma) {
    if (sigma < 1e-6f) { // Treat very small sigma as effectively zero
        return (dx == 0 && dy == 0 && dz == 0) ? 1.0f : 0.0f;
    }
    float variance = sigma * sigma;
    // dx, dy, dz are integer offsets. Convert to float for calculation.
    float fdx = static_cast<float>(dx);
    float fdy = static_cast<float>(dy);
    float fdz = static_cast<float>(dz);
    float r_squared = fdx*fdx + fdy*fdy + fdz*fdz;
    // Ensure variance is not zero to prevent division by zero if sigma was extremely small but not < 1e-6f
    if (variance < 1e-9f) { // A very small positive variance
         return (r_squared < 1e-9f) ? 1.0f : 0.0f; // Effectively 1 if at center, 0 otherwise
    }
    return expf(-r_squared / (2.0f * variance));
}

__global__ void addDensityGaussianKernel(float* d_density,
                                       float centerX, float centerY, float centerZ,
                                       float amount, float sigma, int radius,
                                       int width, int height, int depth) {
    // Each thread handles one combination of (dx_offset, dy_offset, dz_offset)
    // where offsets range from -radius to +radius.
    // Total span for iteration is 2*radius + 1.
    int span = 2 * radius + 1;

    int flat_idx_x = blockIdx.x * blockDim.x + threadIdx.x;
    int flat_idx_y = blockIdx.y * blockDim.y + threadIdx.y;
    int flat_idx_z = blockIdx.z * blockDim.z + threadIdx.z;

    if (flat_idx_x >= span || flat_idx_y >= span || flat_idx_z >= span) {
        return; // This thread is outside the required (0 to span-1) range
    }

    int dx_offset = flat_idx_x - radius; // Converts flat_idx (0 to span-1) to offset (-radius to +radius)
    int dy_offset = flat_idx_y - radius;
    int dz_offset = flat_idx_z - radius;

    // Calculate target cell coordinates based on float center and integer offset, truncating like static_cast<int>
    int nx = static_cast<int>(centerX + static_cast<float>(dx_offset));
    int ny = static_cast<int>(centerY + static_cast<float>(dy_offset));
    int nz = static_cast<int>(centerZ + static_cast<float>(dz_offset));

    // Check if the target cell is within the simulation grid
    if (nx < 0 || nx >= width || ny < 0 || ny >= height || nz < 0 || nz >= depth) {
        return;
    }

    // Calculate Gaussian weight using the integer offsets
    float weight = computeGaussianWeightForOffset(dx_offset, dy_offset, dz_offset, sigma);

    if (weight > 1e-6f) { // Add only if there's a non-negligible contribution
        float val_to_add = amount * weight * 0.7f; // As per user's example
        atomicAdd(&d_density[IX_device(nx, ny, nz, width, height)], val_to_add);
    }
}

// Get density field as host vector
std::vector<float> FluidSimulation::getDensityFieldHost() {
    size_t size = params.width * params.height * params.depth;
    std::vector<float> host_data(size);
    CUDA_CHECK(cudaMemcpy(host_data.data(), d_density, size * sizeof(float), cudaMemcpyDeviceToHost));
    return host_data;
}

// Get velocity field as host vector
std::vector<float> FluidSimulation::getVelocityFieldHost(int component) {
    size_t size = params.width * params.height * params.depth;
    std::vector<float> host_data(size);
    
    float* d_field;
    switch(component) {
        case 0: d_field = d_Vx; break;
        case 1: d_field = d_Vy; break;
        case 2: d_field = d_Vz; break;
        default: return std::vector<float>(); // Return empty vector for invalid component
    }
    
    CUDA_CHECK(cudaMemcpy(host_data.data(), d_field, size * sizeof(float), cudaMemcpyDeviceToHost));
    return host_data;
}

// Add velocity at a point (using CUDA kernel)
void FluidSimulation::addVelocity(int x, int y, int z, float vx, float vy, float vz) {
    if (x < 0 || x >= params.width || y < 0 || y >= params.height || z < 0 || z >= params.depth) return;
    
    dim3 block(1);
    dim3 grid(1);
    addVelocityKernel<<<grid, block>>>(d_Vx, d_Vy, d_Vz, x, y, z, vx, vy, vz,
                                      params.width, params.height, params.depth);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

// Add density at a point (using CUDA kernel)
void FluidSimulation::addDensity(int x, int y, int z, float amount) {
    if (x < 0 || x >= params.width || y < 0 || y >= params.height || z < 0 || z >= params.depth) return;
    
    dim3 block(1);
    dim3 grid(1);
    addDensityKernel<<<grid, block>>>(d_density, x, y, z, amount,
                                     params.width, params.height, params.depth);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // Debug print
    printf("Added density %f at (%d,%d,%d)\n", amount, x, y, z);
}

// Add density with a Gaussian distribution
void FluidSimulation::addDensityGaussian(float centerX, float centerY, float centerZ, float amount, float sigma) {
    if (sigma < 0.0f) { // Sigma should be non-negative
        sigma = 0.0f;
    }

    // Determine the radius of influence for the Gaussian, as per user's example snippet
    // The loop iterates from -radius to +radius for offsets.
    int radius = static_cast<int>(sigma * 2.0f);
    if (radius < 0) radius = 0; 

    // The kernel iterates for offsets from -radius to +radius.
    // The total number of iterations (span) in each dimension is 2*radius + 1.
    int span = 2 * radius + 1;
    // If sigma is 0, radius is 0, span is 1. This is correct for a single point.

    dim3 threadsPerBlock(8, 8, 4); // Example block size
    dim3 numBlocks(
        (span + threadsPerBlock.x - 1) / threadsPerBlock.x,
        (span + threadsPerBlock.y - 1) / threadsPerBlock.y,
        (span + threadsPerBlock.z - 1) / threadsPerBlock.z
    );
    
    // Ensure numBlocks components are at least 1 if span is 1.
    if (numBlocks.x == 0) numBlocks.x = 1;
    if (numBlocks.y == 0) numBlocks.y = 1;
    if (numBlocks.z == 0) numBlocks.z = 1;


    addDensityGaussianKernel<<<numBlocks, threadsPerBlock>>>(
        d_density, // Target the main density buffer
        centerX, centerY, centerZ, amount, sigma, radius,
        params.width, params.height, params.depth
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize()); // Ensure kernel completes and effect is visible
}

// Get current density field
float* FluidSimulation::getDensityField() {
    return d_density;
}

// Save grid data to files
void FluidSimulation::saveGridData(int step) {
    // Create plots directory if it doesn't exist
    system("mkdir -p ../visualizer/renderData");
    
    // Get dimensions
    int N = params.width;
    size_t size = N * N * N * sizeof(float);
    
    // Create temporary arrays for CPU
    std::vector<float> density(N * N * N);
    std::vector<float> velocityX(N * N * N);
    std::vector<float> velocityY(N * N * N);
    std::vector<float> velocityZ(N * N * N);
    std::vector<float> velocityMag(N * N * N);
    
    // Copy data from GPU to CPU
    CUDA_CHECK(cudaMemcpy(density.data(), d_density, size, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(velocityX.data(), d_Vx, size, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(velocityY.data(), d_Vy, size, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(velocityZ.data(), d_Vz, size, cudaMemcpyDeviceToHost));
    // Create filename with step number
    char filename[256];
    
    // Save density
    snprintf(filename, sizeof(filename), "../visualizer/renderData/density_%04d.bin", step);
    FILE* fp = fopen(filename, "wb");
    if (fp) {
        // Write dimensions first
        // fwrite(&N, sizeof(int), 1, fp);
        // fwrite(&N, sizeof(int), 1, fp);
        // fwrite(&N, sizeof(int), 1, fp);
        // Write density data
        fwrite(density.data(), sizeof(float), N * N * N, fp);
        fclose(fp);
    }

    // Save density
    snprintf(filename, sizeof(filename), "../visualizer/renderData/color_%04d.bin", step);
    fp = fopen(filename, "wb");
    if (fp) {
        // Write dimensions first
        // fwrite(&N, sizeof(int), 1, fp);
        // fwrite(&N, sizeof(int), 1, fp);
        // fwrite(&N, sizeof(int), 1, fp);
        // Write density data
        fwrite(density.data(), sizeof(float), N * N * N, fp);
        fclose(fp);
    }
    
    // // Save velocity X
    // snprintf(filename, sizeof(filename), "../visualizer/renderData/velocityX_%04d.bin", step);
    // fp = fopen(filename, "wb");
    // if (fp) {
    //     // Write dimensions first
    //     fwrite(&N, sizeof(int), 1, fp);
    //     fwrite(&N, sizeof(int), 1, fp);
    //     fwrite(&N, sizeof(int), 1, fp);
    //     // Write velocity data
    //     fwrite(velocityX.data(), sizeof(float), N * N * N, fp);
    //     fclose(fp);
    // }
    
    // // Save velocity Y
    // snprintf(filename, sizeof(filename), "../visualizer/renderData/velocityY_%04d.bin", step);
    // fp = fopen(filename, "wb");
    // if (fp) {
    //     // Write dimensions first
    //     fwrite(&N, sizeof(int), 1, fp);
    //     fwrite(&N, sizeof(int), 1, fp);
    //     fwrite(&N, sizeof(int), 1, fp);
    //     // Write velocity data
    //     fwrite(velocityY.data(), sizeof(float), N * N * N, fp);
    //     fclose(fp);
    // }
    
    // // Save velocity Z
    // snprintf(filename, sizeof(filename), "../visualizer/renderData/velocityZ_%04d.bin", step);
    // fp = fopen(filename, "wb");
    // if (fp) {
    //     // Write dimensions first
    //     fwrite(&N, sizeof(int), 1, fp);
    //     fwrite(&N, sizeof(int), 1, fp);
    //     fwrite(&N, sizeof(int), 1, fp);
    //     // Write velocity data
    //     fwrite(velocityZ.data(), sizeof(float), N * N * N, fp);
    //     fclose(fp);
    // }
    // Save velocity magnitude
    snprintf(filename, sizeof(filename), "../visualizer/renderData/velocity_%04d.bin", step);
    fp = fopen(filename, "wb");
    if (fp) {
        // Write dimensions first
        // fwrite(&N, sizeof(int), 1, fp); 
        // fwrite(&N, sizeof(int), 1, fp);
        // fwrite(&N, sizeof(int), 1, fp);

        // Calculate velocity magnitude
        for (int i = 0; i < N * N * N; i++) {
            velocityMag[i] = sqrt(velocityX[i] * velocityX[i] + velocityY[i] * velocityY[i] + velocityZ[i] * velocityZ[i]);
        }

        // Write velocity data
        fwrite(velocityMag.data(), sizeof(float), N * N * N, fp);
        fclose(fp);
    }   
    

    printf("Saved grid data to file \"../visualizer/renderData/dataZ_%04d.bin\"\n", step);
}

// Simulation Methods:

// Base Methods:

// CUDA kernel for linear solver
// __global__ void linSolveKernel(float* x, float* x0, float a, float c, int N) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
//     int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
//     int k = blockIdx.z * blockDim.z + threadIdx.z + 1;
    
//     if (i < N-1 && j < N-1 && k < N-1) {
//         float cRecip = 1.0f / c;
        
//         // Get current cell value
//         float current = x0[IX_device(i, j, k, N, N)];
        
//         // Get neighboring cell values
//         float x_plus = x[IX_device(i+1, j, k, N, N)];
//         float x_minus = x[IX_device(i-1, j, k, N, N)];
//         float y_plus = x[IX_device(i, j+1, k, N, N)];
//         float y_minus = x[IX_device(i, j-1, k, N, N)];
//         float z_plus = x[IX_device(i, j, k+1, N, N)];
//         float z_minus = x[IX_device(i, j, k-1, N, N)];
        
//         // Calculate new value
//         float new_val = (current + a * (x_plus + x_minus + y_plus + y_minus + z_plus + z_minus)) * cRecip;
        
//         // Store result
//         x[IX_device(i, j, k, N, N)] = new_val;
//     }
// }

// // Linear solver for diffusion and pressure
// void FluidSimulation::linSolve(int b, float* x, float* x0, float a, float c, int iter, int N) {
//     // Calculate grid and block dimensions
//     dim3 blockDim(8, 8, 8);  // 8x8x8 threads per block
//     dim3 gridDim(
//         (N + blockDim.x - 1) / blockDim.x,
//         (N + blockDim.y - 1) / blockDim.y,
//         (N + blockDim.z - 1) / blockDim.z
//     );
    
//     // Perform iterations
//     for (int k = 0; k < iter; k++) {
//         // Launch kernel
//         linSolveKernel<<<gridDim, blockDim>>>(x, x0, a, c, N);
//         CUDA_CHECK(cudaGetLastError());
//         CUDA_CHECK(cudaDeviceSynchronize());
        
//         // Set boundary conditions
//         set_bnd(b, x, N);
//     }
// }

#define RADIUS 1                     // Stencil radius (1 for 7-point Laplace)
#define BLOCK_X 8                    // Threads per block in X
#define BLOCK_Y 8                    // Threads per block in Y
#define BLOCK_Z 8                    // Threads per block in Z
//
// Shared-memory footprint = (BLOCK_Z+2R) * (BLOCK_Y+2R) * (BLOCK_X+2R) floats
// → (8+2)³ = 1000 floats ≈ 4 kB  → plenty of occupancy headroom on most GPUs
// -----------------------------------------------------------------------------

__global__ void jacobiSweep3D_shared(float*       __restrict__ dst,
                                     const float* __restrict__ src,
                                     const float* __restrict__ rhs,
                                     float                     a,
                                     float                     cRecip,
                                     int                       N)
{
    // ----------------------------- Global coordinates ------------------------
    int gi = blockIdx.x * BLOCK_X + threadIdx.x;           // 0 … N-1
    int gj = blockIdx.y * BLOCK_Y + threadIdx.y;
    int gk = blockIdx.z * BLOCK_Z + threadIdx.z;

    // Skip the outer one-cell frame; guards later assume 1 ≤ gi<N-1, etc.
    if (gi >= N || gj >= N || gk >= N) return;

    // ----------------------------- Shared memory tile ------------------------
    extern __shared__ float sh[];   // 3-D slab, flattened
    // Dimensions inside shared memory
    const int shX = BLOCK_X + 2*RADIUS;
    const int shY = BLOCK_Y + 2*RADIUS;
    const int shZ = BLOCK_Z + 2*RADIUS;

    // Lambda for flattening (z,y,x) → 1-D
    auto sidx = [=] __device__ (int z,int y,int x) {
        return (z*shY + y)*shX + x;
    };

    // Local coords *inside* the shared tile (including halo offset)
    int li = threadIdx.x + RADIUS;
    int lj = threadIdx.y + RADIUS;
    int lk = threadIdx.z + RADIUS;

    // ----------------------------- Load centre cell --------------------------
    sh[sidx(lk, lj, li)] = src[IX_device(gi, gj, gk, N, N)];

    // ----------------------------- Load halo cells ---------------------------
    // Every thread cooperatively pulls in at most 6 neighbours (faces only).
    // We keep it branch-free; if a neighbour lies outside the true domain we
    // re-use the centre value (Dirichlet zero-gradient boundary).

    // Offsets for 6 face neighbours
#pragma unroll
    for (int face = 0; face < 6; ++face) {
        int di = (face == 0) - (face == 1);   // +1 x, −1 x
        int dj = (face == 2) - (face == 3);   // +1 y, −1 y
        int dk = (face == 4) - (face == 5);   // +1 z, −1 z

        int gni = gi + di;
        int gnj = gj + dj;
        int gnk = gk + dk;

        // Position in shared memory for that neighbour
        int lni = li + di;
        int lnj = lj + dj;
        int lnk = lk + dk;

        // Bounds check once (all threads do identical comparisons → no warp diverge)
        bool inside = (gni >= 0 && gni < N &&
                       gnj >= 0 && gnj < N &&
                       gnk >= 0 && gnk < N);

        sh[sidx(lnk, lnj, lni)] =
            inside ? src[IX_device(gni, gnj, gnk, N, N)]
                   : sh[sidx(lk, lj, li)];  // replicate centre for out-of-domain
    }
    __syncthreads();

    // ----------------------------- Compute Jacobi update ---------------------
    // Skip the global outer frame to avoid reading unallocated memory
    if (gi > 0 && gj > 0 && gk > 0 && gi < N-1 && gj < N-1 && gk < N-1) {

        float nbrSum =
              sh[sidx(lk,     lj,     li+1)] + sh[sidx(lk,     lj,     li-1)]
            + sh[sidx(lk,     lj+1,   li  )] + sh[sidx(lk,     lj-1,   li  )]
            + sh[sidx(lk+1,   lj,     li  )] + sh[sidx(lk-1,   lj,     li  )];

        int gIdx = IX_device(gi, gj, gk, N, N);
        dst[gIdx] = (rhs[gIdx] + a * nbrSum) * cRecip;
    }
}

__global__ void jacobiSweep3D(float*       dst,
                              const float* src,
                              const float* rhs,
                              float        a,
                              float        cRecip,
                              int          N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.z * blockDim.z + threadIdx.z + 1;
    if (i >= N-1 || j >= N-1 || k >= N-1) return;

    int idx = IX_device(i,j,k,N,N);

    // Sum neighbors from the source array (previous iteration)
    float nbrSum =  src[IX_device(i+1,j  ,k  ,N,N)] + src[IX_device(i-1,j  ,k  ,N,N)]
                  + src[IX_device(i  ,j+1,k  ,N,N)] + src[IX_device(i  ,j-1,k  ,N,N)]
                  + src[IX_device(i  ,j  ,k+1,N,N)] + src[IX_device(i  ,j  ,k-1,N,N)];

    // Calculate new value and write to destination array
    dst[idx] = (rhs[idx] + a * nbrSum) * cRecip;
}

// // host wrapper for Jacobi solver
void FluidSimulation::linSolve(int    b,        // field type 0,1,2,3 for boundaries
                                       float* x,        // IN: initial guess, OUT: solution
                                       float* rhs,      // immutable rhs (like divergence)
                                       float  a, float c,
                                       int    iter, int N)
{
    // Use the pre-allocated temporary buffer
    float* tmp = d_tmp;          // same size as x

    float* src = x;              // start reading from x
    float* dst = tmp;            // write into tmp first

    dim3 blk(8,8,8);
    // Grid dimensions cover the interior cells (1 to N-2)
    dim3 grd( (N-2 + blk.x-1) / blk.x,
              (N-2 + blk.y-1) / blk.y,
              (N-2 + blk.z-1) / blk.z );

    float cRecip = 1.0f / c;

    for (int n = 0; n < iter; ++n) {
        // Perform one Jacobi sweep
        // dim3 block(BLOCK_X, BLOCK_Y, BLOCK_Z);

        // // Domain interior is (N-2)³ points (we skip the boundaries),
        // // but we can simply launch enough blocks to cover the full N³
        // // cube and rely on the i>0 && i<N-1 guard inside the kernel.
        // dim3 grid( (N + BLOCK_X - 1) / BLOCK_X,
        //         (N + BLOCK_Y - 1) / BLOCK_Y,
        //         (N + BLOCK_Z - 1) / BLOCK_Z );

        // size_t shMemBytes = (BLOCK_X + 2*RADIUS) *
        //                     (BLOCK_Y + 2*RADIUS) *
        //                     (BLOCK_Z + 2*RADIUS) * sizeof(float);

        // jacobiSweep3D_shared<<<grid, block, shMemBytes>>>(
        //         dst, src, rhs, a, cRecip, N);
        jacobiSweep3D<<<grd, blk>>>(dst, src, rhs, a, cRecip, N);
        CUDA_CHECK(cudaGetLastError());
        // No synchronize needed here, boundary kernel launch will sync.

        // Enforce boundary conditions on the *newly written* data in 'dst'
        set_bnd(b, dst, N);

        // Swap pointers for next iteration:
        // the destination of this step becomes the source for the next
        std::swap(src, dst);
    }

    // After the loop, 'src' points to the array holding the final result
    // (because of the last swap). If 'src' is not the original 'x' array,
    // we need to copy the result back into 'x'.
    if (src != x) {
        size_t total_elements = (size_t)N * N * N;
        size_t size_bytes = total_elements * sizeof(float);
        CUDA_CHECK(cudaMemcpy(x, src, size_bytes, cudaMemcpyDeviceToDevice));
    }
}

// CUDA kernel for setting face boundaries
__global__ void set_bnd_kernel(float* x, int b, int N) {
    // Handle X-Face boundaries (i=0 and i=N-1)
    if (blockIdx.z == 0) {
        int j = blockIdx.x * blockDim.x + threadIdx.x + 1;
        int k = blockIdx.y * blockDim.y + threadIdx.y + 1;
        
        if (j < N-1 && k < N-1) {
            if (b == 1) {
                x[IX_device(0, j, k, N, N)] = -x[IX_device(1, j, k, N, N)];
                x[IX_device(N-1, j, k, N, N)] = -x[IX_device(N-2, j, k, N, N)];
            } else {
                x[IX_device(0, j, k, N, N)] = x[IX_device(1, j, k, N, N)];
                x[IX_device(N-1, j, k, N, N)] = x[IX_device(N-2, j, k, N, N)];
            }
        }
    }
    
    // Handle Y-Face boundaries (j=0 and j=N-1)
    if (blockIdx.z == 1) {
        int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
        int k = blockIdx.y * blockDim.y + threadIdx.y + 1;
        
        if (i < N-1 && k < N-1) {
            if (b == 2) {
                x[IX_device(i, 0, k, N, N)] = -x[IX_device(i, 1, k, N, N)];
                x[IX_device(i, N-1, k, N, N)] = -x[IX_device(i, N-2, k, N, N)];
            } else {
                x[IX_device(i, 0, k, N, N)] = x[IX_device(i, 1, k, N, N)];
                x[IX_device(i, N-1, k, N, N)] = x[IX_device(i, N-2, k, N, N)];
            }
        }
    }
    
    // Handle Z-Face boundaries (k=0 and k=N-1)
    if (blockIdx.z == 2) {
        int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
        int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
        
        if (i < N-1 && j < N-1) {
            if (b == 3) {
                x[IX_device(i, j, 0, N, N)] = -x[IX_device(i, j, 1, N, N)];
                x[IX_device(i, j, N-1, N, N)] = -x[IX_device(i, j, N-2, N, N)];
            } else {
                x[IX_device(i, j, 0, N, N)] = x[IX_device(i, j, 1, N, N)];
                x[IX_device(i, j, N-1, N, N)] = x[IX_device(i, j, N-2, N, N)];
            }
        }
    }
}

// CUDA kernel for setting corner boundaries
__global__ void set_bnd_corner_kernel(float* x, int N) {
    // Only one thread needed for corners
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        // Corner (0,0,0)
        x[IX_device(0, 0, 0, N, N)] = 0.33f * (x[IX_device(1, 0, 0, N, N)] + x[IX_device(0, 1, 0, N, N)] + x[IX_device(0, 0, 1, N, N)]);
        
        // Corner (0,N-1,0)
        x[IX_device(0, N-1, 0, N, N)] = 0.33f * (x[IX_device(1, N-1, 0, N, N)] + x[IX_device(0, N-2, 0, N, N)] + x[IX_device(0, N-1, 1, N, N)]);
        
        // Corner (0,0,N-1)
        x[IX_device(0, 0, N-1, N, N)] = 0.33f * (x[IX_device(1, 0, N-1, N, N)] + x[IX_device(0, 1, N-1, N, N)] + x[IX_device(0, 0, N-2, N, N)]);
        
        // Corner (0,N-1,N-1)
        x[IX_device(0, N-1, N-1, N, N)] = 0.33f * (x[IX_device(1, N-1, N-1, N, N)] + x[IX_device(0, N-2, N-1, N, N)] + x[IX_device(0, N-1, N-2, N, N)]);
        
        // Corner (N-1,0,0)
        x[IX_device(N-1, 0, 0, N, N)] = 0.33f * (x[IX_device(N-2, 0, 0, N, N)] + x[IX_device(N-1, 1, 0, N, N)] + x[IX_device(N-1, 0, 1, N, N)]);
        
        // Corner (N-1,N-1,0)
        x[IX_device(N-1, N-1, 0, N, N)] = 0.33f * (x[IX_device(N-2, N-1, 0, N, N)] + x[IX_device(N-1, N-2, 0, N, N)] + x[IX_device(N-1, N-1, 1, N, N)]);
        
        // Corner (N-1,0,N-1)
        x[IX_device(N-1, 0, N-1, N, N)] = 0.33f * (x[IX_device(N-2, 0, N-1, N, N)] + x[IX_device(N-1, 1, N-1, N, N)] + x[IX_device(N-1, 0, N-2, N, N)]);
        
        // Corner (N-1,N-1,N-1)
        x[IX_device(N-1, N-1, N-1, N, N)] = 0.33f * (x[IX_device(N-2, N-1, N-1, N, N)] + x[IX_device(N-1, N-2, N-1, N, N)] + x[IX_device(N-1, N-1, N-2, N, N)]);
    }
}

// Set boundary conditions
void FluidSimulation::set_bnd(int b, float* x, int N) {
    // Calculate grid and block dimensions for face boundaries
    dim3 blockDim(16, 16);
    dim3 gridDim(
        (N + blockDim.x - 1) / blockDim.x,
        (N + blockDim.y - 1) / blockDim.y,
        3  // 3 types of faces: X, Y, Z
    );
    
    // Set face boundaries
    set_bnd_kernel<<<gridDim, blockDim>>>(x, b, N);
    CUDA_CHECK(cudaGetLastError());

    // Synchronize to ensure all operations are complete
    CUDA_CHECK(cudaDeviceSynchronize());

    // Set corner boundaries
    set_bnd_corner_kernel<<<1, 1>>>(x, N);
    CUDA_CHECK(cudaGetLastError());
    
    // Synchronize to ensure all operations are complete
    CUDA_CHECK(cudaDeviceSynchronize());
}

// Set corner boundaries
void FluidSimulation::set_bnd_corner(float* x, int N) {
    set_bnd_corner_kernel<<<1, 1>>>(x, N);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

// Diffusion step
void FluidSimulation::diffuse(int b, float* x, float* x0, float diff, float dt, int iter, int N) {
    float a = dt * diff * (N - 2) * (N - 2);
    linSolve(b, x, x0, a, 1 + 6 * a, iter, N);
}

// CUDA kernel for advection
__global__ void advectKernel(float* d, float* d0, float* velocX, float* velocY, float* velocZ, 
                            float dt, int N) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.z * blockDim.z + threadIdx.z + 1;
    
    if (i < N-1 && j < N-1 && k < N-1) {
        float dtx = dt * (N - 2);
        float dty = dt * (N - 2);
        float dtz = dt * (N - 2);
        
        // Get velocity at current position
        float tmp1 = dtx * velocX[IX_device(i, j, k, N, N)];
        float tmp2 = dty * velocY[IX_device(i, j, k, N, N)];
        float tmp3 = dtz * velocZ[IX_device(i, j, k, N, N)];
        
        // Calculate back-traced position
        float x = (float)i - tmp1;
        float y = (float)j - tmp2;
        float z = (float)k - tmp3;
        
        // Clamp positions to valid range for interpolation
        float Nfloat = (float)N;
        if(x < 0.5f) x = 0.5f; 
        if(x > Nfloat - 1.5f) x = Nfloat - 1.5f; 
        if(y < 0.5f) y = 0.5f; 
        if(y > Nfloat - 1.5f) y = Nfloat - 1.5f; 
        if(z < 0.5f) z = 0.5f;
        if(z > Nfloat - 1.5f) z = Nfloat - 1.5f;
        
        // Get integer positions for interpolation
        int i0 = floorf(x);
        int i1 = i0 + 1;
        int j0 = floorf(y);
        int j1 = j0 + 1;
        int k0 = floorf(z);
        int k1 = k0 + 1;
        
        // Calculate interpolation weights
        float s1 = x - i0;
        float s0 = 1.0f - s1;
        float t1 = y - j0;
        float t0 = 1.0f - t1;
        float u1 = z - k0;
        float u0 = 1.0f - u1;
        
        // Perform trilinear interpolation
        d[IX_device(i, j, k, N, N)] = 
            s0 * (t0 * (u0 * d0[IX_device(i0, j0, k0, N, N)] +
                        u1 * d0[IX_device(i0, j0, k1, N, N)]) +
                  t1 * (u0 * d0[IX_device(i0, j1, k0, N, N)] +
                        u1 * d0[IX_device(i0, j1, k1, N, N)])) +
            s1 * (t0 * (u0 * d0[IX_device(i1, j0, k0, N, N)] +
                        u1 * d0[IX_device(i1, j0, k1, N, N)]) +
                  t1 * (u0 * d0[IX_device(i1, j1, k0, N, N)] +
                        u1 * d0[IX_device(i1, j1, k1, N, N)]));
    }
}

// Advection step
void FluidSimulation::advect(int b, float* d, float* d0, float* velocX, float* velocY, float* velocZ, float dt, int N) {
    // Calculate grid and block dimensions
    dim3 blockDim(8, 8, 8);  // 8x8x8 threads per block
    dim3 gridDim(
        (N + blockDim.x - 1) / blockDim.x,
        (N + blockDim.y - 1) / blockDim.y,
        (N + blockDim.z - 1) / blockDim.z
    );
    
    // Launch kernel
    advectKernel<<<gridDim, blockDim>>>(d, d0, velocX, velocY, velocZ, dt, N);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // Set boundary conditions
    set_bnd(b, d, N);
}

// Helper function to sum field values on the host (for debugging)
float FluidSimulation::sumFieldHost(float* d_field) {
    size_t num_elements = params.width * params.height * params.depth;
    size_t size_bytes = num_elements * sizeof(float);
    std::vector<float> h_field(num_elements);

    CUDA_CHECK(cudaMemcpy(h_field.data(), d_field, size_bytes, cudaMemcpyDeviceToHost));

    float sum = 0.0f;
    for (size_t i = 0; i < num_elements; ++i) {
        sum += h_field[i];
    }
    return sum;
}

__global__ void computeDivergenceKernel(float* d_divergence, float* d_pressure,
                                      float* d_Vx, float* d_Vy, float* d_Vz, int N) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

    if (i < N - 1 && j < N - 1 && k < N - 1) {
        int idx = IX_device(i, j, k, N, N);
        
        float div_val = -0.5f * (
            (d_Vx[IX_device(i + 1, j, k, N, N)] - d_Vx[IX_device(i - 1, j, k, N, N)]) +
            (d_Vy[IX_device(i, j + 1, k, N, N)] - d_Vy[IX_device(i, j - 1, k, N, N)]) +
            (d_Vz[IX_device(i, j, k + 1, N, N)] - d_Vz[IX_device(i, j, k - 1, N, N)])
        ) / N; // As per user's provided formula

        d_divergence[idx] = div_val;
        d_pressure[idx] = 0.0f; // Initialize pressure to zero
    }
}

__global__ void subtractPressureGradientKernel(float* d_Vx, float* d_Vy, float* d_Vz,
                                             float* d_pressure, int N) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

    if (i < N - 1 && j < N - 1 && k < N - 1) {
        int idx = IX_device(i, j, k, N, N);

        d_Vx[idx] -= 0.5f * (d_pressure[IX_device(i + 1, j, k, N, N)] - d_pressure[IX_device(i - 1, j, k, N, N)]) * N;
        d_Vy[idx] -= 0.5f * (d_pressure[IX_device(i, j + 1, k, N, N)] - d_pressure[IX_device(i, j - 1, k, N, N)]) * N;
        d_Vz[idx] -= 0.5f * (d_pressure[IX_device(i, j, k + 1, N, N)] - d_pressure[IX_device(i, j, k - 1, N, N)]) * N;
    }
}

// Projection step: Enforces divergence-free condition on velocity field
void FluidSimulation::project(float *velocX, float *velocY, float *velocZ, float *p_buffer, float *div_buffer, int iter, int N) {
    // Calculate grid and block dimensions
    dim3 blockDim(8, 8, 8);
    dim3 gridDim(
        (N + blockDim.x - 1) / blockDim.x,
        (N + blockDim.y - 1) / blockDim.y,
        (N + blockDim.z - 1) / blockDim.z
    );

    // Step 1: Compute divergence and initialize pressure buffer to 0
    computeDivergenceKernel<<<gridDim, blockDim>>>(div_buffer, p_buffer, velocX, velocY, velocZ, N);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    set_bnd(0, div_buffer, N); // Boundary condition for divergence
    set_bnd(0, p_buffer, N);   // Boundary condition for pressure (initial guess)
    CUDA_CHECK(cudaDeviceSynchronize()); // Ensure BCs are set before copying
    // Note: linSolve modifies the 'x' buffer in place (p_buffer here)
    linSolve(0, p_buffer, div_buffer, 1.0f, 6.0f, iter, N);

    // --- Continue with original project logic ---

    // Step 3: Subtract the pressure gradient from the velocity field
    // Note: This uses the p_buffer which now contains the GPU result
    subtractPressureGradientKernel<<<gridDim, blockDim>>>(velocX, velocY, velocZ, p_buffer, N);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    set_bnd(1, velocX, N); // Boundary conditions for velocity components
    set_bnd(2, velocY, N);
    set_bnd(3, velocZ, N);
}

__global__ void addVelocityGaussianKernel(float* d_Vx, float* d_Vy, float* d_Vz,
                                       float centerX, float centerY, float centerZ,
                                       float amountVx, float amountVy, float amountVz, float sigma, int radius,
                                       int width, int height, int depth) {
    // Each thread handles one combination of (dx_offset, dy_offset, dz_offset)
    int span = 2 * radius + 1;

    int flat_idx_x = blockIdx.x * blockDim.x + threadIdx.x;
    int flat_idx_y = blockIdx.y * blockDim.y + threadIdx.y;
    int flat_idx_z = blockIdx.z * blockDim.z + threadIdx.z;

    if (flat_idx_x >= span || flat_idx_y >= span || flat_idx_z >= span) {
        return;
    }

    int dx_offset = flat_idx_x - radius;
    int dy_offset = flat_idx_y - radius;
    int dz_offset = flat_idx_z - radius;

    int nx = static_cast<int>(centerX + static_cast<float>(dx_offset));
    int ny = static_cast<int>(centerY + static_cast<float>(dy_offset));
    int nz = static_cast<int>(centerZ + static_cast<float>(dz_offset));

    if (nx < 0 || nx >= width || ny < 0 || ny >= height || nz < 0 || nz >= depth) {
        return;
    }

    float weight = computeGaussianWeightForOffset(dx_offset, dy_offset, dz_offset, sigma);

    if (weight > 1e-6f) { // Add only if there's a non-negligible contribution
        int target_idx = IX_device(nx, ny, nz, width, height);
        atomicAdd(&d_Vx[target_idx], amountVx * weight);
        atomicAdd(&d_Vy[target_idx], amountVy * weight);
        atomicAdd(&d_Vz[target_idx], amountVz * weight);
    }
}

// Add velocity with a Gaussian distribution
void FluidSimulation::addVelocityGaussian(float centerX, float centerY, float centerZ, float vx, float vy, float vz, float sigma) {
    if (sigma < 0.0f) {
        sigma = 0.0f;
    }

    int radius = static_cast<int>(sigma * 2.0f);
    if (radius < 0) radius = 0;

    int span = 2 * radius + 1;

    dim3 threadsPerBlock(8, 8, 4); // Example block size, can be tuned
    dim3 numBlocks(
        (span + threadsPerBlock.x - 1) / threadsPerBlock.x,
        (span + threadsPerBlock.y - 1) / threadsPerBlock.y,
        (span + threadsPerBlock.z - 1) / threadsPerBlock.z
    );

    if (numBlocks.x == 0) numBlocks.x = 1;
    if (numBlocks.y == 0) numBlocks.y = 1;
    if (numBlocks.z == 0) numBlocks.z = 1;

    addVelocityGaussianKernel<<<numBlocks, threadsPerBlock>>>(
        d_Vx, d_Vy, d_Vz,
        centerX, centerY, centerZ,
        vx, vy, vz, sigma, radius,
        params.width, params.height, params.depth
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

