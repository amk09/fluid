#pragma once

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <vector>
#include <string>

// Simulation parameters
struct SimulationParams {
    int width;
    int height;
    int depth;      // Add depth for 3D
    float dt;        // Time step
    float visc; // Viscosity coefficient
    float diff; // Diffusion coefficient
};

// Fluid simulation class
class FluidSimulation {
public:
    FluidSimulation(int width, int height, int depth);
    ~FluidSimulation();

    int IX(int x, int y, int z);
    // Initialize the simulation
    void initialize();

    // Main simulation step
    void step();

    // Add velocity at a point (3D)
    void addVelocity(int x, int y, int z, float vx, float vy, float vz);

    // Add density at a point (3D)
    void addDensity(int x, int y, int z, float amount);

    // Get current density field
    float* getDensityField();

    // Get density field as host vector
    std::vector<float> getDensityFieldHost();

    // Get velocity field as host vector (component: 0=x, 1=y, 2=z)
    std::vector<float> getVelocityFieldHost(int component);

    // Save all grid data (density and velocity) to binary files
    void saveGridData(int step);

    // Linear solver for diffusion and pressure
    void linSolve(int b, float* x, float* x0, float a, float c, int iter, int N);

    // Diffusion step
    void diffuse(int b, float* x, float* x0, float diff, float dt, int iter, int N);

    // Advection step
    void advect(int b, float* d, float* d0, float* velocX, float* velocY, float* velocZ, float dt, int N);

    // Projection step
    void project(float *velocX, float *velocY, float *velocZ, float *p, float *div_field, int iter, int N);
    // Set boundary conditions
    void set_bnd(int b, float* x, int N);
    void set_bnd_corner(float* x, int N);

    // Helper function to sum field values on the host (for debugging)
    float sumFieldHost(float* d_field);
    std::vector<float> getFieldHost(float* d_field);  

    // Interaction Methods:
    void addDensityGaussian(float centerX, float centerY, float centerZ, float amount, float sigma);
    void addVelocityGaussian(float centerX, float centerY, float centerZ, float vx, float vy, float vz, float sigma);
private:

    SimulationParams params;
    
    // Device memory pointers (3D)
    float *d_density0;
    float *d_density;
    
    float *d_Vx;
    float *d_Vy;
    float *d_Vz;

    float *d_Vx0;
    float *d_Vy0;
    float *d_Vz0;

    float *d_tmp;

    // Helper functions
    void allocateMemory();
    void freeMemory();
    void resetFields();
}; 