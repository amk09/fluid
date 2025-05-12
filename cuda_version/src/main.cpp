#include "fluid_simulation.h"
#include <iostream>
#include <fstream>
#include <chrono>
int main() {
    // Create a 64x64x64 fluid simulation
    int w = 64, h = 64, d = 64;
    FluidSimulation sim(w, h, d);
    
    // Initialize the simulation
    sim.initialize();
    
    // Open CSV file for writing timing data
    std::ofstream timing_file("save_times_" + std::to_string(w) + "x" + std::to_string(h) + "x" + std::to_string(d) + "_cuda.csv");
    timing_file << "frame,elapsed_ms\n";  // Write header
    
    // Add velocity

    // for (int i = 0; i < w/3; i++) {
    //     sim.addVelocity(i, h/2, d/2, 0.1, 0, 0);
    // }
    // // Add density
    // sim.addDensity(w*0.7, h/2, d/2, 0.5f);
    // sim.addDensity(w*0.7-1, h/2, d/2, 0.5f);

    // for (int i = 0; i < w/3; i++) {
    //     sim.addVelocity(w-1-i, h/2, d/2, -0.1, 0, 0);
    // }

    // sim.addDensity(w*0.3+1, h/2, d/2, 0.5f);
    // sim.addDensity(w*0.3, h/2, d/2, 0.5f);


    // for (int i = 0; i <= w*0.2; i++) {
    //     sim.addVelocity(i, d/2, i, -0.1, -0.1, 0);
    // }
    // sim.addDensity(w*0.2, d/2, w*0.2, 0.5f);
    // sim.addDensity(w*0.2, d/2, w*0.2+1, 0.5f);

    sim.addDensityGaussian(w*0.1, h/2, d/2, 1, 5);
    sim.addVelocityGaussian(w*0.1, h/2, d/2, 5, 0, 0, 5);

    sim.addDensityGaussian(w*0.9, h/2, d/2, 1, 5);
    sim.addVelocityGaussian(w*0.9, h/2, d/2, -5, 0, 0, 5);

    // sim.addDensityGaussian(w*0.8, h/2, d/2, 1, 5);
    // sim.addVelocityGaussian(w*0.8, h/2, d/2, -5, 0, 0, 5);

    
    for (int i = 0; i < 100; i++) {
        // Time each step, save to csv
        auto start = std::chrono::high_resolution_clock::now();
        sim.step();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        // Convert to milliseconds
        timing_file << i << "," << duration.count() << "\n";
        std::cout << "Step " << i << " took " << duration.count() << " milliseconds" << std::endl;
        sim.saveGridData(i);
    }

    timing_file.close();

    // // Get data back to host
    // std::vector<float> density_field = sim.getDensityFieldHost();
    // std::vector<float> velocity_x = sim.getVelocityFieldHost(0);  // x component
    // std::vector<float> velocity_y = sim.getVelocityFieldHost(1);  // y component
    // std::vector<float> velocity_z = sim.getVelocityFieldHost(2);  // z component

    // for (int i = 0; i < w * h * d; i++) {
    //     std::cout << "Density field: " << density_field[i] << std::endl;
    //     std::cout << "Velocity x: " << velocity_x[i] << std::endl;
    //     std::cout << "Velocity y: " << velocity_y[i] << std::endl;
    //     std::cout << "Velocity z: " << velocity_z[i] << std::endl;
    // }
    // // Add some initial velocity and density at the center
    // int cx = w / 2, cy = h / 2, cz = d / 2;
    // sim.addVelocity(cx, cy, cz, 1.0f, 0.0f, 0.0f);  // Add velocity at center
    // sim.addDensity(cx, cy, cz, 1.0f);               // Add density at center
    
    // // Main simulation loop
    // for (int i = 0; i < 100; i++) {
    //     sim.step();
        
    //     // Print progress every 10 steps
    //     if (i % 10 == 0) {
    //         std::cout << "Step " << i << " completed" << std::endl;
    //     }
    // }
    
    return 0;
} 