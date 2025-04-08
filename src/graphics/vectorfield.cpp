#include "vectorfield.h"

void VectorFields::addVortexField(
    std::vector<float> &vX, std::vector<float> &vY, std::vector<float> &vZ, 
    std::vector<float> &density, int size, float strength) 
{
    //this function creates a vortex field where everything spining around the axis of y.
    int center = size / 2;
    float radius = size / 3.0f;
    
    for (int z = 0; z < size; z++) {
        for (int y = 0; y < size; y++) {
            for (int x = 0; x < size; x++) {
                
                // Vector from center to current point
                float dx = x - center;
                float dy = y - center;
                float dz = z - center;
                
                //we will do the distance-based scaling:
                float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if (dist < radius && dist > 0.1f) {
                    
                    float factor = strength * (1.0f - dist/radius);
                    
                    // Cross product with Y axis to create rotation
                    float vx = -dz * factor / dist;
                    float vz = dx * factor / dist;
                    
                    int idx = index(x, y, z, size);
                    vX[idx] += vx;
                    vZ[idx] += vz;
                    
                    // Add some density in the vortex to see it
                    if (y > center - radius/3 && y < center + radius/3) {
                        density[idx] = std::min(density[idx] + factor * 0.3f, 1.0f);
                    }
                }
            }
        }
    }
}

