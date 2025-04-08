#pragma once

#include <vector>
#include <cmath>
#include <random>
#include <algorithm>

class VectorFields {
public:
    //https://www.cambridge.org/core/journals/flow/article/applications-of-the-vortexsurface-field-to-flow-visualization-modelling-and-simulation/C4E0529215FA712FFBD0C2C3DB84ECD8
    static void addVortexField(
        std::vector<float> &vX, std::vector<float> &vY, std::vector<float> &vZ, 
        std::vector<float> &density, int size, float strength = 2.0f);
        
        
    // Helper function to calculate array index from 3D coordinates
    static inline int index(int x, int y, int z, int N) {
        return x + N * (y + N * z);
    }
};