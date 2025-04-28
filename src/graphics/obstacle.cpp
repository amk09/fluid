#include "obstacle.h"

// Global arrays definition
std::vector<int> g_obstacle;
std::vector<int> g_obstacleCells;

// Calculate 1D array index from 3D coordinates (x, y, z)
int index(int x, int y, int z, int N) {
    return x + y*N + z*N*N;
}

// Initialize obstacle map with "totalCells" entries, all set to 0 (no obstacle)
void initObstacleMap(int totalCells) {
    g_obstacle.resize(totalCells, 0);
    g_obstacleCells.clear();
}

// Clear all obstacles by resetting to 0 and clearing recorded obstacle cells
void clearObstacles(std::vector<float>& density) {

    for (int idx : g_obstacleCells) {
        if (idx >= 0 && idx < density.size()) {
            density[idx] = 0.0f;
        }
    }

    std::fill(g_obstacle.begin(), g_obstacle.end(), 0);
    g_obstacleCells.clear();
}

// Add a cube obstacle to the map, from starting (x0, y0, z0) with width w, height h, depth d
void addObstacleCube(int x0, int y0, int z0, int w, int h, int d, int size) {
    for (int z = z0; z < z0+d; ++z) {
        for (int y = y0; y < y0+h; ++y) {
            for (int x = x0; x < x0+w; ++x) {
                int idx = index(x, y, z, size);
                // Only add valid indices within bounds
                if (idx >= 0 && idx < size*size*size) {
                    // Mark obstacle if the cell is not already an obstacle
                    if (g_obstacle[idx] == 0) {
                        g_obstacle[idx] = 1;
                        g_obstacleCells.push_back(idx);
                    }
                }
            }
        }
    }
}

// Keep the density of the Obstacle
void densityInsideObstacles(std::vector<float>& density) {
    for (int idx : g_obstacleCells) {
        if (idx >= 0 && idx < density.size()) {
            density[idx] = 1.0f;
        }
    }
}

// Set special color values for obstacles
void setObstacleColors(std::vector<int>& colorField) {
    for (int idx : g_obstacleCells) {
        if (idx >= 0 && idx < colorField.size()) {
            // Use special value 999 to indicate obstacle
            colorField[idx] = 999;
        }
    }
}

// Check if a cell is an obstacle
bool isObstacle(int idx) {
    if (idx < 0 || idx >= g_obstacle.size()) {
        return false;
    }
    return g_obstacle[idx] == 1;
}
