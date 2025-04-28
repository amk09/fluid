#ifndef OBSTACLE_H
#define OBSTACLE_H

#include <vector>

// Global arrays for obstacle status and obstacle cell indices
extern std::vector<int> g_obstacle;
extern std::vector<int> g_obstacleCells;

// Initialize the obstacle map with total number of cells
void initObstacleMap(int totalCells);

// Clear all obstacles from the map
void clearObstacles();

// Add a cube-shaped obstacle to the map
void addObstacleCube(int x0, int y0, int z0, int w, int h, int d, int size);

// Clear the density inside obstacle cells to zero
void clearDensityInsideObstacles(std::vector<float>& density);

#endif // OBSTACLE_H
