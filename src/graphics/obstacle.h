#ifndef OBSTACLE_H
#define OBSTACLE_H

#include <vector>

// Global arrays for obstacle status and obstacle cell indices
extern std::vector<int> g_obstacle;
extern std::vector<int> g_obstacleCells;

// Initialize the obstacle map with total number of cells
void initObstacleMap(int totalCells);

// Clear all obstacles from the map
void clearObstacles(std::vector<float>& density);

// Add a cube-shaped obstacle to the map
void addObstacleCube(int x0, int y0, int z0, int w, int h, int d, int size);

// Keep the density of the obstacle
void densityInsideObstacles(std::vector<float>& density);

// Set special color values for obstacles
void setObstacleColors(std::vector<int>& colorField);


// Check if a cell is an obstacle
bool isObstacle(int idx);

#endif // OBSTACLE_H
