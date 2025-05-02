#ifndef OBSTACLE_H
#define OBSTACLE_H

#include <vector>

// Global arrays for obstacle status and obstacle cell indices
extern std::vector<int> g_obstacle;
extern std::vector<int> g_obstacleCells;

// Fruit type enum with color codes
enum FruitType {
    PINEAPPLE = 1000,       // Yellow pineapple body
    PINEAPPLE_PATTERN = 1001, // Darker pattern on pineapple
    PINEAPPLE_CROWN = 1002,  // Green pineapple crown/leaves

    ORANGE = 1010,      // Orange fruit
    ORANGE_STEM = 1011, // Brown orange stem
    ORANGE_LEAF = 1012, // Green orange leaf

    GRAPE = 1020,       // Dark purple grape
    GRAPE_STEM = 1021,  // Brown grape stem
    GRAPE_LEAF = 1022,  // Green grape leaf
    GRAPE_HIGHLIGHT = 1023, // Light purple highlight

    PEAR = 1040,        // Yellow-green pear
    PEAR_STEM = 1041,   // Brown pear stem

    // Cherry types
    CHERRY = 1050,       // Red cherry fruit
    CHERRY_STEM = 1051,  // Green cherry stem
    CHERRY_LEAF = 1052,  // Green cherry leaf

};

// Initialize the obstacle map with total number of cells
void initObstacleMap(int totalCells);

// Clear all obstacles from the map
void clearObstacles(std::vector<float>& density);

// Add fruit-shaped obstacles
void addApple(int x0, int y0, int z0, int size);
void addOrange(int x0, int y0, int z0, int size);
void addGrapes(int x0, int y0, int z0, int size);
void addCherry(int x0, int y0, int z0, int size);
void addPear(int x0, int y0, int z0, int size);


// Original function - now creates multiple fruits
void addObstacleCube(int x0, int y0, int z0, int w, int h, int d, int size);

// Keep the density of the obstacle
void densityInsideObstacles(std::vector<float>& density);

// Set special color values for obstacles
void setObstacleColors(std::vector<int>& colorField);

// Check if a cell is an obstacle
bool isObstacle(int idx);

// Helper function to check if a point is inside a sphere
bool isInsideSphere(int x, int y, int z, int centerX, int centerY, int centerZ, float radius);

#endif // OBSTACLE_H
