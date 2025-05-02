#include "obstacle.h"
#include <cmath>
#include <iostream>
#include <algorithm>

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

// Helper function to check if a point is inside a sphere
bool isInsideSphere(int x, int y, int z, int centerX, int centerY, int centerZ, float radius) {
    float dx = x - centerX;
    float dy = y - centerY;
    float dz = z - centerZ;
    return (dx*dx + dy*dy + dz*dz) <= (radius*radius);
}

// Helper function to add a cell as an obstacle with a specific type
void addObstacleCell(int x, int y, int z, int size, FruitType type) {
    // Basic boundary check
    if (x < 0 || x >= size || y < 0 || y >= size || z < 0 || z >= size) {
        return;
    }

    int idx = index(x, y, z, size);
    if (idx >= 0 && idx < size*size*size) {
        g_obstacle[idx] = type;
        g_obstacleCells.push_back(idx);
    }
}

void addPineapple(int x0, int y0, int z0, int size) {
    int gridSize = g_obstacle.size();
    int n = static_cast<int>(std::cbrt(gridSize));

    float baseRadius = std::min(size / 2.5f, n/8.0f);
    float bodyRadius = baseRadius * 1.1f;
    float bodyHeight = baseRadius * 2.0f;

    // body
    for (int z = z0 - bodyRadius; z <= z0 + bodyRadius; z++) {
        for (int y = y0 - bodyHeight/2; y <= y0 + bodyHeight/2; y++) {
            for (int x = x0 - bodyRadius; x <= x0 + bodyRadius; x++) {

                float dx = (x - x0) / bodyRadius;
                float dy = (y - y0) / bodyHeight;
                float dz = (z - z0) / bodyRadius;

                float distSquared = dx*dx + dy*dy*4 + dz*dz;

                if (distSquared <= 1.0f) {
                    int patternX = static_cast<int>((x - x0 + 100) / (bodyRadius/4)) % 2;
                    int patternY = static_cast<int>((y - y0 + 100) / (bodyHeight/6)) % 2;

                    if ((patternX + patternY) % 2 == 0) {
                        addObstacleCell(x, y, z, n, PINEAPPLE);
                    } else {
                        addObstacleCell(x, y, z, n, PINEAPPLE_PATTERN);
                    }
                }
            }
        }
    }


    int leafBaseY = y0 + bodyHeight/2 + 4.5;

    int leafPattern[4][9] = {
        {0, 0, 0, 0, 2, 1, 2, 0, 0},
        {0, 0, 0, 2, 1, 1, 1, 2, 0},
        {0, 0, 2, 1, 0, 1, 0, 1, 2},
        {0, 0, 0, 0, 0, 1, 0, 0, 0}
    };

    float leafScale = baseRadius / 5.0f;


    for (int row = 3; row >= 0; row--) {
        for (int col = 0; col < 9; col++) {
            if (leafPattern[row][col] > 0) {

                int startX = x0 + (col - 4.5) * ceil(leafScale);
                int endX = x0 + (col - 4.5 + 1) * ceil(leafScale);


                int startY = leafBaseY - (3-row) * ceil(leafScale);
                int endY = leafBaseY - (3-row-1) * ceil(leafScale);


                for (int y = startY; y < endY; y++) {
                    for (int x = startX; x < endX; x++) {
                        for (int z = z0 - ceil(leafScale/2); z <= z0 + ceil(leafScale/2); z++) {
                            if (leafPattern[row][col] == 1) {
                                addObstacleCell(x, y, z, n, PINEAPPLE_CROWN);
                            } else {
                                addObstacleCell(x, y, z, n, GRAPE_LEAF);
                            }
                        }
                    }
                }
            }
        }
    }
}


void addGrapes(int x0, int y0, int z0, int size) {
    int gridSize = g_obstacle.size();
    int n = static_cast<int>(std::cbrt(gridSize));

    // Scale factor - adjust for appropriate size
    float scale = 1.0f;

    // Define the grape pattern with slightly more rounded bottom
    // 0=empty, 1=dark purple grape, 2=light purple highlight, 3=brown stem
    int grapePattern[16][12] = {
        {0,0,0,0,0,3,0,0,0,0,0,0}, // Row 1
        {0,0,0,0,3,3,3,3,0,0,0,0}, // Row 2
        {0,0,3,3,0,3,0,0,3,0,0,0}, // Row 3
        {0,0,0,0,0,3,0,0,0,0,0,0}, // Row 4
        {0,0,0,0,1,1,1,0,0,0,0,0}, // Row 5
        {0,0,0,1,1,2,1,1,1,0,0,0}, // Row 6
        {0,1,1,2,1,1,1,1,1,1,0,0}, // Row 7
        {1,1,1,1,1,2,1,1,2,1,1,0}, // Row 8
        {1,1,1,1,1,1,1,1,1,1,1,0}, // Row 9
        {0,1,1,2,1,1,1,1,2,1,0,0}, // Row 10
        {0,1,1,1,1,1,2,1,1,1,0,0}, // Row 11
        {0,0,1,1,1,2,1,1,1,0,0,0}, // Row 12
        {0,0,0,1,1,1,1,1,0,0,0,0}, // Row 13
        {0,0,0,0,1,1,1,0,0,0,0,0}, // Row 14 - Kept original
        {0,0,0,0,0,1,0,0,0,0,0,0}, // Row 15 - Just made slightly rounder
        {0,0,0,0,0,0,0,0,0,0,0,0}  // Row 16
    };

    // Create grape from pattern - FLIPPED VERTICALLY
    for (int row = 0; row < 16; row++) {
        for (int col = 0; col < 12; col++) {
            if (grapePattern[row][col] > 0) {
                // Calculate position - FLIP by using (15-row) instead of row
                int px = x0 + (col - 6) * scale;
                int py = y0 + (15-row - 8) * scale; // Flipped vertically

                // Add depth for 3D appearance
                for (int z = z0 - 1; z <= z0 + 1; z++) {
                    // Select color type based on pattern
                    if (grapePattern[row][col] == 1) {
                        // Dark purple grape
                        addObstacleCell(px, py, z, n, GRAPE);
                    }
                    else if (grapePattern[row][col] == 2) {
                        // Light purple highlight
                        addObstacleCell(px, py, z, n, GRAPE_HIGHLIGHT);
                    }
                    else if (grapePattern[row][col] == 3) {
                        // Brown stem
                        addObstacleCell(px, py, z, n, GRAPE_STEM);
                    }
                }
            }
        }
    }
}

void addPear(int x0, int y0, int z0, int size) {
    int gridSize = g_obstacle.size();
    int n = static_cast<int>(std::cbrt(gridSize));
    // Use appropriate radius for pear
    float baseRadius = std::min(size / 2.5f, n/8.0f);
    float centerDist = baseRadius * 0.6f;
    float topRadius = baseRadius * 0.6f;
    float bottomRadius = baseRadius * 1.1f;

    // Top part of pear (smaller ball) - NOW AT TOP (positive Y)
    for (int z = z0 - topRadius; z <= z0 + topRadius; z++) {
        for (int y = y0 + centerDist - topRadius; y <= y0 + centerDist + topRadius; y++) {
            for (int x = x0 - topRadius; x <= x0 + topRadius; x++) {
                float dx = x - x0;
                float dy = y - (y0 + centerDist);
                float dz = z - z0;

                // Calculate distance for sphere
                float distance = dx*dx + dy*dy + dz*dz;

                // Create top part of pear
                if (distance <= topRadius*topRadius) {
                    addObstacleCell(x, y, z, n, PEAR);
                }
            }
        }
    }

    // Bottom part of pear (larger ball) - NOW AT BOTTOM (negative Y)
    for (int z = z0 - bottomRadius; z <= z0 + bottomRadius; z++) {
        for (int y = y0 - centerDist - bottomRadius; y <= y0 - centerDist + bottomRadius; y++) {
            for (int x = x0 - bottomRadius; x <= x0 + bottomRadius; x++) {
                float dx = x - x0;
                float dy = y - (y0 - centerDist);
                float dz = z - z0;

                // Calculate distance for sphere
                float distance = dx*dx + dy*dy + dz*dz;

                // Create bottom part of pear
                if (distance <= bottomRadius*bottomRadius) {
                    addObstacleCell(x, y, z, n, PEAR);
                }
            }
        }
    }

    // Connect the two parts with a smooth neck - NOW CORRECTLY ORIENTED
    for (int z = z0 - topRadius; z <= z0 + topRadius; z++) {
        for (int y = y0 - centerDist + bottomRadius; y <= y0 + centerDist - topRadius; y++) {
            for (int x = x0 - topRadius; x <= x0 + topRadius; x++) {
                // Calculate how far we are between the two centers (0 = bottom sphere, 1 = top sphere)
                float t = (y - (y0 - centerDist)) / (2 * centerDist);

                // Interpolate radius based on position
                float radius = bottomRadius * (1-t) + topRadius * t;

                float dx = x - x0;
                float dz = z - z0;

                // Create connecting neck
                if (dx*dx + dz*dz <= radius*radius) {
                    addObstacleCell(x, y, z, n, PEAR);
                }
            }
        }
    }

    // Add stem at top of pear
    int stemHeight = std::max(2, static_cast<int>(baseRadius / 3.0f));
    int stemWidth = std::max(1, static_cast<int>(baseRadius / 5.0f));

    // Position stem at the top of the pear (positive Y values)
    for (int y = y0 + centerDist + topRadius; y < y0 + centerDist + topRadius + stemHeight; y++) {
        for (int x = x0 - stemWidth/2; x <= x0 + stemWidth/2; x++) {
            for (int z = z0 - stemWidth/2; z <= z0 + stemWidth/2; z++) {
                addObstacleCell(x, y, z, n, PEAR_STEM);
            }
        }
    }

}

void addCherry(int x0, int y0, int z0, int size) {
    int gridSize = g_obstacle.size();
    int n = static_cast<int>(std::cbrt(gridSize));

    // Cherry pattern with rounder right cherry
    // 0=empty, 1=red cherry, 2=dark green stem, 3=light green leaf
    int cherryPattern[15][12] = {
        {0,0,0,0,0,2,0,0,0,0,0,0}, // Row 1
        {0,0,0,0,2,2,0,0,0,0,0,0}, // Row 2
        {0,0,0,0,2,0,3,3,3,0,0,0}, // Row 3
        {0,0,0,2,0,3,3,3,3,3,0,0}, // Row 4
        {0,0,0,2,0,3,3,3,3,0,0,0}, // Row 5
        {0,0,2,0,0,0,2,0,0,0,0,0}, // Row 6
        {0,0,2,0,0,2,0,0,0,0,0,0}, // Row 7
        {0,0,2,0,0,0,0,0,0,0,0,0}, // Row 8
        {0,0,2,0,1,0,0,0,1,0,0,0}, // Row 9
        {0,1,0,1,1,1,0,1,1,1,1,0}, // Row 10 - Extended right side to make cherry rounder
        {1,1,1,1,1,1,1,1,1,1,1,1}, // Row 11 - Extended furthest right to make it rounder
        {1,1,1,1,1,1,0,0,1,1,1,1}, // Row 12
        {0,1,1,1,1,1,0,0,1,1,1,1}, // Row 13 - Extended right cherry
        {0,0,1,1,1,0,0,0,0,1,1,0}, // Row 14 - Extended right cherry
        {0,0,0,0,0,0,0,0,0,0,0,0}  // Row 15
    };

    // Scale factor - keep cherry small
    float scale = 0.8f;

    // Create cherry from pattern - FLIPPED VERTICALLY to fix upside-down issue
    for (int row = 0; row < 15; row++) {
        for (int col = 0; col < 12; col++) {
            if (cherryPattern[row][col] > 0) {
                // Calculate position - FLIP by using (14-row) instead of row
                int px = x0 + (col - 6) * scale;
                int py = y0 + (14-row - 7) * scale; // Flipped vertically

                // Add minimal depth for a small cherry
                for (int z = z0 - 1; z <= z0; z++) {
                    // Select color type based on pattern
                    if (cherryPattern[row][col] == 1) {
                        // Red cherry
                        addObstacleCell(px, py, z, n, CHERRY);
                    }
                    else if (cherryPattern[row][col] == 2) {
                        // Green stem
                        addObstacleCell(px, py, z, n, CHERRY_STEM);
                    }
                    else if (cherryPattern[row][col] == 3) {
                        // Green leaf
                        addObstacleCell(px, py, z, n, CHERRY_LEAF);
                    }
                }
            }
        }
    }
}

void addOrange(int x0, int y0, int z0, int size) {
    int gridSize = g_obstacle.size();
    int n = static_cast<int>(std::cbrt(gridSize));

    // Scale factor - adjust for appropriate size
    float scale = 1.0f;

    // Define orange pattern based on reference image without black outline
    // 0=empty, 1=orange, 3=brown stem, 4=green leaf
    int orangePattern[16][16] = {
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // Row 1
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // Row 2
        {0,0,0,0,0,0,0,4,4,4,0,0,0,0,0,0}, // Row 3
        {0,0,0,0,0,0,4,4,4,4,4,0,0,0,0,0}, // Row 4
        {0,0,0,0,0,0,4,4,4,4,0,0,0,0,0,0}, // Row 5
        {0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0}, // Row 6
        {0,0,0,0,0,0,3,3,0,0,0,0,0,0,0,0}, // Row 7
        {0,0,0,0,0,1,1,3,1,1,0,0,0,0,0,0}, // Row 8 - Changed from outline to orange
        {0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0}, // Row 9
        {0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0}, // Row 10
        {0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0}, // Row 11
        {0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0}, // Row 12
        {0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0}, // Row 13
        {0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0}, // Row 14
        {0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0}, // Row 15
        {0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0}  // Row 16
    };

    // Create orange from pattern - FLIPPED VERTICALLY to fix upside-down issue
    for (int row = 0; row < 16; row++) {
        for (int col = 0; col < 16; col++) {
            if (orangePattern[row][col] > 0) {
                // Calculate position - FLIP by using (15-row) instead of row
                int px = x0 + (col - 8) * scale;
                int py = y0 + (15-row - 8) * scale; // Flipped vertically

                // Add depth for 3D appearance
                for (int z = z0 - 1; z <= z0 + 1; z++) {
                    // Select color type based on pattern
                    if (orangePattern[row][col] == 1) {
                        // Orange fruit
                        addObstacleCell(px, py, z, n, ORANGE);
                    }
                    else if (orangePattern[row][col] == 3) {
                        // Brown stem
                        addObstacleCell(px, py, z, n, ORANGE_STEM);
                    }
                    else if (orangePattern[row][col] == 4) {
                        // Green leaf
                        addObstacleCell(px, py, z, n, ORANGE_LEAF);
                    }
                }
            }
        }
    }
}

void addObstacleCube(int x0, int y0, int z0, int w, int h, int d, int size) {
    // Clear previous obstacles first
    std::fill(g_obstacle.begin(), g_obstacle.end(), 0);
    g_obstacleCells.clear();

    // Calculate center of the cube
    int centerX = x0 + w/2;
    int centerY = y0 + h/2;
    int centerZ = z0 + d/2;

    // Use smaller size for fruits to ensure they fit
    int fruitSize = std::min(w, std::min(h, d));

    // Calculate distances for fruit placement
    int horizDist = size/3; //horizontal distance
    int vertDist = size/2.5;  // vertical distance

    // Adjust vertical spacing for grape
    int grapeOffset = vertDist * 0.7; // Reduced vertical offset for grape


    // Pineapple
    int pineappleX = centerX - horizDist;
    int pineappleY = centerY * 0.8;
    int pineappleZ = centerZ;
    addPineapple(pineappleX, pineappleY, pineappleZ, fruitSize);

    // Orange
    int orangeX = centerX + size*0.3;
    int orangeY = centerY + size*0.1;
    int orangeZ = centerZ;
    addOrange(orangeX, orangeY, orangeZ, fruitSize);

    // Grape
    int grapeX = centerX;
    int grapeY = centerY + grapeOffset;
    int grapeZ = centerZ;
    addGrapes(grapeX, grapeY, grapeZ, fruitSize);

    // Cherry
    int cherryX = centerX - size*0.15;
    int cherryY = centerY - vertDist * 0.8;
    int cherryZ = centerZ  - vertDist/2;
    addCherry(cherryX, cherryY, cherryZ, fruitSize);

    // Pear
    int pearX = centerX + vertDist/2;
    int pearY = centerY - vertDist * 0.8;
    int pearZ = centerZ - vertDist/2;
    addPear(pearX, pearY, pearZ, fruitSize);

    std::cout << "Added fruits with adjusted spacing" << std::endl;
    std::cout << "Total obstacle cells: " << g_obstacleCells.size() << std::endl;
}










// Keep the density of the obstacle
void densityInsideObstacles(std::vector<float>& density) {
    for (int idx : g_obstacleCells) {
        if (idx >= 0 && idx < density.size()) {
            density[idx] = 1.0f;
        }
    }
}

// Set special color values for obstacles based on fruit type
void setObstacleColors(std::vector<int>& colorField) {
    for (int idx : g_obstacleCells) {
        if (idx >= 0 && idx < colorField.size()) {
            // Use fruit type as color identifier
            colorField[idx] = g_obstacle[idx];
        }
    }
}

// Check if a cell is an obstacle
bool isObstacle(int idx) {
    if (idx < 0 || idx >= g_obstacle.size()) {
        return false;
    }
    return g_obstacle[idx] != 0;
}
