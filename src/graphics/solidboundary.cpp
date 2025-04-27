#include "solidboundary.h"
#include <cmath>
#include <algorithm>

void addCube(const Eigen::Vector3f& position, const Eigen::Vector3f& size, std::vector<bool>& solidCells, std::vector<SolidObject>& objects, int N) {
    SolidObject obj;
    obj.type = SolidObject::CUBE;
    obj.position = position;
    obj.size = size;
    
    objects.push_back(obj);
    updateSolidCells(solidCells, objects, N);
}

void addSphere(const Eigen::Vector3f& position, float radius, std::vector<bool>& solidCells, std::vector<SolidObject>& objects, int N) {
    SolidObject obj;
    obj.type = SolidObject::SPHERE;
    obj.position = position;
    obj.radius = radius;
    objects.push_back(obj);
    
    // Update solid cells mapping
    updateSolidCells(solidCells, objects, N);
}

void clear(std::vector<bool>& solidCells, std::vector<SolidObject>& objects) {
    objects.clear();
    std::fill(solidCells.begin(), solidCells.end(), false);
}

//below is still under development

void updateSolidCells(std::vector<bool>& solidCells, const std::vector<SolidObject>& objects, int size) {
    // First, clear all solid cells
    std::fill(solidCells.begin(), solidCells.end(), false);
    
    // Then, mark cells for each object
    for (const auto& obj : objects) {
        if (obj.type == SolidObject::CUBE) {
            // Handle cube
            Eigen::Vector3f halfSize = obj.size * 0.5f;
            int minX = std::max(0, static_cast<int>(obj.position.x() - halfSize.x()));
            int maxX = std::min(size-1, static_cast<int>(obj.position.x() + halfSize.x()));
            int minY = std::max(0, static_cast<int>(obj.position.y() - halfSize.y()));
            int maxY = std::min(size-1, static_cast<int>(obj.position.y() + halfSize.y()));
            int minZ = std::max(0, static_cast<int>(obj.position.z() - halfSize.z()));
            int maxZ = std::min(size-1, static_cast<int>(obj.position.z() + halfSize.z()));
            
            for (int k = minZ; k <= maxZ; k++) {
                for (int j = minY; j <= maxY; j++) {
                    for (int i = minX; i <= maxX; i++) {
                        solidCells[IX(i, j, k, size)] = true;
                    }
                }
            }
        } 
        else if (obj.type == SolidObject::SPHERE) {
            // Handle sphere
            float radiusSq = obj.radius * obj.radius;
            int radius = static_cast<int>(std::ceil(obj.radius));
            
            int minX = std::max(0, static_cast<int>(obj.position.x() - radius));
            int maxX = std::min(size-1, static_cast<int>(obj.position.x() + radius));
            int minY = std::max(0, static_cast<int>(obj.position.y() - radius));
            int maxY = std::min(size-1, static_cast<int>(obj.position.y() + radius));
            int minZ = std::max(0, static_cast<int>(obj.position.z() - radius));
            int maxZ = std::min(size-1, static_cast<int>(obj.position.z() + radius));
            
            for (int k = minZ; k <= maxZ; k++) {
                for (int j = minY; j <= maxY; j++) {
                    for (int i = minX; i <= maxX; i++) {
                        float dx = i - obj.position.x();
                        float dy = j - obj.position.y();
                        float dz = k - obj.position.z();
                        float distSq = dx*dx + dy*dy + dz*dz;
                        
                        if (distSq <= radiusSq) {
                            solidCells[IX(i, j, k, size)] = true;
                        }
                    }
                }
            }
        }
    }
}




bool isInside(int i, int j, int k, int m_size, std::vector<bool>& solidCells) {
    // Boundary check
    if (!inBounds(i, j, k, m_size))
        return true;  // Treat outside domain as solid
    return solidCells[IX(i, j, k, m_size)];
}

std::vector<BoundaryOffset> createVelocityLUT() {
    std::vector<BoundaryOffset> lut;
    lut.resize(32);
    
    // First 16 entries: fluid cells (DIR_CENTER not set)
    // Cell is fluid with various boundary configurations
    for (int i = 0; i < 16; i++) {
        lut[i] = {0, 0, 0, 1.0f, 1.0f, 1.0f};  // Default handling for fluid cells
    }
    
    // Special case for fluid cell completely surrounded by boundaries
    lut[15] = {0, 0, 0, 0.0f, 0.0f, 0.0f};  // All directions are boundaries
    
    // Second 16 entries: solid cells (DIR_CENTER is set)
    // These are the cases where the center cell is solid
    
    lut[16] = {0, 0, 0, -1.0f, -1.0f, -1.0f};  // Lone solid cell (should never happen)
    
    // Single fluid neighbor cases
    lut[17] = {-1, 0, 0, -1.0f, -1.0f, -1.0f};  // East neighbor is fluid
    lut[18] = {0, -1, 0, -1.0f, -1.0f, -1.0f};  // North neighbor is fluid
    lut[20] = {1, 0, 0, -1.0f, -1.0f, -1.0f};   // West neighbor is fluid
    lut[24] = {0, 1, 0, -1.0f, -1.0f, -1.0f};   // South neighbor is fluid

 
    // Two fluid neighbors
    lut[19] = {-1, -1, 0, -1.0f, -1.0f, -1.0f}; // East and North are fluid
    lut[22] = {1, -1, 0, -1.0f, -1.0f, -1.0f};  // West and North are fluid
    lut[25] = {-1, 1, 0, -1.0f, -1.0f, -1.0f};  // East and South are fluid
    lut[28] = {1, 1, 0, -1.0f, -1.0f, -1.0f};   // West and South are fluid
    
    // Handle cases where opposing sides are fluid
    lut[21] = {0, -1, 0, -1.0f, -1.0f, -1.0f};  // East and West are fluid
    lut[26] = {-1, 0, 0, -1.0f, -1.0f, -1.0f};  // North and South are fluid
    
    // Three fluid neighbors
    lut[23] = {0, -1, 0, -1.0f, -1.0f, -1.0f};  // East, North, and West are fluid
    lut[27] = {0, 1, 0, -1.0f, -1.0f, -1.0f};   // East, West, and South are fluid
    lut[29] = {1, 0, 0, -1.0f, -1.0f, -1.0f};   // North, West, South are fluid
    lut[30] = {-1, 0, 0, -1.0f, -1.0f, -1.0f};  // East, North, South are fluid
    
    // All sides fluid
    lut[31] = {0, 0, 0, 0.0f, 0.0f, 0.0f};      // All neighbors are fluid
    
    return lut;
}


void updateBoundaryMask(std::vector<bool>& solidCells, std::vector<int>& boundaryTypes, int size) {
    boundaryTypes.resize(size * size * size, 0);
    
    for (int k = 0; k < size; k++) {
        for (int j = 0; j < size; j++) {
            for (int i = 0; i < size; i++) {
                int idx = IX(i, j, k, size);
                
                // Skip non-solid cells
                if (!solidCells[idx]) continue;
                
                int boundaryType = DIR_CENTER;  // This is a solid cell
                
                // Check which neighbors are fluid
                if (i > 0 && !solidCells[IX(i-1, j, k, size)]) boundaryType |= DIR_WEST;
                if (i < size-1 && !solidCells[IX(i+1, j, k, size)]) boundaryType |= DIR_EAST;
                if (j > 0 && !solidCells[IX(i, j-1, k, size)]) boundaryType |= DIR_SOUTH;
                if (j < size-1 && !solidCells[IX(i, j+1, k, size)]) boundaryType |= DIR_NORTH;
                // Add front/back directions for 3D (modify DIR constants accordingly)
                
                boundaryTypes[idx] = boundaryType;
            }
        }
    }
}

std::vector<BoundaryOffset> createDensityLUT() {
    std::vector<BoundaryOffset> lut;
    lut.resize(32);
    
    // Initialize with defaults
    for (int i = 0; i < 32; i++) {
        lut[i] = {0, 0, 0, 1.0f, 1.0f, 1.0f};
    }
    
    // First 16 entries: fluid cells (DIR_CENTER not set)
    // For fluid cells, we generally keep the density as is
    
    // Second 16 entries: solid cells (DIR_CENTER is set)
    // For solid cells, we generally want to zero out density inside obstacles
    for (int i = 16; i < 32; i++) {
        lut[i] = {0, 0, 0, 0.0f, 0.0f, 0.0f};  // Zero density inside solids
    }
    
    // For solid cells adjacent to fluid, you might want to sample from fluid neighbors
    // This creates a smoother visualization near boundaries
    lut[DIR_CENTER | DIR_EAST] = {1, 0, 0, 0.5f, 0.0f, 0.0f};    // Sample from east, at lower intensity
    lut[DIR_CENTER | DIR_NORTH] = {0, 1, 0, 0.0f, 0.5f, 0.0f};   // Sample from north, at lower intensity
    lut[DIR_CENTER | DIR_WEST] = {-1, 0, 0, 0.5f, 0.0f, 0.0f};   // Sample from west, at lower intensity
    lut[DIR_CENTER | DIR_SOUTH] = {0, -1, 0, 0.0f, 0.5f, 0.0f};  // Sample from south, at lower intensity
    
    return lut;
}

void applyVelocityBoundaries(std::vector<float>& vX, std::vector<float>& vY, std::vector<float>& vZ,
    const std::vector<bool>& solidCells, const std::vector<int>& boundaryTypes,
    const std::vector<BoundaryOffset>& velocityLUT, int size) {

    for (int k = 1; k < size - 1; k++) {
        for (int j = 1; j < size - 1; j++) {
            for (int i = 1; i < size - 1; i++) {
                int idx = IX(i, j, k, size);

                // Skip non-boundary cells
                if (!solidCells[idx]) continue;

                int boundaryType = boundaryTypes[idx];
                const BoundaryOffset& offset = velocityLUT[boundaryType];

                // Calculate neighbor position to sample from
                int ni = i + offset.dx;
                int nj = j + offset.dy;
                int nk = k + offset.dz;

                // Ensure we don't go out of bounds
                ni = std::max(0, std::min(ni, size - 1));
                nj = std::max(0, std::min(nj, size - 1));
                nk = std::max(0, std::min(nk, size - 1));

                int neighborIdx = IX(ni, nj, nk, size);

                // Apply velocity transfer with reflection
                vX[idx] = vX[neighborIdx] * offset.xMultiplier;
                vY[idx] = vY[neighborIdx] * offset.yMultiplier;
                vZ[idx] = vZ[neighborIdx] * offset.zMultiplier;
            }
        }
    }
}

void applyDensityBoundaries(std::vector<float>& density,
    const std::vector<bool>& solidCells, 
    const std::vector<int>& boundaryTypes,
    const std::vector<BoundaryOffset>& densityLUT, 
    int size) {

    for (int k = 1; k < size - 1; k++) {
        for (int j = 1; j < size - 1; j++) {
            for (int i = 1; i < size - 1; i++) {
            int idx = IX(i, j, k, size);

            // Skip non-boundary cells
            if (!solidCells[idx]) continue;

            int boundaryType = boundaryTypes[idx];
            const BoundaryOffset& offset = densityLUT[boundaryType];

            // Calculate neighbor position
            int ni = i + offset.dx;
            int nj = j + offset.dy;
            int nk = k + offset.dz;

            // Ensure we don't go out of bounds
            ni = std::max(0, std::min(ni, size - 1));
            nj = std::max(0, std::min(nj, size - 1));
            nk = std::max(0, std::min(nk, size - 1));

            int neighborIdx = IX(ni, nj, nk, size);

            // Apply density with multipliers
            // This allows for partial density near boundaries
            density[idx] = density[neighborIdx] * 
                    std::max(offset.xMultiplier, 
                            std::max(offset.yMultiplier, offset.zMultiplier));
            }
        }
    }
}


