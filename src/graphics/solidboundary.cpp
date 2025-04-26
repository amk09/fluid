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
    std::fill(solidCells.begin(), solidCells.end(), false);
    for(int k =0; k<size; k++){
        for(int j =0; j<size; j++){
            for(int i =0; i<size; i++){
                int idx = index(i, j, k, size);
                if(isInside(i, j, k)){
                    solidCells[idx] = true;
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

Eigen::Vector3f getNormal(int i, int j, int k, int m_size) {
    // Calculate normal by checking neighboring cells
    Eigen::Vector3f normal(0.0f, 0.0f, 0.0f);
    
    // Check six neighbors
    if (i > 0 && !isInside(i-1, j, k)) normal[0] -= 1.0f;
    if (i < m_size-1 && !isInside(i+1, j, k)) normal[0] += 1.0f;
    
    if (j > 0 && !isInside(i, j-1, k)) normal[1] -= 1.0f;
    if (j < m_size-1 && !isInside(i, j+1, k)) normal[1] += 1.0f;
    
    if (k > 0 && !isInside(i, j, k-1)) normal[2] -= 1.0f;
    if (k < m_size-1 && !isInside(i, j, k+1)) normal[2] += 1.0f;
    
    // Normalize if not zero
    float length = normal.norm();
    if (length > 0.001f) {
        normal /= length;
    }
    
    return normal;
}
