#ifndef SOLIDBOUNDARY_H
#define SOLIDBOUNDARY_H
#include "fluidcube.h"


void addCube(const Eigen::Vector3f& position, const Eigen::Vector3f& size, 
                std::vector<bool>& solidCells, std::vector<SolidObject>& objects, int N);
void addSphere(const Eigen::Vector3f& position, float radius, 
                std::vector<bool>& solidCells, std::vector<SolidObject>& objects, int N);
void clear(std::vector<bool>& solidCells, std::vector<SolidObject>& objects);

void updateSolidCells(std::vector<bool>& solidCells, const std::vector<SolidObject>& objects, int size);

bool isInside(int i, int j, int k);
Eigen::Vector3f getNormal(int i, int j, int k, int m_size);

void applyBoundaries(int b, std::vector<float>& field);

inline int IX(int i, int j, int k, int m_size)  {
    return i + j * m_size + k * m_size * m_size;
}

// Helper function to check if index is in bounds
inline bool inBounds(int i, int j, int k, int m_size)  {
    return (i >= 0 && i < m_size && 
            j >= 0 && j < m_size && 
            k >= 0 && k < m_size);
}
inline int index(int x, int y, int z, int size) {
    return (x) + (y) * size + (z) * size * size;
}

#endif // SOLIDBOUNDARY_H

