#include "solidboundary.h"
#include <cmath>
#include <algorithm>

void SolidBoundary::addCube(const Eigen::Vector3f& position, const Eigen::Vector3f& size, std::vector<bool>& solidCells, std::vector<SolidObject>& objects) {
    SolidObject obj;
    obj.type = SolidObject::CUBE;
    obj.position = position;
    obj.size = size;
    
    objects.push_back(obj);
    updateSolidCells(solidCells, objects);
}

void SolidBoundary::addSphere(const Eigen::Vector3f& position, float radius, td::vector<bool>& solidCells) {
    SolidObject obj;
    obj.type = SolidObject::SPHERE;
    obj.position = position;
    obj.radius = radius;
    objects.push_back(obj);
    
    // Update solid cells mapping
    updateSolidCells(solidCells, objects);
}

void SolidBoundary::clear(std::vector<bool>& solidCells, std::vector<SolidObject>& objects) {
    objects.clear();
    std::fill(solidCells.begin(), solidCells.end(), false);
}

void SolidBoundary::updateSolidCells(std::vector<bool>& solidCells, const std::vector<SolidObject>& objects, int size) {
    std::fill(m_solidCells.begin(), m_solidCells.end(), false);
    for(int k =0; k<size; k++){
        for(int j =0; j<size; j++){
            for(int i =0; i<size; i++){
                int idx = index(i, j, k);
                if(isInside(i, j, k)){
                    m_solidCells[idx] = true;
                }
            }
        }
    }



}
