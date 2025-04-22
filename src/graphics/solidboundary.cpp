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
                // int idx = index(i, j, k);
                if(isInside(i, j, k)){
                    solidCells[idx] = true;
                }
            }
        }
    }



}
