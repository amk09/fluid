#ifndef SOLIDBOUNDARY_H
#define SOLIDBOUNDARY_H
#include "fluidcube.h"


void addCube(const Eigen::Vector3f& position, const Eigen::Vector3f& size, 
                std::vector<bool>& solidCells, std::vector<SolidObject>& objects, int N);
void addSphere(const Eigen::Vector3f& position, float radius, 
                std::vector<bool>& solidCells, std::vector<SolidObject>& objects, int N);
void clear(std::vector<bool>& solidCells, std::vector<SolidObject>& objects);

void updateSolidCells(std::vector<bool>& solidCells, const std::vector<SolidObject>& objects, int size);
bool isInside(int i, int j, int k) const;


Eigen::Vector3f getNormal(int i, int j, int k) const;
void applyBoundaries(int b, std::vector<float>& field);




#endif // SOLIDBOUNDARY_H

