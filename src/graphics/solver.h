#ifndef SOLVER_H
#define SOLVER_H
#include "fluidcube.h"

// b here is coord number (setting up boundry!), where 1 for x-coord, 2 for y-coord, 3 for z-coord (0 for nothing)
// diffuse will be (visc for velocity diffuse, diff for density diffuse)
// iter and N will not be used here because we already have this in this class
void diffuse_velocity(int b, vector<float> &dataWrittenTo, vector<float> &dataReadFrom, float viscosity, float dt, int iter, int N);

void diffuse_density(int b, vector<float> &dataWrittenTo, vector<float> &dataReadFrom, float diffuse, float dt, int iter, int N);

// a & c are the parameters required in the paper.
// iter and N will not be used here because we already have this in this class
void lin_solve(int b, vector<float> &dataWrittenTo, vector<float> &dataReadFrom, float a, float c, int iter, int N);

// iter and N will not be used here because we already have this in this class
void project(vector<float> &velocityX, vector<float> &velocityY, vector<float> &velocityZ, vector<float> &p, vector<float> div, int iter, int N);

// b is the same meaning shown above
// N will not be used here because we already have this in this class
void advect(int b, vector<float> &dataWrittenTo, vector<float> &dataReadFrom, vector<float> &velocityX, vector<float> &velocityY, vector<float> &velocityZ, float dt, int N);


// b is the same meaning shown above
// N will not be used here because we already have this in this class
void set_bnd(int b, vector<float> &dataWrittenTo, int N);


//handle the density conservation
void conserveDensity(std::vector<float> &density, float targetDensity, int N);


void addFountainForce(std::vector<float> &velocityX, std::vector<float> &velocityY, 
    std::vector<float> &velocityZ, std::vector<float> &density, 
    int N, float strength, float dt);

void circulateDensity(std::vector<float> &density, int N, float dt, float rate);

void addVorticityConfinement(std::vector<float> &vX, std::vector<float> &vY,
                             std::vector<float> &vZ, float dt, int N, float vorticityStrength);

#endif // SOLVER_H
