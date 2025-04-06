#ifndef SOLVER_H
#define SOLVER_H
#include "fluidcube.h"

// b here is coord number (setting up boundry!), where 1 for x-coord, 2 for y-coord, 3 for z-coord (0 for nothing)
// diffuse will be (visc for velocity diffuse, diff for density diffuse)
// iter and N will not be used here because we already have this in this class
void diffuse_velocity(int b, vector<float> &dataWrittenTo, vector<float> &dataReadFrom, float diffuse, float dt);

void diffuse_density(int b, vector<float> &dataWrittenTo, vector<float> &dataReadFrom, float diffuse, float dt);



#endif // SOLVER_H