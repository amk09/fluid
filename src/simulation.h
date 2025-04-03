#pragma once

#include "graphics/fluidcube.h"
#include "graphics/shape.h"

class Shader;

class Simulation
{
public:
    Simulation();

    void init();

    void update(double seconds);

    void draw(Shader *shader);

    void toggleWire();
private:

    // Below are the original FEM parameters
    Shape m_shape;

    Shape m_ground;
    void initGround();

    // Below are our project parameters
    FluidCube fluidCube;
};
