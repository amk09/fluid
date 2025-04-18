#pragma once

#include "graphics/fluidcube.h"
#include "graphics/shape.h"
#include <QMouseEvent>

class Shader;

class Simulation
{
public:



    Simulation();

    void init();

    void update(double seconds);

    void draw(Shader *shader);

    void toggleWire();

    // Add mouse event handlers
    void handleMousePress(int x, int y, int width, int height);
    void handleMouseMove(int x, int y, int width, int height);


    // Below are our project parameters
    FluidCube fluidCube;


private:

    // Below are the original FEM parameters
    Shape m_shape;

    Shape m_ground;
    void initGround();



    // Last mouse position
    bool m_mousePressed;
    int m_lastMouseX = -1;
    int m_lastMouseY = -1;

    void addDensityWithGaussian(float centerX, float centerY, float centerZ, float amount, float sigma);
    void addVelocityWithGaussian(float centerX, float centerY, float centerZ,
                                 float velX, float velY, float velZ, float sigma);
};
