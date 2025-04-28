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
    void clearAllFluids();

    // Below are our project parameters
    FluidCube fluidCube;


    // Obstacle Stuff
    void addObstacle();
    void clearObstacle();


private:

    // Below are the original FEM parameters
    Shape m_shape;

    Shape m_ground;
    void initGround();

    // Grid Size
    int gridSize;

    // Last mouse position
    bool m_mousePressed;
    int m_lastMouseX = -1;
    int m_lastMouseY = -1;

    void addDensityWithGaussian(float centerX, float centerY, float centerZ, float amount, float sigma);
    void addVelocityWithGaussian(float centerX, float centerY, float centerZ,
                                 float velX, float velY, float velZ, float sigma);

    // Helper function to compute Gaussian weight
    float computeGaussianWeight(int dx, int dy, int dz, float sigma);

    // Helper function to check cell validity
    bool isValidGridCell(int nx, int ny, int nz);

};
