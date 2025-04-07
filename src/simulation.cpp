#include "simulation.h"
#include "graphics/meshloader.h"

#include <iostream>

using namespace Eigen;

Simulation::Simulation() {}

void Simulation::init()
{

    fluidCube.init(256, 0.0, 0.0);
}

void Simulation::update(double seconds)
{
    // Fluid Cube will handle its own update
    fluidCube.update(seconds);
}

void Simulation::draw(Shader *shader)
{
    // m_shape.draw(shader);
    // m_ground.draw(shader);

    // Fluid Cube will handle its own draw
    fluidCube.draw(shader);
}

void Simulation::toggleWire()
{
    // m_shape.toggleWireframe();

    fluidCube.toggleWireframe();
}

void Simulation::initGround()
{
    std::vector<Vector3d> groundVerts;
    std::vector<Vector3i> groundFaces;
    groundVerts.emplace_back(-5, 0, -5);
    groundVerts.emplace_back(-5, 0, 5);
    groundVerts.emplace_back(5, 0, 5);
    groundVerts.emplace_back(5, 0, -5);
    groundFaces.emplace_back(0, 1, 2);
    groundFaces.emplace_back(0, 2, 3);
    m_ground.init(groundVerts, groundFaces);
}
