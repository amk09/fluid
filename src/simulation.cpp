#include "simulation.h"
#include "graphics/meshloader.h"

#include <iostream>

using namespace Eigen;

Simulation::Simulation() {}

void Simulation::init()
{

    fluidCube.init(32, 0.001, 0.001);
}


void Simulation::handleMousePress(int x, int y, int width, int height)
{
    m_lastMouseX = x;
    m_lastMouseY = y;

    // Convert screen coordinates to normalized coordinates (0-1)
    float normalizedX = 1.0f - static_cast<float>(x) / width;  // Flip X to fix mirroring
    float normalizedY = 1.0f - static_cast<float>(y) / height;  // Also flip Y to match OpenGL coordinates

    // Map to fluid grid coordinates
    int gridSize = 32; // Match the size used in init()
    int gridX = static_cast<int>(normalizedX * gridSize);
    int gridY = static_cast<int>(normalizedY * gridSize);
    int gridZ = gridSize / 2; // Middle of the grid

    // Add density at this position with surrounding area for better visibility
    if (gridX >= 1 && gridX < gridSize-1 &&
        gridY >= 1 && gridY < gridSize-1 &&
        gridZ >= 1 && gridZ < gridSize-1) {
        // Add density in a small area around the click point
        for(int dx = -2; dx <= 2; dx++) {
            for(int dy = -2; dy <= 2; dy++) {
                int nx = gridX + dx;
                int ny = gridY + dy;
                if(nx >= 1 && nx < gridSize-1 && ny >= 1 && ny < gridSize-1) {
                    fluidCube.addDensity(nx, ny, gridZ, 15.0f); // Even higher density value for visibility
                }
            }
        }
    }
}


void Simulation::handleMouseMove(int x, int y, int width, int height)
{
    if (x == m_lastMouseX && y == m_lastMouseY) {
        return;
    }

    // Calculate velocity based on mouse movement - REVERSE BOTH directions
    float velocityX = static_cast<float>(m_lastMouseX - x) * 0.5f; // REVERSED direction
    float velocityY = static_cast<float>(m_lastMouseY - y) * 0.5f; // REVERSED direction

    // Convert screen coordinates to normalized coordinates
    float normalizedX = 1.0f - static_cast<float>(x) / width;  // Flip X to fix mirroring
    float normalizedY = 1.0f - static_cast<float>(y) / height;  // Also flip Y to match OpenGL coordinates

    // Map to fluid grid coordinates
    int gridSize = 32; // Same as in init()
    int gridX = static_cast<int>(normalizedX * gridSize);
    int gridY = static_cast<int>(normalizedY * gridSize);
    int gridZ = gridSize / 2; // Middle of the grid

    // Add velocity and density for visualization with enhanced effect
    if (gridX >= 1 && gridX < gridSize-1 &&
        gridY >= 1 && gridY < gridSize-1 &&
        gridZ >= 1 && gridZ < gridSize-1) {
        // Add velocity in a small area with STRONGER effect
        for(int dx = -1; dx <= 1; dx++) {
            for(int dy = -1; dy <= 1; dy++) {
                int nx = gridX + dx;
                int ny = gridY + dy;
                if(nx >= 1 && nx < gridSize-1 && ny >= 1 && ny < gridSize-1) {
                    // Use much larger velocity values for more dramatic effect
                    fluidCube.addVelocity(nx, ny, gridZ, velocityX * 3.0f, velocityY * 3.0f, 0.0f);
                    fluidCube.addDensity(nx, ny, gridZ, 10.0f); // Higher density for better visibility
                }
            }
        }
    }

    // Update last position
    m_lastMouseX = x;
    m_lastMouseY = y;
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
