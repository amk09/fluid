#include "simulation.h"
#include "graphics/meshloader.h"

#include <iostream>

using namespace Eigen;

Simulation::Simulation() {}

void Simulation::init()
{
    // Initialize the fluid with much lower diffusion and viscosity parameters
    // This will make the fluid less likely to spread out horizontally
    fluidCube.init(32, 0.00001, 0.00001);
}


void Simulation::handleMousePress(int x, int y, int width, int height)
{
    m_lastMouseX = x;
    m_lastMouseY = y;

    // Convert screen coordinates to normalized coordinates
    float normalizedX = 1.0f - static_cast<float>(x) / width;
    float normalizedY = 1.0f - static_cast<float>(y) / height;

    // Map to fluid grid coordinates
    int gridSize = 32; // Match the size used in init()
    float gridX = normalizedX * gridSize;
    float gridY = normalizedY * gridSize;
    int gridZ = gridSize / 2; // Middle of the grid

    // Add density with Gaussian distribution for smoother effect
    addDensityWithGaussian(gridX, gridY, gridZ, 8.0f, 1.5f);
}


void Simulation::addDensityWithGaussian(float centerX, float centerY, float centerZ, float amount, float sigma)
{
    int gridSize = 32;
    int radius = static_cast<int>(sigma * 2.0f);

    for (int dz = -radius; dz <= radius; dz++) {
        for (int dy = -radius; dy <= radius; dy++) {
            for (int dx = -radius; dx <= radius; dx++) {
                int nx = static_cast<int>(centerX + dx);
                int ny = static_cast<int>(centerY + dy);
                int nz = static_cast<int>(centerZ + dz);

                if (nx >= 1 && nx < gridSize-1 &&
                    ny >= 1 && ny < gridSize-1 &&
                    nz >= 1 && nz < gridSize-1) {

                    // Calculate Gaussian weight
                    float distSq = dx*dx + dy*dy + dz*dz;
                    float weight = exp(-distSq / (1.5f * sigma * sigma));

                    // Add weighted density
                    fluidCube.addDensity(nx, ny, nz, amount * weight* 0.7f);
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
    float normalizedX = 1.0f - static_cast<float>(x) / width;
    float normalizedY = 1.0f - static_cast<float>(y) / height;

    // Map to fluid grid coordinates
    int gridSize = 32; // Same as in init()
    float gridX = normalizedX * gridSize;
    float gridY = normalizedY * gridSize;
    int gridZ = gridSize / 2; // Middle of the grid

    // Add velocity and density with Gaussian distribution for smoother effect
    addVelocityWithGaussian(gridX, gridY, gridZ, velocityX * 3.5f, velocityY * 3.5f, 0.0f, 1.3f);
    addDensityWithGaussian(gridX, gridY, gridZ, 6.0f, 1.3f);

    // Update last position
    m_lastMouseX = x;
    m_lastMouseY = y;
}


void Simulation::addVelocityWithGaussian(float centerX, float centerY, float centerZ,
                                         float velX, float velY, float velZ, float sigma)
{
    int gridSize = 32;
    int radius = static_cast<int>(sigma * 1.8f);

    for (int dz = -radius; dz <= radius; dz++) {
        for (int dy = -radius; dy <= radius; dy++) {
            for (int dx = -radius; dx <= radius; dx++) {
                int nx = static_cast<int>(centerX + dx);
                int ny = static_cast<int>(centerY + dy);
                int nz = static_cast<int>(centerZ + dz);

                if (nx >= 1 && nx < gridSize-1 &&
                    ny >= 1 && ny < gridSize-1 &&
                    nz >= 1 && nz < gridSize-1) {

                    // Calculate Gaussian weight
                    float distSq = dx*dx + dy*dy + dz*dz;
                    float weight = exp(-distSq / (1.2f * sigma * sigma));

                    // Add weighted velocity
                    fluidCube.addVelocity(nx, ny, nz,
                                          velX * weight * 0.8f,
                                          velY * weight * 0.8f,
                                          velZ * weight * 0.8f);
                }
            }
        }
    }
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
