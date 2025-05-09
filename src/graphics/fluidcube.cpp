#include "fluidcube.h"
#include "graphics/shader.h"
#include "graphics/solver.h"
#include "graphics/obstacle.h"
#include "graphics/vectorfield.h"
#include <iomanip>
#include <iostream>
#include <filesystem>  // C++17
#include <fstream>

namespace fs = std::filesystem;

#define Log(x) std::cout << x << std::endl;

// For offline rendering
static bool offlineRendering = false; // Just change this
static bool offlineFileLoaded = false;
static bool startRenderingOnce = false;
static bool offlineRenderingFF = true; // Just change this
static bool offlineFileLoadedFF = false;
static bool startRenderingOnceFF = false;
const int totalFrames = 500;
static int frameIndex = 0;
std::string densityPath = "offline_data/densityFrames.bin";
std::string colorPath   = "offline_data/colorFrames.bin";
std::string velocityPath = "offline_data/velocityFrames.bin";

FluidCube::FluidCube()
    : m_vao(0),
    m_vbo(0),
    m_ibo(0),
    m_densityTexture(0),
    m_shellDensityTexture(0),
    m_voxelVao(0),
    m_voxelVbo(0),
    m_voxelIbo(0),
    m_fullscreen_vao(0),
    m_fullscreen_vbo(0),
    m_modelMatrix(Eigen::Matrix4f::Identity()),
    m_wireframe(false),
    size(0),
    diff(0.0f),
    visc(0.0f)
{
}

void FluidCube::fountainGeneration() {
    int center = size / 2;
    float radius = size / 6.0f; // Diameter = size / 3

    int range = static_cast<int>(std::ceil(radius)); // Ensure full coverage

    int colorType = m_colorMapType;

    for (int dx = -range; dx <= range; dx++) {
        for (int dz = -range; dz <= range; dz++) {
            float dist2 = dx * dx + dz * dz;
            if (dist2 > radius * radius) continue;

            int x = center + dx;
            int z = center + dz;

            for (int y = 1; y <= 2; y++) {
                float densityAmount = 0.7f;
                float velocityAmount = 10.0f + (rand() / (float)RAND_MAX) * 2.0f;

                addDensityWithColor(x, y, z, densityAmount, colorType);
                addVelocity(x, y, z, 0.0f, velocityAmount, 0.0f);
            }
        }
    }
}

void FluidCube::fountainGenerationTopDown() {
    int center = size / 2;
    float radius = size / 6.0f;

    int range = static_cast<int>(std::ceil(radius));
    int colorType = m_colorMapType;

    for (int dx = -range; dx <= range; dx++) {
        for (int dz = -range; dz <= range; dz++) {
            float dist2 = dx * dx + dz * dz;
            if (dist2 > radius * radius) continue;

            int x = center + dx;
            int z = center + dz;

            int y = size - 3; // Topmost usable layer

            float densityAmount = 0.7f;
            float velocityAmount = - 10.0f - (rand() / (float)RAND_MAX) * 2.0f; // spray downward

            addDensityWithColor(x, y, z, densityAmount, colorType);
            addVelocity(x, y, z, 0.0f, velocityAmount, 0.0f);
        }
    }
}


void FluidCube::draw(Shader *shader)
{
    // // For Voxel
    // drawVoxel(shader);

    // For Ray Marchign
    drawVolume(shader);
}

void FluidCube::update(float dt)
{
    if ((offlineRenderingFF && offlineFileLoadedFF) || startRenderingOnceFF) {
        // Load Density and Color for next stage (just move on to just one place for each update call, start from zero. and stop in the end.)
        renderNextOfflineFrame();
        return;
    }

    if ((offlineRendering && offlineFileLoaded) || startRenderingOnce) {
        // Load Density and Color for next stage (just move on to just one place for each update call, start from zero. and stop in the end.)
        if (frameIndex < densityFrames.size()) {
            density = densityFrames[frameIndex];
            fluidColors = colorFrames[frameIndex];
            velocity = velocityFrames[frameIndex];
            uploadDensityToGPU();
            uploadColorFieldToGPU();

            if (m_velocityTexture == 0)
                glGenTextures(1, &m_velocityTexture);

            setupTextureParameters(m_velocityTexture);
            glTexImage3D(
                GL_TEXTURE_3D, 0, GL_R32F,
                size, size, size, 0, GL_RED, GL_FLOAT,
                velocity.data()
                );

            if (m_renderMode == 1) {
                detectShell();
            }

            frameIndex++;
            // Log(frameIndex);
        }
        return;  // Skip simulation logic
    }


#pragma omp parallel sections
    {
#pragma omp section
        {
            addSource(vX, vX0);
        }

#pragma omp section
        {
            addSource(vY, vY0);
        }

#pragma omp section
        {
            addSource(vZ, vZ0);
        }
    }

    // 1. Add vorticity confinement
    fountainGeneration();
    // fountainGenerationTopDown();
    addVorticityConfinement(vX0, vY0, vZ0, dt, size, m_vorticityStrength);

    // 2. Diffuse velocity
#pragma omp parallel sections
    {
#pragma omp section
        {
            diffuse_velocity(1, vX0, vX, visc, dt, iter, size);
        }

#pragma omp section
        {
            diffuse_velocity(2, vY0, vY, visc, dt, iter, size);
        }

#pragma omp section
        {
            diffuse_velocity(3, vZ0, vZ, visc, dt, iter, size);
        }
    }

    // 3. Project velocity
    project(vX0, vY0, vZ0, vX, vY, iter, size);

    // 4. Advect velocity
#pragma omp parallel sections
    {
#pragma omp section
        {
            advect(1, vX, vX0, vX0, vY0, vZ0, dt, size);
        }

#pragma omp section
        {
            advect(2, vY, vY0, vX0, vY0, vZ0, dt, size);
        }

#pragma omp section
        {
            advect(3, vZ, vZ0, vX0, vY0, vZ0, dt, size);
        }
    }

    // 5. Project again
    project(vX, vY, vZ, vX0, vY0, iter, size);

    // Clean the velocity value in v0
    empty_vel();

    addSource(density, density0);

    // 6. Diffuse density
    diffuse_density(0, density0, density, diff, dt, iter, size);

    // 7. Advect density
    advect(0, density, density0, vX, vY, vZ, dt, size);

    // Fade the density value in density: To avoid the inaccurate calculation's increasing density
    densityFade(dt);

    // Make Sure Obstacle's Density is 1.0
    densityInsideObstacles(density);

    // Clean the density value in density0
    empty_den();

    // Handle boundary cells - set density and color to 0
    for (int k = 0; k < size; k++) {
        for (int j = 0; j < size; j++) {
            // X boundaries
            density[index(0, j, k)] = 0.0f;
            density[index(size-1, j, k)] = 0.0f;
            fluidColors[index(0, j, k)] = 0;
            fluidColors[index(size-1, j, k)] = 0;

            // Y boundaries
            density[index(j, 0, k)] = 0.0f;
            density[index(j, size-1, k)] = 0.0f;
            fluidColors[index(j, 0, k)] = 0;
            fluidColors[index(j, size-1, k)] = 0;

            // Z boundaries
            density[index(j, k, 0)] = 0.0f;
            density[index(j, k, size-1)] = 0.0f;
            fluidColors[index(j, k, 0)] = 0;
            fluidColors[index(j, k, size-1)] = 0;
        }
    }

    // First mark all obstacle cells with special color
    setObstacleColors(fluidColors);

    // Perform extra color enhancement to fix possible inconsistencies
    for (int i = 0; i < totalCells; i++) {
        // Skip obstacle and boundary cells
        if (isObstacle(i)) continue;

        // Skip boundary cells - extra check
        int x = i % size;
        int y = (i / size) % size;
        int z = i / (size * size);
        if (x == 0 || x == size-1 || y == 0 || y == size-1 || z == 0 || z == size-1) continue;

        // For cells with significant density but no color, force default color
        if (density[i] > 0.03f && fluidColors[i] == 0) {
            fluidColors[i] = m_colorMapType;
        }
    }

    // Propagate colors to ensure they follow density properly
    propagateColors();

    // Re-mark all obstacle cells with special color to ensure they weren't overwritten
    setObstacleColors(fluidColors);

    if (m_renderMode == 1) {
        detectShell();
    }

    // Upload both density and color to GPU
    uploadDensityToGPU();
    uploadVelocityToGPU();
    uploadColorFieldToGPU();
}

// TODO::Change the logic here to reduce the density fade speed
void FluidCube::densityFade(float dt)
{
    for (int i = 0; i < totalCells; i++)
    {
        density[i] *= (1 - 1.5*dt);
    }
}





void FluidCube::detectShell() {
    // Create a temporary array to store the shell data
    std::vector<float> shellDensity(size * size * size, 0.0f);

    // First, identify shell cells - cells that are fluid but adjacent to non-fluid
    for (int k = 1; k < size - 1; k++) {
        for (int j = 1; j < size - 1; j++) {
            for (int i = 1; i < size - 1; i++) {
                int idx = index(i, j, k);

                // Skip processing cells with negligible density
                if (density[idx] < 0.01f) continue;

                // Check if this is a shell cell (has at least one neighbor with significantly lower density)
                bool isShell = false;

                // Check 6-neighborhood (faces)
                if (density[index(i+1, j, k)] < 0.01f ||
                    density[index(i-1, j, k)] < 0.01f ||
                    density[index(i, j+1, k)] < 0.01f ||
                    density[index(i, j-1, k)] < 0.01f ||
                    density[index(i, j, k+1)] < 0.01f ||
                    density[index(i, j, k-1)] < 0.01f) {
                    isShell = true;
                }

                if (isShell) {
                    // This is a shell cell - boost its density for rendering
                    shellDensity[idx] = density[idx] * 2.0f; // Boost opacity
                }
            }
        }
    }

    // Store the shell density in a separate texture
    // If you don't have m_shellDensityTexture defined yet, you'll need to add it
    if (m_shellDensityTexture == 0) {
        glGenTextures(1, &m_shellDensityTexture);
    }

    glBindTexture(GL_TEXTURE_3D, m_shellDensityTexture);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_R32F, size, size, size, 0, GL_RED, GL_FLOAT, shellDensity.data());
    glBindTexture(GL_TEXTURE_3D, 0);
}






// Below are the functions that do not need to change











// Add initialization of color field in init method
void FluidCube::init(int size, float diffuse, float viscosity){
    this->size = size;
    this->diff = diffuse;
    this->visc= viscosity;
    this->totalCells = size*size*size;

    offlineRendering = size > 32;

    // We resize the density and velocity vectors with 0.0f value
    density0.resize(totalCells, 0.0f);
    density.resize(totalCells, 0.0f);

    vX0.resize(totalCells, 0.0f);
    vX.resize(totalCells, 0.0f);

    vY0.resize(totalCells, 0.0f);
    vY.resize(totalCells, 0.0f);

    vZ0.resize(totalCells, 0.0f);
    vZ.resize(totalCells, 0.0f);

    // Initialize color field with zeros
    fluidColors.resize(totalCells, 0);
    m_colorTexture = 0;

    initFullscreenQuad();

    // Obstacle Init
    initObstacleMap(totalCells);

    // Inject Density and then upload it to GPU
    uploadDensityToGPU();
    // Initialize color texture
    uploadColorFieldToGPU();

    if(offlineRendering){
        // Offline Calculation
        offRenderingCheck();
    }

    if (offlineRenderingFF) {
        offRenderingCheckFbyF();
    }
}

void FluidCube::addSource(vector<float>& x, vector<float>& x0)
{
    for (int i = 0; i < totalCells; i++)
    {
        x[i] += dt * x0[i];
    }
}

void FluidCube::empty_vel(){
    for(int i = 0; i < vX0.size(); i++){
        vX0[i] = 0.0f;
        vY0[i] = 0.0f;
        vZ0[i] = 0.0f;
    }
}

void FluidCube::empty_den(){
    for(int i = 0; i < density0.size(); i++){
        density0[i] = 0.0f;
    }
}


// Modified addDensity to use addDensityWithColor
void FluidCube::addDensity(int x, int y, int z, float amount) {
    // Call addDensityWithColor with current global color type
    addDensityWithColor(x, y, z, amount, m_colorMapType);
}


void FluidCube::addDensityWithColor(int x, int y, int z, float amount, int colorType) {
    // Clamp coordinates to valid range
    x = std::max(1, std::min(x, size-2));
    y = std::max(1, std::min(y, size-2));
    z = std::max(1, std::min(z, size-2));

    int idx = index(x, y, z);

    // Add to center cell
    float oldDensity = density[idx];
    addDensityToCell(x, y, z, amount, colorType, oldDensity);

    // Add small amounts to neighboring cells for smoothness
    float neighborAmount = amount * 0.2f;

    // Add to 6 direct neighbors
    for (int i = -1; i <= 1; i += 2) {
        // X neighbors
        if (isValidCell(x+i, y, z)) {
            oldDensity = density[index(x+i, y, z)];
            addDensityToCell(x+i, y, z, neighborAmount, colorType, oldDensity);
        }

        // Y neighbors
        if (isValidCell(x, y+i, z)) {
            oldDensity = density[index(x, y+i, z)];
            addDensityToCell(x, y+i, z, neighborAmount, colorType, oldDensity);
        }

        // Z neighbors
        if (isValidCell(x, y, z+i)) {
            oldDensity = density[index(x, y, z+i)];
            addDensityToCell(x, y, z+i, neighborAmount, colorType, oldDensity);
        }
    }
}


// Function to upload color field to GPU
void FluidCube::uploadColorFieldToGPU() {
    if (m_colorTexture == 0)
        glGenTextures(1, &m_colorTexture);

    setupTextureParameters(m_colorTexture);

    // Convert int vector to float for OpenGL
    std::vector<float> colorData(totalCells);
    for (int i = 0; i < totalCells; i++) {
        colorData[i] = static_cast<float>(fluidColors[i]);
    }

    glTexImage3D(
        GL_TEXTURE_3D,
        0,
        GL_R32F,
        size, size, size,
        0,
        GL_RED,
        GL_FLOAT,
        colorData.data()
        );

    glBindTexture(GL_TEXTURE_3D, 0);
}

void FluidCube::uploadVelocityToGPU() {
    if (m_velocityTexture == 0) {
        glGenTextures(1, &m_velocityTexture);
    }

    glBindTexture(GL_TEXTURE_3D, m_velocityTexture);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    // Calculate velocity magnitude
    std::vector<float> velocityMagnitude(size * size * size);

    for (int k = 0; k < size; k++) {
        for (int j = 0; j < size; j++) {
            for (int i = 0; i < size; i++) {
                int idx = index(i, j, k);
                float vx = vX[idx];
                float vy = vY[idx];
                float vz = vZ[idx];

                // Calculate magnitude
                velocityMagnitude[idx] = sqrt(vx*vx + vy*vy + vz*vz);
            }
        }
    }

    glTexImage3D(GL_TEXTURE_3D, 0, GL_R32F, size, size, size, 0,
                 GL_RED, GL_FLOAT, velocityMagnitude.data());
    glBindTexture(GL_TEXTURE_3D, 0);
}

// Upload custom color field for shell rendering
void FluidCube::uploadCustomColorToGPU(const std::vector<float>& customColors) {
    if (m_colorTexture == 0)
        glGenTextures(1, &m_colorTexture);

    glBindTexture(GL_TEXTURE_3D, m_colorTexture);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    glTexImage3D(
        GL_TEXTURE_3D,
        0,
        GL_R32F,
        size, size, size,
        0,
        GL_RED,
        GL_FLOAT,
        customColors.data()
        );

    glBindTexture(GL_TEXTURE_3D, 0);
}

void FluidCube::drawVolume(Shader* shader) {
    // For shell-only rendering (as an option)
    if (m_renderMode == 1) {
        drawShellOnly(shader);
    } else {
        // Make sure density is uploaded normally
        uploadDensityToGPU();
        uploadColorFieldToGPU();
    }

    // Pass both textures to the shader
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, m_densityTexture);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_3D, m_colorTexture);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_3D, m_velocityTexture);

    shader->setUniform("densityTex", 0);
    shader->setUniform("colorTex", 1);
    shader->setUniform("velocityTex", 2);

    shader->setUniform("useVelocityColor", m_useVelocityColor);
    shader->setUniform("velocityScale", m_velocityScale);
    shader->setUniform("velocityBlend", m_velocityBlend);

    shader->setUniform("colorMapType", m_colorMapType);
    shader->setUniform("renderMode", m_renderMode);
    shader->setUniform("size", size);

    // Pass current time for animation effects
    static float timeValue = 0.0f;
    timeValue += 0.016f; // Increment by roughly 1/60 second
    shader->setUniform("time", timeValue);

    glBindVertexArray(m_fullscreen_vao);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glBindVertexArray(0);
}


void FluidCube::drawShellOnly(Shader *shader) {
    // Create temporary arrays for shell information
    if (m_shellDensityTexture == 0) {
        // Generate the shell texture if not already done
        detectShell();
    }

    // Pass the shell texture to the shader
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, m_shellDensityTexture);

    // Also pass color texture
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_3D, m_colorTexture);

    shader->setUniform("densityTex", 0);
    shader->setUniform("colorTex", 1);
    shader->setUniform("colorMapType", m_colorMapType);
    shader->setUniform("renderMode", m_renderMode);
    shader->setUniform("size", size);

    // Pass current time for animation effects
    static float timeValue = 0.0f;
    timeValue += 0.016f; // Increment by roughly 1/60 second
    shader->setUniform("time", timeValue);

    // Draw the fullscreen quad
    glBindVertexArray(m_fullscreen_vao);
    glDrawArrays(GL_TRIANGLES, 0, 6);
}




// Add a method to clear all fluids
void FluidCube::clearAllFluids() {
    // Reset all densities and velocities to zero
    std::fill(density.begin(), density.end(), 0.0f);
    std::fill(density0.begin(), density0.end(), 0.0f);
    std::fill(vX.begin(), vX.end(), 0.0f);
    std::fill(vX0.begin(), vX0.end(), 0.0f);
    std::fill(vY.begin(), vY.end(), 0.0f);
    std::fill(vY0.begin(), vY0.end(), 0.0f);
    std::fill(vZ.begin(), vZ.end(), 0.0f);
    std::fill(vZ0.begin(), vZ0.end(), 0.0f);
    std::fill(fluidColors.begin(), fluidColors.end(), 0);

    // Update GPU textures
    uploadDensityToGPU();
    uploadVelocityToGPU();
    uploadColorFieldToGPU();
}

void FluidCube::propagateColors() {
    // For each cell, if it has significant density but no color assigned yet,
    // get color from neighboring cells with highest density
    for (int k = 1; k < size-1; k++) {
        for (int j = 1; j < size-1; j++) {
            for (int i = 1; i < size-1; i++) {
                int idx = index(i, j, k);

                // Skip obstacle cells
                if (isObstacle(idx)) continue;

                // Skip boundary cells (should be handled already, this is extra protection)
                if (i == 0 || i == size-1 || j == 0 || j == size-1 || k == 0 || k == size-1) continue;

                // Only process cells with density but no color - even lower threshold
                if (density[idx] > 0.005f && fluidColors[idx] == 0) {
                    findBetterColorFromNeighbors(idx, i, j, k);
                }
            }
        }
    }
}



void FluidCube::addVelocity(int x, int y, int z, float amountX, float amountY, float amountZ){
    // Do we have some limitations on Velocity Amount?
    vX[index(x, y, z)] += amountX;
    vY[index(x, y, z)] += amountY;
    vZ[index(x, y, z)] += amountZ;
}

// Below is one of the methods to visualize the density
void FluidCube::initFullscreenQuad() {
    std::vector<GLfloat> fullscreen_quad_data = {
        // POSITIONS        // UV COORDINATES
        -1.0f,  1.0f, 0.0f,   0.0f, 1.0f,  // Upper Left
        -1.0f, -1.0f, 0.0f,   0.0f, 0.0f,  // Bottom Left
        1.0f, -1.0f, 0.0f,   1.0f, 0.0f,  // Bottom Right

        -1.0f,  1.0f, 0.0f,   0.0f, 1.0f,  // Upper Left
        1.0f, -1.0f, 0.0f,   1.0f, 0.0f,  // Bottom Right
        1.0f,  1.0f, 0.0f,   1.0f, 1.0f   // Upper Right
    };

    glGenBuffers(1, &m_fullscreen_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, m_fullscreen_vbo);
    glBufferData(GL_ARRAY_BUFFER, fullscreen_quad_data.size() * sizeof(GLfloat),
                 fullscreen_quad_data.data(), GL_STATIC_DRAW);

    glGenVertexArrays(1, &m_fullscreen_vao);
    glBindVertexArray(m_fullscreen_vao);

    // layout(location = 0): vec3 aPos
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE,
                          5 * sizeof(GLfloat), reinterpret_cast<void*>(0));

    // layout(location = 1): vec2 aUV
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE,
                          5 * sizeof(GLfloat), reinterpret_cast<void*>(3 * sizeof(GLfloat)));

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void FluidCube::uploadDensityToGPU() {
    if (m_densityTexture == 0)
        glGenTextures(1, &m_densityTexture);

    // Set up texture parameters with our helper function
    glBindTexture(GL_TEXTURE_3D, m_densityTexture);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    glTexImage3D(
        GL_TEXTURE_3D,
        0,
        GL_R32F,
        size, size, size,
        0,
        GL_RED,
        GL_FLOAT,
        density.data()
        );

    glBindTexture(GL_TEXTURE_3D, 0);
}


// Below is one of the methods to visualize the density
void FluidCube::initVoxelGeometry(){
    float cubeVertices[] = {
        // 8 corner vertices
        -0.5f, -0.5f, -0.5f, // 0
        0.5f, -0.5f, -0.5f, // 1
        0.5f,  0.5f, -0.5f, // 2
        -0.5f,  0.5f, -0.5f, // 3
        -0.5f, -0.5f,  0.5f, // 4
        0.5f, -0.5f,  0.5f, // 5
        0.5f,  0.5f,  0.5f, // 6
        -0.5f,  0.5f,  0.5f  // 7
    };

    unsigned int cubeIndices[] = {
        // Front face (z = 0.5)
        7, 6, 4, 6, 5, 4,
        // Back face (z = -0.5)
        0, 1, 3, 1, 2, 3,
        // Left face (x = -0.5)
        3, 7, 0, 7, 4, 0,
        // Right face (x = 0.5)
        6, 2, 5, 2, 1, 5,
        // Top face (y = 0.5)
        3, 2, 7, 2, 6, 7,
        // Bottom face (y = -0.5)
        4, 5, 0, 5, 1, 0
    };

    glGenVertexArrays(1, &m_voxelVao);
    glGenBuffers(1, &m_voxelVbo);
    glGenBuffers(1, &m_voxelIbo);

    glBindVertexArray(m_voxelVao);
    glBindBuffer(GL_ARRAY_BUFFER, m_voxelVbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cubeVertices), cubeVertices, GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_voxelIbo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(cubeIndices), cubeIndices, GL_STATIC_DRAW);

    glBindVertexArray(0);
}

void FluidCube::drawVoxel(Shader* shader){
    glBindVertexArray(m_voxelVao);

    float voxelSize = 1.0f / size;

#pragma omp parallel for
    for (int z = 0; z < size; ++z)
        for (int y = 0; y < size; ++y)
            for (int x = 0; x < size; ++x) {
                float d = density[index(x, y, z)];
                if (d < 0.01f) continue;

                float px = (x + 0.5f) / size - 0.5f;
                float py = (y + 0.5f) / size - 0.5f;
                float pz = (z + 0.5f) / size - 0.5f;

                Eigen::Affine3f model =
                    Eigen::Translation3f(px, py, pz) *
                    Eigen::Scaling(voxelSize);

                shader->setUniform("model", model.matrix());
                shader->setUniform("gray", d);
                glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
            }

    glBindVertexArray(0);
}

void FluidCube::setModelMatrix(const Eigen::Affine3f &model)
{
    m_modelMatrix = model.matrix();
}

void FluidCube::toggleWireframe()
{
    m_wireframe = !m_wireframe;
}

void FluidCube::toggleShellRendering() {
    m_renderMode = 1 - m_renderMode; // Toggle between 0 (volume) and 1 (shell)
    // Force shell update if shell mode
    if (m_renderMode == 1) {
        detectShell();
    }
}

void FluidCube::visualizeVelocity() {
    float maxVelocity = 0.0f;

    // Find maximum velocity magnitude
    for (int i = 0; i < vX.size(); i++) {
        float magnitude = std::sqrt(vX[i]*vX[i] + vY[i]*vY[i] + vZ[i]*vZ[i]);
        maxVelocity = std::max(maxVelocity, magnitude);
    }

    std::cout << "Maximum velocity: " << maxVelocity << std::endl;

    // Convert velocity to density for visualization
    if (maxVelocity > 0.001f) {
        for (int i = 0; i < density.size(); i++) {
            float magnitude = std::sqrt(vX[i]*vX[i] + vY[i]*vY[i] + vZ[i]*vZ[i]);
            // Add a portion of velocity magnitude to density
            density[i] = std::min(1.0f, density[i] + magnitude / maxVelocity * 0.2f);
        }
    }
}




void FluidCube::setColorMap(int colorType) {
    m_colorMapType = colorType;
    uploadColorFieldToGPU();
}


void FluidCube::uploadCustomDensityToGPU(const std::vector<float>& customDensity) {
    if (m_densityTexture == 0)
        glGenTextures(1, &m_densityTexture);

    glBindTexture(GL_TEXTURE_3D, m_densityTexture);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    glTexImage3D(
        GL_TEXTURE_3D,
        0,
        GL_R32F,
        size, size, size,
        0,
        GL_RED,
        GL_FLOAT,
        customDensity.data()
        );

    glBindTexture(GL_TEXTURE_3D, 0);
}




// New helper method to set texture parameters
void FluidCube::setupTextureParameters(GLuint texture) {
    glBindTexture(GL_TEXTURE_3D, texture);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
}

// Helper method for adding density to a cell with color consideration
void FluidCube::addDensityToCell(int x, int y, int z, float amount, int colorType, float oldDensity) {
    int idx = index(x, y, z);

    // Add density, clamping to valid range
    density[idx] = std::clamp(density[idx] + amount, 0.0f, 1.0f);

    // Apply color if significant change or low previous density
    if (oldDensity < 0.1f || amount > 0.05f) {
        fluidColors[idx] = colorType;
    }
}

// Method to find a better color from neighbors during propagation
bool FluidCube::findBetterColorFromNeighbors(int idx, int i, int j, int k) {
    float maxDensity = 0.0f;
    int bestColor = 0;
    bool foundBetterColor = false;

    // Check direct neighbors first (more accurate color propagation)
    const int directNeighbors[6][3] = {
        {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}
    };

    for (int n = 0; n < 6; n++) {
        int nx = i + directNeighbors[n][0];
        int ny = j + directNeighbors[n][1];
        int nz = k + directNeighbors[n][2];

        if (isValidCell(nx, ny, nz)) {
            int nidx = index(nx, ny, nz);

            // Skip obstacle neighbors - don't propagate obstacle color
            if (isObstacle(nidx)) continue;

            // If direct neighbor has color and significant density
            if (density[nidx] > 0.05f && fluidColors[nidx] > 0 && fluidColors[nidx] != 999) {
                fluidColors[idx] = fluidColors[nidx];
                foundBetterColor = true;
                return true;
            }
        }
    }

    // If no direct neighbor with color, check all neighbors
    if (!foundBetterColor) {
        // Check all 26 neighbors (including diagonals)
        for (int nz = -1; nz <= 1; nz++) {
            for (int ny = -1; ny <= 1; ny++) {
                for (int nx = -1; nx <= 1; nx++) {
                    // Skip center cell
                    if (nx == 0 && ny == 0 && nz == 0) continue;

                    int ni = i + nx;
                    int nj = j + ny;
                    int nk = k + nz;

                    // Check if neighbor is valid
                    if (isValidCell(ni, nj, nk)) {
                        int nidx = index(ni, nj, nk);

                        // Skip obstacle neighbors
                        if (isObstacle(nidx)) continue;

                        // If neighbor has higher density and a valid color (not obstacle)
                        if (density[nidx] > maxDensity &&
                            fluidColors[nidx] > 0 &&
                            fluidColors[nidx] != 999) {
                            maxDensity = density[nidx];
                            bestColor = fluidColors[nidx];
                            foundBetterColor = true;
                        }
                    }
                }
            }
        }

        // Apply the best color found
        if (bestColor > 0) {
            fluidColors[idx] = bestColor;
            return true;
        }
    }

    return foundBetterColor;
}

// Implement the clearObstacles method (that was previously added)
void FluidCube::clearObstacles() {
    // Clear obstacle markers
    for (int idx : g_obstacleCells) {
        if (idx >= 0 && idx < density.size()) {
            density[idx] = 0.0f;
            fluidColors[idx] = 0; // Also clear colors
        }
    }

    std::fill(g_obstacle.begin(), g_obstacle.end(), 0);
    g_obstacleCells.clear();

    // Update GPU textures immediately
    uploadDensityToGPU();
    uploadColorFieldToGPU();
}


void FluidCube::offRenderingCheck() {
    namespace fs = std::filesystem;

    std::string densityPath = "offline_data/densityFrames.bin";
    std::string colorPath   = "offline_data/colorFrames.bin";
    std::string velocityPath = "offline_data/velocityFrames.bin";

    fs::create_directory("offline_data");

    if (fs::exists(densityPath) && fs::exists(colorPath)) {
        offlineFileLoaded = true;
        // Load from file
        std::ifstream dFile(densityPath, std::ios::binary);
        std::ifstream cFile(colorPath, std::ios::binary);
        std::ifstream vFile(velocityPath, std::ios::binary);

        int numFrames = 0;
        dFile.read(reinterpret_cast<char*>(&numFrames), sizeof(int));

        for (int f = 0; f < numFrames; ++f) {
            std::vector<float> densityFrame(totalCells);
            std::vector<int> colorFrame(totalCells);
            std::vector<float> velocityFrame(totalCells);
            dFile.read(reinterpret_cast<char*>(densityFrame.data()), totalCells * sizeof(float));
            cFile.read(reinterpret_cast<char*>(colorFrame.data()), totalCells * sizeof(int));
            vFile.read(reinterpret_cast<char*>(velocityFrame.data()), totalCells * sizeof(float));
            densityFrames.push_back(std::move(densityFrame));
            colorFrames.push_back(std::move(colorFrame));
            velocityFrames.push_back(std::move(velocityFrame));
        }

        Log("Loaded offline frames from disk.");
    } else {
        offlineFileLoaded = false;
        // Run simulation and save
        for (int i = 0; i < totalFrames; ++i) {
            update(0.00833f); // 120 FPS
            Log(i);
            densityFrames.push_back(density);
            colorFrames.push_back(fluidColors);
            std::vector<float> velocityMagnitudes(totalCells);
            for (int j = 0; j < totalCells; j++) {
                float vx = vX[j];
                float vy = vY[j];
                float vz = vZ[j];
                velocityMagnitudes[j] = sqrt(vx*vx + vy*vy + vz*vz);
            }
            velocityFrames.push_back(std::move(velocityMagnitudes));

        }

        std::ofstream dFile(densityPath, std::ios::binary);
        std::ofstream cFile(colorPath, std::ios::binary);
        std::ofstream vFile(velocityPath, std::ios::binary);

        int numFrames = densityFrames.size();
        dFile.write(reinterpret_cast<const char*>(&numFrames), sizeof(int));

        for (int f = 0; f < numFrames; ++f) {
            dFile.write(reinterpret_cast<const char*>(densityFrames[f].data()), totalCells * sizeof(float));
            cFile.write(reinterpret_cast<const char*>(colorFrames[f].data()), totalCells * sizeof(int));
            vFile.write(reinterpret_cast<const char*>(velocityFrames[f].data()), totalCells * sizeof(float));

        }

        Log("Saved offline frames to disk.");
        startRenderingOnce = true;
    }
}

void FluidCube::renderNextOfflineFrame()
{
    frameIndex = frameIndex % totalFrames;
    if (!offlineRenderingFF) return;

    auto makePath = [](const std::string& prefix, int idx) {
        std::ostringstream oss;
        oss << "offline_data/" << prefix << "_" << std::setw(4) << std::setfill('0') << idx << ".bin";
        return oss.str();
    };

    const std::string dPath = makePath("density",  frameIndex);
    const std::string cPath = makePath("color",    frameIndex);
    const std::string vPath = makePath("velocity", frameIndex);

    if (!fs::exists(dPath)) return;


    std::vector<float> scratchDensity(totalCells);
    std::vector<int>   scratchColor  (totalCells);
    std::vector<float> scratchVelMag (totalCells);

    {
        std::ifstream f(dPath, std::ios::binary);
        f.read(reinterpret_cast<char*>(scratchDensity.data()),
               totalCells * sizeof(float));
    }
    {
        std::ifstream f(cPath, std::ios::binary);
        f.read(reinterpret_cast<char*>(scratchColor.data()),
               totalCells * sizeof(int));
    }
    {
        std::ifstream f(vPath, std::ios::binary);
        f.read(reinterpret_cast<char*>(scratchVelMag.data()),
               totalCells * sizeof(float));
    }

    // ------------ upload ------------
    density      = scratchDensity;
    fluidColors  = scratchColor;
    velocity     = scratchVelMag;
    uploadDensityToGPU();
    uploadColorFieldToGPU();

    if (m_velocityTexture == 0)
        glGenTextures(1, &m_velocityTexture);

    setupTextureParameters(m_velocityTexture);
    glTexImage3D(
        GL_TEXTURE_3D, 0, GL_R32F,
        size, size, size, 0, GL_RED, GL_FLOAT,
        velocity.data()
        );

    if (m_renderMode == 1) {
        detectShell();
    }

    cout << "Rendered Fram: " << frameIndex;
    frameIndex++;
}

void FluidCube::offRenderingCheckFbyF() {
    namespace fs = std::filesystem;

    auto makePath = [](const std::string& prefix, int idx) {
        std::ostringstream oss;
        oss << "offline_data/" << prefix << "_" << std::setw(4) << std::setfill('0') << idx << ".bin";
        return oss.str();
    };

    const std::string firstDensity = makePath("density", 0);

    fs::create_directory("offline_data");

    if (fs::exists(firstDensity)) {
        offlineFileLoadedFF = true;
        return;
    }

    offlineFileLoadedFF = false;

    for (int i = 0; i < totalFrames; ++i) {
        update(0.0083f);                    // 120 FPS step
        Log("Sim frame " + std::to_string(i));

        densityFrames.push_back(density);
        colorFrames.push_back(fluidColors);

        // store velocity magnitude (same convention as before)
        std::vector<float> velMag(totalCells);
        for (int j = 0; j < totalCells; ++j) {
            float vx = vX[j], vy = vY[j], vz = vZ[j];
            velMag[j] = std::sqrt(vx*vx + vy*vy + vz*vz);
        }
        velocityFrames.push_back(std::move(velMag));

        {
            std::ofstream dFile(makePath("density",  i), std::ios::binary);
            dFile.write(reinterpret_cast<const char*>(densityFrames.back().data()),
                        totalCells * sizeof(float));
        }
        {
            std::ofstream cFile(makePath("color",    i), std::ios::binary);
            cFile.write(reinterpret_cast<const char*>(colorFrames.back().data()),
                        totalCells * sizeof(int));
        }
        {
            std::ofstream vFile(makePath("velocity", i), std::ios::binary);
            vFile.write(reinterpret_cast<const char*>(velocityFrames.back().data()),
                        totalCells * sizeof(float));
        }
    }

    Log("Saved " + std::to_string(totalFrames) + " offline frames to disk.");
    startRenderingOnceFF = true;
}



void FluidCube::restartOfflineRendering(){
    frameIndex = 0;
}
