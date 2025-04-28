#include "fluidcube.h"
#include "graphics/shader.h"
#include "graphics/solver.h"
#include "graphics/obstacle.h"
#include "graphics/vectorfield.h"
#include <iomanip>
#include <iostream>

FluidCube::FluidCube()
    : m_vao(0),
    m_vbo(0),
    m_ibo(0),
    m_densityTexture(0),
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




void FluidCube::draw(Shader *shader)
{
    // // For Voxel
    // drawVoxel(shader);

    // For Ray Marchign
    drawVolume(shader);
}


void FluidCube::addSource(vector<float> x, vector<float> x0)
{
    for (int i = 0; i < totalCells; i++)
    {
        x[i] += dt * x0[i];
    }
}

void FluidCube::densityFade(float dt)
{
    for (int i = 0; i < totalCells; i++)
    {
        density[i] *= (1 - 0.5*dt);
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
    density[idx] = std::clamp(density[idx] + amount, 0.0f, 1.0f);

    // More aggressive color setting - apply color even for small density changes
    if (oldDensity < 0.1f || amount > 0.05f) {
        fluidColors[idx] = colorType;
    }

    // Add small amounts to neighboring cells for smoothness
    float neighborAmount = amount * 0.2f;

    // Add to 6 direct neighbors
    for (int i = -1; i <= 1; i += 2) {
        int nx = x + i;
        int ny = y;
        int nz = z;
        if (nx >= 1 && nx < size-1) {
            int nidx = index(nx, ny, nz);
            oldDensity = density[nidx];
            density[nidx] = std::clamp(density[nidx] + neighborAmount, 0.0f, 1.0f);
            // More aggressive color propagation
            if (oldDensity < 0.1f || neighborAmount > 0.02f) {
                fluidColors[nidx] = colorType;
            }
        }
    }

    for (int i = -1; i <= 1; i += 2) {
        int nx = x;
        int ny = y + i;
        int nz = z;
        if (ny >= 1 && ny < size-1) {
            int nidx = index(nx, ny, nz);
            oldDensity = density[nidx];
            density[nidx] = std::clamp(density[nidx] + neighborAmount, 0.0f, 1.0f);
            // More aggressive color propagation
            if (oldDensity < 0.1f || neighborAmount > 0.02f) {
                fluidColors[nidx] = colorType;
            }
        }
    }

    for (int i = -1; i <= 1; i += 2) {
        int nx = x;
        int ny = y;
        int nz = z + i;
        if (nz >= 1 && nz < size-1) {
            int nidx = index(nx, ny, nz);
            oldDensity = density[nidx];
            density[nidx] = std::clamp(density[nidx] + neighborAmount, 0.0f, 1.0f);
            // More aggressive color propagation
            if (oldDensity < 0.1f || neighborAmount > 0.02f) {
                fluidColors[nidx] = colorType;
            }
        }
    }
}



// Function to upload color field to GPU
void FluidCube::uploadColorFieldToGPU() {
    if (m_colorTexture == 0)
        glGenTextures(1, &m_colorTexture);

    glBindTexture(GL_TEXTURE_3D, m_colorTexture);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

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

    shader->setUniform("densityTex", 0);
    shader->setUniform("colorTex", 1);
    shader->setUniform("fluidStyle", 1);
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
    std::vector<float> shellDensity(size * size * size, 0.0f);
    std::vector<float> shellColors(size * size * size, 0.0f);

    // Detect shell cells (cells that have density and at least one neighbor without density)
    // Using precise threshold for better shell definition
    for (int k = 1; k < size-1; k++) {
        for (int j = 1; j < size-1; j++) {
            for (int i = 1; i < size-1; i++) {
                int idx = index(i, j, k);

                // Skip cells with very low density
                if (density[idx] < 0.01f) continue;

                // Check if any of the six adjacent cells has low density
                bool isShell = false;
                if (density[index(i-1, j, k)] < 0.01f || density[index(i+1, j, k)] < 0.01f ||
                    density[index(i, j-1, k)] < 0.01f || density[index(i, j+1, k)] < 0.01f ||
                    density[index(i, j, k-1)] < 0.01f || density[index(i, j, k+1)] < 0.01f) {
                    isShell = true;
                }

                // If it's a shell cell, enhance its density and preserve color
                if (isShell) {
                    // Use moderate density enhancement for clarity
                    shellDensity[idx] = density[idx] * 3.5f;
                    shellColors[idx] = static_cast<float>(fluidColors[idx]);
                }
            }
        }
    }

    // Upload this modified density field and color field to GPU
    uploadCustomDensityToGPU(shellDensity);
    uploadCustomColorToGPU(shellColors);
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
    uploadColorFieldToGPU();
}

// Add initialization of color field in init method
void FluidCube::init(int size, float diffuse, float viscosity){
    this->size = size;
    this->diff = diffuse;
    this->visc= viscosity;
    this->totalCells = size*size*size;

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
}



void FluidCube::update(float dt)
{
    // 1. Add vorticity confinement
    addVorticityConfinement(vX0, vY0, vZ0, dt, size, m_vorticityStrength);

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

    // Make Sure Obstacle's Density is 0.0
    densityInsideObstacles(density);

    // Clean the density value in density0
    empty_den();

    // Perform extra color enhancement to fix possible inconsistencies
    for (int i = 0; i < totalCells; i++) {
        // For cells with significant density but no color, force default color
        if (density[i] > 0.03f && fluidColors[i] == 0) {
            // Assign a safe default color (0) instead of using the global color
            fluidColors[i] = 0;
        }
    }


    // Propagate colors to ensure they follow density properly
    propagateColors();

    // Upload both density and color to GPU
    uploadDensityToGPU();
    uploadColorFieldToGPU();
}


// Add this new function to propagate colors with fluid movement
void FluidCube::propagateColors() {
    // For each cell, if it has significant density but no color assigned yet,
    // get color from neighboring cells with highest density
    for (int k = 1; k < size-1; k++) {
        for (int j = 1; j < size-1; j++) {
            for (int i = 1; i < size-1; i++) {
                int idx = index(i, j, k);

                // Only process cells with density but no color - even lower threshold
                if (density[idx] > 0.005f && fluidColors[idx] == 0) {
                    float maxDensity = 0.0f;
                    int bestColor = 0;

                    // Check immediate neighbors first (6-connected)
                    // This prioritizes direct color transfer over diagonal connections
                    bool foundDirectNeighbor = false;

                    // Check direct neighbors first (more accurate color propagation)
                    const int directNeighbors[6][3] = {
                        {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}
                    };

                    for (int n = 0; n < 6; n++) {
                        int nx = i + directNeighbors[n][0];
                        int ny = j + directNeighbors[n][1];
                        int nz = k + directNeighbors[n][2];

                        if (nx >= 1 && nx < size-1 &&
                            ny >= 1 && ny < size-1 &&
                            nz >= 1 && nz < size-1) {

                            int nidx = index(nx, ny, nz);

                            // If direct neighbor has color and significant density
                            if (density[nidx] > 0.05f && fluidColors[nidx] > 0) {
                                fluidColors[idx] = fluidColors[nidx];
                                foundDirectNeighbor = true;
                                break;
                            }
                        }
                    }

                    // If no direct neighbor with color, check all neighbors
                    if (!foundDirectNeighbor) {
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
                                    if (ni >= 1 && ni < size-1 &&
                                        nj >= 1 && nj < size-1 &&
                                        nk >= 1 && nk < size-1) {

                                        int nidx = index(ni, nj, nk);

                                        // If neighbor has higher density and a color
                                        if (density[nidx] > maxDensity && fluidColors[nidx] > 0) {
                                            maxDensity = density[nidx];
                                            bestColor = fluidColors[nidx];
                                        }
                                    }
                                }
                            }
                        }

                        // Apply the best color found
                        if (bestColor > 0) {
                            fluidColors[idx] = bestColor;
                        }
                    }
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


    std::vector<float> enhancedDensity(density.size());

    float minDensity = 1.0f;
    float maxDensity = 0.0f;
    for (int i = 0; i < density.size(); i++) {
        minDensity = std::min(minDensity, density[i]);
        maxDensity = std::max(maxDensity, density[i]);
    }

    // Apply contrast enhancement
    float range = maxDensity - minDensity;
    if (range > 0.01f) {  // Only if there's actual variation
        for (int i = 0; i < density.size(); i++) {
            // Normalize to [0,1] range
            float normalized = (density[i] - minDensity) / range;
            // Apply power curve for more contrast
            enhancedDensity[i] = std::pow(normalized, 1.5f);
        }
    } else {
        enhancedDensity = density;  // No enhancement if uniform
    }

    glBindTexture(GL_TEXTURE_3D, m_densityTexture);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    // glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    // glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
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


