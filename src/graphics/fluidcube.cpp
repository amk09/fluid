#include "fluidcube.h"
#include "graphics/shader.h"
#include "graphics/solver.h"
#include "graphics/vectorfield.h" 
#include <iomanip>
#include <iostream>

float FluidCube::getTotalDensity() {
    float total = 0.0f;
    for (int i = 0; i < density.size(); i++) {
        total += density[i];
    }
    return total;
}

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

void FluidCube::update(float dt)
{   

    if (initialDensity < 0.0f) {
        initialDensity = getTotalDensity();
    }

    addFountainForce(vX, vY, vZ, density, size, 5.f, dt);
    circulateDensity(density, size, dt, 0.1f);
    
    // // Debug output
    // static int frameCount = 0;
    // frameCount++;
    // if (frameCount % 30 == 0) {
    //     std::cout << "Current density: " << getTotalDensity() 
    //               << " Initial: " << initialDensity << std::endl;
    // }

//                        WHAT WE HAD
    // 1. Diffuse velocity
    diffuse_velocity(1, vX0, vX, visc, dt, iter, size);
    diffuse_velocity(2, vY0, vY, visc, dt, iter, size);
    diffuse_velocity(3, vZ0, vZ, visc, dt, iter, size);

    // 2. Project velocity
    project(vX0, vY0, vZ0, vX, vY, iter, size);

    // 3. Advect velocity
    advect(1, vX, vX0, vX0, vY0, vZ0, dt, size);
    advect(2, vY, vY0, vX0, vY0, vZ0, dt, size);
    advect(3, vZ, vZ0, vX0, vY0, vZ0, dt, size);

    // 4. Project again
    project(vX, vY, vZ, vX0, vY0, iter, size);

    // 5. Diffuse density
    diffuse_density(0, density0, density, diff, dt, iter, size);

    // 6. Advect density
    advect(0, density, density0, vX, vY, vZ, dt, size);

    conserveDensity(density, initialDensity, size);

    uploadDensityToGPU();
}

void FluidCube::draw(Shader *shader)
{
    // // For Voxel
    // drawVoxel(shader);

    // For Ray Marchign
    drawVolume(shader);
}

void FluidCube::init(int size, float diffuse, float viscosity){
    this->size = size;
    this->diff = diffuse;
    this->visc= viscosity;

    // Below is the number of cells in this cube
    int totalCells = size * size * size;

    // We resize the density and velocity vectors with 0.0f value
    // Not just make sure they are pre-allocated, but also there is no density and velocity in the fluid cube
    density0.resize(totalCells, 0.0f);
    density.resize(totalCells, 0.0f);

    vX0.resize(totalCells, 0.0f);
    vX.resize(totalCells, 0.0f);

    vY0.resize(totalCells, 0.0f);
    vY.resize(totalCells, 0.0f);

    vZ0.resize(totalCells, 0.0f);
    vZ.resize(totalCells, 0.0f);

    // // Call Voxel to Build
    // initVoxelGeometry();
    initFullscreenQuad();

    // Inject Density and then upload it to GPU openGL
    test();
    uploadDensityToGPU();
}

void FluidCube::addDensity(int x, int y, int z, float amount){
    // The reason why we clamp here is just to make sure density is always between 0 to 1
    // Clamp after density
    int idx = index(x, y, z);
    density[idx] = std::clamp(density[idx] + amount, 0.0f, 1.0f);
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

void FluidCube::drawVolume(Shader* shader){
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, m_densityTexture);
    shader->setUniform("fluidStyle", 1);
    // shader->setUniform("densityTex", 0);
    shader->setUniform("size" , size);

    glBindVertexArray(m_fullscreen_vao);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glBindVertexArray(0);
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


// Below is for the test methods
void FluidCube::test() {
    // int center = size / 2;

    // // A Sphere
    // int radius = size / 8;
    // int radius2 = radius * radius;

    // for (int z = 0; z < size; ++z){
    //     for (int y = 0; y < size; ++y){
    //         for (int x = 0; x < size; ++x) {
    //             int dx = x - center;
    //             int dy = y - center;
    //             int dz = z - center;
    //             int dist2 = dx * dx + dy * dy + dz * dz;

    //             if (dist2 < radius2)
    //                 addDensity(x, y, z, 1.0f);
    //         }
    //     }
    // }

    // // VectorFields::addVortexField(vX, vY, vZ, density, size, 2.0f);

    // // // A Cube
    // // for (int z = 0; z < size; ++z)
    // //     for (int y = 0; y < size; ++y)
    // //         for (int x = 0; x < size; ++x)
    // //             addDensity(x, y, z, 0.5f);

    int center = size / 2;
    
    // Create a pool of water at the bottom third of the container
    int waterHeight = size / 3;
    
    for (int z = 0; z < size; ++z) {
        for (int y = 0; y < waterHeight; ++y) {
            for (int x = 0; x < size; ++x) {
                // Add density that decreases with height for natural water look
                float heightFactor = 1.0f - (float)y / waterHeight;
                addDensity(x, y, z, 0.8f * heightFactor);
            }
        }
    }

    if (1) {
        int fountainRadius = size / 10;
        for (int z = center - fountainRadius; z <= center + fountainRadius; z++) {
            for (int x = center - fountainRadius; x <= center + fountainRadius; x++) {
                float dx = x - center;
                float dz = z - center;
                
                if ((dx*dx + dz*dz) <= fountainRadius*fountainRadius) {
                    for (int y = 0; y < 3; y++) {
                        addVelocity(x, y, z, 0.0f, 5.0f, 0.0f);
                    }
                }
            }
        }
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