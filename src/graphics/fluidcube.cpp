#include "fluidcube.h"
#include "graphics/shader.h"
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

void FluidCube::update(float dt)
{
    uploadDensityToGPU();
}

void FluidCube::draw(Shader *shader)
{
    // // For Voxel
    // drawVoxel(shader);

    // For Ray Marchign
    drawVolume(shader);
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

void FluidCube::drawVolume(Shader* shader){
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, m_densityTexture);
    shader->setUniform("densityTex", 0);
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
    int center = size / 2;

    // A Sphere
    int radius = size / 2;
    int radius2 = radius * radius;

    for (int z = 0; z < size; ++z)
        for (int y = 0; y < size; ++y)
            for (int x = 0; x < size; ++x) {
                int dx = x - center;
                int dy = y - center;
                int dz = z - center;
                int dist2 = dx * dx + dy * dy + dz * dz;

                if (dist2 < radius2)
                    addDensity(x, y, z, 0.5f);
            }

    // // A Cube
    // for (int z = 0; z < size; ++z)
    //     for (int y = 0; y < size; ++y)
    //         for (int x = 0; x < size; ++x)
    //             addDensity(x, y, z, 0.5f);

    // std::cout << "=== Density slice at z = " << center << " ===\n";
    for (int y = 0; y < size; ++y) {
        for (int x = 0; x < size; ++x) {
            float d = density[index(x, y, center)];
            // if (d > 0.01f)
            //     std::cout << std::fixed << std::setprecision(2) << d << " ";
            // else
            //     std::cout << " .   ";
        }
        // std::cout << '\n';
    }
}
