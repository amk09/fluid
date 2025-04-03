#ifndef FLUIDCUBE_H
#define FLUIDCUBE_H

#include <GL/glew.h>
#include <vector>

#include <Eigen/Dense>

using namespace std;

class Shader;

class FluidCube
{
public:
    // Initilize the Fluid Cub
    FluidCube();

    // Initilize the Fluid Cub
    void init(int size, float diffuse, float viscosity);

    // Set the modelMatrix, define where we can see the cube
    void setModelMatrix(const Eigen::Affine3f &model);

    // Toggle the Wire Frame to see the cube
    void toggleWireframe();

    // Update using a dt
    void update(float dt);

    // Render the Stable Fluids
    void draw(Shader *shader);

    // Add Density and Velocity
    void addDensity(int x, int y, int z, float amount);
    void addVelocity(int x, int y, int z, float amountX, float amountY, float amountZ);

private:
    // Below are the units for OpenGL rendering
    GLuint m_vao;
    GLuint m_vbo;
    GLuint m_ibo;
    GLuint m_densityTexture;
    Eigen::Matrix4f m_modelMatrix;
    bool m_wireframe;
    // Plus the unit - voxel or the cell
    GLuint m_voxelVao, m_voxelVbo, m_voxelIbo;
    void initVoxelGeometry();

    // Below are the parameters for Stable Fluids
    int size; // N in the paper, how many grids you want to divide the cube edges
    float dt; // Time Step, do we really need this? Since we update it in our tick with a fixed time step? I decieded not to use it at start.
    float diffuse; // The parameter controls the diffuse
    float viscosity; // Thickness of the fluids

    vector<float> density0; // The previous density
    vector<float> density; // The current density
    vector<float> vX0; // Previous x-direction velocity
    vector<float> vX; // Current x-direction velocity
    vector<float> vY0; // Previous y-direction velocity
    vector<float> vY; // Current y-direction velocity
    vector<float> vZ0; // Previous z-direction velocity
    vector<float> vZ; // Current z-direction velocity

    // Some getters here only to reduce the difficulty of extracting data
    int index(int x, int y, int z){
        return (x) + (y) * size + (z) * size * size;
    }
};

#endif // FLUIDCUBE_H
