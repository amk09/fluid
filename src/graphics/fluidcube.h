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
    // These will be added to 'simulation' to call to enable interaction.
    // To test them right now, just add them in the 'init' function
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
    GLuint m_fullscreen_vao;
    GLuint m_fullscreen_vbo;
    void initVoxelGeometry();
    void drawVoxel(Shader* shader);
    void initFullscreenQuad();
    void drawVolume(Shader* shader);

    // Below are the parameters for Stable Fluids
    int size; // 'N' in the paper and the rest of codes, how many grids you want to divide the cube edges
    float dt; // Time Step, do we really need this? Since we update it in our tick with a fixed time step? I decieded not to use it at start.
    float diff; // The parameter controls the diffuse // For Density to Spread Out in the Fluid Cube
    float visc; // Thickness of the fluids // For Velocity to Spread Out in the Fluid Cube

    // Although they are written in Previous and Current versions.
    // Actually they are just stated to swap to use instead of a really previous or current state. （eg. temp? value)
    vector<float> density0; // The previous density
    vector<float> density; // The current density
    vector<float> vX0; // Previous x-direction velocity
    vector<float> vX; // Current x-direction velocity
    vector<float> vY0; // Previous y-direction velocity
    vector<float> vY; // Current y-direction velocity
    vector<float> vZ0; // Previous z-direction velocity
    vector<float> vZ; // Current z-direction velocity

    // Below are the main functions should be done in Stable Fluids
    int iter = 4; // This is the key parameter to solve stable formula (20 in the paper; 10 in the video; 4 in the website)

    // b here is coord number, where 1 for x-coord, 2 for y-coord, 3 for z-coord (0 for nothing)
    // diffuse will be (visc for velocity diffuse, diff for density diffuse)
    // iter and N will not be used here because we already have this in this class
    void diffuse(int b, vector<float> &dataWrittenTo, vector<float> &dataReadFrom, float diffuse, float dt);

    // a & c are the parameters required in the paper.
    // iter and N will not be used here because we already have this in this class
    void lin_solve(int b, vector<float> &dataWrittenTo, vector<float> &dataReadFrom, float a, float c);

    // iter and N will not be used here because we already have this in this class
    void project(vector<float> &velocityX, vector<float> &velocityY, vector<float> &velocityZ, vector<float> &p, vector<float> div);

    // b is the same meaning shown above
    // N will not be used here because we already have this in this class
    void set_bnd(int b, vector<float> &dataWrittenTo);

    // b is the same meaning shown above
    // N will not be used here because we already have this in this class
    void advect(int b, vector<float> &dataWrittenTo, vector<float> &dataReadFrom, vector<float> &velocityX, vector<float> &velocityY, vector<float> &velocityZ, float dt);

    // Some getters here only to reduce the difficulty of extracting data
    int index(int x, int y, int z){
        return (x) + (y) * size + (z) * size * size;
    }

    // Some updates
    void uploadDensityToGPU();

    // Some tests
    void test();
};

#endif // FLUIDCUBE_H
