#ifndef FLUIDCUBE_H
#define FLUIDCUBE_H

#include <GL/glew.h>
#include <vector>

#include <Eigen/Dense>

using namespace std;

class Shader;

struct SolidObject {
    enum Type { CUBE, SPHERE };
    Type type;
    Eigen::Vector3f position;
    Eigen::Vector3f size;      // For cube: half-extents
    float radius;              // For sphere
};

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

    
    //debugger
    void visualizeVelocity();
    float getTotalDensity();


    void drawShellOnly(Shader *shader);
    void uploadCustomDensityToGPU(const std::vector<float>& customDensity);
    void setColorMap(int colorType);
    int getColorMap() const;
    void setRenderMode(int mode);
    int getRenderMode() const;

    void setVorticityStrength(float strength);
    float getVorticityStrength() const;




private:

    float initialDensity = -1.0f;
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



    std::vector<SolidObject> m_solidObjects;
    std::vector<bool> m_solidCells; // Flag for solid cells

    void updateSolidCells();
    bool isInsideSolid(int i, int j, int k) const;

    Eigen::Vector3f getSolidNormal(int i, int j, int k) const; //same as get_norm



    // Below are the parameters for Stable Fluids
    int size; // 'N' in the paper and the rest of codes, how many grids you want to divide the cube edges
    float dt; // Time Step, do we really need this? Since we update it in our tick with a fixed time step? I decieded not to use it at start.
    float diff; // The parameter controls the diffuse // For Density to Spread Out in the Fluid Cube
    float visc; // Thickness of the fluids // For Velocity to Spread Out in the Fluid Cube
    int totalCells;


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
    float m_vorticityStrength = 0.5f;
    // Some getters here only to reduce the difficulty of extracting data
    int index(int x, int y, int z);

    int m_colorMapType = 0;  // 0: default, 1: blue, 2: purple, 3: cyan-yellow, 4: orange-grey
    int m_renderMode = 0;

    // Some updates
    void uploadDensityToGPU();
    
     // Handle the increasing density bug caused by inaccuate calculation
  
    void addSource(vector<float> x, vector<float> x0);
    void empty_vel();
    void empty_den();
    void densityFade(float dt);
    // Some tests
    void test();
};
inline int FluidCube::getColorMap() const { return m_colorMapType; }

inline void FluidCube::setRenderMode(int mode) { m_renderMode = mode; }

inline int FluidCube::getRenderMode() const { return m_renderMode; }

inline void FluidCube::setVorticityStrength(float strength) { m_vorticityStrength = strength; }

inline float FluidCube::getVorticityStrength() const { return m_vorticityStrength; }

inline int FluidCube::index(int x, int y, int z){
    return (x) + (y) * size + (z) * size * size;
}
#endif // FLUIDCUBE_H
