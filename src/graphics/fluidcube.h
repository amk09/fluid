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
    // Initialize the Fluid Cube
    FluidCube();

    // Initialize the Fluid Cube
    void init(int size, float diffuse, float viscosity);

    // Set the modelMatrix, define where we can see the cube
    void setModelMatrix(const Eigen::Affine3f &model);

    // Toggle the Wire Frame to see the cube
    void toggleWireframe();
    void toggleShellRendering();
    // Update using a dt
    void update(float dt);

    // Render the Stable Fluids
    void draw(Shader *shader);

    // Add Density and Velocity
    // These will be added to 'simulation' to call to enable interaction.
    void addDensity(int x, int y, int z, float amount);
    // New method to add density with specific color
    void addDensityWithColor(int x, int y, int z, float amount, int colorType);
    void addVelocity(int x, int y, int z, float amountX, float amountY, float amountZ);

    // Clear all fluids from the cube
    void clearAllFluids();

    // New method to clear obstacles
    void clearObstacles();

    //debugger
    void visualizeVelocity();
    float getTotalDensity();

    void drawShellOnly(Shader *shader);
    void uploadCustomDensityToGPU(const std::vector<float>& customDensity);
    void uploadCustomColorToGPU(const std::vector<float>& customColors);
    void setColorMap(int colorType);
    int getColorMap() const;
    void setRenderMode(int mode);
    int getRenderMode() const;

    void setVorticityStrength(float strength);
    float getVorticityStrength() const;

    // Get current color type
    int getCurrentColorType() const { return m_colorMapType; }
    int getSize() const {return size;}
    std::vector<float>& getDensity() { return density; }

    void detectShell();

    // For Offline Rendering
    void restartOfflineRendering();
    void uploadVelocityToGPU();
    void toggleVelocityColoring() { m_useVelocityColor = !m_useVelocityColor; }
    void setVelocityScale(float scale) { m_velocityScale = scale; }
    void setVelocityBlend(float blend) { m_velocityBlend = blend; }
    bool isUsingVelocityColor() const { return m_useVelocityColor; }
    float getVelocityScale() const { return m_velocityScale; }
    float getVelocityBlend() const { return m_velocityBlend; }

    void offRenderingCheckFbyF();
    void renderNextOfflineFrame();
    float frameRate = 120.0f;
private:
    // Below are the units for OpenGL rendering
    GLuint m_vao;
    GLuint m_vbo;
    GLuint m_ibo;
    GLuint m_densityTexture;
    GLuint m_colorTexture;  // New texture for color field
    GLuint m_shellDensityTexture;
    GLuint m_velocityTexture = 0;
    bool m_useVelocityColor = true;
    float m_velocityScale = 3.0f;  // Starting scale - adjust based on your simulation
    float m_velocityBlend = 0.5f;  // 50% blend initially
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
    float dt; // Time Step, do we really need this? Since we update it in our tick with a fixed time step? I decided not to use it at start.
    float diff; // The parameter controls the diffuse // For Density to Spread Out in the Fluid Cube
    float visc; // Thickness of the fluids // For Velocity to Spread Out in the Fluid Cube
    int totalCells;

    // Although they are written in Previous and Current versions.
    // Actually they are just stated to swap to use instead of a really previous or current state. (eg. temp? value)
    vector<float> density0; // The previous density
    vector<float> density; // The current density
    vector<float> vX0; // Previous x-direction velocity
    vector<float> vX; // Current x-direction velocity
    vector<float> vY0; // Previous y-direction velocity
    vector<float> vY; // Current y-direction velocity
    vector<float> vZ0; // Previous z-direction velocity
    vector<float> vZ; // Current z-direction velocity
    vector<int> fluidColors; // Color map type for each cell
    vector<float> velocity; // Velocity field for the fluid

    // Below are the main functions should be done in Stable Fluids
    int iter = 4; // This is the key parameter to solve stable formula (20 in the paper; 10 in the video; 4 in the website)

    float m_vorticityStrength = 1000.0f;  // Default vorticity strength

    // Some getters here only to reduce the difficulty of extracting data
    int index(int x, int y, int z);

    int m_colorMapType = 0;  // 0: default, 1: blue, 2: purple, 3: cyan-yellow, 4: orange-grey, 5: fire, 6: ocean, 7: plasma
    int m_renderMode = 0;    // 0: volume, 1: shell

    // Some updates
    void uploadDensityToGPU();
    void uploadColorFieldToGPU();

    // Helper method for texture uploads
    void setupTextureParameters(GLuint texture);

    // Helper method for adding density to adjacent cells
    void addDensityToCell(int x, int y, int z, float amount, int colorType, float oldDensity);

    // Helper method for checking neighbor validity
    bool isValidCell(int x, int y, int z);

    // Helper method for color propagation
    bool findBetterColorFromNeighbors(int idx, int i, int j, int k);

    // Handle the increasing density bug caused by inaccurate calculation
    void addSource(vector<float>& x, vector<float>& x0);
    void empty_vel();
    void empty_den();
    void densityFade(float dt);

    // Some tests
    void test();
    void propagateColors();

    void fountainGeneration();
    void fountainGenerationTopDown();

    // Offline Rendering
    std::vector<std::vector<float>> densityFrames;
    std::vector<std::vector<int>> colorFrames;
    std::vector<std::vector<float>> velocityFrames;
    void offRenderingCheck();
};

inline int FluidCube::getColorMap() const { return m_colorMapType; }

inline void FluidCube::setRenderMode(int mode) { m_renderMode = mode; }

inline int FluidCube::getRenderMode() const { return m_renderMode; }

inline void FluidCube::setVorticityStrength(float strength) { m_vorticityStrength = strength; }

inline float FluidCube::getVorticityStrength() const { return m_vorticityStrength; }

inline int FluidCube::index(int x, int y, int z){
    return (x) + (y) * size + (z) * size * size;
}

// Helper method to check if a cell is valid
inline bool FluidCube::isValidCell(int x, int y, int z) {
    return x >= 1 && x < size-1 &&
           y >= 1 && y < size-1 &&
           z >= 1 && z < size-1;
}

#endif // FLUIDCUBE_H
