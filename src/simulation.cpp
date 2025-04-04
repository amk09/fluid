#include "simulation.h"
#include "graphics/meshloader.h"

#include <iostream>

using namespace Eigen;

Simulation::Simulation() {}

void Simulation::init()
{
    // STUDENTS: This code loads up the tetrahedral mesh in 'example-meshes/single-tet.mesh'
    //    (note: your working directory must be set to the root directory of the starter code
    //    repo for this file to load correctly). You'll probably want to instead have this code
    //    load up a tet mesh based on e.g. a file path specified with a command line argument.
    // std::vector<Vector3d> vertices;
    // std::vector<Vector4i> tets;
    // if (MeshLoader::loadTetMesh(":/example-meshes/single-tet.mesh", vertices, tets)) {
    //     // STUDENTS: This code computes the surface mesh of the loaded tet mesh, i.e. the faces
    //     //    of tetrahedra which are on the exterior surface of the object. Right now, this is
    //     //    hard-coded for the single-tet mesh. You'll need to implement surface mesh extraction
    //     //    for arbitrary tet meshes. Think about how you can identify which tetrahedron faces
    //     //    are surface faces...
    //     std::vector<Vector3i> faces;
    //     faces.emplace_back(1, 0, 2);
    //     faces.emplace_back(2, 0, 3);
    //     faces.emplace_back(3, 1, 2);
    //     faces.emplace_back(3, 0, 1);
    //     m_shape.init(vertices, faces, tets);
    // }
    // m_shape.setModelMatrix(Affine3f(Eigen::Translation3f(0, 2, 0)));

    // initGround();

    fluidCube.init(256, 0.0, 0.0);
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
