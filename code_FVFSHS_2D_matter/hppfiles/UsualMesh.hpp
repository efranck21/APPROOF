#ifndef USUALMESH_HPP
#define USUALMESH_HPP
#include "ClassMesh.hpp"
#include "Variables.hpp"
#include "Data.hpp"



Mesh ChoiceMesh(Data & d);

Mesh CartesianMeshT(Data & d);

Mesh RandomCartesianMeshT(Data & d);

Mesh KershawMesh(Data & d);

Mesh CartesianMesh(Data & d);

Mesh RandomCartesianMesh(Data & d);

Mesh SmoothCartesianMesh(Data & d);

Mesh CheckMesh(Data & d);

Mesh DoubleSquareMeshQ(Data & d);

#endif



 
