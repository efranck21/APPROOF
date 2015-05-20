#ifndef INITIALISATION_HPP
#define INITIALISATION_HPP
#include "UsualMesh.hpp"
#include "Data.hpp"
#include "ClassMesh.hpp"
#include "ParamPhys.hpp"

using namespace std;


void ChoiceInit(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param, double & InitTime); 

void InitP1(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param);

void InitE(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param);

void InitDiff(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param);

void Initadv(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param);

void InitM1(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param);

void InitP1Matter(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param);

void InitM1Matter(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param);

void InitP1Compton(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param);

#endif
