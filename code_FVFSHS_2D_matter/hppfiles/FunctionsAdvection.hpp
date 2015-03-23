#ifndef FUNCTIONSADVECTION_HPP
#define FUNCTIONSADVECTION_HPP
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Variables.hpp"
#include "ClassMesh.hpp"
#include "Functions.hpp"
#include "Initialisation.hpp"
#include "Functionsinit.hpp"
#include "Tensor.hpp"
#include "R2.hpp"
#include "ParamPhys.hpp"


int EnsembleR(Data & d,Mesh & Mh, variable & v,ParamPhysic & Param,int r,int numCell, R2 a);

double variablekr(Data & d,Mesh & Mh, variable & v,int var,TabConnecInv & tab,ParamPhysic & Param,int numCell,int numGr, R2 a,int ordre);

double VertexUpwind(Data & d,Mesh & Mh, variable & v,int var,TabConnecInv & tab,ParamPhysic & Param,int numCell,int r, R2 a);

double EdgeUpwind(Data & d,Mesh & Mh, variable & v,int var,TabConnecInv & tab,ParamPhysic & Param,int numCell,int k,int r, R2 a);

R2 Gradiant(Data & d,Mesh & Mh,variable & v,int var,TabConnecInv & tab,ParamPhysic & Param,int numGr);

double VertexSlope(Data & d,Mesh & Mh,variable & v,int var,TabConnecInv & tab,ParamPhysic & Param,int numCell,R2 Grad);

double EdgeSlope(Data & d,Mesh & Mh,variable & v,int var,TabConnecInv & tab,ParamPhysic & Param,int numCell,R2 Grad);

double Limiter(ParamPhysic & Param,double m);

double InterLineVertex(Data d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamPhysic & Param,int var,int ordre,int j,int r);

#endif
