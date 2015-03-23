#ifndef FUNCTIONSGEO_HPP
#define FUNCTIONSGEO_HPP

#include <iostream>
#include <cmath>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include "Variables.hpp"
#include "Data.hpp"
#include "FunctionsConnect.hpp" 
#include "ParamPhys.hpp"
R2 GravityCenter(cell cc,int nbnode);

int centerdomain(Data & d,Mesh & Mh);

double AverageQuantity(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamPhysic & Param,int numGr,int var);

double Vjr(Mesh & Mh,int j,int r);

double Vr(Mesh & Mh,int r, TabConnecInv & tab);

double geomCorrection(Data & d,Mesh & Mh,int j);

cell ControlVolume(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,int numGr);

cell SubControlVolume(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,int numGr,int j);

void OneLayerStencilNodeQ(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,int *& Stencil,int NStencil,int numGr);

void TwoLayerStencilNodeQ(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,int *& Stencil,int NStencil,int numGr);
#endif
