#ifndef TIMEDISCRET_HPP
#define TIMEDISCRET_HPP
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "Flux.hpp"
#include "Source.hpp"
#include "Initialisation.hpp"
#include <fstream>
#include "Functions.hpp"
#include "ClassMesh.hpp"
#include "Variables.hpp"
#include "WriteData.hpp"
#include "ParamPhys.hpp"

using namespace std;

void ExplicitDiscret(Data & d,Mesh & Mh, variable & v,TabConnecInv & TConnectInv, ParamPhysic & Param,double dt,int & ntAnim, double & time);


void SemiImplicit(Data & d,Mesh & Mh, variable & v,TabConnecInv & TConnectInv, ParamPhysic & Param,double dt,int & ntAnim, double & time);

int SIschemeImplemented(Data & d,Mesh & Mh, variable & v,TabConnecInv & TConnectInv, ParamPhysic & Param);

#endif
