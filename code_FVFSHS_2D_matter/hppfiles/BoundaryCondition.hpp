#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

#include <iostream>
#include <cmath>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include <cstdlib>
#include "Variables.hpp"
#include "Data.hpp"
#include "Functions.hpp"
#include "AnalyticSolutions.hpp"
#include "Diagnostics.hpp"
#include "Flux.hpp"
#include "Source.hpp"
#include "ParamPhys.hpp"

//////// General functions ///////////

R2 Norext(Data & d,Mesh & Mh, int numGr);

void BoundaryCondition(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

//////// Advection functions ////////

void BoundaryConditionAdvection(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BCDirichlet_SolExact_Advection(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

/////// Diffusion functions /////////

void BoundaryConditionDiff(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BCDirichlet_SolExact_diff(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BC_Periodic_diff(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

/////////// P1 functions ////////////

void BoundaryConditionP1(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BCNeumann_P1(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BC_Periodic_P1(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BCDirichlet_P1(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BC_2DRiemann_P1(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

/////////// P1 Matter functions ////////////

void BoundaryConditionP1Matter(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BC_Neumann_P1Matter(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BC_BlackBody_Left(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

///////// Euler functions //////////

void BoundaryConditionEuler(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BCNeumannEuler(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BCDirichlet_ST_Euler(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BCDirichlet_ST_EulerAverage(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BC_Periodic_Euler(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BC_NeumanY_PeriodicX_Euler(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);



/////////// M1 functions ////////////

void BoundaryConditionM1(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BC_Neumann_M1(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BC_Periodic_M1(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BC_Dirichlet_SolExact_M1(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);


/////////// M1 functions Matter ////////////

void BoundaryConditionM1Matter(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BC_Neumann_M1M(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);

void BC_Periodic_M1M(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time);



#endif
