#ifndef DIAGNOSTICS_HPP
#define DIAGNOSTICS_HPP

#include <iostream>
#include <cmath>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include "Functions.hpp"
#include "Variables.hpp"
#include "Data.hpp"
#include "Flux.hpp"
#include "Source.hpp"
#include "ParamPhys.hpp"

double NormLP_Error_One_variable(Data & d, Mesh & Mh,variable & v,ParamPhysic & Param,double p,double t,int var);

double NormLP_Error_All_variables_Euler(Data & d, Mesh & Mh,variable & v,ParamPhysic & Param,double p,double t);

double NormLP_Error_All_variables(Data & d, Mesh & Mh,variable & v,ParamPhysic & Param,double p,double t);

double NormLP_Stability(Data & d, Mesh & Mh,variable & v,ParamPhysic & Param,double p);

double Mass(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,double time); 

void Diagnostics_Quantities(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,int nt,double time); 

double Diagnostics_Error(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,double time, double p);

//////////// Diagnostics Euler ///////

double NormLP_Velocity(Data & d, Mesh & Mh,variable & v,ParamPhysic & Param,double p,double t);

#endif
