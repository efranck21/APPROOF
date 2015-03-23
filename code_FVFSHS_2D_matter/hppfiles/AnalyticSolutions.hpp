#ifndef ANALYTICSOLUTIONS_HPP
#define ANALYTICSOLUTIONS_HPP

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

///////// general fucntions //////

double SolFon(Data & d,ParamPhysic & Param,Vertex p,double t,Vertex center,int var);

double FinalTime(Data & d,ParamPhysic & Param);

int TraceSolfond(Data & d,ParamPhysic & Param);

double Heaviside(Vertex y);

////////// Euler functions ///////

double SolFonEuler(Data & d,ParamPhysic & Param,Vertex p,double t,Vertex center,int var);

double SolFonEulerAverage(Data & d,Mesh & Mh, ParamPhysic & Param,Vertex p,double t,Vertex center,int numCell,int var);


////////// Euler functions ///////

double SolFonM1(Data & d,ParamPhysic & Param,Vertex p,double t,Vertex center,int var);

double CasTestM1Diff(Data & d, ParamPhysic & Param,Vertex p, double t,Vertex center, int var );

double CasTestM1DiffPeriodic(Data & d, ParamPhysic & Param,Vertex p, double t,Vertex center,int var);

double CasTestM1TransportSmooth(Data & d, ParamPhysic & Param,Vertex p, double t,Vertex center,int var);

double CasTestM1Transport(Data & d, ParamPhysic & Param,Vertex p, double t,Vertex center,int var);


///////// Diffusion functions //////

double SolFonDiff(Data & d,ParamPhysic & Param,Vertex p,double t,Vertex center);

double SolFonChaleur(Data & d,ParamPhysic & Param,Vertex p,double t,Vertex center);

double SolNeumannChaleur(Data & d,ParamPhysic & Param,Vertex p,double t);


//////// Advection functions ////////

double SolFonAdvection(Data & d,ParamPhysic & Param,Vertex p,double t,Vertex center);

double TransportSol(Data & d,ParamPhysic & Param,Vertex p,double t,Vertex center);



/////// P1 functions /////////

double SolFonP1(Data & d,ParamPhysic & Param,Vertex p,double t,Vertex center,int var);

double CasTestTelegraph(Data & d,ParamPhysic & Param,Vertex p,double t,int var);

double CasTestOnde(Data & d,ParamPhysic & Param,Vertex p,double t,int var);

double UniformSolution(Data & d,ParamPhysic & Param,Vertex p,double t,int var);

#endif
