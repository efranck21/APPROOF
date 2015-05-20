#ifndef FLUX_HPP
#define FLUX_HPP
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Variables.hpp"
#include "ClassMesh.hpp"
#include "Functions.hpp"
#include "Initialisation.hpp"
#include "Functionsinit.hpp"
#include "AnalyticSolutions.hpp"
#include "Diagnostics.hpp"
#include "ParamPhys.hpp"



vectorflux ChoiceFluxE(Data & d ,int numCell,Mesh &  Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param);

vectorflux FluxVF4Diff(Data & d,int numCell, Mesh & Mh,  variable & v,  TabConnecInv & tab,ParamPhysic & Param);

vectorflux FluxEdgeAdv(Data & d,int numCell,Mesh &  Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param);

vectorflux ChoiceFluxN(Data & d ,int numCell,Mesh & Mh, variable & v,TabConnecInv & tab,ParamPhysic & Param,R2 ** ur);

R2 SolveurNodal(Data & d,int numGr,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param, int group);


///////// flux d'advection/////////

vectorflux FluxVertexadv(Data & d,int numCell, Mesh &  Mh, variable & v,TabConnecInv & tab,ParamPhysic & Param);


/////////////// flux de diffusion/////////


vectorflux FluxVertexDiff(Data & d,int numCell,Mesh &  Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 ** ur);

vectorflux FluxVertexNlDiff(Data & d,int numCell,Mesh &  Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 ** ur);

double FluxadvectionDiff(Data & d,int numCell,int r,Mesh &  Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 a,R2** ur);

void MatrixDiff(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr,double &a11,double &a12,double &a21,double &a22,double &b1,double &b2); 

////////////// flux du modèle P1//////////////////


vectorflux FluxVertexP1(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 ** ur);

vectorflux FluxVertexP1Gosse(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 ** ur);

void MrConstructP1(Data & d,int numGr,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,double  Mr[2][2]);

void MatrixP1(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr,double &a11,double &a12,double &a21,double &a22,double &b1,double &b2); 


////////////// flux du modèle P1 Matter//////////////////


vectorflux FluxVertexP1Matter(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 ** ur);

vectorflux FluxVertexP1MatterGosse(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 ** ur);

void MrConstructP1Matter(Data & d,int numGr,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,double  Mr[2][2]);

void MatrixP1Matter(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr,double &a11,double &a12,double &a21,double &a22,double &b1,double &b2); 


/////////////// Euler  Flux ////////

vectorflux FluxVertexEulerGosse(Data & d,int numCell,Mesh & Mh,variable & v,TabConnecInv & tab,ParamPhysic & Param,R2 ** ur);

vectorflux FluxVertexEuler(Data & d,int numCell,Mesh & Mh,variable & v,TabConnecInv & tab,ParamPhysic & Param,R2 ** ur);

void MrConstructE(Data & d,int numGr,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,double  Mr[2][2]);

R2 ClassicalSolver(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr);

void MatrixEuler(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr,double &a11,double &a12,double &a21,double &a22,double &b1,double &b2); 


/////////////// Euler functions //////

double PressureLaw(Data & d,Mesh & Mh, variable & v,ParamEuler & Euler, int numCell);

double EnergyIn(Data & d,Mesh & Mh, variable & v,ParamEuler & Euler, int numCell);

double rhor(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamEuler & Euler,int numGr);

double c(Data & d,Mesh & Mh, variable & v,ParamEuler & Euler, int numCell);

double WaveSpeed(Data & d,Mesh & Mh, variable & v, int numCell,int r,TabConnecInv & tab,ParamEuler & Euler);

double remap(Data & d,Mesh & Mh, variable & v,int var,TabConnecInv & tab,ParamPhysic & Param,int numCell,int r, R2 ** ur, R2 a);

R2 IntGradiantPressureInterpolation(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamEuler & Euler,int numGr);

double IntDensityInterpolation(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamEuler & Euler,int numGr);

R2 RhoGravityVector(Data & d,Mesh &Mh, variable & v, TabConnecInv & tab,ParamEuler & Euler,int numGr);

R2 GravityVector(Data & d,ParamEuler & Euler);

void InitTab_rhoGravity(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamEuler & Euler,int numGr);

#endif
 



//////////////// M1 functions ////////////
double q(Data & d, ParamM1 & M1,double E,double F1,double F2);

double coefk(Data & d, ParamM1 & M1,double E,double F1,double F2);

void Calcul_u(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab, ParamM1 & M1,R2 & u,int numCell);

void MatrixM1(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr,double &a11,double &a12,double &a21,double &a22,double &b1,double &b2);

void MatrixAPM1(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr,double &a11,double &a12,double &a21,double &a22,double &b1,double &b2);

void MatrixClassicM1(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr,double &a11,double &a12,double &a21,double &a22,double &b1,double &b2); 

vectorflux FluxVertexAPM1(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab, ParamPhysic & Param,R2 ** ur);

vectorflux FluxVertexClassicM1(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab, ParamPhysic & Param,R2 ** ur);

vectorflux FluxVertexM1Gosse(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab, ParamPhysic & Param,R2 ** ur);

double remap(Data & d,Mesh & Mh, variable & v,int var,TabConnecInv & tab, ParamPhysic & Param,int numCell,int r,R2 ** ur, R2 a);

void MrConstructM1(Data & d,int numGr,Mesh & Mh,variable & v, TabConnecInv & tab, ParamPhysic & Param,double Mr[2][2]);



//////////////// M1 Matter functions ////////////

void MatrixM1Matter(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr,double &a11,double &a12,double &a21,double &a22,double &b1,double &b2);

vectorflux FluxVertexClassicM1Matter(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab, ParamPhysic & Param,R2 ** ur);


double q_M1M(Data & d, ParamM1Matter & M1Matter,double E,double F1, double F2);

void  Calcul_u_M1M(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab, ParamM1Matter & M1Matter,R2 & u,int numCell);

double coefk_M1M(Data & d, ParamM1Matter & M1Matter,double E,double F1, double F2);



////////////// flux du modèle P1 Compton//////////////////


vectorflux FluxVertexP1Compton(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 ** ur);

vectorflux FluxVertexP1ComptonGosse(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 ** ur);

void MrConstructP1Compton(Data & d,int numGr,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,double  Mr[2][2]);

void MatrixP1Compton(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr,double &a11,double &a12,double &a21,double &a22,double &b1,double &b2,int group); 
