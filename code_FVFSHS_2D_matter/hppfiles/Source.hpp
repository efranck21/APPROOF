#ifndef SOURCE_HPP
#define SOURCE_HPP
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

//////////////// General functions ////////////

vectorflux ChoiceSource(Data & d,int numCell,Mesh & Mh,variable & v,TabConnecInv & tab,ParamPhysic & Param, R2 * ur);

R2 SemiImplicitPart(Data & d ,int numCell,Mesh & Mh, variable & v,variable & temp, TabConnecInv & tab,ParamPhysic & Param,double dt);

//////////// Advection functions /////////

vectorflux Sourceadv(Data & d,int numCell,Mesh & Mh, variable & v,ParamPhysic & Param);


//////////// Diffusion functions /////////

vectorflux SourceDiff(Data & d,int numCell,Mesh & Mh,variable & v,ParamPhysic & Param);




//////////// P1 explicit functions /////////

vectorflux SourceP1(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param);

//////////// P1 semi-implicit functions /////////

void MatrixSourceP1(Data & d,int numCell,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,double  Mrsource[2][2]);

R2 SolveurImpliciteP1(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,double dt,double b1,double b2);



//////////// P1 explicit functions /////////

vectorflux SourceP1Matter(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param);

//////////// P1 semi-implicit functions /////////

void MatrixSourceP1Matter(Data & d,int numCell,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,double  Mrsource[2][2]);

R2 SolveurImpliciteP1Matter(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,double dt,double b1,double b2);

R2 SplittingCalculMatter(Data & d ,ParamPhysic & Param, int numCell,Mesh &  Mh,variable & vnp,variable  & vn, TabConnecInv & tab,double dt,double dT);


//////////// Euler explicit functions /////////

double FrictionUcarre(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 * ur);

R2 FrictionU(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 * ur);

vectorflux SourceE(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 * ur);

R2 SourceGravityEuler(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param);

double SourceGravityEulerU(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 * ur);




//////////// Euler semi implicit functions /////////

R2 SolveurImpliciteE(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,double dt,double b1,double b2, double rho);

void MatrixSourceWBE(Data & d,int numGr,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,double  Mrsource[2][2]);



/////////// M1 functions ////////

R2 SolveurImpliciteM1(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param, double dt,double b1, double b2);

void MatrixSourceM1(Data & d,int numCell,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param, double  Mrsource[2][2]);


vectorflux SourceM1(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param);



/////////// M1 with Matter functions ////////


vectorflux SourceM1Matter(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param);








#endif
