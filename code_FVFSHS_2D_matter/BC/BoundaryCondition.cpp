#include <iostream>
#include <cmath>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include <cstdlib>
#include "Variables.hpp"
#include "Data.hpp"
#include "BoundaryCondition.hpp"
#include "Functions.hpp"
#include "Tensor.hpp"
#include "Functionsgeo.hpp"
#include "AnalyticSolutions.hpp"
#include "Diagnostics.hpp"
#include "Functionsinit.hpp"
#include "ParamPhys.hpp"

 /** File which contains the function for the boundary conditions common to the models**/


R2 Norext(Data & d,Mesh & Mh, int numGr){
   /** Function which compute the external normal**/
  R2 Next(0,0);
  if(Mh.xr(numGr).lab==2) {
    Next.x=-1;
    Next.y=0; 
  }
  if(Mh.xr(numGr).lab==3) {
    Next.x=0;
    Next.y=1;
  }
  if(Mh.xr(numGr).lab==4) {
    Next.x=1;
    Next.y=0;   
  }
  if(Mh.xr(numGr).lab==5) {
    Next.x=0;
    Next.y=-1;
  } 
  if(Mh.xr(numGr).x==0 && Mh.xr(numGr).y==0) {
    Next.x=-1;
    Next.y=-1;
  }
  if(Mh.xr(numGr).x==d.Tx && Mh.xr(numGr).y==0) {
    Next.x=1;
    Next.y=-1;
  }
  if(Mh.xr(numGr).x==d.Tx && Mh.xr(numGr).y==d.Ty) {
    Next.x=1;
    Next.y=1;
  }
  if(Mh.xr(numGr).x==0 && Mh.xr(numGr).y==d.Ty) {
    Next.x=-1;
    Next.y=1;
  }
  return Next;
}


void BoundaryCondition(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,double time){
 /** Function which initialize the phantom cells for the model given by Data and Param **/
  if (Param.Model == 1) {
    BoundaryConditionDiff(d,Mh,v,Param,time);
  }

  if (Param.Model == 2) {
    BoundaryConditionAdvection(d,Mh,v,Param,time);
  }

   if (Param.Model == 3) {
     BoundaryConditionP1(d,Mh,v,Param,time);
  }

   if (Param.Model == 4) {
     BoundaryConditionP1Matter(d,Mh,v,Param,time);
  }
  
   if (Param.Model == 5) {
     BoundaryConditionEuler(d,Mh,v,Param,time);
  }

  if (Param.Model == 6) {
    BoundaryConditionM1(d,Mh,v,Param,time);
  }
  if (Param.Model == 7) {
    BoundaryConditionM1Matter(d,Mh,v,Param,time);
  }
   
}
  





