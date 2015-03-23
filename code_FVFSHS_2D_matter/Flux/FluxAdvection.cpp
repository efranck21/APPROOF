#include <iostream>
#include <cmath>
#include <cassert>
#include "Flux.hpp"
#include "Variables.hpp"
#include "Functions.hpp"
#include <stdio.h>
#include "Tensor.hpp"
#include "BoundaryCondition.hpp"
#include "Functionsinit.hpp"
#include "FunctionsAdvection.hpp"


vectorflux FluxVertexadv(Data & d,int numCell,Mesh &  Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param){
  /**  Nodal advection fluxes  **/
 vectorflux res(1);
 R2 a(0,0); 
  int numGr;
  double s=0;
  
  for(int r=0;r<Mh.nbnodelocal;r++){
    numGr=Mh(numCell,r);
    a=AverageVelocityNode(d,Mh,v,tab,Param.Adv,numGr,numCell);
    s=s+VertexUpwind(d,Mh,v,0,tab,Param,numCell,r,a);

  }
  
  res.vflux[0]=-s;
  return res;	


}


