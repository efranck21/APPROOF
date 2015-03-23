
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

vectorflux FluxEdgeAdv(Data & d,int numCell,Mesh &  Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param){
   /**  Edge advection fluxes  **/
 vectorflux res(1);
 R2 a(0.,0.); 
 int k=0;
 int r1=0,r2=0,numGr1=0,numGr2=0;
 int numGr=0;
 double s=0.;
  
  for(int r=0;r<Mh.nbnodelocal;r++){
    numGr=Mh(numCell,r); 
    r1=r;
    if(r==Mh.nbnodelocal-1){
      r2=0;
    }
    else{
      r2=r+1;
    }
    numGr1=Mh(numCell,r1);
    numGr2=Mh(numCell,r2);      

    k=InverseEdge(Mh,numGr1,numGr2,numCell,tab);
    a.x=0.5*(Param.Adv.a[numCell].x+Param.Adv.a[k].x);
    a.y=0.5*(Param.Adv.a[numCell].y+Param.Adv.a[k].y);
    s=s+EdgeUpwind(d,Mh,v,0,tab,Param,numCell,k,r,a);
  }
  
  res.vflux[0]=-s;
  return res;	


}
