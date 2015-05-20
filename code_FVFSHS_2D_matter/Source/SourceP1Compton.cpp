#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Variables.hpp"
#include "Flux.hpp"
#include "Source.hpp"
#include "BoundaryCondition.hpp"
#include "Tensor.hpp"
#include "FunctionsAdvection.hpp"

vectorflux SourceP1Compton(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param) {
  /** Function wich compute the explicit source terms for P1 with matter**/
  vectorflux res(d);
  R2 sol(0,0);

   for(int g=0;g<d.ngroup;g++){
 
   res.vflux[0+g*Param.P1C.nb_moment]=0;
   res.vflux[1+g*Param.P1C.nb_moment]=0;
   res.vflux[2+g*Param.P1C.nb_moment]=0;
   
   }
   
   res.vflux[Param.P1C.nbvar-1]=0;
  
   return res;
}
