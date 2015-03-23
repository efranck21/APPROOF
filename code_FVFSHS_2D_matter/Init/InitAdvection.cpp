#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include "Initialisation.hpp"
#include <math.h>
#include "AnalyticSolutions.hpp"
#include "Diagnostics.hpp"
#include "Functions.hpp"
#include "Flux.hpp"
#include "ParamPhys.hpp"

void Initadv(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param){
    /** function which construct the initial datas 
      for the advection model
  **/
  int in=0;

  switch(d.nTest)
    {
    case 1 :
      
      in=centerdomain(d,Mh);
      for(int j=0;j<Mh.nc;j++){
	v.var[0][j]=TransportSol(d,Param,Mh.xj(j),0,Mh.xj(in));
      }
      break;
      
    }

}
