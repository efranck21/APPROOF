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

void InitDiff(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param){
   /** function which construct the initial datas 
      for the diffusion model
  **/
  int in=0;

  switch(d.nTest)
    {
    case 1 : 
   	       
      in=centerdomain(d,Mh);
      for(int j=0;j<Mh.nc;j++){
	v.var[0][j]=0;
      }
      v.var[0][in]=1/(Mh.area(in));
      break;
      
      case 2 :

	 in=centerdomain(d,Mh);

	 for(int j=0;j<Mh.nc;j++){
	   v.var[0][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),0);
	 }
      break;

      case 3 :
	 
	 for(int j=0;j<Mh.nc;j++){
	   v.var[0][j]=cos(M_PI*Mh.xj(j).y)*cos(M_PI*Mh.xj(j).x);
	 }
	 break;
    }
}
  
