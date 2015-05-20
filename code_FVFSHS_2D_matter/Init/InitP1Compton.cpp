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


void InitP1Compton(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param){
   /** function which construct the initial datas 
      for the p1 model with matter
  **/
  double T=0;
  
  switch(d.nTest)
    {
      
    case 1 :
      
      for(int j=0;j<Mh.nc;j++){
	for(int k=0; k<Param.P1C.nb_group;k++){
	  v.var[0+k*Param.P1C.nb_moment][j]=1;
	  v.var[1+k*Param.P1C.nb_moment][j]=0;	
	  v.var[2+k*Param.P1C.nb_moment][j]=0;
	   }
	v.var[Param.P1C.nbvar-1][j]=1;
     
      }
      break;
      
 
      
    }
}
