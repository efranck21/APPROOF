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

  void InitM1(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param){
     /** function which construct the initial datas 
      for the M1 model
  **/
    int in=0;
    in=centerdomain(d,Mh);
    
       switch(d.nTest)
       { 

      case 1 : 	
      in=centerdomain(d,Mh);
      for(int j=0;j<Mh.nc;j++){
	v.var[0][j]=SolFonM1(d,Param,Mh.xj(j),0,Mh.xj(in),0);
	v.var[1][j]=SolFonM1(d,Param,Mh.xj(j),0,Mh.xj(in),1);
	v.var[2][j]=SolFonM1(d,Param,Mh.xj(j),0,Mh.xj(in),2);
      }

      break;
     
       case 2: 

      for(int j=0;j<Mh.nc;j++){
	v.var[0][j]=SolFonM1(d,Param,Mh.xj(j),0,Mh.xj(in),0);
	v.var[1][j]=SolFonM1(d,Param,Mh.xj(j),0,Mh.xj(in),1);
	v.var[2][j]=SolFonM1(d,Param,Mh.xj(j),0,Mh.xj(in),2);
      }
            

      break;
	case 3 :
	 in=centerdomain(d,Mh);
	 for(int j=0;j<Mh.nc;j++){
	   v.var[0][j]=SolFonM1(d,Param,Mh.xj(j),0,Mh.xj(in),0);
	   v.var[1][j]=SolFonM1(d,Param,Mh.xj(j),0,Mh.xj(in),1);
	   v.var[2][j]=SolFonM1(d,Param,Mh.xj(j),0,Mh.xj(in),2);
	   
	 }
	
	 break; 
       case 4:	 
	 for(int j=0;j<Mh.nc;j++){
	   
	   v.var[0][j]=SolFonM1(d,Param,Mh.xj(j),0,Mh.xj(in),0);
	   v.var[1][j]=SolFonM1(d,Param,Mh.xj(j),0,Mh.xj(in),1);
	   v.var[2][j]=SolFonM1(d,Param,Mh.xj(j),0,Mh.xj(in),2);
	 }
	 break;
    
  }
}
