
#include <iostream>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include "Functions.hpp"
#include "AnalyticSolutions.hpp"
#include "ParamPhys.hpp"

/** This file contains the functions generic for all the models associated with the analytical solutions **/




double SolFon(Data & d,ParamPhysic & Param,Vertex p,double t,Vertex center,int var){
 /** The function Solfond five the analytical solution for a vertex "p", 
      the time "t", the variable "var" with a center of the scheme given by "center **/
  double res=0;
  
  
 if(Param.Model == 1){   
   res=SolFonDiff(d,Param,p,t,center);

  }

  if(Param.Model == 2){   
    res=SolFonAdvection(d,Param,p,t,center);
  }

  if(Param.Model == 3){   
    res=SolFonP1(d,Param,p,t,center,var);
  }
  
  if(Param.Model == 5){   
    res=SolFonEuler(d,Param,p,t,center,var);
  }
  if(Param.Model == 6){   
    res=SolFonM1(d,Param,p,t,center,var);
  }
  return res;
}




int TraceSolfond(Data & d, ParamPhysic & Param){
  /** This fucntion return 1 if we have a analytical solution 
and 0 if there are not analytical solution **/
  int p=0;

  if(Param.Model == 1){
    p=1;
  }
  if(Param.Model == 2){
    p=1;
  }
  if(Param.Model == 3 && d.nTest < 5){
    p=1;
  }
  if(Param.Model == 5 && d.nTest < 8){
    p=1;
  }
  if(Param.Model == 6){
    p=1;
  }
  	

 return p;
}
