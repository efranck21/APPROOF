
#include <iostream>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include "Functions.hpp"
#include "AnalyticSolutions.hpp"
#include "ParamPhys.hpp"

double SolFonAdvection(Data & d,ParamPhysic & Param,Vertex p,double t,Vertex center){
  /** fucntion which gives the advection analytical solution for p, t and center  **/ 
  double res;


  if(d.nTest==1){   
    res=TransportSol(d,Param,p,t,center);
    return res;
  }

    
    else {
      cout << "Pas de solution fondamentale associé a ce cas test"<<endl;exit(1);
    }
}

double TransportSol(Data & d,ParamPhysic & Param,Vertex p,double t,Vertex center){
   /** Analytical solution: advection of gaussian function  **/ 
  R2 pp(0,0);
  double res=0;
  R2 a(0,0);
  a.x=Param.Adv.velocity.x;
  a.y=Param.Adv.velocity.y;
  if(!strcmp(d.Typemodel,"M1")){ a.x=1.; a.y=0;}
  pp.x=p.x-a.x*t;
  pp.y=p.y-a.y*t;
  res=(1/(4*(1.)*M_PI*0.01))*exp(-((pp-center).norme2()/(4*(1.)*0.01)));

 
  return res;
		  
}
