
#include <iostream>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include "Functions.hpp"
#include "AnalyticSolutions.hpp"
#include "ParamPhys.hpp"


double SolFonM1(Data & d, ParamPhysic & Param,Vertex p,double t,Vertex center,int var){
 /** function which gives the M1 analytical solution for p, t, center and the variable "var"  **/  
  double res=0;
  
  if(d.nTest==1){
    res=CasTestM1Diff(d,Param,p,t,center,var);
    return res;
  }
  if(d.nTest==2){
    res=CasTestM1DiffPeriodic(d,Param,p,t,center,var);
    return res;
  }
   if(d.nTest==3){
    res=CasTestM1TransportSmooth(d,Param,p,t,center,var);
    return res;
  }
  if(d.nTest==4){
    res=CasTestM1Transport(d,Param,p,t,center,var);
    return res;
  }
   

    
    else {
      cout << "Pas de solution fondamentale associé a ce cas test"<<endl;exit(1);
    }
}

double CasTestM1Diff(Data & d, ParamPhysic & Param,Vertex p, double t,Vertex center, int var ){
   /** Analytical solution: M1 solution in the diffusion regime for eps<<1**/
  double res;
 
  if(var == 0) {
     res=SolFonChaleur(d,Param,p,t,center);
   }
   if(var == 1) {
     res=0.0000;
   }
   if(var == 2) {
     res=0.0000;
   }
  
  
  return res;
}

double CasTestM1DiffPeriodic(Data & d, ParamPhysic & Param,Vertex p, double t,Vertex center,int var){
    /** Analytical solution: periodic M1 solution in the diffusion regime for eps<<1**/
  double res;
  double co=0;

  co=Param.M1.sig_value;
   if(var == 0) {
     res=SolNeumannChaleur(d,Param,p,t)+2;
   }
   if(var == 1) {
     res=0.0000;
   }
   if(var == 2) {
     res=0.0000;
   }
  
  return res;
}

double CasTestM1TransportSmooth(Data & d, ParamPhysic & Param,Vertex p, double t,Vertex center,int var){
    /** Analytical solution: Free transport smooth solution for M1 model for sigma=0 and eps=1**/
  double res=0;
  double co=0;

  co=Param.M1.sig_value;
   if(var == 0) {
     res=TransportSol(d,Param,p,t,center)+0.0001;
   }
   if(var == 1) {
     res=TransportSol(d,Param,p,t,center)+0.0001;
   }
   if(var == 2) {
     res=0.0000;
   }
  
  return res;
}

double CasTestM1Transport(Data & d, ParamPhysic & Param,Vertex p, double t,Vertex center,int var){
   /** Analytical solution: Free transport discontinuous solution for M1 model for sigma=0 and eps=1**/
  double res;
  double cc=0;
  cc=0.0001;
  if(p.x-t>=(d.Tx/2-0.1) && p.x-t<=(d.Tx/2+0.1)){
     if(p.y>=(d.Ty/2-0.1) && p.y<=(d.Ty/2+0.1)){
	  cc=1.;
     }    
  }
   if(var == 0) {
     res = cc;
   }
   if(var == 1) {
     res = 0.0000 ;
   }
   if(var == 2) {
     res = 0.0000 ;
   }
  
  return res;
}
