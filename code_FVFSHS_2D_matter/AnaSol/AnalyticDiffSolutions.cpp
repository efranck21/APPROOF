#include <iostream>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include "Functions.hpp"
#include "AnalyticSolutions.hpp"
#include "ParamPhys.hpp"


double SolFonDiff(Data & d, ParamPhysic & Param,Vertex p,double t,Vertex center){
   /** fucntion which gives the diffusion analytical solution for p, t and center  **/
  double res;

     
  if(d.nTest==1){   
    res=SolFonChaleur(d,Param,p,t,center);
    return res;
  }
  
  if(d.nTest==2){   
    res=SolNeumannChaleur(d,Param,p,t);
    return res;
  }
    
    else {
      cout << "Pas de solution fondamentale associé a ce cas test"<<endl;exit(1);
    }
}




double SolFonChaleur(Data & d, ParamPhysic & Param,Vertex p,double t,Vertex center){
   /** Analytical solution: fundamental solution of the heat equation for different diffusion coefficient
 and with initial time =0.001  **/
  double res=0;
  double cD=0;
  double t0=0.01;
  if(Param.Model == 1 && Param.Diff.Coefdiff_value !=-1) {
    cD=Param.Diff.Coefdiff_value;
  }
   if(Param.Model == 3 && Param.P1.sig_value !=-1) {
     cD=1./Param.P1.sig_value;
  }

   if(Param.Model ==6 && Param.M1.sig_value !=-1) {
     cD=1./(3.*Param.M1.sig_value);
   }
   res=(1./(4*cD*M_PI*(t+t0)))*exp(-((p-center).norme2()/(4*cD*(t+t0))))+1.;
  return res;
}


double SolNeumannChaleur(Data & d, ParamPhysic & Param,Vertex p,double t){
     /** Analytical solution: periodic solution of the heat equation for different diffusion coefficient**/
  double res=0;
   double cD=0;

  if(Param.Model == 1 && Param.Diff.Coefdiff_value !=-1) {
    cD=Param.Diff.Coefdiff_value;
  }
   if(Param.Model == 1 && Param.P1.sig_value !=-1) {
    cD=1./Param.P1.sig_value;
  }
  if(Param.Model == 1 && Param.P1.sig_value !=-1) {
     cD=1./(3.*Param.M1.sig_value);
  }

  res=exp(-cD*(2*M_PI*M_PI*t))*cos(M_PI*p.x)*cos(M_PI*p.y);
 
  
  return res;
  }
