
#include <iostream>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include "Functions.hpp"
#include "AnalyticSolutions.hpp"
#include "ParamPhys.hpp"


double SolFonP1(Data & d, ParamPhysic & Param,Vertex p,double t,Vertex center,int var){
    /** function which gives the P1 analytical solution for p, t and center  **/ 
  double res=0;
  
   if(d.nTest==2 && Param.P1.sig_value == 1){
     res=CasTestTelegraph(d,Param,p,t,var);
     return res;
   }

   
    if(d.nTest==1 && Param.P1.sig_value == 0){
      res=CasTestOnde(d,Param,p,t,var);
     return res;
   }
  
   if(d.nTest==3){
      if(var == 0) {
        res=SolFonChaleur(d,Param,p,t,center);
     }
     if(var == 1) {
       res=0.; 
      }
     if(var == 2) {
       res=0.; 
     }
     return res;
   }

   if(d.nTest==4){
     res=UniformSolution(d,Param,p,t,var);
     return res;
   }


    
    else {
      cout << "Pas de solution fondamentale associé a ce cas test"<<endl;exit(1);
    }
}

double CasTestOnde(Data & d, ParamPhysic & Param,Vertex p, double t,int var){
      /** Analytical solution: periodic solution of the P1 equation for sigma=0 and eps=1**/
  double res;
   if(var == 0) {
     res=-M_PI*sqrt(2)*sin(M_PI*sqrt(2)*t)*cos(M_PI*p.x)*cos(M_PI*p.y);
     }
     if(var == 1) {
       res=M_PI*cos(M_PI*sqrt(2)*t)*sin(M_PI*p.x)*cos(M_PI*p.y);
       }
     if(var == 2) {
       res=M_PI*cos(M_PI*sqrt(2)*t)*cos(M_PI*p.x)*sin(M_PI*p.y); 
     }
  
  return res;
}


double CasTestTelegraph(Data & d, ParamPhysic & Param,Vertex p, double t,int var){
      /** Analytical solution: periodic solution of the P1 equation for sigma=1 and eps=1**/
  double a,res1, res2, res;
  a=sqrt(8.*pow(M_PI,2.)-1.);
  if(var == 0) {
    res1=-(1/2.)*exp(-(t/2.))*(cos((a/2.)*t)+(1./a)*sin((a/2.)*t));
    res2=-exp(-(t/2.))*((a/2.)*sin((a/2.)*t)-(1./2.)*cos((a/2.)*t));
    res=(res1+res2)*cos(M_PI*p.x)*cos(M_PI*p.y);
  }
  if(var == 1) {
    res=M_PI*exp(-(t/2.))*(cos((a/2.)*t)+(1./a)*sin((a/2.)*t))*sin(M_PI*p.x)*cos(M_PI*p.y);
  }
  if(var == 2) {
    res=M_PI*exp(-(t/2.))*(cos((a/2.)*t)+(1./a)*sin((a/2.)*t))*cos(M_PI*p.x)*sin(M_PI*p.y); 
  }
  
  return res;
}


double UniformSolution(Data & d, ParamPhysic & Param,Vertex p, double t,int var){
/** Analytical solution: periodic solution of the P1 equation for  eps<1/sqrt(8 pi)**/
  
  double lambda_1,lambda_2,alpha,alpha_prime,res;
  lambda_1=-Param.P1.sig_value*(sqrt(1.-(pow(Param.P1.eps_value,2)/pow(Param.P1.sig_value,2))*8.*pow(M_PI,2))+1);
  lambda_1=lambda_1/(2*pow(Param.P1.eps_value,2.));
  lambda_2=Param.P1.sig_value*(sqrt(1-(pow(Param.P1.eps_value,2.)/pow(Param.P1.sig_value,2))*8.*pow(M_PI,2))-1);
  lambda_2=lambda_2/(2*pow(Param.P1.eps_value,2.));

  
  alpha= (lambda_2/(lambda_2-lambda_1))*exp(lambda_1*t)-(lambda_1/(lambda_2-lambda_1))*exp(lambda_2*t);
  alpha_prime=((lambda_2*lambda_1)/(lambda_2-lambda_1))*exp(lambda_1*t)-((lambda_1*lambda_2)/(lambda_2-lambda_1))*exp(lambda_2*t);
  
  if(var == 0) { 
    res=(alpha+(pow(Param.P1.eps_value,2.)/pow(Param.P1.sig_value,2.))*alpha_prime)*cos(M_PI*p.x)*cos(M_PI*p.y);
  }
  if(var == 1) { 
    res=(Param.P1.eps_value/Param.P1.sig_value)*alpha*M_PI*sin(M_PI*p.x)*cos(M_PI*p.y);
  }
  if(var == 2) { 
    res=(Param.P1.eps_value/Param.P1.sig_value)*alpha*M_PI*cos(M_PI*p.x)*sin(M_PI*p.y);
  }
  return res;
}


