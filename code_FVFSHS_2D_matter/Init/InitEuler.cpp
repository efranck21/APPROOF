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

 void InitE(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param){
 /** function which construct the initial datas 
      for the Euler model
  **/
   int in=0;
   double lambda=0.001;
   double norm2=0;
   R2 g(0,0);
  g=GravityVector(d,Param.Euler);

     switch(d.nTest)
    {
    case 1 : 		
      //condition initiale
      for(int j=0;j<Mh.nc;j++){
        if(Param.Euler.LHO <2) {
	  v.var[0][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),0);
	  v.var[1][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),1);
	  v.var[2][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),2);
	  v.var[3][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),3);
        }
	else {
	  v.var[0][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),0,Mh.xj(in),j,0);
          v.var[1][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),0,Mh.xj(in),j,1);
	  v.var[2][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),0,Mh.xj(in),j,2);
	  v.var[3][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),0,Mh.xj(in),j,3);
        }

      }
    
      break;
 
    case 2 : 
      for(int j=0;j<Mh.nc;j++){
	v.var[0][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),0);
        v.var[1][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),1)+lambda*cos(M_PI*Mh.xj(j).x)*cos(M_PI*Mh.xj(j).y);
	v.var[2][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),2)+lambda*cos(M_PI*Mh.xj(j).x)*cos(M_PI*Mh.xj(j).y);
	v.var[3][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),3);
	
      }     
      break;
   
      case 3 : 
       for(int j=0;j<Mh.nc;j++){
       if(Param.Euler.LHO <2) {
	  v.var[0][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),0);
	  v.var[1][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),1);
	  v.var[2][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),2);
	  v.var[3][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),3);
        }
	else {
	  v.var[0][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),0,Mh.xj(in),j,0);
          v.var[1][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),0,Mh.xj(in),j,1);
	  v.var[2][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),0,Mh.xj(in),j,2);
	  v.var[3][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),0,Mh.xj(in),j,3);
        }
      }
      break;
 
    case 4 : 
      for(int j=0;j<Mh.nc;j++){
	v.var[0][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),0);
        v.var[1][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),1)+lambda*sin(M_PI*Mh.xj(j).x)*sin(M_PI*Mh.xj(j).y);
	v.var[2][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),2)+lambda*sin(M_PI*Mh.xj(j).x)*sin(M_PI*Mh.xj(j).y);
	v.var[3][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),3);
	
      }     
      break;


      case 7 : 

       in=centerdomain(d,Mh);
      for(int j=0;j<Mh.nc;j++){
	v.var[0][j]= SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),0);
	v.var[1][j]= SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),1);
	v.var[2][j]= SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),2);
	v.var[3][j]= SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),3);

      }

     case 8 : 	
      //condition initiale

      for(int j=0;j<Mh.nc;j++){
	norm2=(Mh.xj(j),Mh.xj(j));
	if(sqrt(norm2)<0.5){
	  v.var[0][j]=1.;
	  v.var[3][j]=1./(Param.Euler.gamma-1.);	
	  v.var[1][j]=0;
	  v.var[2][j]=0;
	}
	else{
	  v.var[0][j]=0.125;
	  v.var[3][j]=0.1/(Param.Euler.gamma-1.);	
	  v.var[1][j]=0;
	  v.var[2][j]=0;
	}
      }
      break;

     case 5 : 
		
      //condition initiale 
       in=centerdomain(d,Mh);
      for(int j=0;j<Mh.nc;j++){
       if(Param.Euler.LHO <2) {
	  v.var[0][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),0);
	  v.var[1][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),1);
	  v.var[2][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),2);
	  v.var[3][j]=SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),3);
        }
	else {
	  v.var[0][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),0,Mh.xj(in),j,0);
          v.var[1][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),0,Mh.xj(in),j,1);
	  v.var[2][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),0,Mh.xj(in),j,2);
	  v.var[3][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),0,Mh.xj(in),j,3);
        }
       }
      break;

       case 6 : 
		
      //condition initiale
       in=centerdomain(d,Mh);
      for(int j=0;j<Mh.nc;j++){
	v.var[0][j]= SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),0);
	v.var[1][j]= SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),1)+lambda*sin(M_PI*Mh.xj(j).x)*sin(M_PI*Mh.xj(j).y);
	v.var[2][j]= SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),2)+lambda*sin(M_PI*Mh.xj(j).x)*sin(M_PI*Mh.xj(j).y);
	v.var[3][j]= SolFon(d,Param,Mh.xj(j),0,Mh.xj(in),3);

      }
      break;

      case 9 : 	
      //condition initiale

      for(int j=0;j<Mh.nc;j++){
	if(Mh.xj(j).y<(0.5*d.Ty-0.05*cos(2*M_PI*Mh.xj(j).x/d.Tx)*d.Ty)){
	  v.var[0][j]=1.;
	  v.var[1][j]=0.;
	  v.var[2][j]=0.;
	  v.var[3][j]=(1./(Param.Euler.gamma-1.))*(2*(g.y*d.Ty+g.y*d.Tx+1)-(g.y*Mh.xj(j).y+g.x*Mh.xj(j).x)*v.var[0][j]);	
	}
	else{
	  v.var[0][j]=2.;
	  v.var[3][j]=(1./(Param.Euler.gamma-1.))*(2*(g.y*d.Ty+g.y*d.Tx+1)-(g.y*Mh.xj(j).y+g.x*Mh.xj(j).x)*v.var[0][j]);	
	  v.var[1][j]=0;
	  v.var[2][j]=0;
	}
	v.var[2][j]=0;//v.var[0][j]*lambda*0.25*(1+cos(2*M_PI*Mh.xj(j).x/0.5))*(1+cos(2*M_PI*Mh.xj(j).y/1.5));
}	
      break;

    }


 }





