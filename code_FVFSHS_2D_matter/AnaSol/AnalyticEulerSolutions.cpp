
#include <iostream>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include "Functions.hpp"
#include "AnalyticSolutions.hpp"
#include "ParamPhys.hpp"
#include "HOfunctions.hpp"

/////////////////////////Fonctions liées aux solutions fondamental/////////////



double SolFonEuler(Data & d, ParamPhysic & Param,Vertex p,double t,Vertex center,int var){
    /** function which gives the Euler analytical solution for p, t, center and the variable "var"  **/ 
  double res=0;
  R2 g(0,0);
  g=GravityVector(d,Param.Euler);
  R2 u(2.,0.);
  double rho1=0,rho2=0,X=0,Y=0,p0=0,u0=0,v0=0,r=0,rho=0;
 
  if(d.nTest==1){
    if(var==0){
      res=1.;
    }
    if(var==1){
      res=0.;
    }
    if(var==2){
      res=0.;
    }
    if(var==3){
      res=-(1./(Param.Euler.gamma-1.))*(p,g)+6;
    }
  }

  if(d.nTest==2){
    if(var==0){
      res=1.;
    }
    if(var==1){
      res=0.;
    }
    if(var==2){
      res=0.;
    }
    if(var==3){
      res=-(1./(Param.Euler.gamma-1.))*(p,g)+6;
    }
  }
  if(d.nTest==3){
    if(var==0){
      res=p.y+0.5;
    }
    if(var==1){
      res=0.;
    }
    if(var==2){
      res=0.;
    }
    if(var==3){
      res=(1./(Param.Euler.gamma-1.))*(-(0.5*pow(p.y,2.)+0.5*p.y)*g.y+10);
    }
  }

   if(d.nTest==4){
    if(var==0){
      res=p.y+0.5;
    }
    if(var==1){
      res=0;
    }
    if(var==2){
      res=0;
    }
    if(var==3){
      res=(1./(Param.Euler.gamma-1.))*(-(0.5*pow(p.y,2.)+0.5*p.y)*g.y+10);
    }
  }

   if(d.nTest==7){
     u0=1;
     v0=1.125;
     X=0.75+u0*t-0.5*(g.x*g.x)*t*t;
     Y=0.75+v0*t-0.5*(g.y*g.y)*t*t;
     p0=1;
     rho1=1.;
     rho2 =0.05;
     r=sqrt(pow((p.x-X),2.)+pow((p.y-Y),2.));
    if(var==0){
      if(r <= 0.5){
       	res=rho1+rho2*(sin((r/0.5)*M_PI+0.5*M_PI)+1);
        }
	else {res=rho1; }
    }
    if(var==1){
       if(r <= 0.5){
	 rho=rho1+rho2*(sin((r/0.5)*M_PI+0.5*M_PI)+1); }
       else {
	 rho=rho1;}
       res=rho*(u0+g.x*t);
    
    }
    if(var==2){
       if(r <= 0.5){
	 rho=rho1+rho2*(sin((r/0.5)*M_PI+0.5*M_PI)+1); }
       else {
	 rho=rho1;}
       res=rho*(v0+g.y*t);
    }
    if(var==3){
       if(r <= 0.5){
	 rho=rho1+rho2*(sin((r/0.5)*M_PI+0.5*M_PI)+1); }
       else {
	 rho=rho1;}
       res=rho*((p0/((Param.Euler.gamma-1)*rho))+0.5*((u0+g.x*t)*(u0+g.x*t)+(v0+g.y*t)*(v0+g.y*t)));
    }
   
   }


  if(d.nTest==5){
    if(var==0){
      res=exp(-(g.x*p.x+g.y*p.y));
    }
    if(var==1){
      res=0.;
    }
    if(var==2){
      res=0.;
    }
    if(var==3){
      res=(1./(Param.Euler.gamma-1.))*exp(-(g.x*p.x+g.y*p.y));
    }
  }

   if(d.nTest==6){
    if(var==0){
      res=exp(-(g.x*p.x+g.y*p.y));
    }
    if(var==1){
      res=0;
    }
    if(var==2){
      res=0;
    }
    if(var==3){
      res=(1./(Param.Euler.gamma-1.))*exp(-(g.x*p.x+g.y*p.y));
    }
  }
  return res;
   
}



double SolFonEulerAverage(Data & d,Mesh & Mh, ParamPhysic & Param,Vertex p,double t,Vertex center,int numCell,int var){
    /** function which gives the Euler analytical solution for p, t, center and the variable "var"  **/ 
  double res=0;
  R2 g(0,0);
  g=GravityVector(d,Param.Euler);
  R2 q(0.,0.);
   double J=0;
  double Qp1=0;
  double Qp2=0;
  double w1=0,w2=0;
  int Order=0;
  Order=Param.Euler.LHO;
 
  if(d.nTest==1){
     for(int i=0;i<Order+3;i++){
       for(int j=0;j<Order+3;j++){
	 QuadratureUnitSquare(Qp1,w1,i,Order+3);
	 QuadratureUnitSquare(Qp2,w2,j,Order+3);

	 Mapping(d,Mh,q,J,Qp1,Qp2,Mh.cells[numCell]);
	 if(var==0){
	   res=res+w1*w2*Valabs(J)*1.;
	 }
	 if(var==1){
	   res=0.;
	 }
	 if(var==2){
	   res=0.;
	 }
	 if(var==3){
	   res=res+w1*w2*Valabs(J)*(-(1./(Param.Euler.gamma-1.))*(q,g)+6);
	 }
       }
     }
 
  }

  if(d.nTest==3){
    res=0; 
     for(int i=0;i<Order+3;i++){
       for(int j=0;j<Order+3;j++){
	 QuadratureUnitSquare(Qp1,w1,i,Order+3);
	 QuadratureUnitSquare(Qp2,w2,j,Order+3);;
	 Mapping(d,Mh,q,J,Qp1,Qp2,Mh.cells[numCell]);
	 if(var==0){
	   res=res+w1*w2*Valabs(J)*(q.y+0.5);
	 }
	 if(var==1){
	   res=0.;
	 }
	 if(var==2){
	   res=0.;
	 }
	 if(var==3){
	   res=res+w1*w2*Valabs(J)*((1./(Param.Euler.gamma-1.))*(-(0.5*pow(q.y,2.)+0.5*q.y)*g.y+10));
	 }
       }
     }
     
  }

  if(d.nTest==5){
     for(int i=0;i<Order+3;i++){
       for(int j=0;j<Order+3;j++){
	 QuadratureUnitSquare(Qp1,w1,i,Order+3);
	 QuadratureUnitSquare(Qp2,w2,j,Order+3);

	 Mapping(d,Mh,q,J,Qp1,Qp2,Mh.cells[numCell]);
	 if(var==0){
	   res=res+w1*w2*Valabs(J)*(exp(-(g.x*q.x+g.y*q.y)));
	  
	 }
	 if(var==1){
	   res=0.;
	 }
	 if(var==2){
	   res=0.;
	 }
	 if(var==3){
	   res=res+w1*w2*Valabs(J)*((1./(Param.Euler.gamma-1.))*exp(-(g.x*q.x+g.y*q.y)));
	 }
       }
     }
  }



  res=res/Mh.area(numCell);

  return res;
   
}
