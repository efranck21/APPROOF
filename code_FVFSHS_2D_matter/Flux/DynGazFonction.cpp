
#include <iostream>
#include <cmath>
#include <cassert>
#include "Flux.hpp"
#include "Variables.hpp"
#include "Functions.hpp"
#include <stdio.h>
#include "Tensor.hpp"
#include "BoundaryCondition.hpp"
#include "Functionsinit.hpp"
#include "FunctionsAdvection.hpp"
#include "ParamPhys.hpp"
#include "HOfunctions.hpp"


/** This file contains different functions useful for the discretization of the Euler equations **/

double PressureLaw(Data & d,Mesh & Mh, variable & v,ParamEuler & Euler, int numCell){
  /** This function compute the pressure law. Actuzlly only the perfect gas law is implemented **/
  double res=0;

    if(Euler.lp==1){
      res=(Euler.gamma-1.)*v.var[0][numCell]*EnergyIn(d,Mh,v,Euler,numCell); 
    }


  return res;
}

double rhor(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamEuler & Euler,int numGr){
  /** This function compute a value of rho at the node numGr. 
This value is given by the value in the cell NumCell or by an average around the node numGr **/
  int nbmeans=0,typerho=0;
  int jG=0;
  double m=0;
  nbmeans=2;
  typerho=1;
  if(typerho==1){

    if(nbmeans==1){
      for(int j=0;j<tab.TabInv[numGr].taille;j++){
	jG=tab.TabInv[numGr].TabCell[j];
	m=m+Vjr(Mh,jG,numGr)*v.var[0][jG];
      }

      m=m/Vr(Mh,numGr,tab);
    }
    if(nbmeans==2){
      for(int j=0;j<tab.TabInv[numGr].taille;j++){
	jG=tab.TabInv[numGr].TabCell[j];
	m=m+v.var[0][jG];
      }
    m=m/tab.TabInv[numGr].taille;
    
    }
  }
  if(typerho==2){
    m=v.var[0][jG];
  } 

  return m;
}


double c(Data & d,Mesh & Mh, variable & v,ParamEuler & Euler, int numCell){
  /** This function compute the sound velocity using the internal energy. **/
  double res=0;

  res=sqrt((Euler.gamma-1)*Euler.gamma*EnergyIn(d,Mh,v,Euler,numCell)); 
  
  return res;
}


double EnergyIn(Data & d,Mesh & Mh, variable & v,ParamEuler & Euler, int numCell){
  /** This function compute the internal energy. **/
  double res=0.;
  double temp=0.;
  
  temp=pow(v.var[0][numCell],2.);
  res=(v.var[3][numCell]/v.var[0][numCell])-0.5*((pow(v.var[1][numCell],2.)+pow(v.var[2][numCell],2.))/temp);
  if(res<0){
    cout<<"internal energy negative "<<numCell<<" "<<Mh.xj(numCell).lab<<endl;
  }

  return res;
}


double WaveSpeed(Data & d,Mesh & Mh, variable & v, int numCell,int r,TabConnecInv & tab,ParamEuler & Euler){
  /** This function compute the discrete wave speed using an average around the node numGr. **/
  double res=0;
  int numGr=0,jG=0;
 
      numGr=r;//Mh(numCell,r);
      for(int j=0;j<tab.TabInv[numGr].taille;j++){
	jG=tab.TabInv[numGr].TabCell[j];  
	res=res+(1./tab.TabInv[numGr].taille)*v.var[0][jG]*c(d,Mh,v,Euler,jG); 
      }


  return res;
}


double remap(Data & d,Mesh & Mh, variable & v,int var,TabConnecInv & tab,ParamPhysic & Param,int numCell,int r, R2 * ur, R2 a){
  /** This function compute the remap part of the fluxes using the advction fluxes at the velocity a. **/
  double res=0;

    res=VertexUpwind(d,Mh,v,var,tab,Param,numCell,r,a);

  return res;

}


R2 RhoGravityVector(Data & d,Mesh &Mh, variable & v, TabConnecInv & tab,ParamEuler & Euler, int numGr){
   /** Function wich compute rho g at the node r  **/
  R2 res(0,0);
  R2 g(0,0);
  double rhonode=0.;
  double rhonodeHO=0.;
  R2 GradpnodeHO(0,0);
  double a11=0,a12=0,a21=0,a22=0,b1=0,b2=0,Det=0;
  int numLrj,jG;
  R2 sol(0,0);
  tensor beta(d,'s');
  
  
  g=GravityVector(d,Euler);
  
  if(Euler.LHO ==1) {
   
    rhonode=rhor(d,Mh,v,tab,Euler,numGr);
    res.x= rhonode*g.x;
    res.y= rhonode*g.y;

  }
  else {
    
    SolverLSProblem(d,Mh,v,tab,Euler,Euler.TabLS[numGr],Euler.LHO,numGr);
    GradpnodeHO=IntGradiantPressureInterpolation(d,Mh,v,tab,Euler,numGr);
    rhonodeHO=IntDensityInterpolation(d,Mh,v,tab,Euler,numGr);
      
    for(int j=0;j<tab.TabInv[numGr].taille;j++){
      jG=tab.TabInv[numGr].TabCell[j];   
      numLrj=NodeGtoL(Mh,jG,numGr);  
      
      a11=a11-Mh.ljr(jG,numLrj)*Mh.njr(jG,numLrj).x*(Mh.xr(jG,numLrj).x-Mh.xj(jG).x);
      a12=a12-Mh.ljr(jG,numLrj)*Mh.njr(jG,numLrj).x*(Mh.xr(jG,numLrj).y-Mh.xj(jG).y);
      a21=a21-Mh.ljr(jG,numLrj)*Mh.njr(jG,numLrj).y*(Mh.xr(jG,numLrj).x-Mh.xj(jG).x);
      a22=a22-Mh.ljr(jG,numLrj)*Mh.njr(jG,numLrj).y*(Mh.xr(jG,numLrj).y-Mh.xj(jG).y); 
      
      b1=b1+Mh.ljr(jG,numLrj)*PressureLaw(d,Mh,v,Euler,jG)*Mh.njr(jG,numLrj).x;
      b2=b2+Mh.ljr(jG,numLrj)*PressureLaw(d,Mh,v,Euler,jG)*Mh.njr(jG,numLrj).y;
    }
    
    Det=a11*a22-a21*a12;
     
     if(Mh.vertices[numGr].lab >-1){
       
       if(Det==0){cout<<"The nodal matrix for the gradiant WB coef is not invertible"<<numGr<<endl;exit(1);}
       else { 
	 sol.x=(b1*a22-b2*a12)/Det;
	 sol.y=(b2*a11-b1*a21)/Det;
       }
     }
     else{
       sol.x=b1/Vr(Mh,numGr,tab);
       sol.y=b2/Vr(Mh,numGr,tab);
     }


     
    res.x=GradpnodeHO.x+rhonodeHO*g.x-sol.x;
    res.y=GradpnodeHO.y+rhonodeHO*g.y-sol.y;
  }
 
  return res;
}
  

R2 IntGradiantPressureInterpolation(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamEuler & Euler,int numGr){ /** NOTTT FINISH **/
  /** function which compute the integral on the volume control to the derivative of the polynomial reconstruction of the pressure **/
  R2 res(0,0);
  R2 subres(0,0); 
  cell CVlr=0;
  R2 p(0,0);
  int d1=0,d2=0;
  double w1=0,w2=0,J=0;
  double Qp1=0;
  double Qp2=0;
  double area=0;

  
  
  for(int l=0;l<Mh.nbnodelocal;l++)   {
    CVlr=SubControlVolume(d,Mh,v,tab,numGr,l);
    subres.x=0;
    subres.y=0;
      for(int i=0;i<Euler.LHO+3;i++){
	for(int j=0;j<Euler.LHO+3;j++){
	  QuadratureUnitSquare(Qp1,w1,i,Euler.LHO+3);
	  QuadratureUnitSquare(Qp2,w2,j,Euler.LHO+3);
  
	  Mapping(d,Mh,p,J,Qp1,Qp2,CVlr);
          
	  for (int k = 0; k < Euler.TabLS[numGr].NColum; k++ )
	    {
	      LocalDegreePolynomial(d1,d2,k,Euler.LHO);
              if(d1!=0)  {
		subres.x=subres.x+w1*w2*J*d1*Euler.TabLS[numGr].Sol_pressure[k]*pow(p.x,d1-1)*pow(p.y,d2);
	      }
	    }
           for (int k = 0; k < Euler.TabLS[numGr].NColum; k++ )
	    {
	      LocalDegreePolynomial(d1,d2,k,Euler.LHO);
              if(d2!=0)  {
		subres.y=subres.y+w1*w2*J*d2*Euler.TabLS[numGr].Sol_pressure[k]*pow(p.x,d1)*pow(p.y,d2-1);
	      }
	    }

	 }
      }
      res.x=res.x+(subres.x);
      res.y=res.y+(subres.y);
      area=area+CVlr.area;
}
  
  res.x=res.x/area;
  res.y=res.y/area;
 
  return res;
}

double IntDensityInterpolation(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamEuler & Euler,int numGr){
   /** function which compute the integral on the volume control to the polynomial reconstruction of the density **/
  double res=0,subres=0;
  cell CVlr=0;
  R2 p(0,0);
  int d1=0,d2=0;
  double w1=0,w2=0,J=0;
  double Qp1=0;
  double Qp2=0;
  double area=0;

  for(int l=0;l<Mh.nbnodelocal;l++)   {
    CVlr=SubControlVolume(d,Mh,v,tab,numGr,l);
    subres=0;
      for(int i=0;i<Euler.LHO+3;i++){
	for(int j=0;j<Euler.LHO+3;j++){
	  QuadratureUnitSquare(Qp1,w1,i,Euler.LHO+3);
	  QuadratureUnitSquare(Qp2,w2,j,Euler.LHO+3);
      
	  Mapping(d,Mh,p,J,Qp1,Qp2,CVlr);
	  for (int k = 0; k < Euler.TabLS[numGr].NColum; k++ )
	  {
	   LocalDegreePolynomial(d1,d2,k,Euler.LHO);
	      subres=subres+w1*w2* Valabs(J)*Euler.TabLS[numGr].Sol_rho[k]*pow(p.x,d1)*pow(p.y,d2);
	  }

	 }
      }
      res=res+subres;
      area=area+CVlr.area;
  }

  res=res/area;
  return res;
}


void InitTab_rhoGravity(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamEuler & Euler,int numGr){
/** Function which initialize the table of Euler quantities using the given value (paramEuler.dat)
or the test case if the eps, sigma or phi are variables **/
  R2 R2ext(-10000,-10000);
  
  if(Mh.xr(numGr).lab == -1) {
    Euler.TabRhoG[numGr]=R2ext;
  }  
  else
    {
      Euler.TabRhoG[numGr]=RhoGravityVector(d,Mh,v,tab,Euler,numGr);
    }
  
}  

R2 GravityVector(Data & d,ParamEuler & Euler){
  /** Function wich compute g  **/
  R2 res(0,0);

  res.x=Euler.gx;
  res.y=Euler.gy;

  return res;
}
