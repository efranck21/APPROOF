#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Variables.hpp"
#include "Flux.hpp"
#include "Source.hpp"
#include "BoundaryCondition.hpp"
#include "Tensor.hpp"



R2 SolveurImpliciteP1(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param, double dt,double b1, double b2){
     /** Function wich compute the implicit contribution to the source term for P1 **/
  R2 res(0,0);
  double a11=0,a12=0,a21=0,a22=0,Det=0;
  double Mr[2][2];
  tensor alpha(d,'h');
  int numGr=0;

  if(d.scheme == 1 ){
    
    a11=1.+dt*((Param.P1.sigma[numCell])/(Param.P1.eps[numCell]*Param.P1.eps[numCell]));
    a12=0;
    a21=0;
    a22=1.+dt*((Param.P1.sigma[numCell])/(Param.P1.eps[numCell]*Param.P1.eps[numCell]));
  }

  if(d.scheme == 3){
    
    for(int r=0;r<Mh.nbnodelocal;r++){
      
      numGr=Mh(numCell,r);
      alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
      MrConstructP1(d,numGr,Mh,v,tab,Param,Mr);
      a11=a11+Mh.ljr(numCell,r)*(alpha.ten[0][0]*(1.-Mr[0][0])+alpha.ten[0][1]*(-Mr[1][0]));
      a12=a12+Mh.ljr(numCell,r)*(alpha.ten[0][0]*(-Mr[0][1])+alpha.ten[0][1]*(1.-Mr[1][1]));
      a21=a21+Mh.ljr(numCell,r)*(alpha.ten[1][0]*(1.-Mr[0][0])+alpha.ten[1][1]*(-Mr[1][0]));
      a22=a22+Mh.ljr(numCell,r)*(alpha.ten[1][0]*(-Mr[0][1])+alpha.ten[1][1]*(1.-Mr[1][1]));

    }
    
    
    a11=(dt/(Mh.area(numCell)*Param.P1.eps[numCell]))*a11+1.;
    a12=(dt/(Mh.area(numCell)*Param.P1.eps[numCell]))*a12;
    a21=(dt/(Mh.area(numCell)*Param.P1.eps[numCell]))*a21;
    a22=(dt/(Mh.area(numCell)*Param.P1.eps[numCell]))*a22+1.;
    
  }
  
  Det=a11*a22-a21*a12;
  if(Det==0){ cout<<"The matrix for implicit P1 part is not invertible"<<endl;exit(1);}
  else { 
    res.x=(b1*a22-b2*a12)/Det;
    res.y=(b2*a11-b1*a21)/Det;
  }

  return res;

}


vectorflux SourceP1(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param) {
  /** Function wich compute the explicit source terms for P1 **/
  vectorflux res(d);
  R2 sol(0,0);
  double Q=0;


  if(d.Typetime == 'E'){
    if( d.scheme==1 ){
      res.vflux[0]=0;
      res.vflux[1]=-((Param.P1.sigma[numCell])/(Param.P1.eps[numCell]*Param.P1.eps[numCell]))*Mh.area(numCell)*v.var[1][numCell];
      res.vflux[2]=-((Param.P1.sigma[numCell])/(Param.P1.eps[numCell]*Param.P1.eps[numCell]))*Mh.area(numCell)*v.var[2][numCell];
      
    }
    
    
    if( d.scheme==2 ){
      res.vflux[0]=0;
      res.vflux[1]=0;
      res.vflux[2]=0;  
    }
    
    if(d.scheme==3){
      double Mr[2][2];
      MatrixSourceP1(d,numCell,Mh,v,tab,Param,Mr);
      res.vflux[0]=Q;
      res.vflux[1]=Mr[0][0]*v.var[1][numCell]+Mr[0][1]*v.var[2][numCell];
      res.vflux[2]=Mr[1][0]*v.var[1][numCell]+Mr[1][1]*v.var[2][numCell];
    }
  }

  if(d.Typetime == 'S'){
    res.vflux[0]=0;
    res.vflux[1]=0;
    res.vflux[2]=0;  
    
  }

  if(d.nTest ==6 ){
    int sigmaA=0;
    double Q=0;
    if(Mh.cells[numCell].lab==2) {sigmaA=0; Q=1;}
    if(Mh.cells[numCell].lab==0) {sigmaA=0; Q=0;}
    if(Mh.cells[numCell].lab==1) {sigmaA=30; Q=0;}
    res.vflux[0]=res.vflux[0]+(Q-sigmaA*v.var[1][numCell])*Mh.area(numCell);
  }
  
  return res;
}
  


void MatrixSourceP1(Data & d,int numCell,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,double Mrsource[2][2]){
   /** Function wich compute the nodal matrix for explicit source term **/
  R2 res(0,0);
  double a11=0,a12=0,a21=0,a22=0;
  double Mr[2][2];
  tensor alpha(d,'h');
  int numGr=0;

  for(int r=0;r<Mh.nbnodelocal;r++){
    
    numGr=Mh(numCell,r);
    alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
    MrConstructP1(d,numGr,Mh,v,tab,Param,Mr);
    a11=a11+Mh.ljr(numCell,r)*(alpha.ten[0][0]*(1.-Mr[0][0])+alpha.ten[0][1]*(-Mr[1][0]));
    a12=a12+Mh.ljr(numCell,r)*(alpha.ten[0][0]*(-Mr[0][1])+alpha.ten[0][1]*(1.-Mr[1][1]));
    a21=a21+Mh.ljr(numCell,r)*(alpha.ten[1][0]*(1.-Mr[0][0])+alpha.ten[1][1]*(-Mr[1][0]));
    a22=a22+Mh.ljr(numCell,r)*(alpha.ten[1][0]*(-Mr[0][1])+alpha.ten[1][1]*(1.-Mr[1][1]));

  }
  Mrsource[0][0]=-(1./Param.P1.eps[numCell])*a11;
  Mrsource[0][1]=-(1./Param.P1.eps[numCell])*a12;
  Mrsource[1][0]=-(1./Param.P1.eps[numCell])*a21;
  Mrsource[1][1]=-(1./Param.P1.eps[numCell])*a22;


}
