#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Variables.hpp"
#include "Flux.hpp"
#include "Source.hpp"
#include "BoundaryCondition.hpp"
#include "Tensor.hpp"
#include "FunctionsAdvection.hpp"

vectorflux SourceP1Matter(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param) {
  /** Function wich compute the explicit source terms for P1 with matter**/
  vectorflux res(d);
  R2 sol(0,0);

  if( d.Typetime=='E' ){
    if( d.scheme==1 ){
      res.vflux[0]=0;
      res.vflux[1]=-((Param.P1M.sigmaA[numCell])/(Param.P1M.eps[numCell]*Param.P1M.eps[numCell]))*Mh.area(numCell)*v.var[1][numCell];
      res.vflux[2]=-((Param.P1M.sigmaA[numCell])/(Param.P1M.eps[numCell]*Param.P1M.eps[numCell]))*Mh.area(numCell)*v.var[2][numCell];
    }
    
    if( d.scheme==2 ){
      res.vflux[0]=0;
      res.vflux[1]=0;
      res.vflux[2]=0;  
    }
    
    if(d.scheme==3){
      double Mr[2][2];
      MatrixSourceP1Matter(d,numCell,Mh,v,tab,Param,Mr);
      res.vflux[0]=0;
      res.vflux[1]=Mr[0][0]*v.var[1][numCell]+Mr[0][1]*v.var[2][numCell];
      res.vflux[2]=Mr[1][0]*v.var[1][numCell]+Mr[1][1]*v.var[2][numCell];
    }
    
    res.vflux[0]=((Param.P1M.sigmaA[numCell])/(Param.P1M.eps[numCell]*Param.P1M.eps[numCell]))*(Param.P1M.a_value*pow(v.var[v.nbvar-1][numCell],4.)-v.var[0][numCell]);
    
    if(Param.P1M.Cv_value!=0){
      res.vflux[3]=((Param.P1M.sigmaA[numCell])/(Param.P1M.eps[numCell]*Param.P1M.eps[numCell]))*(Param.P1M.a_value*pow(v.var[v.nbvar-1][numCell],4.)-v.var[0][numCell])/Param.P1M.Cv_value;
    }
    else {
      res.vflux[3]=0;
    }
  }

   if( d.Typetime=='S' ){
      res.vflux[0]=0;
      res.vflux[1]=0;
      res.vflux[2]=0;
      res.vflux[3]=0;
    }
    
  
  return res;
}


void MatrixSourceP1Matter(Data & d,int numCell,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,double Mrsource[2][2]){
    /** Function wich compute the nodal matrix for explicit source term **/
  R2 res(0,0);
  double a11=0,a12=0,a21=0,a22=0;
  double Mr[2][2];
  tensor alpha(d,'h');
  int numGr=0;

  for(int r=0;r<Mh.nbnodelocal;r++){
    
    numGr=Mh(numCell,r);
    alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
    MrConstructP1Matter(d,numGr,Mh,v,tab,Param,Mr);
    a11=a11+Mh.ljr(numCell,r)*(alpha.ten[0][0]*(1.-Mr[0][0])+alpha.ten[0][1]*(-Mr[1][0]));
    a12=a12+Mh.ljr(numCell,r)*(alpha.ten[0][0]*(-Mr[0][1])+alpha.ten[0][1]*(1.-Mr[1][1]));
    a21=a21+Mh.ljr(numCell,r)*(alpha.ten[1][0]*(1.-Mr[0][0])+alpha.ten[1][1]*(-Mr[1][0]));
    a22=a22+Mh.ljr(numCell,r)*(alpha.ten[1][0]*(-Mr[0][1])+alpha.ten[1][1]*(1.-Mr[1][1]));

  }
  Mrsource[0][0]=-(1./Param.P1M.eps[numCell])*a11;
  Mrsource[0][1]=-(1./Param.P1M.eps[numCell])*a12;
  Mrsource[1][0]=-(1./Param.P1M.eps[numCell])*a21;
  Mrsource[1][1]=-(1./Param.P1M.eps[numCell])*a22;


}


R2 SolveurImpliciteP1Matter(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param, double dt,double b1, double b2){
   /** Function wich compute the implicit contribution to the source term for P1 with matter (no relaxation with matter) **/
  R2 res(0,0);
  double a11=0,a12=0,a21=0,a22=0,Det=0;
  double Mr[2][2];
  tensor alpha(d,'h');
  int numGr=0;

  if(d.scheme == 1 ){
    
    a11=1.+dt*((Param.P1M.sigmaA[numCell])/(Param.P1M.eps[numCell]*Param.P1M.eps[numCell]));
    a12=0;
    a21=0;
    a22=1.+dt*((Param.P1M.sigmaA[numCell])/(Param.P1M.eps[numCell]*Param.P1M.eps[numCell]));
  }

  if(d.scheme == 3){
    
    for(int r=0;r<Mh.nbnodelocal;r++){
      
      numGr=Mh(numCell,r);
      alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
      MrConstructP1Matter(d,numGr,Mh,v,tab,Param,Mr);
      a11=a11+Mh.ljr(numCell,r)*(alpha.ten[0][0]*(1.-Mr[0][0])+alpha.ten[0][1]*(-Mr[1][0]));
      a12=a12+Mh.ljr(numCell,r)*(alpha.ten[0][0]*(-Mr[0][1])+alpha.ten[0][1]*(1.-Mr[1][1]));
      a21=a21+Mh.ljr(numCell,r)*(alpha.ten[1][0]*(1.-Mr[0][0])+alpha.ten[1][1]*(-Mr[1][0]));
      a22=a22+Mh.ljr(numCell,r)*(alpha.ten[1][0]*(-Mr[0][1])+alpha.ten[1][1]*(1.-Mr[1][1]));

    }
    
    
    a11=(dt/(Mh.area(numCell)*Param.P1M.eps[numCell]))*a11+1.;
    a12=(dt/(Mh.area(numCell)*Param.P1M.eps[numCell]))*a12;
    a21=(dt/(Mh.area(numCell)*Param.P1M.eps[numCell]))*a21;
    a22=(dt/(Mh.area(numCell)*Param.P1M.eps[numCell]))*a22+1.;
    
  }
  
  Det=a11*a22-a21*a12;
  if(Det==0){ cout<<"la matrice du solveur nodale P1 matter n'est pas inversible"<<endl;exit(1);}
  else { 
    res.x=(b1*a22-b2*a12)/Det;
    res.y=(b2*a11-b1*a21)/Det;
  }

  return res;

}


R2 SplittingCalculMatter(Data & d,ParamPhysic & Param ,int numCell,Mesh &  Mh,variable & vnp,variable  & vn, TabConnecInv & tab,double dt,double dT){
 /** Function wich compute the implicit contribution to the matter relaxation **/
  double Ej,Tj;
  double Ep=0,Tp=0,Tpm=0;
  double Kp=0;
  double pmax=25;
  double dteps=0.;
  R2 res(0,0);
 
  double a11,a12,a21,a22,Det,b1,b2;
  double mup, thetap;


  Ej=vnp.var[0][numCell];
  Tj=vnp.var[vnp.nbvar-1][numCell];
 
  Ep=Ej;
  Tp=Tj;
   
   Kp=(Param.P1M.Cv_value*Tp)/(Param.P1M.a_value*pow(Tj,4.));
   mup=1/(4*pow(Tp,3.));
   
   int p=0;
   for(p=0;p<pmax;p++){
     
     dteps=(1./(pow(Param.P1M.eps[numCell],2.)))*Param.P1M.sigmaA[numCell]*dt;
     
     a11=1+dteps;
     a12=-dteps;
     a21=-dteps;
     a22=Param.P1M.Cv_value*mup+dteps;
     
     b1=Ej;
     b2=Param.P1M.Cv_value*mup*Param.P1M.a_value*pow(Tj,4.);
     Tpm=Tp;
     Det=a11*a22-a21*a12;
     
     if(Det==0){ cout<<"the matrix associated with the matter source term is not invertible"<<endl;exit(1);}
     else { 
       Ep=(b1*a22-b2*a12)/Det;
       thetap=(b2*a11-b1*a21)/Det;  
     }
     
     Tp=pow(thetap/Param.P1M.a_value,1./4.);     
     mup=1/(Param.P1M.a_value*(pow(Tp,3.)+Tj*pow(Tp,2.)+Tp*pow(Tj,2.)+pow(Tj,3.)));
     
     if(Valabs(Tp-Tpm)/Tpm<dT){
       break;
     }
     
   }
   
   
   res.x=Ep;
   res.y=Tp;
      
   return res;
   
}
