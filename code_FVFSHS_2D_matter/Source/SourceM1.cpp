#include <iostream>
#include <cmath>
#include <cassert>
#include "Variables.hpp"
#include "Flux.hpp"
#include "Source.hpp"
#include "BoundaryCondition.hpp"
#include "Tensor.hpp"
#include "FunctionsAdvection.hpp"



vectorflux SourceM1(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab, ParamPhysic & Param){
 /** Function wich compute the explicit source terms for M1 **/

  
  vectorflux res(d);
  double kj=0;
  double Mr[2][2];
  R2 uj(0,0);
    
  int compton_scattering = 1 ; //(1 if taken into account, 0 if not)

  kj=coefk(d,Param.M1,v.var[0][numCell],v.var[1][numCell],v.var[2][numCell]);
  Calcul_u(d,Mh,v,tab,Param.M1,uj,numCell);
  
  if(d.Typetime == 'E')
  {
    if( d.scheme==1 )
    {
      res.vflux[0] = 0. ;
      res.vflux[1] = -((Param.M1.sigma[numCell])/(Param.M1.eps[numCell]*Param.M1.eps[numCell]))*Mh.area(numCell)*v.var[1][numCell];
      res.vflux[2] = -((Param.M1.sigma[numCell])/(Param.M1.eps[numCell]*Param.M1.eps[numCell]))*Mh.area(numCell)*v.var[2][numCell];
      
      if (compton_scattering == 1)
      {
          
          double T = 1. ; // temperature electron
          double sigma_a = Param.M1.sigma[numCell] ;
          
          double temp = 0. ;
          double source_E = 0. ;
          double source_F = 0. ;
          R2 g ;
          
          g.x = v.var[1][numCell] / v.var[0][numCell] ;
          g.y = v.var[2][numCell] / v.var[0][numCell] ;
                    
          temp = 2. + sqrt(4.-3.*(g,g)) ;
          
          source_E = sqrt(sqrt( v.var[0][numCell] )) * temp*sqrt(sqrt(temp)) * ( (4.+(g,g))/(2.*temp-(g,g)) ) * sqrt(sqrt( (temp-(g,g))/(1.-(g,g)) )) ;
          
          source_F = sqrt(sqrt( v.var[0][numCell] )) * sqrt(sqrt(temp)) * ( (20.+(g,g))/(2.*temp+(g,g)) ) * sqrt(sqrt( (temp-(g,g))/(1.-(g,g)) )) ;

          res.vflux[0] += Mh.area(numCell)*sigma_a*(T*T*T*T - v.var[0][numCell]) + Mh.area(numCell)*Param.M1.sigma[numCell]*v.var[0][numCell]*( 4.*T - source_E ) ;
          res.vflux[1] += Mh.area(numCell)*( -sigma_a -Param.M1.sigma[numCell]/3. +Param.M1.sigma[numCell]*( 4.*T - source_F ) ) * v.var[1][numCell] ;
          res.vflux[2] += Mh.area(numCell)*( -sigma_a -Param.M1.sigma[numCell]/3. +Param.M1.sigma[numCell]*( 4.*T - source_F ) ) * v.var[2][numCell] ;
          
      }
        
    }
    
    
    if( d.scheme==2 ){
      res.vflux[0]=0;
      res.vflux[1]=0;
      res.vflux[2]=0;  
    }
    
    if(d.scheme==3){
      MatrixSourceM1(d,numCell,Mh,v,tab,Param,Mr);
      res.vflux[0]=0;
      res.vflux[1]=Mr[0][0]*v.var[1][numCell]+Mr[0][1]*v.var[2][numCell];
      res.vflux[2]=Mr[1][0]*v.var[1][numCell]+Mr[1][1]*v.var[2][numCell];
    }
  }

  if(d.Typetime == 'S')
  {
    res.vflux[0]=0;
    res.vflux[1]=0;
    res.vflux[2]=0;  
    
  }
  
  return res;
}




R2 SolveurImpliciteM1(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab, ParamPhysic & Param, double dt,double b1, double b2){
   /** Function wich compute the implicit contribution to the source term for M1 **/
  
  R2 res(0,0);
  double a11=0,a12=0,a21=0,a22=0,Det=0,rjr=0;
  double Mr[2][2];
  tensor alpha(d,'h');
  int numGr=0;
  double kj=0;
  R2 uj(0,0);

  kj=coefk(d,Param.M1,v.var[0][numCell],v.var[1][numCell],v.var[2][numCell]);
  
  if(d.scheme == 1 ){
    
    a11=1.+(1./kj)*dt*((Param.M1.sigma[numCell])/(Param.M1.eps[numCell]*Param.M1.eps[numCell]));
    a12=0;
    a21=0;
    a22=1.+(1./kj)*dt*((Param.M1.sigma[numCell])/(Param.M1.eps[numCell]*Param.M1.eps[numCell]));
  }

  if(d.scheme == 3){
    for(int r=0;r<Mh.nbnodelocal;r++){
      
      numGr=Mh(numCell,r);
      alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
      MrConstructM1(d,numGr,Mh,v,tab,Param,Mr);
      Calcul_u(d,Mh,v,tab,Param.M1,uj,numCell);
      rjr=(4./sqrt(3.))*(v.var[0][numCell]/(3+(uj,uj)));

      a11=a11+Mh.ljr(numCell,r)*rjr*(alpha.ten[0][0]*(1.-Mr[0][0])+alpha.ten[0][1]*(-Mr[1][0]));
      a12=a12+Mh.ljr(numCell,r)*rjr*(alpha.ten[0][0]*(-Mr[0][1])+alpha.ten[0][1]*(1.-Mr[1][1]));
      a21=a21+Mh.ljr(numCell,r)*rjr*(alpha.ten[1][0]*(1.-Mr[0][0])+alpha.ten[1][1]*(-Mr[1][0]));
      a22=a22+Mh.ljr(numCell,r)*rjr*(alpha.ten[1][0]*(-Mr[0][1])+alpha.ten[1][1]*(1.-Mr[1][1]));
      
    }
    a11=(1./kj)*(dt/(Mh.area(numCell)*Param.M1.eps[numCell]))*a11+1.;
    a12=(1./kj)*(dt/(Mh.area(numCell)*Param.M1.eps[numCell]))*a12;
    a21=(1./kj)*(dt/(Mh.area(numCell)*Param.M1.eps[numCell]))*a21;
    a22=(1./kj)*(dt/(Mh.area(numCell)*Param.M1.eps[numCell]))*a22+1.;
    
  }
  
  Det=a11*a22-a21*a12;
  if(Det==0){ cout<<"The matrix for implicit M1 part is not invertible"<<endl;exit(1);}
  else { 
    res.x=(b1*a22-b2*a12)/Det;
    res.y=(b2*a11-b1*a21)/Det;
  }


  return res;

}


void MatrixSourceM1(Data & d,int numCell,Mesh & Mh,variable & v, TabConnecInv & tab, ParamPhysic & Param,double Mrsource[2][2]){
   /** Function wich compute the nodal matrix for explicit source term **/
  R2 res(0,0);
  double a11=0,a12=0,a21=0,a22=0,rjr=0;
  double Mr[2][2];
  tensor alpha(d,'h');
  int numGr=0;
  R2 uj(0,0);

  for(int r=0;r<Mh.nbnodelocal;r++){
    
    numGr=Mh(numCell,r);
    alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
    MrConstructM1(d,numGr,Mh,v,tab,Param,Mr);
    Calcul_u(d,Mh,v,tab,Param.M1,uj,numCell);
    rjr=(4./sqrt(3.))*(v.var[0][numCell]/(3.+(uj,uj)));
    a11=a11+Mh.ljr(numCell,r)*rjr*(alpha.ten[0][0]*(1.-Mr[0][0])+alpha.ten[0][1]*(-Mr[1][0]));
    a12=a12+Mh.ljr(numCell,r)*rjr*(alpha.ten[0][0]*(-Mr[0][1])+alpha.ten[0][1]*(1.-Mr[1][1]));
    a21=a21+Mh.ljr(numCell,r)*rjr*(alpha.ten[1][0]*(1.-Mr[0][0])+alpha.ten[1][1]*(-Mr[1][0]));
    a22=a22+Mh.ljr(numCell,r)*rjr*(alpha.ten[1][0]*(-Mr[0][1])+alpha.ten[1][1]*(1.-Mr[1][1]));

  }
  Mrsource[0][0]=-(1./Param.M1.eps[numCell])*a11;
  Mrsource[0][1]=-(1./Param.M1.eps[numCell])*a12;
  Mrsource[1][0]=-(1./Param.M1.eps[numCell])*a21;
  Mrsource[1][1]=-(1./Param.M1.eps[numCell])*a22;


}







