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
#include "Source.hpp"



vectorflux SourceE(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 * ur){
/** Function wich compute the explicit source terms for Euler **/
  
  vectorflux res(4);
  R2 g(0,0);
  g=GravityVector(d,Param.Euler);

  if((d.scheme==1 || d.scheme==4) && d.Typetime=='E'){
    res.vflux[0]=0.;
    res.vflux[1]=-(Mh.area(numCell)/Param.Euler.eps[numCell])*(v.var[0][numCell]*g.x+(Param.Euler.sigma[numCell]/Param.Euler.eps[numCell])*v.var[1][numCell]);
    res.vflux[2]=-(Mh.area(numCell)/Param.Euler.eps[numCell])*(v.var[0][numCell]*g.y+(Param.Euler.sigma[numCell]/Param.Euler.eps[numCell])*v.var[2][numCell]);
    res.vflux[3]=-(Mh.area(numCell)/Param.Euler.eps[numCell])*((g.x*v.var[1][numCell]+g.y*v.var[2][numCell])+(Param.Euler.sigma[numCell]/Param.Euler.eps[numCell])*((pow(v.var[1][numCell],2.)+pow(v.var[2][numCell],2.))/v.var[0][numCell]));
  }

  if((d.scheme==1 || d.scheme==4) && d.Typetime=='S'){
    res.vflux[0]=0.;
    res.vflux[1]=-(Mh.area(numCell)/Param.Euler.eps[numCell])*v.var[0][numCell]*g.x;
    res.vflux[2]=-(Mh.area(numCell)/Param.Euler.eps[numCell])*v.var[0][numCell]*g.y;
    res.vflux[3]=-(Mh.area(numCell)/Param.Euler.eps[numCell])*(g.x*v.var[1][numCell]+g.y*v.var[2][numCell]);
  }
  if(d.scheme==2){
    res.vflux[0]=0.;
    res.vflux[1]=0.;
    res.vflux[2]=0.;
    res.vflux[3]=0.;
  }
  if(d.scheme==3 && d.Typetime=='E'){
    R2 vectorg(0,0);
    double gu=0.;
    double uu=0.;
    R2 friction(0,0);
    vectorg=SourceGravityEuler(d,numCell,Mh,v,tab,Param);
    gu=SourceGravityEulerU(d,numCell,Mh,v,tab,Param,ur);
    uu=FrictionUcarre(d,numCell,Mh,v,tab,Param,ur);
    friction=FrictionU(d,numCell,Mh,v,tab,Param,ur);
    res.vflux[0]=0.;
    res.vflux[1]=vectorg.x+friction.x;
    res.vflux[2]=vectorg.y+friction.y;
    res.vflux[3]=gu+uu;
  }

  if(d.scheme==3 && d.Typetime=='S'){
    R2 vectorg(0,0);
    double gu=0.;
    double uu=0.;
    vectorg=SourceGravityEuler(d,numCell,Mh,v,tab,Param);
    gu=SourceGravityEulerU(d,numCell,Mh,v,tab,Param,ur);
    uu=FrictionUcarre(d,numCell,Mh,v,tab,Param,ur);
    res.vflux[0]=0.;
    res.vflux[1]=vectorg.x;
    res.vflux[2]=vectorg.y;
    res.vflux[3]=gu+uu;
  }

  return res;
}

void MrConstructE(Data & d,int numGr,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,double  Mr[2][2]){
/** Function wich compute the Magical Matrix Mr for the Euler schemes **/
  double a11=0.,a12=0.,a21=0.,a22=0.,Det=0.,temp=0.;
 double c11=0.,c12=0.,c21=0.,c22=0.;
 int numLrj=0,jG=0;
 double signode=0., epsnode=0.;
 tensor alpha(d,'h');
 tensor beta(d,'s');
 double rjr=0.,rhonode=0.;


 for(int j=0;j<tab.TabInv[numGr].taille;j++){
      jG=tab.TabInv[numGr].TabCell[j];
      numLrj=NodeGtoL(Mh,jG,numGr);  
      
      beta=inittensor(d,Mh,v,tab,'s',numGr,jG);
      alpha=inittensor(d,Mh,v,tab,'h',numGr,jG);   
      signode=AverageSigNodeEuler(d,Mh,v,tab,Param.Euler,numGr,jG);
      epsnode=AverageEpsNodeEuler(d,Mh,v,tab,Param.Euler,numGr,jG);
       
      rjr=WaveSpeed(d,Mh,v,jG,numGr,tab,Param.Euler);
      rhonode=rhor(d,Mh,v,tab,Param.Euler,numGr);
	
       c11=c11+Mh.ljr(jG,numLrj)*rjr*alpha.ten[0][0];
       a11=a11+rhonode*(signode/epsnode)*Mh.ljr(jG,numLrj)*beta.ten[0][0]+Mh.ljr(jG,numLrj)*rjr*alpha.ten[0][0];
       c12=c12+Mh.ljr(jG,numLrj)*rjr*alpha.ten[0][1];
       a12=a12+rhonode*(signode/epsnode)*Mh.ljr(jG,numLrj)*beta.ten[0][1]+Mh.ljr(jG,numLrj)*rjr*alpha.ten[0][1];
       c21=c21+Mh.ljr(jG,numLrj)*rjr*alpha.ten[1][0];
       a21=a21+rhonode*(signode/epsnode)*Mh.ljr(jG,numLrj)*beta.ten[1][0]+Mh.ljr(jG,numLrj)*rjr*alpha.ten[1][0];
       c22=c22+Mh.ljr(jG,numLrj)*rjr*alpha.ten[1][1];
       a22=a22+rhonode*(signode/epsnode)*Mh.ljr(jG,numLrj)*beta.ten[1][1]+Mh.ljr(jG,numLrj)*rjr*alpha.ten[1][1];      
 }
 Det=a11*a22-a21*a12;
 temp=a11;

 if(Mh.xr(numGr).lab <0){
   a11=1.;
   a12=0.;
   a21=0.;
   a22=1.;
 }
 else{
   a11=a22/Det;
   a12=-a12/Det;
   a21=-a21/Det;
   a22=temp/Det;
 }
 Mr[0][0]=a11*c11+a12*c21;
 Mr[0][1]=a11*c12+a12*c22;
 Mr[1][0]=a21*c11+a22*c21;
 Mr[1][1]=a21*c12+a22*c22;
 
}




R2 SolveurImpliciteE(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,double dt,double b1,double b2, double rho){
  /** Function wich compute the implicit contribution to the source term for Euler **/
  R2 res(0,0);

  double a11=0.,a12=0.,a21=0.,a22=0.,Det=0.;
  double Mr[2][2];
  tensor alpha(d,'h');
  int numGr=0;
  double rjr=0;


  if(d.scheme == 1 || d.scheme==4){

    a11=(Param.Euler.sigma[numCell])/(Mh.area(numCell)*Param.Euler.eps[numCell]);
    a12=0;
    a21=0;
    a22=(Param.Euler.sigma[numCell])/(Mh.area(numCell)*Param.Euler.eps[numCell]);
  }

   if(d.scheme == 3){

     for(int r=0;r<Mh.nbnodelocal;r++){
    
       numGr=Mh(numCell,r);
       alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
       MrConstructE(d,numGr,Mh,v,tab,Param,Mr);
       rjr=WaveSpeed(d,Mh,v,numCell,numGr,tab,Param.Euler);
       a11=a11+Mh.ljr(numCell,r)*rjr*(alpha.ten[0][0]*(1.-Mr[0][0])+alpha.ten[0][1]*(-Mr[1][0]));
       a12=a12+Mh.ljr(numCell,r)*rjr*(alpha.ten[0][0]*(-Mr[0][1])+alpha.ten[0][1]*(1.-Mr[1][1]));
       a21=a21+Mh.ljr(numCell,r)*rjr*(alpha.ten[1][0]*(1.-Mr[0][0])+alpha.ten[1][1]*(-Mr[1][0]));
       a22=a22+Mh.ljr(numCell,r)*rjr*(alpha.ten[1][0]*(-Mr[0][1])+alpha.ten[1][1]*(1.-Mr[1][1]));

     }
  }
  a11=(dt/(Mh.area(numCell)*Param.Euler.eps[numCell]*rho))*a11+1.;
  a12=(dt/(Mh.area(numCell)*Param.Euler.eps[numCell]*rho))*a12;
  a21=(dt/(Mh.area(numCell)*Param.Euler.eps[numCell]*rho))*a21;
  a22=(dt/(Mh.area(numCell)*Param.Euler.eps[numCell]*rho))*a22+1.;

  
  
  Det=a11*a22-a21*a12;
  if(Det==0){ cout<<"The implicit matrix for the euler scheme is not invertible"<<endl;exit(1);}
  else { 
    res.x=(b1*a22-b2*a12)/Det;
    res.y=(b2*a11-b1*a21)/Det;
  }
 

 
  return res;
}

R2 FrictionU(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 * ur){
/** Function wich compute the nodal friction term for momemtum equation **/
  R2 res(0,0);
  tensor alpha(d,'h');
  R2 sol(0,0);
  int numGr=0;
  double rjr=0;
  double Mr[2][2];
  double a11=0,a12=0,a21=0,a22=0;
  double b1=0,b2=0;
  R2 uj(0,0);

  uj.x=v.var[1][numCell]/v.var[0][numCell];
  uj.y=v.var[2][numCell]/v.var[0][numCell];
  a11=0; a12=0; a21=0; a22=0;

    for(int r=0;r<Mh.nbnodelocal;r++){
       numGr=Mh(numCell,r);
       MrConstructE(d,numGr,Mh,v,tab,Param,Mr);
      sol=ur[numGr];
      rjr=WaveSpeed(d,Mh,v,numCell,numGr,tab,Param.Euler);
      alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
      a11=a11+Mh.ljr(numCell,r)*rjr*(alpha.ten[0][0]*(1-Mr[0][0])+alpha.ten[0][1]*(-Mr[1][0]));
      a12=a12+Mh.ljr(numCell,r)*rjr*(alpha.ten[0][0]*(-Mr[0][1])+alpha.ten[0][1]*(1-Mr[1][1]));
      a21=a21+Mh.ljr(numCell,r)*rjr*(alpha.ten[1][0]*(1-Mr[0][0])+alpha.ten[1][1]*(-Mr[1][0]));
      a22=a22+Mh.ljr(numCell,r)*rjr*(alpha.ten[1][0]*(-Mr[0][1])+alpha.ten[1][1]*(1-Mr[1][1]));
    }  
      b1=a11*uj.x+a12*uj.y;
      b2=a21*uj.x+a22*uj.y;
        
      res.x=-b1/Param.Euler.eps[numCell];
      res.y=-b2/Param.Euler.eps[numCell];
  return res;
 
}

double FrictionUcarre(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 * ur){
/** Function wich compute the nodal friction term for energy equation **/
double res=0;
  tensor alpha(d,'h');
  R2 g(0,0);
  R2 sol(0,0);
  g=GravityVector(d,Param.Euler);
  int numGr=0;
  double rjr=0;
  double Mr[2][2];
  double a11=0,a12=0,a21=0,a22=0;
  double c1=0,c2=0;
  double b1=0,b2=0;
  R2 uj(0,0);

  uj.x=v.var[1][numCell]/v.var[0][numCell];
  uj.y=v.var[2][numCell]/v.var[0][numCell];

    for(int r=0;r<Mh.nbnodelocal;r++){
       numGr=Mh(numCell,r);
       MrConstructE(d,numGr,Mh,v,tab,Param,Mr);
      sol=ur[numGr];
      rjr=WaveSpeed(d,Mh,v,numCell,numGr,tab,Param.Euler);
      alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
      a11=Mh.ljr(numCell,r)*rjr*(alpha.ten[0][0]*(1-Mr[0][0])+alpha.ten[0][1]*(-Mr[1][0]));
      a12=Mh.ljr(numCell,r)*rjr*(alpha.ten[0][0]*(-Mr[0][1])+alpha.ten[0][1]*(1-Mr[1][1]));
      a21=Mh.ljr(numCell,r)*rjr*(alpha.ten[1][0]*(1-Mr[0][0])+alpha.ten[1][1]*(-Mr[1][0]));
      a22=Mh.ljr(numCell,r)*rjr*(alpha.ten[1][0]*(-Mr[0][1])+alpha.ten[1][1]*(1-Mr[1][1]));
    
      b1=a11*uj.x+a12*uj.y;
      b2=a21*uj.x+a22*uj.y;
      
      c1=Mr[0][0]*sol.x+Mr[0][1]*sol.y;
      c2=Mr[1][0]*sol.x+Mr[1][1]*sol.y;
        
      res=res-(b1*c1+b2*c2)/Param.Euler.eps[numCell];
    }

  return res;
 
}

R2 SourceGravityEuler(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param){
  /** Function wich compute the nodal gravity source term (momemtum equation) for Euler equations  **/
  R2 res(0,0);
  tensor alpha(d,'h');
  R2 rhog(0,0);
  R2 g1(0,0);
  int numGr=0;
  double rjr=0;
  double Mr[2][2];
  double Nr[2][2];
  double a11=0,a12=0,a21=0,a22=0;
  double c11=0,c12=0,c21=0,c22=0;


    for(int r=0;r<Mh.nbnodelocal;r++){
      numGr=Mh(numCell,r);
      MrConstructE(d,numGr,Mh,v,tab,Param,Mr);
      MatrixSourceWBE(d,numGr,Mh,v,tab,Param,Nr);

      alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
      rhog=Param.Euler.TabRhoG[numGr];
      rjr=WaveSpeed(d,Mh,v,numCell,numGr,tab,Param.Euler);
      a11=Mh.ljr(numCell,r)*rjr*(alpha.ten[0][0]*Mr[0][0]+alpha.ten[0][1]*Mr[1][0]);
      a12=Mh.ljr(numCell,r)*rjr*(alpha.ten[0][0]*Mr[0][1]+alpha.ten[0][1]*Mr[1][1]);
      a21=Mh.ljr(numCell,r)*rjr*(alpha.ten[1][0]*Mr[0][0]+alpha.ten[1][1]*Mr[1][0]);
      a22=Mh.ljr(numCell,r)*rjr*(alpha.ten[1][0]*Mr[0][1]+alpha.ten[1][1]*Mr[1][1]);

      c11=(a11*Nr[0][0]+a12*Nr[1][0]);
      c12=(a11*Nr[0][1]+a12*Mr[1][1]);
      c21=(a21*Nr[0][0]+a22*Nr[1][0]);
      c22=(a21*Nr[0][1]+a22*Nr[1][1]);

      g1.x=c11*rhog.x+c12*rhog.y;
      g1.y=c21*rhog.x+c22*rhog.y;

      res.x=res.x-g1.x;
      res.y=res.y-g1.y;

    }

  return res;
}

double SourceGravityEulerU(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 * ur){
   /** Function wich compute the nodal gravity source term (energy equation) for Euler equations  **/
  double res=0;
  tensor alpha(d,'h');
  R2 rhog(0,0);
  R2 sol(0,0);
  R2 g1(0,0);
  R2 g2(0,0);
  int numGr=0;
  double rjr=0;
  double Mr[2][2];
  double Nr[2][2];
  double a11=0,a12=0,a21=0,a22=0;
  double c11=0,c12=0,c21=0,c22=0;


    for(int r=0;r<Mh.nbnodelocal;r++){
      numGr=Mh(numCell,r);
      sol=ur[numGr];
      MrConstructE(d,numGr,Mh,v,tab,Param,Mr);
      MatrixSourceWBE(d,numGr,Mh,v,tab,Param,Nr);
      alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);

      rhog=Param.Euler.TabRhoG[numGr];
      rjr=WaveSpeed(d,Mh,v,numCell,numGr,tab,Param.Euler);
      a11=Mh.ljr(numCell,r)*rjr*(alpha.ten[0][0]*Mr[0][0]+alpha.ten[0][1]*Mr[1][0]);
      a12=Mh.ljr(numCell,r)*rjr*(alpha.ten[0][0]*Mr[0][1]+alpha.ten[0][1]*Mr[1][1]);
      a21=Mh.ljr(numCell,r)*rjr*(alpha.ten[1][0]*Mr[0][0]+alpha.ten[1][1]*Mr[1][0]);
      a22=Mh.ljr(numCell,r)*rjr*(alpha.ten[1][0]*Mr[0][1]+alpha.ten[1][1]*Mr[1][1]);

      c11=(a11*Nr[0][0]+a12*Nr[1][0]);
      c12=(a11*Nr[0][1]+a12*Mr[1][1]);
      c21=(a21*Nr[0][0]+a22*Nr[1][0]);
      c22=(a21*Nr[0][1]+a22*Nr[1][1]);

      g1.x=c11*rhog.x+c12*rhog.y;
      g1.y=c21*rhog.x+c22*rhog.y;

      g2.x=Mr[0][0]*sol.x+Mr[0][1]*sol.y;
      g2.y=Mr[1][0]*sol.x+Mr[1][1]*sol.y;


      res=res-(g1,g2);

    }

  return res;
}

void MatrixSourceWBE(Data & d,int numGr,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,double  Mrsource[2][2]){
   /** Function wich compute the nodal matrix Nr useful for wb scheme for Euler equations  **/
  double a11=0,a12=0,a21=0,a22=0,Det=0,temp;
 double c11=0,c12=0,c21=0,c22=0;
 int numLrj=0,jG=0;
 tensor alpha(d,'h');
 tensor beta(d,'s');
 double rjr=0;


 for(int j=0;j<tab.TabInv[numGr].taille;j++){
      jG=tab.TabInv[numGr].TabCell[j];
      numLrj=NodeGtoL(Mh,jG,numGr);  
      

      beta=inittensor(d,Mh,v,tab,'s',numGr,jG);
      alpha=inittensor(d,Mh,v,tab,'h',numGr,jG);

      rjr=WaveSpeed(d,Mh,v,jG,numGr,tab,Param.Euler);
    
       c11=c11+Mh.ljr(jG,numLrj)*rjr*alpha.ten[0][0];
       a11=a11+Mh.ljr(jG,numLrj)*beta.ten[0][0];
       c12=c12+Mh.ljr(jG,numLrj)*rjr*alpha.ten[0][1];
       a12=a12+Mh.ljr(jG,numLrj)*beta.ten[0][1];
       c21=c21+Mh.ljr(jG,numLrj)*rjr*alpha.ten[1][0];
       a21=a21+Mh.ljr(jG,numLrj)*beta.ten[1][0];
       c22=c22+Mh.ljr(jG,numLrj)*rjr*alpha.ten[1][1];
       a22=a22+Mh.ljr(jG,numLrj)*beta.ten[1][1];      
 }
 Det=c11*c22-c21*c12;
 temp=c11;

 if(Mh.xr(numGr).lab < 0){
   c11=1.;
   c12=0;
   c21=0;
   c22=1.;
 }
 else{
   c11=c22/Det;
   c12=-c12/Det;
   c21=-c21/Det;
   c22=temp/Det;
 }
 Mrsource[0][0]=c11*a11+c12*a21;
 Mrsource[0][1]=c11*a12+c12*a22;
 Mrsource[1][0]=c21*a11+c22*a21;
 Mrsource[1][1]=c21*a12+c22*a22;

}







