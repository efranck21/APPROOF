
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
#include "Source.hpp"





vectorflux FluxVertexEulerGosse(Data & d,int numCell,Mesh & Mh,variable & v,TabConnecInv & tab,ParamPhysic & Param,R2 * ur){
  /**  Nodal Lagrange+remap JL-(b) fluxes with local source term for the Euler model **/
  vectorflux res(4);
  R2 sol(0.,0.);
  R2 sol2(0.,0.);
  double s[4]={0.,0.,0.,0.};
  tensor alpha(d,'h');
  tensor beta(d,'s');
  int numGr;
  double rjr=0;
  R2 Gjr(0,0);
  double pj=0;
  R2 uj(0,0);
  double Mr[2][2];
  double temp=0;

  Mr[0][0]=0;
  Mr[0][1]=0;
  Mr[1][0]=0;
  Mr[1][1]=0;

  uj.x=v.var[1][numCell]/v.var[0][numCell];
  uj.y=v.var[2][numCell]/v.var[0][numCell];

    for(int r=0;r<Mh.nbnodelocal;r++){
    // calcul du flux nodale vr

     numGr=Mh(numCell,r);
     MrConstructE(d,numGr,Mh,v,tab,Param,Mr);
     sol=ur[numGr];
     sol2=ClassicalSolver(d,Mh,v,tab,Param,numGr);
  
     beta=inittensor(d,Mh,v,tab,'s',numGr,numCell);
     alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
     rjr=WaveSpeed(d,Mh,v,numCell,numGr,tab,Param.Euler);
      
     pj=PressureLaw(d,Mh,v,Param.Euler,numCell);


     Gjr.x=Mh.ljr(numCell,r)*(pj*Mh.njr(numCell,r).x+rjr*(alpha.ten[0][0]*Mr[0][0]+alpha.ten[0][1]*Mr[1][0])*(uj.x-sol2.x)+rjr*(alpha.ten[0][0]*Mr[0][1]+alpha.ten[0][1]*Mr[1][1])*(uj.y-sol2.y));
    
     Gjr.y=Mh.ljr(numCell,r)*(pj*Mh.njr(numCell,r).y+rjr*(alpha.ten[1][0]*Mr[0][0]+alpha.ten[1][1]*Mr[1][0])*(uj.x-sol2.x)+rjr*(alpha.ten[1][0]*Mr[0][1]+alpha.ten[1][1]*Mr[1][1])*(uj.y-sol2.y));

  

     temp=sol.x;
     sol.x=Mr[0][0]*sol.x+Mr[0][1]*sol.y;
     sol.y=Mr[1][0]*temp+Mr[1][1]*sol.y;
          


     s[0]=s[0]+remap(d,Mh,v,0,tab,Param,numCell,r,ur,sol);
     s[1]=s[1]+remap(d,Mh,v,1,tab,Param,numCell,r,ur,sol)+Gjr.x;                    
     s[2]=s[2]+remap(d,Mh,v,2,tab,Param,numCell,r,ur,sol)+Gjr.y;
     s[3]=s[3]+remap(d,Mh,v,3,tab,Param,numCell,r,ur,sol)+(Gjr,sol);
    
  }
   
  res.vflux[0]=-s[0]/Param.Euler.eps[numCell];
  res.vflux[1]=-s[1]/Param.Euler.eps[numCell];
  res.vflux[2]=-s[2]/Param.Euler.eps[numCell];
  res.vflux[3]=-s[3]/Param.Euler.eps[numCell];


  return res;
}

vectorflux FluxVertexEuler(Data & d,int numCell,Mesh & Mh,variable & v,TabConnecInv & tab,ParamPhysic & Param,R2 * ur){
  /**  Nodal Lagrange+remap JL-(b) fluxes for the Euler model **/
  vectorflux res(4);
    R2 sol(0,0);
    double s[4]={0.,0.,0.,0.};
  tensor alpha(d,'h');
  tensor beta(d,'s');
  int numGr=0;
  double rjr=0.;
  R2 Gjr(0,0);
  R2 Gjr2(0,0);
  double rhonode=0.,pj=0.;
  double epsnode=0,signode=0; 
  double sign=0.;
  R2 uj(0,0);
  R2 rhog(0,0);


  uj.x=v.var[1][numCell]/v.var[0][numCell];
  uj.y=v.var[2][numCell]/v.var[0][numCell];

  if(d.scheme==1){sign=1.;}
  if(d.scheme==2){sign=0.;}
  if(d.scheme==4){sign=0.;}

  for(int r=0;r<Mh.nbnodelocal;r++){
    // calcul du flux nodale vr

     numGr=Mh(numCell,r);
     sol=ur[numGr];

     beta=inittensor(d,Mh,v,tab,'s',numGr,numCell);
     alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
     rjr=WaveSpeed(d,Mh,v,numCell,numGr,tab,Param.Euler);
     signode=AverageSigNodeEuler(d,Mh,v,tab,Param.Euler,numGr,numCell);
      epsnode=AverageEpsNodeEuler(d,Mh,v,tab,Param.Euler,numGr,numCell);
      
     pj=PressureLaw(d,Mh,v,Param.Euler,numCell);
     rhonode=rhor(d,Mh,v,tab,Param.Euler,numGr);
     rhog=Param.Euler.TabRhoG[numGr];

     Gjr.x=Mh.ljr(numCell,r)*(pj*Mh.njr(numCell,r).x+rjr*(alpha.ten[0][0]*(uj.x-sol.x)+alpha.ten[0][1]*(uj.y-sol.y))-sign*(rhog.x*beta.ten[0][0]+rhog.y*beta.ten[0][1]+rhonode*(signode/epsnode)*(sol.x*beta.ten[0][0]+sol.y*beta.ten[0][1])));
     
     Gjr.y=Mh.ljr(numCell,r)*(pj*Mh.njr(numCell,r).y+rjr*(alpha.ten[1][0]*(uj.x-sol.x)+alpha.ten[1][1]*(uj.y-sol.y))-sign*(rhog.x*beta.ten[1][0]+rhog.y*beta.ten[1][1]+rhonode*(signode/epsnode)*(sol.x*beta.ten[1][0]+sol.y*beta.ten[1][1])));

     Gjr2.x=Mh.ljr(numCell,r)*(pj*Mh.njr(numCell,r).x+rjr*(alpha.ten[0][0]*(uj.x-sol.x)+alpha.ten[0][1]*(uj.y-sol.y))-sign*(rhog.x*beta.ten[0][0]+rhog.y*beta.ten[0][1]+rhonode*(signode/epsnode)*(sol.x*beta.ten[0][0]+sol.y*beta.ten[0][1])));
     
     Gjr2.y=Mh.ljr(numCell,r)*(pj*Mh.njr(numCell,r).y+rjr*(alpha.ten[1][0]*(uj.x-sol.x)+alpha.ten[1][1]*(uj.y-sol.y))-sign*(rhog.x*beta.ten[1][0]+rhog.y*beta.ten[1][1]+rhonode*(signode/epsnode)*(sol.x*beta.ten[1][0]+sol.y*beta.ten[1][1])));
          
    
     s[0]=s[0]+remap(d,Mh,v,0,tab,Param,numCell,r,ur,sol);
     s[1]=s[1]+remap(d,Mh,v,1,tab,Param,numCell,r,ur,sol)+Gjr.x;                    
     s[2]=s[2]+remap(d,Mh,v,2,tab,Param,numCell,r,ur,sol)+Gjr.y;
     s[3]=s[3]+remap(d,Mh,v,3,tab,Param,numCell,r,ur,sol)+(Gjr,sol);

    
  }

       
  res.vflux[0]=-s[0]/Param.Euler.eps[numCell];
  res.vflux[1]=-s[1]/Param.Euler.eps[numCell];
  res.vflux[2]=-s[2]/Param.Euler.eps[numCell];
  res.vflux[3]=-s[3]/Param.Euler.eps[numCell];



  return res;
}






void MatrixEuler(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr,double &a11,double &a12,double &a21,double &a22,double &b1,double &b2){
/**  Nodal matrix for the Lagrange+remap JL-(b) scheme for the Euler model **/
  
  double pj=0.;
  R2 uj(0.,0.);
  int jG=0,numLrj=0;
  tensor alpha(d,'h');
  tensor beta(d,'s');
  double rjr=0.,rhonode=0;
  R2 rhog(0.,0.);
  R2 source(0.,0.);
  double signode=0., epsnode=0.;
 
  for(int j=0;j<tab.TabInv[numGr].taille;j++){
      jG=tab.TabInv[numGr].TabCell[j];
      numLrj=NodeGtoL(Mh,jG,numGr);  
      
      beta=inittensor(d,Mh,v,tab,'s',numGr,jG);
      alpha=inittensor(d,Mh,v,tab,'h',numGr,jG);//numero local de r dans j    
      signode=AverageSigNodeEuler(d,Mh,v,tab,Param.Euler,numGr,jG);
      epsnode=AverageEpsNodeEuler(d,Mh,v,tab,Param.Euler,numGr,jG);
      if(d.scheme==3 || d.scheme==4){
	signode=0;
	}
   
      rjr=WaveSpeed(d,Mh,v,jG,numGr,tab,Param.Euler);
      rhonode=rhor(d,Mh,v,tab,Param.Euler,numGr);
      rhog=Param.Euler.TabRhoG[numGr];
       if(d.scheme==4){
	 rhonode=0.;
	 rhog.x=0.;
	 rhog.y=0.;
       }
      

      a11=a11+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[0][0]+rhonode*(signode/epsnode)*beta.ten[0][0]);
      a12=a12+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[0][1]+rhonode*(signode/epsnode)*beta.ten[0][1]);
      a21=a21+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[1][0]+rhonode*(signode/epsnode)*beta.ten[1][0]);
      a22=a22+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[1][1]+rhonode*(signode/epsnode)*beta.ten[1][1]);

    
      pj=PressureLaw(d,Mh,v,Param.Euler,jG);
      uj.x=v.var[1][jG]/v.var[0][jG];
      uj.y=v.var[2][jG]/v.var[0][jG];

      source.x=source.x-Mh.ljr(jG,numLrj)*(beta.ten[0][0]*rhog.x+beta.ten[0][1]*rhog.y);
      source.y=source.y-Mh.ljr(jG,numLrj)*(beta.ten[1][0]*rhog.x+beta.ten[1][1]*rhog.y);
	
      b1=b1+Mh.ljr(jG,numLrj)*(pj*Mh.njr(jG,numLrj).x+rjr*alpha.ten[0][0]*uj.x+rjr*alpha.ten[0][1]*uj.y);
      b2=b2+Mh.ljr(jG,numLrj)*(pj*Mh.njr(jG,numLrj).y+rjr*alpha.ten[1][0]*uj.x+rjr*alpha.ten[1][1]*uj.y);
       
 
  }
  
      b1=b1+source.x;
      b2=b2+source.y;

}

R2 ClassicalSolver(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamPhysic & Param,int numGr){
/**  Classical Nodal solver for the classical  Lagrange+remap scheme for the Euler model **/
  double a11=0,a12=0,a21=0,a22=0,b1=0,b2=0,Det=0;
  double pj=0;
  R2 uj(0,0);
  int jG=0,numLrj=0;
  tensor alpha(d,'h');
  double rjr=0;
  R2 g(0,0);
  R2 source(0,0);
  R2 sol(0,0);
  

  if(Mh.xr(numGr).lab > -1){

  for(int j=0;j<tab.TabInv[numGr].taille;j++){
      jG=tab.TabInv[numGr].TabCell[j];
      numLrj=NodeGtoL(Mh,jG,numGr);  
      
      alpha=inittensor(d,Mh,v,tab,'h',numGr,jG);    
      rjr=WaveSpeed(d,Mh,v,jG,numGr,tab,Param.Euler);

      a11=a11+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[0][0]);
      a12=a12+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[0][1]);
      a21=a21+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[1][0]);
      a22=a22+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[1][1]);
        
    
      pj=PressureLaw(d, Mh,v,Param.Euler,jG);
      uj.x=v.var[1][jG]/v.var[0][jG];
      uj.y=v.var[2][jG]/v.var[0][jG];

      
      b1=b1+Mh.ljr(jG,numLrj)*(pj*Mh.njr(jG,numLrj).x+rjr*alpha.ten[0][0]*uj.x+rjr*alpha.ten[0][1]*uj.y);
      b2=b2+Mh.ljr(jG,numLrj)*(pj*Mh.njr(jG,numLrj).y+rjr*alpha.ten[1][0]*uj.x+rjr*alpha.ten[1][1]*uj.y);
      
      
  }

  Det=a11*a22-a21*a12;
  if(Det==0){ cout<<"The matrix associated with the classical sovler for euler equation is not invertible"<<endl;exit(1);}
  else { 
    sol.x=(b1*a22-b2*a12)/Det;
    sol.y=(b2*a11-b1*a21)/Det;
  }
  }
  else {
    sol.x=0;
    sol.y=0;
  }
  return sol;
}
