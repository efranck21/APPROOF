
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

/** Flux du schéma de diffusion GLACE linéaire**/
vectorflux FluxVertexDiff(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 ** ur){
   /**  Nodal fluxes for linear diffusion scheme   **/
  vectorflux res(1);
  R2 sol(0,0);
  int numGr;
  double s=0;
  R2 stab(0,0);
    
  for(int r=0;r<Mh.nbnodelocal;r++){
     numGr=Mh(numCell,r);

     sol=ur[0][numGr];  

     s=s+Mh.ljr(numCell,r)*(sol,Mh.njr(numCell,r));
    
  }
    res.vflux[0]=-s;
    return res;	
}

vectorflux FluxVertexNlDiff(Data & d,int numCell,Mesh &  Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 ** ur){
   /**  Nodal fluxes for nonlinear diffusion scheme based on advection   **/
 vectorflux res(1);
 R2 sol(0,0);
  int numGr;
  double s=0;  
  double cd=0;
  tensor beta(d,'s');
  double Er;

    for(int r=0;r<Mh.nbnodelocal;r++){
     numGr=Mh(numCell,r);
     
     cd=AverageCoefDiffNode(d,Mh,v,tab,Param.Diff,numGr,numCell);   
     
     Er=AverageQuantity(d,Mh,v,tab,Param,numGr,0);
     sol=ur[0][numGr];
       
     s=s+VertexUpwind(d,Mh,v,0,tab,Param,numCell,r,sol/Er);
  }

    res.vflux[0]=-s;
    return res;	
}


void MatrixDiff(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr,double  & a11,double &a12,double &a21,double &a22,double &b1,double &b2){
   /**  Nodal matrix for nodal diffusion scheme   **/
  int jG=0,numLrj=0;
  double cd=0;
  tensor beta(d,'s');
  b1 =0;
  b2 =0;

  for(int j=0;j<tab.TabInv[numGr].taille;j++){
    jG=tab.TabInv[numGr].TabCell[j];   
    numLrj=NodeGtoL(Mh,jG,numGr);  
    
    beta=inittensor(d,Mh,v,tab,'s',numGr,jG);
    cd=AverageCoefDiffNode(d,Mh,v,tab,Param.Diff,numGr,jG);
    
    a11=a11+cd*Mh.ljr(jG,numLrj)*beta.ten[0][0];
    a12=a12+cd*Mh.ljr(jG,numLrj)*beta.ten[0][1];
    a21=a21+cd*Mh.ljr(jG,numLrj)*beta.ten[1][0];
    a22=a22+cd*Mh.ljr(jG,numLrj)*beta.ten[1][1]; 
 
    b1=b1+Mh.ljr(jG,numLrj)*(v.var[0][jG])*Mh.njr(jG,numLrj).x;
    b2=b2+Mh.ljr(jG,numLrj)*(v.var[0][jG])*Mh.njr(jG,numLrj).y;
  }

}



