#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Variables.hpp"
#include "Flux.hpp"
#include "Source.hpp"
#include "BoundaryCondition.hpp"
#include "Tensor.hpp"
#include "FunctionsAdvection.hpp"


vectorflux FluxVertexP1Matter(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 * ur){
  /**  Nodal flux for the nodal JL-(b) scheme for the P1 model with matter **/
    vectorflux res(4);
  double Fu1=0,Fu2=0;
  R2 sol(0,0);
  double s[3]={0,0,0};
  tensor alpha(d,'h');
  tensor beta(d,'s');
  int numGr=0;
  int z=0;
  double sigAnode=0;
  double epsnode=0;

  if(d.scheme==1){
    z=1;} /** JL-(a) scheme**/
  if(d.scheme==2){
    z=0;}  /** JL-(b) scheme. The part which depend of the source term is killed by the nodal source term **/

  for(int r=0;r<Mh.nbnodelocal;r++){
    
    numGr=Mh(numCell,r);
    sol=ur[numGr];

    beta=inittensor(d,Mh,v,tab,'s',numGr,numCell);
    alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
    sigAnode=AverageSigANodeP1Matter(d,Mh,v,tab,Param.P1M,numGr,numCell);
    epsnode=AverageEpsNodeP1Matter(d,Mh,v,tab,Param.P1M,numGr,numCell);
    
    Fu1=v.var[0][numCell]*Mh.njr(numCell,r).x+alpha.ten[0][0]*(v.var[1][numCell]-sol.x)+alpha.ten[0][1]*(v.var[2][numCell]-sol.y)
      -z*(sigAnode/epsnode)*(sol.x*beta.ten[0][0]+sol.y*beta.ten[0][1]);

    Fu2=v.var[0][numCell]*Mh.njr(numCell,r).y+alpha.ten[1][0]*(v.var[1][numCell]-sol.x)+alpha.ten[1][1]*(v.var[2][numCell]-sol.y)
      -z*(sigAnode/epsnode)*(sol.x*beta.ten[1][0]+sol.y*beta.ten[1][1]);
    
     s[0]=s[0]+Mh.ljr(numCell,r)*(sol.x*Mh.njr(numCell,r).x+sol.y*Mh.njr(numCell,r).y);
     s[1]=s[1]+Mh.ljr(numCell,r)*Fu1;
     s[2]=s[2]+Mh.ljr(numCell,r)*Fu2;
      
  }
 
  res.vflux[0]=-s[0]/Param.P1M.eps[numCell];
  res.vflux[1]=-s[1]/Param.P1M.eps[numCell];
  res.vflux[2]=-s[2]/Param.P1M.eps[numCell];
  res.vflux[3]=0;
  return res;
}



vectorflux FluxVertexP1MatterGosse(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 * ur){
  /**  Nodal flux for the nodal JL-(b) scheme with local source term for the P1 model with matter **/
    vectorflux res(4);
  double Fu1=0,Fu2=0;
  R2 sol(0,0);
  double s[3]={0,0,0};
  tensor alpha(d,'h');
  tensor beta(d,'s');
  int numGr=0;
    double sigAnode=0;
  double epsnode=0;
  double Mr[2][2];
  double temp=0;

  for(int r=0;r<Mh.nbnodelocal;r++){
      
    numGr=Mh(numCell,r);
    MrConstructP1Matter(d,numGr,Mh,v,tab,Param,Mr);
    
    sol=ur[numGr];   
    
    alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
    sigAnode=AverageSigANodeP1Matter(d,Mh,v,tab,Param.P1M,numGr,numCell);
    epsnode=AverageEpsNodeP1Matter(d,Mh,v,tab,Param.P1M,numGr,numCell);
         
    Fu1=v.var[0][numCell]*Mh.njr(numCell,r).x+(alpha.ten[0][0]*Mr[0][0]+alpha.ten[0][1]*Mr[1][0])*(v.var[1][numCell]-sol.x)+(alpha.ten[0][0]*Mr[0][1]+alpha.ten[0][1]*Mr[1][1])*(v.var[2][numCell]-sol.y);
    
    Fu2=v.var[0][numCell]*Mh.njr(numCell,r).y+(alpha.ten[1][0]*Mr[0][0]+alpha.ten[1][1]*Mr[1][0])*(v.var[1][numCell]-sol.x)+(alpha.ten[1][0]*Mr[0][1]+alpha.ten[1][1]*Mr[1][1])*(v.var[2][numCell]-sol.y);

    temp=sol.x;
    sol.x=Mr[0][0]*sol.x+Mr[0][1]*sol.y;
    sol.y=Mr[1][0]*temp+Mr[1][1]*sol.y;
      
    s[0]=s[0]+Mh.ljr(numCell,r)*(sol.x*Mh.njr(numCell,r).x+sol.y*Mh.njr(numCell,r).y);
    s[1]=s[1]+Mh.ljr(numCell,r)*Fu1;
    s[2]=s[2]+Mh.ljr(numCell,r)*Fu2;

  }  
  
  res.vflux[0]=-s[0]/Param.P1M.eps[numCell];
  res.vflux[1]=-s[1]/Param.P1M.eps[numCell];
  res.vflux[2]=-s[2]/Param.P1M.eps[numCell];
  res.vflux[3]=0;	
  return res;
}


void MatrixP1Matter(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr,double &a11,double &a12,double &a21,double &a22,double &b1,double &b2){
  /**  Nodal matrix for the JL-(b) and JL-(a) schemes for the P1 model with matter **/
  
  tensor alpha(d,'h');
  tensor beta(d,'s');
  int numLrj=0, jG=0; //  numero globale de r
  double sigAnode=0;
  double epsnode=0;

  
  for(int j=0;j<tab.TabInv[numGr].taille;j++){
    jG=tab.TabInv[numGr].TabCell[j];
    numLrj=NodeGtoL(Mh,jG,numGr);  

    beta=inittensor(d,Mh,v,tab,'s',numGr,jG);
    alpha=inittensor(d,Mh,v,tab,'h',numGr,jG);//numero local de r dans j    
    sigAnode=AverageSigANodeP1Matter(d,Mh,v,tab,Param.P1M,numGr,jG);
    epsnode=AverageEpsNodeP1Matter(d,Mh,v,tab,Param.P1M,numGr,jG);
    if(d.scheme==3 ){
      sigAnode=0.;
    }
   
    
    a11=a11+Mh.ljr(jG,numLrj)*(alpha.ten[0][0]+(sigAnode/epsnode)*beta.ten[0][0]);
    a12=a12+Mh.ljr(jG,numLrj)*(alpha.ten[0][1]+(sigAnode/epsnode)*beta.ten[0][1]);
    a21=a21+Mh.ljr(jG,numLrj)*(alpha.ten[1][0]+(sigAnode/epsnode)*beta.ten[1][0]);
    a22=a22+Mh.ljr(jG,numLrj)*(alpha.ten[1][1]+(sigAnode/epsnode)*beta.ten[1][1]);      
    
    b1=b1+Mh.ljr(jG,numLrj)*((v.var[0][jG])*Mh.njr(jG,numLrj).x+alpha.ten[0][0]*v.var[1][jG]+alpha.ten[0][1]*v.var[2][jG]);
    b2=b2+Mh.ljr(jG,numLrj)*((v.var[0][jG])*Mh.njr(jG,numLrj).y+alpha.ten[1][0]*v.var[1][jG]+alpha.ten[1][1]*v.var[2][jG]);   
  }

}


void MrConstructP1Matter(Data & d,int numGr,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,double Mr[2][2]){
/**  Nodal Magical Matrix Mr for the JL-(b) and JL-(a) schemes for the P1 model with matter **/
  
  double a11=0,a12=0,a21=0,a22=0,Det=0,temp;
  double c11=0,c12=0,c21=0,c22=0;
  int numLrj=0,jG=0;
  double sigAnode=0;
  double epsnode=0;
  tensor alpha(d,'h');
  tensor beta(d,'s');


 for(int j=0;j<tab.TabInv[numGr].taille;j++){
      jG=tab.TabInv[numGr].TabCell[j];
      numLrj=NodeGtoL(Mh,jG,numGr);  
      

      beta=inittensor(d,Mh,v,tab,'s',numGr,jG);
      alpha=inittensor(d,Mh,v,tab,'h',numGr,jG);   
      sigAnode=AverageSigANodeP1Matter(d,Mh,v,tab,Param.P1M,numGr,jG);
      epsnode=AverageEpsNodeP1Matter(d,Mh,v,tab,Param.P1M,numGr,jG);
	
       c11=c11+Mh.ljr(jG,numLrj)*alpha.ten[0][0];
       a11=a11+(sigAnode/epsnode)*Mh.ljr(jG,numLrj)*beta.ten[0][0]+Mh.ljr(jG,numLrj)*alpha.ten[0][0];
       c12=c12+Mh.ljr(jG,numLrj)*alpha.ten[0][1];
       a12=a12+(sigAnode/epsnode)*Mh.ljr(jG,numLrj)*beta.ten[0][1]+Mh.ljr(jG,numLrj)*alpha.ten[0][1];
       c21=c21+Mh.ljr(jG,numLrj)*alpha.ten[1][0];
       a21=a21+(sigAnode/epsnode)*Mh.ljr(jG,numLrj)*beta.ten[1][0]+Mh.ljr(jG,numLrj)*alpha.ten[1][0];
       c22=c22+Mh.ljr(jG,numLrj)*alpha.ten[1][1];
       a22=a22+(sigAnode/epsnode)*Mh.ljr(jG,numLrj)*beta.ten[1][1]+Mh.ljr(jG,numLrj)*alpha.ten[1][1];      
 }
 Det=a11*a22-a21*a12;
 temp=a11;

 if(Mh.xr(numGr).lab <0){
   a11=1.;
   a12=0;
   a21=0;
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
