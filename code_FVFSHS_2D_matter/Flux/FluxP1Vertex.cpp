#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Variables.hpp"
#include "Flux.hpp"
#include "Source.hpp"
#include "BoundaryCondition.hpp"
#include "Tensor.hpp"
#include "FunctionsAdvection.hpp"




vectorflux FluxVertexP1(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 * ur){
  /**  Nodal JL-(b) and JL-(a) fluxes for the P1 model **/
    vectorflux res(3);
  double Fu1=0,Fu2=0;
  R2 sol(0,0);
  double s[3]={0,0,0};
  tensor alpha(d,'h');
  tensor beta(d,'s');
  int numGr=0;
  int z=0;
  double signode=0;
  double epsnode=0;
  double var0jr=0, var1jr=0, var2jr=0;

  if(d.scheme==1){
    z=1;} /** JL-(a) scheme**/
  if(d.scheme==2){
    z=0;} /** JL-(b) scheme. The part which depend of the source term is killed by the nodal source term **/

  for(int r=0;r<Mh.nbnodelocal;r++){
   
    numGr=Mh(numCell,r);
    sol=ur[numGr];
    
     beta=inittensor(d,Mh,v,tab,'s',numGr,numCell);
     alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
     
     signode=AverageSigNodeP1(d,Mh,v,tab,Param.P1,numGr,numCell);
     epsnode=AverageEpsNodeP1(d,Mh,v,tab,Param.P1,numGr,numCell);

     var0jr=InterLineVertex(d,Mh,v,tab,Param,0,d.Total_order,numCell,r);
     var1jr=InterLineVertex(d,Mh,v,tab,Param,1,d.Total_order,numCell,r);
     var2jr=InterLineVertex(d,Mh,v,tab,Param,2,d.Total_order,numCell,r);
       
     Fu1=var0jr*Mh.njr(numCell,r).x+alpha.ten[0][0]*(var1jr-sol.x)+alpha.ten[0][1]*(var2jr-sol.y)
       -z*(signode/epsnode)*(sol.x*beta.ten[0][0]+sol.y*beta.ten[0][1]);
     
     Fu2=var0jr*Mh.njr(numCell,r).y+alpha.ten[1][0]*(var1jr-sol.x)+alpha.ten[1][1]*(var2jr-sol.y)
       -z*(signode/epsnode)*(sol.x*beta.ten[1][0]+sol.y*beta.ten[1][1]);

     // Fu1=v.var[0][numCell]*Mh.njr(numCell,r).x+alpha.ten[0][0]*(v.var[1][numCell]-sol.x)+alpha.ten[0][1]*(v.var[2][numCell]-sol.y)
     //  -z*(signode/epsnode)*(sol.x*beta.ten[0][0]+sol.y*beta.ten[0][1]);
     
     //Fu2=v.var[0][numCell]*Mh.njr(numCell,r).y+alpha.ten[1][0]*(v.var[1][numCell]-sol.x)+alpha.ten[1][1]*(v.var[2][numCell]-sol.y)
      // -z*(signode/epsnode)*(sol.x*beta.ten[1][0]+sol.y*beta.ten[1][1]);
    

     s[0]=s[0]+Mh.ljr(numCell,r)*(sol.x*Mh.njr(numCell,r).x+sol.y*Mh.njr(numCell,r).y);
     s[1]=s[1]+Mh.ljr(numCell,r)*Fu1;
     s[2]=s[2]+Mh.ljr(numCell,r)*Fu2;
      
  }
 
  res.vflux[0]=-s[0]/Param.P1.eps[numCell];
  res.vflux[1]=-s[1]/Param.P1.eps[numCell];
  res.vflux[2]=-s[2]/Param.P1.eps[numCell];
  return res;
}



vectorflux FluxVertexP1Gosse(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 * ur){
   /**  Nodal JL-(b) fluxes with local source term for the P1 model **/
    vectorflux res(3);
  double Fu1=0,Fu2=0;
  R2 sol(0,0);
  double s[3]={0,0,0};
  tensor alpha(d,'h');
  tensor beta(d,'s');
  int numGr=0;
    double signode=0;
  double epsnode=0;
  double Mr[2][2];
  double temp=0;
 double var0jr=0, var1jr=0, var2jr=0;
  

  
  for(int r=0;r<Mh.nbnodelocal;r++){
      
    numGr=Mh(numCell,r);
    MrConstructP1(d,numGr,Mh,v,tab,Param,Mr);
    
    sol=ur[numGr];

     
    var0jr=InterLineVertex(d,Mh,v,tab,Param,0,d.Total_order,numCell,r);
     var1jr=InterLineVertex(d,Mh,v,tab,Param,1,d.Total_order,numCell,r);
     var2jr=InterLineVertex(d,Mh,v,tab,Param,2,d.Total_order,numCell,r);
    
    alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
    signode=AverageSigNodeP1(d,Mh,v,tab,Param.P1,numGr,numCell);
    epsnode=AverageEpsNodeP1(d,Mh,v,tab,Param.P1,numGr,numCell);
      
           // calcul du flux ujr      
    Fu1=var0jr*Mh.njr(numCell,r).x+(alpha.ten[0][0]*Mr[0][0]+alpha.ten[0][1]*Mr[1][0])*(var1jr-sol.x)+(alpha.ten[0][0]*Mr[0][1]+alpha.ten[0][1]*Mr[1][1])*(var2jr-sol.y);
    
    Fu2=var0jr*Mh.njr(numCell,r).y+(alpha.ten[1][0]*Mr[0][0]+alpha.ten[1][1]*Mr[1][0])*(var1jr-sol.x)+(alpha.ten[1][0]*Mr[0][1]+alpha.ten[1][1]*Mr[1][1])*(var2jr-sol.y);
           // sommation des flux
    temp=sol.x;
    sol.x=Mr[0][0]*sol.x+Mr[0][1]*sol.y;
    sol.y=Mr[1][0]*temp+Mr[1][1]*sol.y;
      
    s[0]=s[0]+Mh.ljr(numCell,r)*(sol.x*Mh.njr(numCell,r).x+sol.y*Mh.njr(numCell,r).y);
    s[1]=s[1]+Mh.ljr(numCell,r)*Fu1;
    s[2]=s[2]+Mh.ljr(numCell,r)*Fu2;

  }  
  
  res.vflux[0]=-s[0]/Param.P1.eps[numCell];
  res.vflux[1]=-s[1]/Param.P1.eps[numCell];
  res.vflux[2]=-s[2]/Param.P1.eps[numCell];
	
  return res;
}

void MatrixP1(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr,double &a11,double &a12,double &a21,double &a22,double &b1,double &b2){
   /**  Nodal matrix for the JL-(b) and JL-(a) schemes for the P1 model **/
  tensor alpha(d,'h');
  tensor beta(d,'s');
  int numLrj=0, jG=0; //  numero globale de r
  double signode=0;
  double epsnode=0;
  double var0jr=0,var1jr=0,var2jr=0;
  
  for(int j=0;j<tab.TabInv[numGr].taille;j++){
    jG=tab.TabInv[numGr].TabCell[j];
    numLrj=NodeGtoL(Mh,jG,numGr);  

    var0jr=InterLineVertex(d,Mh,v,tab,Param,0,d.Total_order,jG,numLrj);
    var1jr=InterLineVertex(d,Mh,v,tab,Param,1,d.Total_order,jG,numLrj);
    var2jr=InterLineVertex(d,Mh,v,tab,Param,2,d.Total_order,jG,numLrj);
    
    beta=inittensor(d,Mh,v,tab,'s',numGr,jG);
    alpha=inittensor(d,Mh,v,tab,'h',numGr,jG);//numero local de r dans j    
    signode=AverageSigNodeP1(d,Mh,v,tab,Param.P1,numGr,jG);
    epsnode=AverageEpsNodeP1(d,Mh,v,tab,Param.P1,numGr,jG);
    if(d.scheme==3 ){
      signode=0.;
    }
   
    
    a11=a11+Mh.ljr(jG,numLrj)*(alpha.ten[0][0]+(signode/epsnode)*beta.ten[0][0]);
    a12=a12+Mh.ljr(jG,numLrj)*(alpha.ten[0][1]+(signode/epsnode)*beta.ten[0][1]);
    a21=a21+Mh.ljr(jG,numLrj)*(alpha.ten[1][0]+(signode/epsnode)*beta.ten[1][0]);
    a22=a22+Mh.ljr(jG,numLrj)*(alpha.ten[1][1]+(signode/epsnode)*beta.ten[1][1]);      
    
    b1=b1+Mh.ljr(jG,numLrj)*((var0jr)*Mh.njr(jG,numLrj).x+alpha.ten[0][0]*var1jr+alpha.ten[0][1]*var2jr);
    b2=b2+Mh.ljr(jG,numLrj)*((var0jr)*Mh.njr(jG,numLrj).y+alpha.ten[1][0]*var1jr+alpha.ten[1][1]*var2jr);   
  }

}






void MrConstructP1(Data & d,int numGr,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,double Mr[2][2]){
    /**  Nodal Magical Matrix Mr for the JL-(b) and JL-(a) schemes for the P1 model **/
  double a11=0,a12=0,a21=0,a22=0,Det=0,temp;
 double c11=0,c12=0,c21=0,c22=0;
 int numLrj=0,jG=0;
    double signode=0;
  double epsnode=0;
 tensor alpha(d,'h');
 tensor beta(d,'s');



 for(int j=0;j<tab.TabInv[numGr].taille;j++){
      jG=tab.TabInv[numGr].TabCell[j];
      numLrj=NodeGtoL(Mh,jG,numGr);  
     
      
      beta=inittensor(d,Mh,v,tab,'s',numGr,jG);
      alpha=inittensor(d,Mh,v,tab,'h',numGr,jG);//numero local de r dans j    
      signode=AverageSigNodeP1(d,Mh,v,tab,Param.P1,numGr,jG);
      epsnode=AverageEpsNodeP1(d,Mh,v,tab,Param.P1,numGr,jG);
       
      /////////////////////////////////////////////////////////Construction de la matrice //////////////////////////////////////////////////////////////
	
       c11=c11+Mh.ljr(jG,numLrj)*alpha.ten[0][0];
       a11=a11+(signode/epsnode)*Mh.ljr(jG,numLrj)*beta.ten[0][0]+Mh.ljr(jG,numLrj)*alpha.ten[0][0];
       c12=c12+Mh.ljr(jG,numLrj)*alpha.ten[0][1];
       a12=a12+(signode/epsnode)*Mh.ljr(jG,numLrj)*beta.ten[0][1]+Mh.ljr(jG,numLrj)*alpha.ten[0][1];
       c21=c21+Mh.ljr(jG,numLrj)*alpha.ten[1][0];
       a21=a21+(signode/epsnode)*Mh.ljr(jG,numLrj)*beta.ten[1][0]+Mh.ljr(jG,numLrj)*alpha.ten[1][0];
       c22=c22+Mh.ljr(jG,numLrj)*alpha.ten[1][1];
       a22=a22+(signode/epsnode)*Mh.ljr(jG,numLrj)*beta.ten[1][1]+Mh.ljr(jG,numLrj)*alpha.ten[1][1];      
 }
 Det=a11*a22-a21*a12;
 temp=a11;

 if(Mh.xr(numGr).lab < 0){
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










