#include <iostream>
#include <cmath>
#include <cassert>
#include "Variables.hpp"
#include "Flux.hpp"
#include "Source.hpp"
#include "BoundaryCondition.hpp"
#include "Tensor.hpp"
#include "FunctionsAdvection.hpp"
#include "ParamPhys.hpp"


/** Calcul de la pression q**/
double q(Data & d, ParamM1 & M1,double E,double F1, double F2){
    /**  function which compute the isotorpic pressure of the M1 model **/
  double res;
  R2 f;
  if(E==0){
    res=0;
  }
  else{
    f.x=F1/(E);
    f.y=F2/(E);
      double khi=(3.+4.*(f,f))/(5.+2.*sqrt(4.-3.*(f,f)));
      res=0.5*(1.-khi)*E;
  }
  return res;
}

/** Calcul du coefficient kr**/
double coefk(Data & d, ParamM1 & M1,double E,double F1, double F2){
    /**  function which compute the coefficienty k define by F=k U, with F the flux and U the velocity in the euler description **/
  double res=0.;
  R2 f(0,0);
 
  if(E==0){
    res=0;
  }
  else{
    f.x=F1/(E);
    f.y=F2/(E);
    if(f.x==0 && f.y==0){
       res=(4./3.)*E;
     }
     
     else{
       double khi=(3.+4.*(f,f))/(5.+2.*sqrt(4.-3.*(f,f)));
         if(khi < (1./3.)+0.000001 && khi >(1./3.)-0.000001){
        res=(4./3.)*E;
       }
       else{
	res=(2.*E*(f,f))/(3.*khi-1.);
       }     
     }
  }
  return res;

}

void  Calcul_u(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab, ParamM1 & M1,R2 & u,int numCell){
    /** Function which compute the velocity U in cell "numcell" **/
  R2 f;
  double khi;

    if(v.var[0][numCell]==0){
      u.x=0;
      u.y=0;
    }
    else{
      f.x=v.var[1][numCell]/(v.var[0][numCell]);
      f.y=v.var[2][numCell]/(v.var[0][numCell]);
       if((f.x==0 && f.y==0)){
	u.x=0;
	u.y=0;
       }
       
      else{
	khi=(3.+4.*(f,f))/(5.+2.*sqrt(4.-3.*(f,f)));
	
	u.x=((3.*khi-1.)*f.x)/(2.*(f,f));
	u.y=((3.*khi-1.)*f.y)/(2.*(f,f));
 
      }
    
    }

}





void MatrixAPM1(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr,double &a11,double &a12,double &a21,double &a22,double &b1,double &b2){
    /**  Nodal matrix for the AP schemes for the M1 model **/
  tensor alpha(d,'h');
  tensor beta(d,'s');
  int numLrj=0, jG=0;  
  double signode=0., epsnode=0.;
  double rjr=0., kr=0.;
  double qj=0.;
  double Er=0.;
  R2  Fr(0,0);
  R2  uj(0,0);

  for(int j=0;j<tab.TabInv[numGr].taille;j++){
      jG=tab.TabInv[numGr].TabCell[j];
      numLrj=NodeGtoL(Mh,jG,numGr);  
      
      beta=inittensor(d,Mh,v,tab,'s',numGr,jG);
      alpha=inittensor(d,Mh,v,tab,'h',numGr,jG);//numero local de r dans j    
      signode=AverageSigNodeM1(d,Mh,v,tab,Param.M1,numGr,jG);
      epsnode=AverageEpsNodeM1(d,Mh,v,tab,Param.M1,numGr,jG);  
      if(d.scheme==3){
	signode=0;
      }
      
      Calcul_u(d,Mh,v,tab,Param.M1,uj,jG);
      qj=q(d,Param.M1,v.var[0][jG],v.var[1][jG],v.var[2][jG]);
      rjr=(4./sqrt(3.))*(v.var[0][jG]/(3+(uj,uj)));
      
      Er=AverageEnergyNodeM1(d,Mh,v,tab,Param.M1,numGr,jG);
      Fr=AverageFluxNodeM1(d,Mh,v,tab,Param.M1,numGr,jG);
      kr=coefk(d,Param.M1,Er,Fr.x,Fr.y);
      a11=a11+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[0][0]+kr*(signode/epsnode)*beta.ten[0][0]);
      a12=a12+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[0][1]+kr*(signode/epsnode)*beta.ten[0][1]);
      a21=a21+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[1][0]+kr*(signode/epsnode)*beta.ten[1][0]);
      a22=a22+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[1][1]+kr*(signode/epsnode)*beta.ten[1][1]);

    
      
      b1=b1+Mh.ljr(jG,numLrj)*(qj*Mh.njr(jG,numLrj).x+rjr*alpha.ten[0][0]*uj.x+rjr*alpha.ten[0][1]*uj.y);
      b2=b2+Mh.ljr(jG,numLrj)*(qj*Mh.njr(jG,numLrj).y+rjr*alpha.ten[1][0]*uj.x+rjr*alpha.ten[1][1]*uj.y);

  }

 
}

void MatrixClassicM1(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr,double &a11,double &a12,double &a21,double &a22,double &b1,double &b2){
    /**  Nodal matrix for the classical scheme for the M1 model **/
  tensor alpha(d,'h');
  tensor beta(d,'s');
  int numLrj=0, jG=0; 
  double rjr=0.;
  double qj=0.;
  R2  uj(0,0);

  for(int j=0;j<tab.TabInv[numGr].taille;j++){
      jG=tab.TabInv[numGr].TabCell[j];
      numLrj=NodeGtoL(Mh,jG,numGr);  
      
      beta=inittensor(d,Mh,v,tab,'s',numGr,jG);
      alpha=inittensor(d,Mh,v,tab,'h',numGr,jG); 
      
      Calcul_u(d,Mh,v,tab,Param.M1,uj,jG);
      qj=q(d,Param.M1,v.var[0][jG],v.var[1][jG],v.var[2][jG]);
      rjr=(4./sqrt(3.))*(v.var[0][jG]/(3+(uj,uj)));

      a11=a11+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[0][0]);
      a12=a12+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[0][1]);
      a21=a21+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[1][0]);
      a22=a22+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[1][1]);  
      
      b1=b1+Mh.ljr(jG,numLrj)*(qj*Mh.njr(jG,numLrj).x+rjr*alpha.ten[0][0]*uj.x+rjr*alpha.ten[0][1]*uj.y);
      b2=b2+Mh.ljr(jG,numLrj)*(qj*Mh.njr(jG,numLrj).y+rjr*alpha.ten[1][0]*uj.x+rjr*alpha.ten[1][1]*uj.y);

  }

 
}

void MatrixM1(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr,double &a11,double &a12,double &a21,double &a22,double &b1,double &b2){ 
  /**  Nodal matrix for the nodal schemes for the M1 model **/

  switch(d.scheme){
    case 1:
      /** Matrix for the classical scheme**/
      MatrixClassicM1(d,Mh,v,tab,Param,numGr,a11,a12,a21,a22,b1,b2);
      break;
    case 2:
       /** Matrix for the AP scheme**/
      MatrixAPM1(d,Mh,v,tab,Param,numGr,a11,a12,a21,a22,b1,b2);
      break;
    case 3:
       /** Matrix for the AP scheme**/
      MatrixAPM1(d,Mh,v,tab,Param,numGr,a11,a12,a21,a22,b1,b2);
      break;
  default: cout <<" the scheme"<<d.scheme<<" does not exist "<<endl;
    }

 
}

/** Flux pour le schéma aux noeuds associé a M1**/
vectorflux FluxVertexAPM1(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab, ParamPhysic & Param,R2 * ur){
  /**  Nodal flux for the nodal JL-(b) scheme for the M1 model **/
    vectorflux res(3);
  R2 sol;
  double s[3]={0,0,0};
  tensor alpha(d,'h');
  tensor beta(d,'s');
  int numGr;
  double kr=0,rjr=0;
  R2 Gjr(0,0);
  R2 Gjr2(0,0);
  double Er=0,qj=0,signode=0, epsnode=0;
  R2 Fr(0,0);
  R2 uj(0,0);
 
  
  for(int r=0;r<Mh.nbnodelocal;r++){

     numGr=Mh(numCell,r);
     sol=ur[numGr];
     beta=inittensor(d,Mh,v,tab,'s',numGr,numCell);
     alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
     signode=AverageSigNodeM1(d,Mh,v,tab,Param.M1,numGr,numCell);
     epsnode=AverageEpsNodeM1(d,Mh,v,tab,Param.M1,numGr,numCell); 

     Er=AverageEnergyNodeM1(d,Mh,v,tab,Param.M1,numGr,numCell);
     Fr=AverageFluxNodeM1(d,Mh,v,tab,Param.M1,numGr,numCell);
     kr=coefk(d,Param.M1,Er,Fr.x,Fr.y);
     
     qj=q(d,Param.M1,v.var[0][numCell],v.var[1][numCell],v.var[2][numCell]);
     Calcul_u(d,Mh,v,tab,Param.M1,uj,numCell);
     rjr=(4./sqrt(3.))*(v.var[0][numCell]/(3+(uj,uj)));
     
     Gjr.x=Mh.ljr(numCell,r)*(qj*Mh.njr(numCell,r).x+rjr*(alpha.ten[0][0]*(uj.x-sol.x)+alpha.ten[0][1]*(uj.y-sol.y))-kr*(signode/epsnode)*(sol.x*beta.ten[0][0]+sol.y*beta.ten[0][1]));
     
     Gjr.y=Mh.ljr(numCell,r)*(qj*Mh.njr(numCell,r).y+rjr*(alpha.ten[1][0]*(uj.x-sol.x)+alpha.ten[1][1]*(uj.y-sol.y))-kr*(signode/epsnode)*(sol.x*beta.ten[1][0]+sol.y*beta.ten[1][1]));
      
     Gjr2.x=Mh.ljr(numCell,r)*(qj*Mh.njr(numCell,r).x+rjr*(alpha.ten[0][0]*(uj.x-sol.x)+alpha.ten[0][1]*(uj.y-sol.y)));
     
     Gjr2.y=Mh.ljr(numCell,r)*(qj*Mh.njr(numCell,r).y+rjr*(alpha.ten[1][0]*(uj.x-sol.x)+alpha.ten[1][1]*(uj.y-sol.y)));
    
     s[0]=s[0]+remap(d,Mh,v,0,tab,Param,numCell,r,ur,sol)+(sol,Gjr);
     s[1]=s[1]+remap(d,Mh,v,1,tab,Param,numCell,r,ur,sol)+Gjr2.x;                    
     s[2]=s[2]+remap(d,Mh,v,2,tab,Param,numCell,r,ur,sol)+Gjr2.y;
    
  }
 
  res.vflux[0]=-s[0]/Param.M1.eps[numCell];
  res.vflux[1]=-s[1]/Param.M1.eps[numCell];
  res.vflux[2]=-s[2]/Param.M1.eps[numCell];

      return res;
}


vectorflux FluxVertexClassicM1(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab, ParamPhysic & Param,R2 * ur){
  /**  Nodal flux for the classical nodal scheme for the M1 model **/
    vectorflux res(3);
  R2 sol;
  double s[3]={0,0,0};
  tensor alpha(d,'h');
  tensor beta(d,'s');
  int numGr;
  double rjr=0;
  R2 Gjr(0,0);
  double qj=0;
  R2 uj(0,0);
 
  
  for(int r=0;r<Mh.nbnodelocal;r++){

     numGr=Mh(numCell,r);
     sol=ur[numGr];
     beta=inittensor(d,Mh,v,tab,'s',numGr,numCell);
     alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
     
     qj=q(d,Param.M1,v.var[0][numCell],v.var[1][numCell],v.var[2][numCell]);
     Calcul_u(d,Mh,v,tab,Param.M1,uj,numCell);
     rjr=(4./sqrt(3.))*(v.var[0][numCell]/(3+(uj,uj)));
     
     Gjr.x=Mh.ljr(numCell,r)*(qj*Mh.njr(numCell,r).x+rjr*(alpha.ten[0][0]*(uj.x-sol.x)+alpha.ten[0][1]*(uj.y-sol.y)));
			      
     Gjr.y=Mh.ljr(numCell,r)*(qj*Mh.njr(numCell,r).y+rjr*(alpha.ten[1][0]*(uj.x-sol.x)+alpha.ten[1][1]*(uj.y-sol.y)));						          
     s[0]=s[0]+remap(d,Mh,v,0,tab,Param,numCell,r,ur,sol)+(sol,Gjr);
     s[1]=s[1]+remap(d,Mh,v,1,tab,Param,numCell,r,ur,sol)+Gjr.x;                    
     s[2]=s[2]+remap(d,Mh,v,2,tab,Param,numCell,r,ur,sol)+Gjr.y;
    
  }
 
  res.vflux[0]=-s[0]/Param.M1.eps[numCell];
  res.vflux[1]=-s[1]/Param.M1.eps[numCell];
  res.vflux[2]=-s[2]/Param.M1.eps[numCell];

      return res;
}


vectorflux FluxVertexM1Gosse(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab, ParamPhysic & Param,R2 * ur){
  /**  Nodal flux for the nodal JL-(b) with local source term scheme for the M1 model **/
    vectorflux res(3);
  R2 sol;
  double s[3]={0,0,0};
  tensor alpha(d,'h');
  tensor beta(d,'s');
  int numGr=0;
  double rjr=0.;
  R2 Gjr(0,0);
  R2 Gjr2(0,0);
  double qj;
  R2 Fr(0,0);
   R2 uj(0,0);
  double Mr[2][2];
  double temp=0.;
  double Er=0.,kr=0.;
  double signode=0.,epsnode=0.;
 
  
  for(int r=0;r<Mh.nbnodelocal;r++){
     numGr=Mh(numCell,r);
     MrConstructM1(d,numGr,Mh,v,tab,Param,Mr);
     sol=ur[numGr];
     signode=AverageSigNodeM1(d,Mh,v,tab,Param.M1,numGr,numCell);
     epsnode=AverageEpsNodeM1(d,Mh,v,tab,Param.M1,numGr,numCell);

     Er=AverageEnergyNodeM1(d,Mh,v,tab,Param.M1,numGr,numCell);
     Fr=AverageFluxNodeM1(d,Mh,v,tab,Param.M1,numGr,numCell);
     kr=coefk(d,Param.M1,Er,Fr.x,Fr.y);

     alpha=inittensor(d,Mh,v,tab,'h',numGr,numCell);
     beta=inittensor(d,Mh,v,tab,'s',numGr,numCell);
  
     qj=q(d,Param.M1,v.var[0][numCell],v.var[1][numCell],v.var[2][numCell]);
     Calcul_u(d,Mh,v,tab,Param.M1,uj,numCell);
     rjr=(4./sqrt(3.))*(v.var[0][numCell]/(3+(uj,uj)));


    Gjr.x=(qj*Mh.njr(numCell,r).x+rjr*(alpha.ten[0][0]*Mr[0][0]+alpha.ten[0][1]*Mr[1][0])*(uj.x-sol.x)+rjr*(alpha.ten[0][0]*Mr[0][1]+alpha.ten[0][1]*Mr[1][1])*(uj.y-sol.y));
          
     Gjr.y=(qj*Mh.njr(numCell,r).y+rjr*(alpha.ten[1][0]*Mr[0][0]+alpha.ten[1][1]*Mr[1][0])*(uj.x-sol.x)+rjr*(alpha.ten[1][0]*Mr[0][1]+alpha.ten[1][1]*Mr[1][1])*(uj.y-sol.y));


     temp=sol.x;
     sol.x=Mr[0][0]*sol.x+Mr[0][1]*sol.y;
     sol.y=Mr[1][0]*temp+Mr[1][1]*sol.y;

     Gjr2.x=qj*Mh.njr(numCell,r).x+rjr*(alpha.ten[0][0]*(uj.x-sol.x)+alpha.ten[0][1]*(uj.y-sol.y))-kr*(signode/epsnode)*(sol.x*beta.ten[0][0]+sol.y*beta.ten[0][1]);
     
     Gjr2.y=qj*Mh.njr(numCell,r).y+rjr*(alpha.ten[1][0]*(uj.x-sol.x)+alpha.ten[1][1]*(uj.y-sol.y))-kr*(signode/epsnode)*(sol.x*beta.ten[1][0]+sol.y*beta.ten[1][1]);

     s[0]=s[0]+remap(d,Mh,v,0,tab,Param,numCell,r,ur,sol)+Mh.ljr(numCell,r)*(sol,Gjr2);
     s[1]=s[1]+remap(d,Mh,v,1,tab,Param,numCell,r,ur,sol)+Mh.ljr(numCell,r)*Gjr.x;                    
     s[2]=s[2]+remap(d,Mh,v,2,tab,Param,numCell,r,ur,sol)+Mh.ljr(numCell,r)*Gjr.y;

     }
     
  res.vflux[0]=-s[0]/Param.M1.eps[numCell];
  res.vflux[1]=-s[1]/Param.M1.eps[numCell];
  res.vflux[2]=-s[2]/Param.M1.eps[numCell];

      return res;
}


void MrConstructM1(Data & d,int numGr,Mesh & Mh,variable & v, TabConnecInv & tab, ParamPhysic & Param,double Mr[2][2]){
  /**  Nodal magical matrix Mr for the M1 model **/
  double a11=0,a12=0,a21=0,a22=0,Det=0,temp=0;
 double c11=0,c12=0,c21=0,c22=0;
 int numLrj=0,jG=0;
 double signode=0, epsnode=0;
 double kr=0,Er=0,rjr=0;
 R2 Fr(0,0);
 tensor alpha(d,'h');
 tensor beta(d,'s');
 R2 uj(0,0); 



 for(int j=0;j<tab.TabInv[numGr].taille;j++){
      jG=tab.TabInv[numGr].TabCell[j];
      numLrj=NodeGtoL(Mh,jG,numGr);  
      

      beta=inittensor(d,Mh,v,tab,'s',numGr,jG);
      alpha=inittensor(d,Mh,v,tab,'h',numGr,jG);   
      signode=AverageSigNodeM1(d,Mh,v,tab,Param.M1,numGr,jG);
      epsnode=AverageEpsNodeM1(d,Mh,v,tab,Param.M1,numGr,jG);
       Calcul_u(d,Mh,v,tab,Param.M1,uj,jG);
      rjr=(4./sqrt(3.))*(v.var[0][jG]/(3+(uj,uj)));
       
      Er=AverageEnergyNodeM1(d,Mh,v,tab,Param.M1,numGr,jG);
      Fr=AverageFluxNodeM1(d,Mh,v,tab,Param.M1,numGr,jG);
      kr=coefk(d,Param.M1,Er,Fr.x,Fr.y);
     
       c11=c11+Mh.ljr(jG,numLrj)*rjr*alpha.ten[0][0];
       a11=a11+kr*(signode/epsnode)*Mh.ljr(jG,numLrj)*beta.ten[0][0]+Mh.ljr(jG,numLrj)*rjr*alpha.ten[0][0];
       c12=c12+Mh.ljr(jG,numLrj)*rjr*alpha.ten[0][1];
       a12=a12+kr*(signode/epsnode)*Mh.ljr(jG,numLrj)*beta.ten[0][1]+Mh.ljr(jG,numLrj)*rjr*alpha.ten[0][1];
       c21=c21+Mh.ljr(jG,numLrj)*rjr*alpha.ten[1][0];
       a21=a21+kr*(signode/epsnode)*Mh.ljr(jG,numLrj)*beta.ten[1][0]+Mh.ljr(jG,numLrj)*rjr*alpha.ten[1][0];
       c22=c22+Mh.ljr(jG,numLrj)*rjr*alpha.ten[1][1];
       a22=a22+kr*(signode/epsnode)*Mh.ljr(jG,numLrj)*beta.ten[1][1]+Mh.ljr(jG,numLrj)*rjr*alpha.ten[1][1]; 
   
       
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




