#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Variables.hpp"
#include "Flux.hpp"
#include "Source.hpp"
#include "BoundaryCondition.hpp"
#include "Tensor.hpp"
#include "FunctionsAdvection.hpp"


/** Calcul de la pression q**/
double q_M1M(Data & d, ParamM1Matter & M1Matter,double E,double F1, double F2){
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
double coefk_M1M(Data & d, ParamM1Matter & M1Matter,double E,double F1, double F2){
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

void  Calcul_u_M1M(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab, ParamM1Matter & M1Matter,R2 & u,int numCell){
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


void MatrixM1Matter(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,ParamPhysic & Param,int numGr,double &a11,double &a12,double &a21,double &a22,double &b1,double &b2){
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
        
        Calcul_u_M1M(d,Mh,v,tab,Param.M1M,uj,jG);
        qj=q_M1M(d,Param.M1M,v.var[0][jG],v.var[1][jG],v.var[2][jG]);
        rjr=(4./sqrt(3.))*(v.var[0][jG]/(3+(uj,uj)));
        
        a11=a11+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[0][0]);
        a12=a12+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[0][1]);
        a21=a21+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[1][0]);
        a22=a22+Mh.ljr(jG,numLrj)*(rjr*alpha.ten[1][1]);
        
        b1=b1+Mh.ljr(jG,numLrj)*(qj*Mh.njr(jG,numLrj).x+rjr*alpha.ten[0][0]*uj.x+rjr*alpha.ten[0][1]*uj.y);
        b2=b2+Mh.ljr(jG,numLrj)*(qj*Mh.njr(jG,numLrj).y+rjr*alpha.ten[1][0]*uj.x+rjr*alpha.ten[1][1]*uj.y);
        
    }
    
    
}



vectorflux FluxVertexClassicM1Matter(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab, ParamPhysic & Param,R2 * ur){
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
        
        qj=q_M1M(d,Param.M1M,v.var[0][numCell],v.var[1][numCell],v.var[2][numCell]);
        Calcul_u_M1M(d,Mh,v,tab,Param.M1M,uj,numCell);
        rjr=(4./sqrt(3.))*(v.var[0][numCell]/(3+(uj,uj)));
        
        Gjr.x=Mh.ljr(numCell,r)*(qj*Mh.njr(numCell,r).x+rjr*(alpha.ten[0][0]*(uj.x-sol.x)+alpha.ten[0][1]*(uj.y-sol.y)));
        
        Gjr.y=Mh.ljr(numCell,r)*(qj*Mh.njr(numCell,r).y+rjr*(alpha.ten[1][0]*(uj.x-sol.x)+alpha.ten[1][1]*(uj.y-sol.y)));

	
        s[0]=s[0]+remap(d,Mh,v,0,tab,Param,numCell,r,ur,sol)+(sol,Gjr);
        s[1]=s[1]+remap(d,Mh,v,1,tab,Param,numCell,r,ur,sol)+Gjr.x;
        s[2]=s[2]+remap(d,Mh,v,2,tab,Param,numCell,r,ur,sol)+Gjr.y;
        
    }

    res.vflux[0]=-s[0]/Param.M1M.eps[numCell];
    res.vflux[1]=-s[1]/Param.M1M.eps[numCell];
    res.vflux[2]=-s[2]/Param.M1M.eps[numCell];
    res.vflux[3]= 0. ;

    return res;
}
