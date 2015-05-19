
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Variables.hpp"
#include "ClassMesh.hpp"
#include "Functions.hpp"
#include "Initialisation.hpp"
#include "Functionsinit.hpp"
#include "FunctionsAdvection.hpp"
#include <algorithm>
#include "BoundaryCondition.hpp"

using std::min;
using std::max;
using std::abs;

/** File which contains the functions for the advetionf fluxes **/


int EnsembleR(Data & d,Mesh & Mh, variable & v,ParamPhysic & Param,int r,int numCell, R2 a){
  /** function which compute if the normal nodal velocity is positive or negative **/
  int res=0;
  double value;
  int numGr=0;
  
  numGr=Mh(numCell,r);
  value=Mh.ljr(numCell,r)*(a,Mh.njr(numCell,r));
  if(value>0){
    res=1;
  }
  if(value<0){
    res=-1;
  }
  return res;

}


double variablekr(Data & d,Mesh & Mh, variable & v,int var,TabConnecInv & tab,ParamPhysic & Param,int numCell,int numGr, R2  a,int ordre){
  /** nodal average the compute the nodal value of the unknown for the advection scheme **/
  double dem=0;
  double num=0;
  double res=0;
  int jG=0,r=0;
  
  for(int j=0;j<tab.TabInv[numGr].taille;j++){
    jG=tab.TabInv[numGr].TabCell[j];
    r=NodeGtoL(Mh,jG,numGr); 
 
    if(EnsembleR(d,Mh,v,Param,r,jG,a)==1){
      num=num+Mh.ljr(jG,r)*(a,Mh.njr(jG,r))*InterLineVertex(d,Mh,v,tab,Param,var,ordre,jG,r);
      dem=dem+Mh.ljr(jG,r)*(a,Mh.njr(jG,r));
      
    } 
  }

  if(dem==0){
      cout<<" "<<a<<" "<<Mh.njr(jG,r)<<endl;
      cout<<" le denominateur dans l'advection est nul"  <<numGr<<endl;exit(1);
  }
  res=num/dem;
  return res;
}



double VertexUpwind(Data & d,Mesh & Mh, variable & v,int var,TabConnecInv & tab,ParamPhysic & Param,int numCell,int r, R2 a){
  /** Vertex upwind fluxes for one interface **/
  double res=0;
  int numGr;
  int ordre=1;


  if(Param.Model == 1){
    ordre=2;
  } 
  if(Param.Model ==2){
    ordre=Param.Adv.OrderAdv;
  }
  if(Param.Model == 5){
    ordre=Param.Euler.OrderAdv;
  }
  if(Param.Model == 6){
    ordre=Param.M1.OrderAdv;
  }

  numGr=Mh(numCell,r);
 
  if(EnsembleR(d,Mh,v,Param,r,numCell,a)==1){
    res=res+Mh.ljr(numCell,r)*(a,Mh.njr(numCell,r))*InterLineVertex(d,Mh,v,tab,Param,var,ordre,numCell,r);
  }
  if(EnsembleR(d,Mh,v,Param,r,numCell,a)==-1){
    
    res=res+Mh.ljr(numCell,r)*(a,Mh.njr(numCell,r))*variablekr(d,Mh,v,var,tab,Param,numCell,numGr,a,ordre);
  }
  
  return res;
}

double EdgeUpwind(Data & d,Mesh & Mh, variable & v,int var,TabConnecInv & tab,ParamPhysic & Param,int numCell,int k,int r, R2 a){
  /** Edge upwind fluxes for one interface **/
  double res=0;
  double an=0;
   
  an=Mh.ljk(numCell,r)*(Mh.njk(numCell,r),a);
  if(an>0){
      res=an*v.var[0][numCell];
      }
  else{
    res=an*v.var[0][k];
  }
  return res;
}




double VertexSlope(Data & d,Mesh & Mh,variable & v,int var,TabConnecInv & tab,ParamPhysic & Param,int numCell, R2 Grad){
  /** Vertex slope for MUSCL **/
  int taille;
  int * tabnumCell=NULL;
  double * Tabpr=NULL;
  double * Tabpi=NULL;
  double slope=0;
  double  m1,m2;

    taille=tabVF9(Mh,tab,tabnumCell,numCell);
    Tabpr = new double[Mh.nbnodelocal];
    Tabpi = new double[taille];

    for(int i=0;i<taille;i++){
      Tabpi[i]=v.var[var][tabnumCell[i]];
    }
    for(int r=0;r<Mh.nbnodelocal;r++){
      Tabpr[r]=v.var[var][numCell]-(Mh.xr(numCell,r)-Mh.xj(numCell),Grad);
    }

    Tridouble(Tabpr,Mh.nbnodelocal);
    Tridouble(Tabpi,taille);
    m1=(Tabpi[0]-v.var[var][numCell])/(Tabpr[0]-v.var[var][numCell]);
    m2=(Tabpi[taille-1]-v.var[var][numCell])/(Tabpr[Mh.nbnodelocal-1]-v.var[var][numCell]);

    slope=Limiter(Param,min(m1,m2));
    
    delete[] Tabpr;
    delete[] Tabpi;
    delete[] tabnumCell;


    return slope;
}

double EdgeSlope(Data & d,Mesh & Mh,variable & v,int var,TabConnecInv & tab,ParamPhysic & Param,int numCell, R2 Grad){
   /** Edge slope for MUSCL **/
  double * Tabvark;
  double * TabGrad;
  double slope=0,slope2=0;
  R2 xjk(0,0);
  int r2=0,numGr1=0,numGr2=0,k=0;
  
    Tabvark = new double[Mh.nbnodelocal];
    TabGrad = new double[Mh.nbnodelocal];

    for(int r=0;r<Mh.nbnodelocal;r++){
      if(r==Mh.nbnodelocal-1){
	r2=0;
      }
      else{
	r2=r+1;
      }
      numGr1=Mh(numCell,r);
      numGr2=Mh(numCell,r2);
      xjk=0.5*(Mh.xr(numGr1)+Mh.xr(numGr2));
      k=InverseEdge(Mh,numGr1,numGr2,numCell,tab);
      Tabvark[r]=v.var[var][k];
      TabGrad[r]=Valabs((Grad,xjk-Mh.xj(numCell)));
      
    }
    Tridouble(Tabvark,Mh.nbnodelocal);
    Tridouble(TabGrad,Mh.nbnodelocal);

    slope2=min(Tabvark[Mh.nbnodelocal-1]-v.var[var][numCell],v.var[var][numCell]-Tabvark[0]);
    slope=min(1.,0.5*(slope2/TabGrad[Mh.nbnodelocal-1]));

    
    delete[] Tabvark;
    delete[] TabGrad;

 
    return slope;
}





double Limiter(ParamPhysic & Param,double m){
   /** Limiter for MUSCL slop **/
  double ml=0;
  /** minmod **/
  if(Param.Adv.Nlim==1){
    double beta=1;
    ml=max(0.,min(m,beta));
  }
  /** superbee **/
 if(Param.Adv.Nlim==2){
   double temp=0;
   double beta=1;
   temp=max(min(beta*m,1.),min(beta,m));
   ml=max(0.,temp);
 }
 /** Van leer **/
 if(Param.Adv.Nlim==3){
   ml=(m+abs(m))/(1+abs(m));
 }
 if(Param.Model==3){
    double beta=1;
    ml=1;
  }
 
 return ml;
}


double InterLineVertex(Data d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamPhysic & Param,int var,int ordre,int j,int r){
   /** function which compute the value used by the scheme. Cell average for the order 1, linear reconstruction for the order 2 **/
  double res=0;
  int numGr;
  
  if(ordre==1){
    res=v.var[var][j];
  }
  if(ordre==2){
      R2 Grad(0,0);
      double m;
      numGr=Mh(j,r);
      Grad=Gradiant(d,Mh,v,var,tab,Param,numGr);
      m=VertexSlope(d,Mh,v,var,tab,Param,j,Grad);
      res=v.var[var][j]-m*((Mh.xr(j,r)-Mh.xj(j)),Grad);

    }
  
  return res;
}



R2 Gradiant(Data & d,Mesh & Mh,variable & v,int var,TabConnecInv & tab,ParamPhysic & Param,int numGr){
   /** Nodal gradiant for the MUSCL reconstruction **/
  double a11=0,a12=0,a21=0,a22=0,b1=0,b2=0,Det=0;
   int numLrj,jG;
   R2 sol(0,0);
   tensor beta(d,'s');


      a11=0,a12=0,a21=0,a22=0,b1=0,b2=0,Det=0;
     for(int j=0;j<tab.TabInv[numGr].taille;j++){
       jG=tab.TabInv[numGr].TabCell[j];   
       numLrj=NodeGtoL(Mh,jG,numGr);  

       a11=a11+Mh.ljr(jG,numLrj)*Mh.njr(jG,numLrj).x*(Mh.xr(jG,numLrj).x-Mh.xj(jG).x);
       a12=a12+Mh.ljr(jG,numLrj)*Mh.njr(jG,numLrj).x*(Mh.xr(jG,numLrj).y-Mh.xj(jG).y);
       a21=a21+Mh.ljr(jG,numLrj)*Mh.njr(jG,numLrj).y*(Mh.xr(jG,numLrj).x-Mh.xj(jG).x);
       a22=a22+Mh.ljr(jG,numLrj)*Mh.njr(jG,numLrj).y*(Mh.xr(jG,numLrj).y-Mh.xj(jG).y); 

       b1=b1+Mh.ljr(jG,numLrj)*v.var[var][jG]*Mh.njr(jG,numLrj).x;
       b2=b2+Mh.ljr(jG,numLrj)*v.var[var][jG]*Mh.njr(jG,numLrj).y;
     }
     
     Det=a11*a22-a21*a12;
     
     if(Mh.vertices[numGr].lab >-1){
       
       if(Det==0){cout<<"The nodal matrix for the gradiant (MUSCL) is not invertible"<<numGr<<endl;exit(1);}
       else { 
	 sol.x=(b1*a22-b2*a12)/Det;
	 sol.y=(b2*a11-b1*a21)/Det;
       }
     }
     else{
       sol.x=b1/Vr(Mh,numGr,tab);
       sol.y=b2/Vr(Mh,numGr,tab);
     }

        return sol;
}

