
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Flux.hpp"

/**
   Function which gives the flux associated with your model and your shceme (in Data)**/


vectorflux ChoiceFluxN(Data & d ,int numCell,Mesh &  Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2 ** ur){
  /**  Function which gives the Nodal flux associated with your model and your scheme (in Data) **/
vectorflux res(d);
 
  if(Param.Model == 1){
      switch(d.scheme){
      case 1:
	/** Linear Nodal diffusion scheme**/
	res=FluxVertexDiff(d,numCell,Mh,v,tab,Param,ur);
	break;
	/** Nonlinear Nodal diffusion scheme aobtain using advection fluxes**/
      case 2:
	res=FluxVertexNlDiff(d,numCell,Mh,v,tab,Param,ur);
      break;
      default: cout <<" the scheme associated with the number"<<d.scheme<<" does not exist "<<endl;
    }
  
  }

  if(Param.Model == 2){
    
    switch(d.scheme){
    case 1:
      /** Linear Nodal advection scheme (order one and two)**/
      res=FluxVertexadv(d,numCell,Mh,v,tab,Param);
      break;
    default: cout <<" the scheme associated with the number"<<d.scheme<<" does not exist "<<endl;
    }
  }
  
   if(Param.Model == 3){
    
    switch(d.scheme){
    case 1:
      /** Linear JL-(a) nodal scheme for P1 model **/
      res=FluxVertexP1(d,numCell,Mh,v,tab,Param,ur);
      break;
    case 2:
      /** Linear JL-(b) nodal scheme for P1 model **/
      res=FluxVertexP1(d,numCell,Mh,v,tab,Param,ur);
      break;  
    case 3:
      /** Linear JL-(b) nodal scheme with local source term for P1 model **/
      res=FluxVertexP1Gosse(d,numCell,Mh,v,tab,Param,ur);
      break;  
    default: cout <<" the scheme associated with the number"<<d.scheme<<" does not exist "<<endl;
    }
  }

   if(Param.Model == 4){
     
     switch(d.scheme){
     case 1:
       /** Linear JL-(a) nodal scheme for P1 model with matter**/
       res=FluxVertexP1Matter(d,numCell,Mh,v,tab,Param,ur);
       break;
     case 2:
       /** Linear JL-(b) nodal scheme for P1 model with matter **/
       res=FluxVertexP1Matter(d,numCell,Mh,v,tab,Param,ur);
       break;  
     case 3:
       /** Linear JL-(b) nodal scheme with local source term for P1 model with matter **/
       res=FluxVertexP1MatterGosse(d,numCell,Mh,v,tab,Param,ur);
       break;  
     default: cout <<" the scheme associated with the number"<<d.scheme<<" does not exist "<<endl;
     }
   }
   
  if(Param.Model == 5){
    
    switch(d.scheme){
    case 1:
      /** Lagrange+remap JL-(a) nodal scheme for Euler model **/
      res=FluxVertexEuler(d,numCell,Mh,v,tab,Param,ur);
      break;
    case 2:
       /** Lagrange+remap JL-(b) nodal scheme for Euler model **/
      res=FluxVertexEuler(d,numCell,Mh,v,tab,Param,ur);
      break;
    case 3:
       /** Lagrange+remap JL-(b) nodal scheme with local source term for Euler model **/
      res=FluxVertexEulerGosse(d,numCell,Mh,v,tab,Param,ur);
      break;
    case 4:
       /** Classic Lagrange+remap nodal scheme for Euler model **/
      res=FluxVertexEuler(d,numCell,Mh,v,tab,Param,ur);
      break;
    default: cout <<" the scheme associated with the number"<<d.scheme<<" does not exist "<<endl;
    }
  }

   if(Param.Model == 6){
    
    switch(d.scheme){
    case 1:
       /** Classic Lagrange+remap nodal scheme for M1 model **/
      res=FluxVertexClassicM1(d,numCell,Mh,v,tab,Param,ur);
      break;
    case 2:
      /** Lagrange+remap JL-(b) nodal scheme for M1 model **/
      res=FluxVertexAPM1(d,numCell,Mh,v,tab,Param,ur);
      break;
    case 3:
      /** Lagrange+remap JL-(b) nodal scheme with local source term for M1 model **/
      res=FluxVertexM1Gosse(d,numCell,Mh,v,tab,Param,ur);
      break;
    default: cout <<" the scheme associated with the number"<<d.scheme<<" does not exist "<<endl;
    }
  }

    if (Param.Model == 7)
    {
        res=FluxVertexClassicM1Matter(d,numCell,Mh,v,tab,Param,ur);
    }

    if(Param.Model == 8){
     
     switch(d.scheme){
     case 1:
       /** Linear JL-(a) nodal scheme for P1 model with matter**/
       res=FluxVertexP1Compton(d,numCell,Mh,v,tab,Param,ur);
       break;
     case 2:
       /** Linear JL-(b) nodal scheme for P1 model with matter **/
       res=FluxVertexP1Compton(d,numCell,Mh,v,tab,Param,ur);
       break;  
     case 3:
       /** Linear JL-(b) nodal scheme with local source term for P1 model with matter **/
       res=FluxVertexP1ComptonGosse(d,numCell,Mh,v,tab,Param,ur);
       break;  
     default: cout <<" the scheme associated with the number"<<d.scheme<<" does not exist "<<endl;
     }
   }
  
  return res;
}

vectorflux ChoiceFluxE(Data & d ,int numCell,Mesh &  Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param){
   /**  Function which gives the Edge flux associated with your model and your shceme (in Data) **/
vectorflux res(d);
 
  if(Param.Model ==1){
      switch(d.scheme){
      case 1:
      /** Linear VF4 Edge scheme for diffusion **/
	res=FluxVF4Diff(d,numCell,Mh,v,tab,Param);
	break;
    default: cout <<" the scheme associated with the number"<<d.scheme<<" does not exist "<<endl;
    }
  
  }

   if(Param.Model == 2){
    
    switch(d.scheme){
    case 1:
      /** Linear Edge scheme for advection **/
      res=FluxEdgeAdv(d,numCell,Mh,v,tab,Param);
      break;
    default: cout <<" the scheme associated with the number"<<d.scheme<<" does not exist "<<endl;
    }
   }  
  return res;
}

R2 SolveurNodal(Data & d,int numGr,Mesh & Mh, variable & v,TabConnecInv & tab,ParamPhysic & Param,int group){
   /**  Function which gives the nodal sovler associated with your model (in Data) **/
  double a11=0,a12=0,a21=0,a22=0,b1=0,b2=0,Det=0;
  R2 sol(0,0);

  if(Mh.vertices[numGr].lab > -1){
   
    if(Param.Model == 1){
       MatrixDiff(d,Mh,v,tab,Param,numGr,a11,a12,a21,a22,b1,b2);
    }
     if(Param.Model == 3){
      MatrixP1(d,Mh,v,tab,Param,numGr,a11,a12,a21,a22,b1,b2);
    }
     if(Param.Model == 4){
       MatrixP1Matter(d,Mh,v,tab,Param,numGr,a11,a12,a21,a22,b1,b2);
     }
     if(Param.Model == 5){
       MatrixEuler(d,Mh,v,tab,Param,numGr,a11,a12,a21,a22,b1,b2);
     }
     if(Param.Model == 6){
       MatrixM1(d,Mh,v,tab,Param,numGr,a11,a12,a21,a22,b1,b2);
     } 
     if(Param.Model == 7){
       MatrixM1Matter(d,Mh,v,tab,Param,numGr,a11,a12,a21,a22,b1,b2);
     }
     if(Param.Model == 8){
       MatrixP1Compton(d,Mh,v,tab,Param,numGr,a11,a12,a21,a22,b1,b2,group);
     }

     
    Det=a11*a22-a21*a12;
    if(Det==0){ cout<<"The matrix associated with the nodel solver is not invertible at the node "<<numGr<<endl;exit(1);}
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








