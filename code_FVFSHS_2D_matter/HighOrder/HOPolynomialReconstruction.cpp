#include <iostream>
#include <cmath>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include "Functions.hpp"
#include "HOfunctions.hpp"
#include "Functionsgeo.hpp"
#include "Flux.hpp" 

void LSInit(Data & d, Mesh & Mh, variable & v,TabConnecInv & tab,ParamEuler & Euler, LSProblem & LS, int PolynomialOrder,int numGr) {
  /** function which initialize a Least square problem for high order  reconstruction **/ 
  int k=0;

  ConstructionLSStencil(d,Mh,v,tab,Euler,LS,PolynomialOrder,numGr); /** initialize M (number of cell) and the table of the stencil **/
   
  k=PolynomialOrder;
  LS.NColum=0.5*(k+1)*(k+2); /** LHO is the order of convergence for WB scheme consequently the order of the polynomial k=LHO-1 **/
  
  LS.b_rho = new double[LS.NLign];
  LS.Sol_rho = new double[LS.NColum];

  LS.b_pressure = new double[LS.NLign];
  LS.Sol_pressure = new double[LS.NColum];

   for (int i = 0; i < LS.NColum; i++ )
    {
      LS.Sol_rho[i] =0;
      LS.Sol_pressure[i] = 0;
    }
 
  LS.A = new double[LS.NLign*LS.NColum];

  for (int i = 0; i < LS.NColum*LS.NLign; i++ )
    {
      LS.A[i] =0;
    }
  

  ConstructionLSMatrix(d,Mh,v,tab,Euler,LS,PolynomialOrder,numGr); 
}

void ConstructionLSStencil(Data & d, Mesh & Mh, variable & v,TabConnecInv & tab,ParamEuler & Euler, LSProblem & LS,int PolynomialOrder,int numGr) {
  /** function which compute the stencil necessary for the high order reconstruction around the node numGr **/
  /** Actually this function is valid only on quadrangular mesh **/

  if(d.Typemesh == 'Q') {
    if(PolynomialOrder == 1) {
      LS.NLign=4;
      LS.Stencil = new int[LS.NLign];
      OneLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr);
    
    }
    if(PolynomialOrder == 2) {
      LS.NLign=16;
      LS.Stencil = new int[LS.NLign];
      if(Mh.xr(numGr).lab == 0){
	 TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr);
      }
      if(Mh.xr(numGr).lab == 2){
	TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr+1);
      }
      if(Mh.xr(numGr).lab == 3){
	TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr-(d.Nx+3));
      }
      if(Mh.xr(numGr).lab == 4){
	TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr-1);
      }
      if(Mh.xr(numGr).lab == 5){
	TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr+(d.Nx+3));
      }
      if(Mh.xr(numGr).lab == 1){
	if(Mh.xr(numGr).x==0 && Mh.xr(numGr).y ==0){
	  TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr+(d.Nx+4));
	}
	if(Mh.xr(numGr).x==0 && Mh.xr(numGr).y ==d.Ty){
	  TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr-(d.Nx+2));
	}
	if(Mh.xr(numGr).x==d.Tx && Mh.xr(numGr).y ==d.Ty){
	  TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr-(d.Nx+4));
	}
	if(Mh.xr(numGr).x==d.Tx && Mh.xr(numGr).y ==0){
	  TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr+(d.Nx+2));
	}
      }

      
    }
    if(PolynomialOrder == 3) {
      LS.NLign=16;
      LS.Stencil = new int[LS.NLign];
      
      if(Mh.xr(numGr).lab == 0){
	 TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr);
      }
      if(Mh.xr(numGr).lab == 2){
	TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr+1);
      }
      if(Mh.xr(numGr).lab == 3){
	TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr-(d.Nx+3));
      }
      if(Mh.xr(numGr).lab == 4){
	TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr-1);
      }
      if(Mh.xr(numGr).lab == 5){
	TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr+(d.Nx+3));
      }
      if(Mh.xr(numGr).lab == 1){
	if(Mh.xr(numGr).x==0 && Mh.xr(numGr).y ==0){
	  TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr+(d.Nx+4));
	}
	if(Mh.xr(numGr).x==0 && Mh.xr(numGr).y ==d.Ty){
	  TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr-(d.Nx+2));
	}
	if(Mh.xr(numGr).x==d.Tx && Mh.xr(numGr).y ==d.Ty){
	  TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr-(d.Nx+4));
	}
	if(Mh.xr(numGr).x==d.Tx && Mh.xr(numGr).y ==0){
	  TwoLayerStencilNodeQ(d,Mh,v,tab,LS.Stencil,LS.NLign,numGr+(d.Nx+2));
	}
    }
  }
 }
  
}

void SolverLSProblem(Data & d, Mesh & Mh, variable & v,TabConnecInv & tab,ParamEuler & Euler, LSProblem & LS,int PolynomialOrder,int numGr) {
   /** function which solve the two least square problems for the high order WB reconstruction around the node numGr **/
  
  double r2_norm_rho=0;
  double r2_norm_pressure=0;
  double * r2_rho;
  double * r2_pressure;
    

  r2_rho = new double[LS.NLign];
  r2_pressure = new double[LS.NLign];

   for (int i = 0; i < LS.NLign; i++ )
    {
      r2_rho[i] =0;
      r2_pressure[i] = 0;
    }
  
  ConstructionLSRHS(d,Mh,v,tab,Euler,LS,PolynomialOrder,numGr);

  svd_solve (LS.NLign, LS.NColum, LS.A, LS.b_rho, LS.Sol_rho );
  svd_solve (LS.NLign, LS.NColum, LS.A, LS.b_pressure, LS.Sol_pressure );

  r2_rho = mult_mv (LS.NLign, LS.NColum, LS.A, LS.Sol_rho );
  r2_pressure = mult_mv (LS.NLign, LS.NColum, LS.A, LS.Sol_pressure);

  for (int i = 0; i < LS.NLign; i++ )
    {
      r2_rho[i] = r2_rho[i] - LS.b_rho[i];
      r2_pressure[i] = r2_pressure[i] - LS.b_pressure[i];
    }
 

  r2_norm_rho = vec_norm2 ( LS.NLign, r2_rho );
  r2_norm_pressure = vec_norm2 ( LS.NLign, r2_pressure );

  if(r2_norm_rho > 0.0000001 ||r2_norm_pressure > 0.0000001)
  {
  //cout<<" node "<<numGr<<" label :"<<Mh.xr(numGr)<<endl;
  // cout<<" LS rho : "<<r2_norm_rho<<" LS pressure : "<<r2_norm_pressure<<endl;
    //  for (int i = 0; i < LS.NLign*LS.NColum; i++ )
    // {
    //cout<<" Aij "<<LS.A[i]<<endl;
    //   }
    //for (int i = 0; i < LS.NLign; i++ )
    // {
    //	cout<<" rhs rho  "<<i<< "  "<<LS.b_rho[i]<<" "<<LS.Stencil[i]<<endl;
    //  }
    //  for (int i = 0; i < LS.NColum; i++ )
    // {
    //	cout<<" sol rho  "<<i<< "  "<<LS.Sol_rho[i]<<endl;
    //	cout<<" sol p  "<<i<< "  "<<LS.Sol_pressure[i]<<endl;
	//}
     }

  delete [] r2_rho;
  delete [] r2_pressure;
 
}

void ConstructionLSMatrix(Data & d, Mesh & Mh, variable & v,TabConnecInv & tab,ParamEuler & Euler, LSProblem & LS,int PolynomialOrder,int numGr) {
 /** function which compute the matrix associated with the least square problem for the high order WB reconstruction around the node numGr **/
  ;  int d1=0,d2=0;
  int numCell=0;

  
  for (int i = 0; i < LS.NLign; i++ )
  {
    numCell=LS.Stencil[i];
    for (int j = 0; j < LS.NColum; j++ )
    {
      LocalDegreePolynomial(d1,d2,j,PolynomialOrder);
      LS.A[i+j*LS.NLign] = ComputeMatrixcoefficent(d,Mh,d1,d2,numCell,numGr,PolynomialOrder);
    }
  }
  
}

void LocalDegreePolynomial(int & d1, int & d2, int j,int PolynomialOrder){

  if(PolynomialOrder ==1 ) {
    switch(j){ 
      case 0 : 	
	d1=0;
	d2=0;	
      break;    
      case 1:
	d1=1;
	d2=0;	
      break;
      case 2:
	d1=0;
	d2=1;	
      break;
     }     
  }
  if(PolynomialOrder ==2 ) {
    switch(j){ 
      case 0 : 	
	d1=0;
	d2=0;	
      break;    
      case 1:
	d1=1;
	d2=0;	
      break;
      case 2:
	d1=2;
	d2=0;	
      break;
      case 3 : 	
	d1=0;
	d2=1;	
      break;    
      case 4:
	d1=1;
	d2=1;	
      break;
      case 5:
	d1=0;
	d2=2;	
      break;
     }  
    
  }
   if(PolynomialOrder ==3 ) {
    switch(j){ 
      case 0 : 	
	d1=0;
	d2=0;	
      break;    
      case 1:
	d1=1;
	d2=0;	
      break;
      case 2:
	d1=2;
	d2=0;	
      break;
      case 3 : 	
	d1=3;
	d2=0;	
      break;    
      case 4:
	d1=0;
	d2=1;	
      break;
      case 5:
	d1=1;
	d2=1;	
      break;
      case 6:
	d1=2;
	d2=1;	
      break;
      case 7: 	
	d1=0;
	d2=2;	
      break;    
      case 8:
	d1=1;
	d2=2;	
      break;
      case 9:
	d1=0;
	d2=3;	
      break;
     }  
    
  }
  
}

void ConstructionLSRHS(Data & d, Mesh & Mh, variable & v,TabConnecInv & tab,ParamEuler & Euler, LSProblem & LS,int PolynomialOrder,int numGr) {
 /** function which compute the two rhs associated with the least square problem for the high order WB reconstruction around the node numGr **/
  
  int nCell=0;
   for (int i = 0; i < LS.NLign; i++ )
    {
      nCell=LS.Stencil[i];
      LS.b_rho[i] =Mh.area(nCell)*v.var[0][nCell];
    }
    for (int i = 0; i < LS.NLign; i++ )
    {
      nCell=LS.Stencil[i];
      LS.b_pressure[i] =Mh.area(nCell)*PressureLaw(d,Mh,v,Euler,nCell);

      // cout<<" b "<<LS.b_pressure[i]<<endl;
    }
  
}


void QuadratureUnitSquare(double & p, double & w, int nb, int Order){ /** NOTTT FINISH **/
 /** function which compute the two rhs associated with the least square problem for the high order WB reconstruction around the node numGr **/

  if(Order == 1) {
    switch(nb){ 
      case 0 : 	
	w=1;
	p=-0.5773502691896257;	
      break;    
      case 1:
	w=1;
	p=0.5773502691896257;	
      break;
     }   
  }
  
  if(Order == 2) {
    switch(nb){ 
      case 0 : 	
	w=1;
	p=-0.5773502691896257;	
      break;    
      case 1:
	w=1;
	p=0.5773502691896257;	
      break;
     }   
  }
  if(Order == 3) {
    switch(nb){ 
      case 0 : 	
	w=0.5555555555555556;
	p=-0.7745966692414834;	
      break;    
      case 1:
	w=0.8888888888888888;
	p=0;	
      break;
      case 2:
	w=0.5555555555555556;
	p=0.7745966692414834;	
      break;
     }  
  }

  if(Order == 4) {
    switch(nb){ 
      case 0 : 	
	w=0.3478548451374538;
	p=-0.8611363115940526;	
      break;    
      case 1:
	w=0.6521451548625461;
	p=-0.3399810435848563;	
      break;
      case 2:
	w=0.6521451548625461;
	p=0.3399810435848563;	
      break;
      case 3:
	w=0.3478548451374538;
	p=0.8611363115940526;	
      break;
     }
    
  }

  if(Order == 5) {
    switch(nb){ 
      case 0 : 	
	w=0.2369268850561891;
	p=-0.9061798459386640;	
      break;    
      case 1:
	w=0.4786286704993665;
	p=-0.5384693101056831;	
      break;
      case 2:
	w=0.5688888888888889;
	p=0;	
      break;
      case 3:
	w=0.4786286704993665;
	p=0.5384693101056831;	
      break;
      case 4:
	w=0.2369268850561891;
	p=0.9061798459386640;	
      break;
     }
    
  }

  if(Order == 6) {
    switch(nb){ 
      case 0 : 	
	w=0.1713244923791704;
	p=-0.9324695142031521;	
      break;    
      case 1:
	w=0.3607615730481386;
	p=-0.6612093864662645;	
      break;
      case 2:
	w=0.4679139345726910;
	p=-0.2386191860831969;	
      break;
      case 3:
	w=0.4679139345726910;
	p=0.2386191860831969;	
      break;
      case 4:
	w=0.3607615730481386;
	p=0.6612093864662645;	
      break;
       case 5 : 	
	w=0.1713244923791704;
	p=0.9324695142031521;	
      break; 
     }
    
  }
   
  
}

double Mapping(Data & d, Mesh & Mh,R2 & p, double & J,double px, double py,cell & C){
  double ShapeFunctions[Mh.nbnodelocal];
  double Jacobian[2][2];
  p.x=0;
  p.y=0;
  
  if(d.Typemesh == 'Q') {
    ShapeFunctions[0]=0.25*(1-px)*(1-py); // shapre function for the node (-1,-1) in the unit square
    ShapeFunctions[1]=0.25*(1+px)*(1-py); // shapre function for the node (-1,1) in the unit square
    ShapeFunctions[2]=0.25*(1+px)*(1+py); // shapre function for the node (1,1) in the unit square
    ShapeFunctions[3]=0.25*(1-px)*(1+py); // shapre function for the node (-1,1) in the unit square
  }
  
  for(int i=0;i<Mh.nbnodelocal;i++){
    p.x=p.x+C.vertices[i].x*ShapeFunctions[i];
    p.y=p.y+C.vertices[i].y*ShapeFunctions[i];
}


  Jacobian[0][0]=0.25*((1+py)*(C.vertices[2].x-C.vertices[3].x) +(1-py)*(C.vertices[1].x-C.vertices[0].x));
  Jacobian[0][1]=0.25*((1+py)*(C.vertices[2].y-C.vertices[3].y) +(1-py)*(C.vertices[1].y-C.vertices[0].y));
  Jacobian[1][0]=0.25*((1-px)*(C.vertices[3].x-C.vertices[0].x) +(1+px)*(C.vertices[2].x-C.vertices[1].x));
  Jacobian[1][1]=0.25*((1-px)*(C.vertices[3].y-C.vertices[0].y) +(1+px)*(C.vertices[2].y-C.vertices[1].y));
   
  J=Jacobian[0][0]*Jacobian[1][1]-Jacobian[1][0]*Jacobian[0][1];
}

double ComputeMatrixcoefficent(Data & d, Mesh & Mh,int d1, int d2, int numCell,int numGr,int Order){
  double J=0;
  R2 p(0,0);
  double Qp1=0;
  double Qp2=0;
  double w1=0,w2=0;
  double res=0;


  for(int i=0;i<Order+2;i++){
   for(int j=0;j<Order+2;j++){
     QuadratureUnitSquare(Qp1,w1,i,Order+2);
     QuadratureUnitSquare(Qp2,w2,j,Order+2);

     Mapping(d,Mh,p,J,Qp1,Qp2,Mh.cells[numCell]);
     res=res+w1*w2*Valabs(J)*pow(p.x,d1)*pow(p.y,d2);
      }
  }
  
  return res;
}

  
