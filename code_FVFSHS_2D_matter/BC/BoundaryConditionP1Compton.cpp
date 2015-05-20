#include <iostream>
#include <cmath>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include <cstdlib>
#include "Variables.hpp"
#include "Data.hpp"
#include "BoundaryCondition.hpp"
#include "Functions.hpp"
#include "Tensor.hpp"
#include "Functionsgeo.hpp"
#include "AnalyticSolutions.hpp"
#include "Diagnostics.hpp"
#include "Functionsinit.hpp"

void BoundaryConditionP1Compton(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time){
 /** Function which initialize the phantom cells for the P1 model with matter equation.
For triangular mesh only the Boundary conditions based on the black body are written **/
   switch(d.nTest)
    {
    case 1 : BC_Neumann_P1Compton(d,Mh,v,Param,time);
    
      break;}
}

void BC_Neumann_P1Compton(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time){
  /** Function which initialize the phantom cells to have Neumann solution **/
  for(int j=0;j<Mh.nc;j++){
    if(Mh.cells[j].lab==-1){
      
	  
	    //cote gauche

	    if((Mh.xj(j).x < 0) && (Mh.xj(j).y> 0 && Mh.xj(j).y< d.Ty)){
	      for(int k=0; k<Param.P1C.nb_group;k++){
		v.var[0+k*Param.P1C.nb_moment][j]=v.var[0+k*Param.P1C.nb_moment][j+1];
		v.var[1+k*Param.P1C.nb_moment][j]=-v.var[1+k*Param.P1C.nb_moment][j+1];	//-
		v.var[2+k*Param.P1C.nb_moment][j]=v.var[2+k*Param.P1C.nb_moment][j+1];
	      }
	      v.var[Param.P1C.nbvar-1][j]=v.var[Param.P1C.nbvar-1][j+1];

	   }
	    //cote droit
	    if(Mh.xj(j).x> d.Tx && (Mh.xj(j).y> 0 && Mh.xj(j).y< d.Ty)){
	     for(int k=0; k<Param.P1C.nb_group;k++){
		v.var[0+k*Param.P1C.nb_moment][j]=v.var[0+k*Param.P1C.nb_moment][j-1];
		v.var[1+k*Param.P1C.nb_moment][j]=-v.var[1+k*Param.P1C.nb_moment][j-1];	//-
		v.var[2+k*Param.P1C.nb_moment][j]=v.var[2+k*Param.P1C.nb_moment][j-1];
	      }
	      v.var[Param.P1C.nbvar-1][j]=v.var[Param.P1C.nbvar-1][j-1];
	   }
	    //cote bas
	    if(Mh.xj(j).y< 0 && (Mh.xj(j).x> 0  && Mh.xj(j).x< d.Tx)){
	      for(int k=0; k<Param.P1C.nb_group;k++){
		v.var[0+k*Param.P1C.nb_moment][j]=v.var[0+k*Param.P1C.nb_moment][j+(d.Nx+2)];
		v.var[1+k*Param.P1C.nb_moment][j]=v.var[1+k*Param.P1C.nb_moment][j+(d.Nx+2)];	//-
		v.var[2+k*Param.P1C.nb_moment][j]=-v.var[2+k*Param.P1C.nb_moment][j+(d.Nx+2)];
	      }
	      v.var[Param.P1C.nbvar-1][j]=v.var[Param.P1C.nbvar-1][j+(d.Nx+2)];
	   }
	    //cote haut
	    if(Mh.xj(j).y> d.Ty && (Mh.xj(j).x > 0  && Mh.xj(j).x < d.Tx)){
	       for(int k=0; k<Param.P1C.nb_group;k++){
		v.var[0+k*Param.P1C.nb_moment][j]=v.var[0+k*Param.P1C.nb_moment][j-(d.Nx+2)];
		v.var[1+k*Param.P1C.nb_moment][j]=v.var[1+k*Param.P1C.nb_moment][j-(d.Nx+2)];	//-
		v.var[2+k*Param.P1C.nb_moment][j]=-v.var[2+k*Param.P1C.nb_moment][j-(d.Nx+2)];
	      }
	      v.var[Param.P1C.nbvar-1][j]=v.var[Param.P1C.nbvar-1][j-(d.Nx+2)];
	    }

      //coin bas gauche

             if(Mh.xj(j).y< 0 && Mh.xj(j).x<0){
	       //v.var[0][j]=v.var[0][j+d.Nx+1];
		//v.var[1][j]=-v.var[1][j+d.Nx+1]; // -
	       //v.var[2][j]=-v.var[2][j+d.Nx+1]; //-
	       //v.var[3][j]=v.var[3][j+d.Nx+1];
	   }
	   //coin haut gauche
	    if(Mh.xj(j).y> d.Ty && Mh.xj(j).x<0){
	      //v.var[0][j]=v.var[0][j-(d.Nx+1)];
	      //v.var[1][j]=-v.var[1][j-(d.Nx+1)]; //-
	      //v.var[2][j]=-v.var[2][j-(d.Nx+1)];//-
	      //v.var[3][j]=v.var[3][j-(d.Nx+1)];

	   }
	    //coin haut droit
	    if(Mh.xj(j).y>d.Ty && Mh.xj(j).x>d.Tx){
	      //v.var[0][j]=v.var[0][j+d.Nx+1];
	      //v.var[1][j]=-v.var[1][j+d.Nx+1]; //-
	      //v.var[2][j]=-v.var[2][j+d.Nx+1]; //-
	      //v.var[3][j]=v.var[3][j+d.Nx+1];
	   }
	    //coin bas droit
	    if(Mh.xj(j).y < 0 && Mh.xj(j).x>d.Tx){
	      //v.var[0][j]=v.var[0][j+d.Nx-1];
	      //v.var[1][j]=-v.var[1][j+d.Nx-1]; //-
	      //v.var[2][j]=-v.var[2][j+d.Nx-1]; //-
	      //v.var[3][j]=v.var[3][j+d.Nx-1];
	    }

    }
  }
}

/**
void BC_BlackBodyCompton_Left(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time){
    Function which initialize the phantom cells to a black body in the left of the domain 
   int in=0;
  in=centerdomain(d,Mh);
  
  for(int j=0;j<Mh.nc;j++){
    if(Mh.cells[j].lab==-1){
	    //cote gauche
    
	    if((Mh.xj(j).x < 0) && (Mh.xj(j).y> 0 && Mh.xj(j).y< d.Ty)){
	      v.var[0][j]=Param.P1M.a_value*pow(Param.P1M.TempBB,4.);
	     v.var[1][j]=0;	//-
	     v.var[2][j]=0;
	     v.var[3][j]=Param.P1M.TempBB;
	   }
	    //cote droit
	    if(Mh.xj(j).x> d.Tx && (Mh.xj(j).y> 0 && Mh.xj(j).y< d.Ty)){
	     v.var[0][j]=v.var[0][j-1];
	     v.var[1][j]=-v.var[1][j-1];
	     v.var[2][j]=v.var[2][j-1]; //-
	     v.var[3][j]=v.var[3][j-1];
	   }
	    //cote bas
	    if(Mh.xj(j).y< 0 && (Mh.xj(j).x> 0  && Mh.xj(j).x< d.Tx)){
	      v.var[0][j]=v.var[0][j+(d.Nx+2)];
	      v.var[1][j]=v.var[1][j+(d.Nx+2)];
	      v.var[2][j]=-v.var[2][j+(d.Nx+2)]; //-
	      v.var[3][j]=v.var[3][j+(d.Nx+2)];
	   }
	    //cote haut
	    if(Mh.xj(j).y> d.Ty && (Mh.xj(j).x > 0  && Mh.xj(j).x < d.Tx)){
	      v.var[0][j]=v.var[0][j-(d.Nx+2)];
	      v.var[1][j]=v.var[1][j-(d.Nx+2)];
	      v.var[2][j]=-v.var[2][j-(d.Nx+2)]; //-
	      v.var[3][j]=-v.var[3][j-(d.Nx+2)]; //-
	    }

      //coin bas gauche

             if(Mh.xj(j).y< 0 && Mh.xj(j).x<0){
	       v.var[0][j]=Param.P1C.a_value*pow(Param.P1C.TempBB,4.);
		v.var[1][j]=0; // -
	       v.var[2][j]=0; //-
	       v.var[3][j]=Param.P1C.TempBB;
	   }
	   //coin haut gauche
	    if(Mh.xj(j).y> d.Ty && Mh.xj(j).x<0){
	      v.var[0][j]=Param.P1C.a_value*pow(Param.P1C.TempBB,4.);
		v.var[1][j]=0; // -
	       v.var[2][j]=0; //-
	       v.var[3][j]=Param.P1C.TempBB;

	   }
	    //coin haut droit
	    if(Mh.xj(j).y>d.Ty && Mh.xj(j).x>d.Tx){
	      //v.var[0][j]=v.var[0][j+d.Nx+1];
	      //v.var[1][j]=-v.var[1][j+d.Nx+1]; //-
	      //v.var[2][j]=-v.var[2][j+d.Nx+1]; //-
	      //v.var[3][j]=v.var[3][j+d.Nx+1];
	   }
	    //coin bas droit
	    if(Mh.xj(j).y < 0 && Mh.xj(j).x>d.Tx){
	      //v.var[0][j]=v.var[0][j+d.Nx-1];
	      //v.var[1][j]=-v.var[1][j+d.Nx-1]; //-
	      //v.var[2][j]=-v.var[2][j+d.Nx-1]; //-
	      //v.var[3][j]=v.var[3][j+d.Nx-1];
	    }

    }
  }
}
**/
