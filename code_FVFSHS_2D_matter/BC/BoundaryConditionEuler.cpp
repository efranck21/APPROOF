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



void BoundaryConditionEuler(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time){
 /** Function which initialize the phantom cells for the Euler model.
For triangular mesh only the Boundary conditions based on th exact solutions are written **/
  
  if (d.Typemesh=='Q'){
    switch(d.nTest)
    {
    case 1 :
      if(Param.Euler.LHO <2) {
	  BCDirichlet_ST_Euler(d,Mh,v,Param,time);
        }
	else {
	  BCDirichlet_ST_EulerAverage(d,Mh,v,Param,time);
        }
      
      break;
      
    case 2 : BCDirichlet_ST_Euler(d,Mh,v,Param,time);
      
      break;

    case 3 :
        if(Param.Euler.LHO <2) {
	  BCDirichlet_ST_Euler(d,Mh,v,Param,time);
        }
	else {
	  BCDirichlet_ST_EulerAverage(d,Mh,v,Param,time);
        }
    
      break;
 
    case 4 : BCDirichlet_ST_Euler(d,Mh,v,Param,time);

      break;

    case 5 : 
      if(Param.Euler.LHO <2) {
	  BCDirichlet_ST_Euler(d,Mh,v,Param,time);
        }
	else {
	  BCDirichlet_ST_EulerAverage(d,Mh,v,Param,time);
        }
      break;
 
    case 6 : BCDirichlet_ST_Euler(d,Mh,v,Param,time);

      break;

    case 8 : BCNeumannEuler(d,Mh,v,Param,time);

      break;

    case 7 : BC_Periodic_Euler(d,Mh,v,Param,time);

      break;
      
    case 9 : BC_NeumanY_PeriodicX_Euler(d,Mh,v,Param,time);

      break;}
  }   
  
  if (d.Typemesh=='T' && d.nTest < 7){
    BCDirichlet_P1(d,Mh,v,Param,time);
  }
  if (d.Typemesh=='T' && d.nTest > 6){
    cout<<"Condition limit for this test case does not exist on triangular meshes"<<endl;
    exit(0);    
  }

}


void BC_NeumanY_PeriodicX_Euler(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time){
  /** Function which initialize the phantom cells to have Neumann condition **/
    R2 g(0,0);
  g=GravityVector(d,Param.Euler);
  for(int j=0;j<Mh.nc;j++){
    if(Mh.cells[j].lab==-1){
      
	  
	    //cote gauche

	    if((Mh.xj(j).x < 0) && (Mh.xj(j).y> 0 && Mh.xj(j).y< d.Ty)){
	     v.var[0][j]=v.var[0][j+1];
	     v.var[1][j]=-v.var[1][j+1];	//-
	     v.var[2][j]=v.var[2][j+1];
	     v.var[3][j]=v.var[3][j+1];
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
	      v.var[0][j]=1;
	      v.var[1][j]=0;
	      v.var[2][j]=0; //-
	      v.var[3][j]=(1./(Param.Euler.gamma-1.))*(2*(g.y*d.Ty+g.y*d.Tx+1)-(g.y*Mh.xj(j).y+g.x*Mh.xj(j).x)*v.var[0][j]);
	   }
	    //cote haut
	    if(Mh.xj(j).y> d.Ty && (Mh.xj(j).x > 0  && Mh.xj(j).x < d.Tx)){
	      v.var[0][j]=2;
	      v.var[1][j]=0;
	      v.var[2][j]=0; //-
	      v.var[3][j]=(1./(Param.Euler.gamma-1.))*(2*(g.y*d.Ty+g.y*d.Tx+1)-(g.y*Mh.xj(j).y+g.x*Mh.xj(j).x)*v.var[0][j]);
	    }
 
             
      //coin bas gauche

             if(Mh.xj(j).y< 0 && Mh.xj(j).x<0){
	       v.var[0][j]=1;
	      v.var[1][j]=0;
	      v.var[2][j]=0; //-
	      v.var[3][j]=(1./(Param.Euler.gamma-1.))*(2*(g.y*d.Ty+g.y*d.Tx+1)-(g.y*Mh.xj(j).y+g.x*Mh.xj(j).x)*v.var[0][j]);
	   }
	   //coin haut gauche
	    if(Mh.xj(j).y> d.Ty && Mh.xj(j).x<0){
	      v.var[0][j]=2;
	      v.var[1][j]=0;
	      v.var[2][j]=0; //-
	      v.var[3][j]=(1./(Param.Euler.gamma-1.))*(2*(g.y*d.Ty+g.y*d.Tx+1)-(g.y*Mh.xj(j).y+g.x*Mh.xj(j).x)*v.var[0][j]);

	   }
	    //coin haut droit
	    if(Mh.xj(j).y>d.Ty && Mh.xj(j).x>d.Tx){
	      v.var[0][j]=2;
	      v.var[1][j]=0;
	      v.var[2][j]=0; //-
	      v.var[3][j]=(1./(Param.Euler.gamma-1.))*(2*(g.y*d.Ty+g.y*d.Tx+1)-(g.y*Mh.xj(j).y+g.x*Mh.xj(j).x)*v.var[0][j]);
	   }
	    //coin bas droit
	    if(Mh.xj(j).y < 0 && Mh.xj(j).x>d.Tx){
	      v.var[0][j]=1;
	      v.var[1][j]=0;
	      v.var[2][j]=0; //-
	      v.var[3][j]=(1./(Param.Euler.gamma-1.))*(2*(g.y*d.Ty+g.y*d.Tx+1)-(g.y*Mh.xj(j).y+g.x*Mh.xj(j).x)*v.var[0][j]);
	    }



    }
  }
}

    

void BCNeumannEuler(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time){
  /** Function which initialize the phantom cells to have Neumann condition **/
  for(int j=0;j<Mh.nc;j++){
    if(Mh.cells[j].lab==-1){
      
      
      //cote gauche
      if((Mh.xj(j).x < 0) && (Mh.xj(j).y> 0 && Mh.xj(j).y< d.Ty)){      
	v.var[0][j]=v.var[0][j+d.Nx];
	v.var[1][j]=v.var[1][j+d.Nx];	//-
	v.var[2][j]=v.var[2][j+d.Nx];
	v.var[3][j]=v.var[3][j+d.Nx];
      }
      //cote droit
      if(Mh.xj(j).x> d.Tx && (Mh.xj(j).y> 0 && Mh.xj(j).y< d.Ty)){
	v.var[0][j]=v.var[0][j-(d.Nx)];
	v.var[1][j]=v.var[1][j-(d.Nx)];	//-
	v.var[2][j]=v.var[2][j-(d.Nx)];
	v.var[3][j]=v.var[3][j-(d.Nx)];
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
	v.var[3][j]=v.var[3][j-(d.Nx+2)];
      }
      
             
      //coin bas gauche

      if(Mh.xj(j).y< 0 && Mh.xj(j).x<0){
	v.var[0][j]=v.var[0][j+d.Nx+2];
	v.var[1][j]=v.var[1][j+d.Nx+2]; // -
	v.var[2][j]=-v.var[2][j+d.Nx+2]; //-
	v.var[3][j]=v.var[3][j+d.Nx+2];
      }
	   //coin haut gauche
      if(Mh.xj(j).y< 0 && Mh.xj(j).x<0){
	v.var[0][j]=v.var[0][j-d.Nx+2];
	v.var[1][j]=v.var[1][j-d.Nx+2]; // -
	v.var[2][j]=-v.var[2][j-d.Nx+2]; //-
	v.var[3][j]=v.var[3][j-d.Nx+2];
      }
	    //coin haut droit
      if(Mh.xj(j).y>d.Ty && Mh.xj(j).x>d.Tx){
	v.var[0][j]=v.var[0][j-d.Nx+2];
	v.var[1][j]=v.var[1][j-d.Nx+2]; // -
	v.var[2][j]=-v.var[2][j-d.Nx+2]; //-
	v.var[3][j]=v.var[3][j-d.Nx+2];
      }
	    //coin bas droit
      if(Mh.xj(j).y < 0 && Mh.xj(j).x>d.Tx){
	v.var[0][j]=v.var[0][j+d.Nx+2];
	v.var[1][j]=v.var[1][j+d.Nx+2]; // -
	v.var[2][j]=-v.var[2][j+d.Nx+2]; //-
	v.var[3][j]=v.var[3][j+d.Nx+2];
      }



    }
  }
}








void BCDirichlet_ST_Euler(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time){
    /** Function which initialize the phantom cells with the exact solution **/
  int in=0;
  in=centerdomain(d,Mh);
  for(int j=0;j<Mh.nc;j++){
    if(Mh.cells[j].lab==-1){
	  
	    //cote gauche

	    if((Mh.xj(j).x < 0) && (Mh.xj(j).y> 0 && Mh.xj(j).y< d.Ty)){
	      v.var[0][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),0);
	      v.var[1][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),1);	//-
	      v.var[2][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),2);
	      v.var[3][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),3);
	   }
	    //cote droit
	    if(Mh.xj(j).x> d.Tx && (Mh.xj(j).y> 0 && Mh.xj(j).y< d.Ty)){
	      v.var[0][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),0);
	      v.var[1][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),1);	//-
	      v.var[2][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),2);
	      v.var[3][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),3);
	   }
	    //cote bas
	    if(Mh.xj(j).y< 0 && (Mh.xj(j).x> 0  && Mh.xj(j).x< d.Tx)){
	      v.var[0][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),0);
	      v.var[1][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),1);	//-
	      v.var[2][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),2);
	      v.var[3][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),3);
	   }
	    //cote haut
	    if(Mh.xj(j).y> d.Ty && (Mh.xj(j).x > 0  && Mh.xj(j).x < d.Tx)){
	      v.var[0][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),0);
	      v.var[1][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),1);	//-
	      v.var[2][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),2);
	      v.var[3][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),3);
	    }
 
             
      //coin bas gauche

             if(Mh.xj(j).y< 0 && Mh.xj(j).x<0){
	       v.var[0][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),0);
	       v.var[1][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),1);	//-
	       v.var[2][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),2);
	       v.var[3][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),3);
	   }
	   //coin haut gauche
	    if(Mh.xj(j).y> d.Ty && Mh.xj(j).x<0){
	      v.var[0][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),0);
	      v.var[1][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),1);	//-
	      v.var[2][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),2);
	      v.var[3][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),3);

	   }
	    //coin haut droit
	    if(Mh.xj(j).y>d.Ty && Mh.xj(j).x>d.Tx){
	      v.var[0][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),0);
	      v.var[1][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),1);	//-
	      v.var[2][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),2);
	      v.var[3][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),3);
	   }
	    //coin bas droit
	    if(Mh.xj(j).y < 0 && Mh.xj(j).x>d.Tx){
	      v.var[0][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),0);
	      v.var[1][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),1);	//-
	      v.var[2][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),2);
	      v.var[3][j]=SolFonEuler(d,Param,Mh.xj(j),time,Mh.xj(in),3);
	    }



    }
  }
}

void BCDirichlet_ST_EulerAverage(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time){
    /** Function which initialize the phantom cells with the exact solution **/
  int in=0;
  in=centerdomain(d,Mh);
  for(int j=0;j<Mh.nc;j++){
    if(Mh.cells[j].lab==-1){
	  
	    //cote gauche

	    if((Mh.xj(j).x < 0) && (Mh.xj(j).y> 0 && Mh.xj(j).y< d.Ty)){
	      v.var[0][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,0);
	      v.var[1][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,1);	//-
	      v.var[2][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,2);
	      v.var[3][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,3);
	   }
	    //cote droit
	    if(Mh.xj(j).x> d.Tx && (Mh.xj(j).y> 0 && Mh.xj(j).y< d.Ty)){
	      v.var[0][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,0);
	      v.var[1][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,1);	//-
	      v.var[2][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,2);
	      v.var[3][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,3);
	   }
	    //cote bas
	    if(Mh.xj(j).y< 0 && (Mh.xj(j).x> 0  && Mh.xj(j).x< d.Tx)){
	      v.var[0][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,0);
	      v.var[1][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,1);	//-
	      v.var[2][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,2);
	      v.var[3][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,3);
	   }
	    //cote haut
	    if(Mh.xj(j).y> d.Ty && (Mh.xj(j).x > 0  && Mh.xj(j).x < d.Tx)){
	      v.var[0][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,0);
	      v.var[1][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,1);	//-
	      v.var[2][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,2);
	      v.var[3][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,3);
	    }
 
             
      //coin bas gauche

             if(Mh.xj(j).y< 0 && Mh.xj(j).x<0){
	       v.var[0][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,0);
	       v.var[1][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,1);	//-
	       v.var[2][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,2);
	       v.var[3][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,3);
	   }
	   //coin haut gauche
	    if(Mh.xj(j).y> d.Ty && Mh.xj(j).x<0){
	      v.var[0][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,0);
	      v.var[1][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,1);	//-
	      v.var[2][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,2);
	      v.var[3][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,3);

	   }
	    //coin haut droit
	    if(Mh.xj(j).y>d.Ty && Mh.xj(j).x>d.Tx){
	      v.var[0][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,0);
	      v.var[1][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,1);	//-
	      v.var[2][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,2);
	      v.var[3][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,3);
	   }
	    //coin bas droit
	    if(Mh.xj(j).y < 0 && Mh.xj(j).x>d.Tx){
	      v.var[0][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,0);
	      v.var[1][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,1);	//-
	      v.var[2][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,2);
	      v.var[3][j]=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),time,Mh.xj(in),j,3);
	    }



    }
  }
}


void BC_Periodic_Euler(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time){
   /** Function which initialize the phantom cells to have periodic solution **/
  for(int j=0;j<Mh.nc;j++){
  
    if(Mh.cells[j].lab==-1){
    
      
	    //cote gauche
      if((Mh.xj(j).x < 0) && (Mh.xj(j).y> 0 && Mh.xj(j).y< d.Ty)){      
	      v.var[0][j]=v.var[0][j+d.Nx];
	      v.var[1][j]=v.var[1][j+d.Nx];	//-
	      v.var[2][j]=v.var[2][j+d.Nx];
	      v.var[3][j]=v.var[3][j+d.Nx];
	   }
	    //cote droit
      if(Mh.xj(j).x> d.Tx && (Mh.xj(j).y> 0 && Mh.xj(j).y< d.Ty)){
	     v.var[0][j]=v.var[0][j-(d.Nx)];
	     v.var[1][j]=v.var[1][j-(d.Nx)];	//-
	     v.var[2][j]=v.var[2][j-(d.Nx)];
	     v.var[3][j]=v.var[3][j-(d.Nx)];
	   }
	    //cote bas
      if(Mh.xj(j).y< 0 && (Mh.xj(j).x> 0  && Mh.xj(j).x< d.Tx)){
	     v.var[0][j]=v.var[0][j+(d.Ny-1)*(d.Nx+2)+d.Nx+2];
	     v.var[1][j]=v.var[1][j+(d.Ny-1)*(d.Nx+2)+d.Nx+2];	//-
	     v.var[2][j]=v.var[2][j+(d.Ny-1)*(d.Nx+2)+d.Nx+2];
	     v.var[3][j]=v.var[3][j+(d.Ny-1)*(d.Nx+2)+d.Nx+2];
	   }
	    //cote haut
	    if(Mh.xj(j).y> d.Ty && (Mh.xj(j).x > 0  && Mh.xj(j).x < d.Tx)){
              
	      v.var[0][j]=v.var[0][j-((d.Ny-1)*(d.Nx+2)+d.Nx+2)];
	      v.var[1][j]=v.var[1][j-((d.Ny-1)*(d.Nx+2)+d.Nx+2)];	//-
	      v.var[2][j]=v.var[2][j-((d.Ny-1)*(d.Nx+2)+d.Nx+2)];
	      v.var[3][j]=v.var[3][j-((d.Ny-1)*(d.Nx+2)+d.Nx+2)];
	    }
 
             
      //coin bas gauche

             if(Mh.xj(j).y< 0 && Mh.xj(j).x<0){
	       v.var[0][j]=v.var[0][j+(d.Nx+2)*(d.Ny+1)-2];
	       v.var[1][j]=v.var[1][j+(d.Nx+2)*(d.Ny+1)-2];	//-
	       v.var[2][j]=v.var[2][j+(d.Nx+2)*(d.Ny+1)-2];
	       v.var[3][j]=v.var[3][j+(d.Nx+2)*(d.Ny+1)-2];
	   }
	   //coin haut gauche
	       if(Mh.xj(j).y> d.Ty && Mh.xj(j).x<0){
		 v.var[0][j]=v.var[0][j-((d.Nx+2)*(d.Ny-1)+2)];
		 v.var[1][j]=v.var[1][j-((d.Nx+2)*(d.Ny-1)+2)];	//-
	       v.var[2][j]=v.var[2][j-((d.Nx+2)*(d.Ny-1)+2)];
	      v.var[3][j]=v.var[3][j-((d.Nx+2)*(d.Ny-1)+2)];

	   }
	    //coin haut droit
	       if(Mh.xj(j).y>d.Ty && Mh.xj(j).x>d.Tx){
	         v.var[0][j]=v.var[0][j-((d.Nx+2)*(d.Ny+1)-2)];
	      v.var[1][j]=v.var[1][j-(d.Nx+2)*(d.Ny+1)-2];	//-
	      v.var[2][j]=v.var[2][j-(d.Nx+2)*(d.Ny+1)-2];
	      v.var[3][j]=v.var[3][j-(d.Nx+2)*(d.Ny+1)-2];
	   }
	    //coin bas droit
	     if(Mh.xj(j).y < 0 && Mh.xj(j).x>d.Tx){
	       v.var[0][j]=v.var[0][j+(d.Nx+2)*(d.Ny-1)+2];
	      v.var[1][j]=v.var[1][j+(d.Nx+2)*(d.Ny-1)+2];	//-
	      v.var[2][j]=v.var[2][j+(d.Nx+2)*(d.Ny-1)+2];
	      v.var[3][j]=v.var[3][j+(d.Nx+2)*(d.Ny-1)+2];
	    }



    }
  }
}
