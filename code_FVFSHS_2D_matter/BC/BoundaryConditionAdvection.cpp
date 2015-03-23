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

void BoundaryConditionAdvection(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time){
 /** Function which initialize the phantom cells for the advection model **/
      switch(d.nTest)
    {
    case 1 : BCDirichlet_SolExact_Advection(d,Mh,v,Param,time);

      break;}
}





void BCDirichlet_SolExact_Advection(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time){
   /** Function which initialize the phantom cells with the exact solution **/
  int in=0;
  in=centerdomain(d,Mh);

  for(int j=0;j<Mh.nc;j++){
    if(Mh.cells[j].lab==-1){
	  
	    //cote gauche

	    if((Mh.xj(j).x < 0) && (Mh.xj(j).y> 0 && Mh.xj(j).y< d.Ty)){
	      v.var[0][j]=SolFon(d,Param,Mh.xj(j),time,Mh.xj(in),0);
             
	   }
	    //cote droit
	    if(Mh.xj(j).x> d.Tx && (Mh.xj(j).y> 0 && Mh.xj(j).y< d.Ty)){
	      v.var[0][j]=SolFon(d,Param,Mh.xj(j),time,Mh.xj(in),0);
	   }
	    //cote bas
	    if(Mh.xj(j).y< 0 && (Mh.xj(j).x> 0  && Mh.xj(j).x< d.Tx)){
	      v.var[0][j]=SolFon(d,Param,Mh.xj(j),time,Mh.xj(in),0);
	   }
	    //cote haut
	    if(Mh.xj(j).y> d.Ty && (Mh.xj(j).x > 0  && Mh.xj(j).x < d.Tx)){
	      v.var[0][j]=SolFon(d,Param,Mh.xj(j),time,Mh.xj(in),0);
	    }

      //coin bas gauche

             if(Mh.xj(j).y< 0 && Mh.xj(j).x<0){
	       v.var[0][j]=SolFon(d,Param,Mh.xj(j),time,Mh.xj(in),0);
	   }
	   //coin haut gauche
	    if(Mh.xj(j).y> d.Ty && Mh.xj(j).x<0){
	      v.var[0][j]=SolFon(d,Param,Mh.xj(j),time,Mh.xj(in),0);

	   }
	    //coin haut droit
	    if(Mh.xj(j).y>d.Ty && Mh.xj(j).x>d.Tx){
	      v.var[0][j]=SolFon(d,Param,Mh.xj(j),time,Mh.xj(in),0);
	   }
	    //coin bas droit
	    if(Mh.xj(j).y < 0 && Mh.xj(j).x>d.Tx){
	      v.var[0][j]=SolFon(d,Param,Mh.xj(j),time,Mh.xj(in),0);
	    }

    }
  }
}
