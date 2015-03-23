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

void BoundaryConditionM1Matter(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time){
    /** Function which initialize the phantom cells for the M1 model.
     For triangular mesh only the Boundary conditions based on th exact solutions are written **/
    if (d.Typemesh=='Q'){
        switch(d.nTest)
        {
            case 1 : BC_Periodic_M1M(d,Mh,v,Param,time);
                
                break;
                
            case 2 : BC_Periodic_M1M(d,Mh,v,Param,time);
                
                break;
                
            case 3 : BC_Periodic_M1M(d,Mh,v,Param,time);
                
                break;
                
            case 4 :BC_Periodic_M1M(d,Mh,v,Param,time);
                
            break;}
    }
    if (d.Typemesh=='T'){
        BCDirichlet_P1(d,Mh,v,Param,time);
    }
}

void BC_Neumann_M1M(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time){
    /** Function which initialize the phantom cells to have Neumann condition **/
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



void BC_Periodic_M1M(Data & d,Mesh & Mh,variable & v, ParamPhysic & Param,double time){
    /** Function which initialize the phantom cells to obtain periodic condition **/
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
