#ifndef PARAMADVECTION_HPP
#define PARAMADVECTION_HPP

#include <cmath>
#include <cstdlib>
#include <iostream> 
#include "R2.hpp"
#include "Data.hpp"
#include "ClassMesh.hpp"
#include <cassert>
#include "FunctionsConnect.hpp"
#include "Variables.hpp"
#include <fstream>
#include <string.h>
#include <stdlib.h>

using namespace std;

/** Class: physical parameter for advection model 
    - nbvar: number of variable for the model
    - nbcell: number of cells
    - a: table of velocity coefficient for each cell
    - Nlim : number of limiter
    - OrderAdv : order of advection scheme (MUSCL for the second order)
    - velocity : advection velocity (-100000 if the velocity is variable)

**/ 

class ParamAdvection {
public:
  int nbvar;
  R2 *a;
  R2 velocity;
  int Nlim;
  int OrderAdv;
  int nbcell; 
  
  ParamAdvection(){
    nbvar=1;
    nbcell=0;
    velocity.x=0.;
    velocity.y=0.;
    a=NULL;
    OrderAdv=1;
    Nlim=1;
  }
  
  ParamAdvection(Data & d,Mesh & Mh){ 
    nbvar=1;
    nbcell=Mh.nc;
    velocity.x=0.;
    velocity.y=0.;    
    a=new R2[nbcell];
    OrderAdv=1;
    Nlim=1;
  }

  ParamAdvection(const ParamAdvection & pa){
    nbvar=pa.nbvar;
    nbcell=pa.nbcell;
    velocity.x=pa.velocity.x;
    velocity.y=pa.velocity.y;
    OrderAdv=pa.OrderAdv;
    Nlim=pa.Nlim;
    a = new R2[nbcell];
    for(int i=0;i<nbvar;i++){
      a[i]= pa.a[i];
    }    
  }
  
  ~ParamAdvection(){
      if(a!=NULL){
      delete [] a;
      }      
  }
	
 ParamAdvection & operator=( ParamAdvection pa){
    nbvar=pa.nbvar;
    nbcell=pa.nbcell;
    OrderAdv=pa.OrderAdv;
    velocity.x=pa.velocity.x;
    velocity.y=pa.velocity.y;
    Nlim=pa.Nlim;
    if(a==NULL){
      a = new R2[nbcell];
    } 
    for(int i=0;i<nbvar;i++){
      a[i]= pa.a[i];
    }    
    return *this;
 }    
}; 


void ParamAdv_Init(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamAdvection & Adv, const char * filename);

void ParamAdv_Initwithfile(ParamAdvection & Adv, const char * filename);

void ParamAdv_InitTab(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamAdvection & Adv);

R2 AverageVelocityNode(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamAdvection & Adv,int numGr,int j);
#endif

