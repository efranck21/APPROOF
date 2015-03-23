#ifndef PARAMDIFF_HPP
#define PARAMDIFF_HPP

#include <cmath>
#include <cstdlib>
#include <iostream> 
#include "R2.hpp"
#include "Data.hpp"
#include "ClassMesh.hpp"
#include "Variables.hpp"
#include "FunctionsConnect.hpp"
#include <cassert>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;

/** Class: physical parameter for diffusion model 
    - nbvar: number of variable for the model
    - nbcell: number of cells
    - D: table of diffusion coefficient for each cell
    - Coefdiff_value :diffusion coefficient (-100000 if the coefficient is variable)

**/

class ParamDiff {
public:
  int nbvar;
  int nbcell;
  double *D;
  int Coefdiff_value;
  
  ParamDiff(){
    nbvar=1;
    nbcell=0;
    Coefdiff_value=0;
    D=NULL;
  }
  
  ParamDiff(Data & d,Mesh & Mh){ 
    nbvar=1;
    nbcell=Mh.nc;
    Coefdiff_value=0;
    D=new double[nbcell];
  }

  ParamDiff(const ParamDiff & pd){
    nbvar=pd.nbvar;
    nbcell=pd.nbcell;
    Coefdiff_value=pd.Coefdiff_value;
    D = new double[nbcell];
    for(int i=0;i<nbcell;i++){
      D[i]= pd.D[i];
    }    
  }
  
  ~ParamDiff(){
      if(D!=NULL){
      delete [] D;
      }      
  }
	
 ParamDiff & operator=( ParamDiff pd){
    nbvar=pd.nbvar;
    nbcell=pd.nbcell;
    Coefdiff_value=pd.Coefdiff_value;
    if(D==NULL){
      D = new double[nbcell];
    } 
    for(int i=0;i<nbcell;i++){
      D[i]= pd.D[i];
    }    
    return *this;
 }    
}; 


void ParamDiff_Init(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab, ParamDiff & Diff, const char * filename);

void ParamDiff_Initwithfile(ParamDiff & Diff, const char * filename);

void ParamDiff_InitTab(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamDiff & Diff);

void ParamDiff_InitPhysicalTestcase(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamDiff & Diff);

double AverageCoefDiffNode(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamDiff & Diff,int numGr,int j);
#endif
