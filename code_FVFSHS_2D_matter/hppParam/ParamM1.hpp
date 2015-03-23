#ifndef PARAMM1_HPP
#define PARAMM1_HPP

#include <cmath>
#include <cstdlib>
#include <iostream> 
#include "R2.hpp"
#include "Data.hpp"
#include "ClassMesh.hpp"
#include "FunctionsConnect.hpp"
#include "Variables.hpp"
#include <cassert>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;

/** Class: physical parameter for M1 model 
    - nbvar: number of variable for the model
    - nbcell: number of cells
    - eps: table of epsilon coefficient for each cell
    - sigma: table of sigma coefficient for each cell
    - eps_value : epsilon value (-100000 if eps is variable)
    - sig_value : sigma value (-100000 if sig is variable)
    - OrderAdv : order of remap scheme (MUSCL for the second order)
**/


class ParamM1 {
public:
  int nbvar;
  double *eps;
  double *sigma;
  int nbcell;
  double eps_value;
  double sig_value;
  int OrderAdv;
  
  ParamM1(){
    nbvar=3;
    eps_value=0;
    sig_value=0;
    nbcell=0;
    OrderAdv=0; 
    eps=NULL;
    sigma=NULL;
  }
  
  ParamM1(Data & d,Mesh & Mh){ 
    nbvar=3;
    nbcell=Mh.nc;
    OrderAdv=1;
    eps_value=0;
    sig_value=0;
    eps=new double[nbcell];
    sigma=new double[nbcell];
  }

  ParamM1(const ParamM1 & pm){
    nbvar=pm.nbvar;
    nbcell=pm.nbcell;
    OrderAdv=pm.OrderAdv;
    eps_value=pm.eps_value;
    sig_value=pm.sig_value;
    eps = new double[nbcell];
    for(int i=0;i<nbcell;i++){
      eps[i]= pm.eps[i];
    }    
    sigma = new double[nbcell];
    for(int i=0;i<nbcell;i++){
      sigma[i]= pm.sigma[i];
    }  
  }
  
  ~ParamM1(){
      if(eps!=NULL){
	delete [] eps;
      }     
       if(sigma!=NULL){
	 delete [] sigma;
      }  
  }
	
 ParamM1 & operator=( ParamM1 pm){
    nbvar=pm.nbvar;
    nbcell=pm.nbcell;
    OrderAdv=pm.OrderAdv;
    eps_value=pm.eps_value;
    sig_value=pm.sig_value;
    if(eps==NULL){
      eps = new double[nbcell];
    }
    for(int i=0;i<nbcell;i++){
      eps[i]= pm.eps[i];
    }    
    if(sigma==NULL){
      sigma = new double[nbcell];
    }
    for(int i=0;i<nbcell;i++){
      sigma[i]= pm.sigma[i];
    }  
    return *this;
 }    

}; 

void ParamM1_Init(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamM1 & M1, const char * filename);

void ParamM1_Initwithfile(ParamM1 & M1, const char * filename);

void ParamM1_InitTab(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamM1 & M1);

void ParamM1_InitPhysicalTestcase(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamM1 & M1);

double AverageSigNodeM1(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamM1 & M1,int numGr,int j);

double AverageEpsNodeM1(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamM1 & M1,int numGr,int j);

double AverageEnergyNodeM1(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamM1 & M1,int numGr,int j);

R2 AverageFluxNodeM1(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamM1 & M1,int numGr,int j);

#endif
