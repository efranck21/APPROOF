#ifndef PARAMP1_HPP
#define PARAMP1_HPP

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

/** Class: physical parameter for P1 model 
    - nbvar: number of variable for the model
    - nbcell: number of cells
    - eps: table of epsilon coefficient for each cell
    - sigma: table of sigma coefficient for each cell
    - eps_value : epsilon value (-100000 if eps is variable)
    - sig_value : sigma value (-100000 if sig is variable)  
**/


class ParamP1 {
public:
  int nbvar;
  double *eps;
  double *sigma;
  int nbcell;
  double eps_value;
  double sig_value;
  
  ParamP1(){
    nbvar=3;
    eps_value=0;
    sig_value=0;
    nbcell=0;
    eps=NULL;
    sigma=NULL;
  }
  
  ParamP1(Data & d,Mesh & Mh){ 
    nbvar=3;
    nbcell=Mh.nc;
    eps_value=0;
    sig_value=0;
    eps=new double[nbcell];
    sigma=new double[nbcell];
  }

  ParamP1(const ParamP1 & pp){
    nbvar=pp.nbvar;
    nbcell=pp.nbcell;
    eps_value=pp.eps_value;
    sig_value=pp.sig_value;
    eps = new double[nbcell];
    for(int i=0;i<nbcell;i++){
      eps[i]= pp.eps[i];
    }    
    sigma = new double[nbcell];
    for(int i=0;i<nbcell;i++){
      sigma[i]= pp.sigma[i];
    }  
  }
  
  ~ParamP1(){
      if(eps!=NULL){
	delete [] eps;
      }     
       if(sigma!=NULL){
	 delete [] sigma;
      }  
  }
	
 ParamP1 & operator=( ParamP1 pp){
    nbvar=pp.nbvar;
    nbcell=pp.nbcell;
    eps_value=pp.eps_value;
    sig_value=pp.sig_value;
    if(eps==NULL){
      eps = new double[nbcell];
    }
    for(int i=0;i<nbcell;i++){
      eps[i]= pp.eps[i];
    }    
    if(sigma==NULL){
      sigma = new double[nbcell];
    }
    for(int i=0;i<nbcell;i++){
      sigma[i]= pp.sigma[i];
    }  
    return *this;
 }    

}; 

void ParamP1_Init(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1 & P1, const char * filename);

void ParamP1_Initwithfile(ParamP1 & P1, const char * filename);

void ParamP1_InitTab(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1 & P1);

void ParamP1_InitPhysicalTestcase(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1 & P1);

double AverageSigNodeP1(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1 & P1,int numGr,int j);

double AverageEpsNodeP1(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1 & P1,int numGr,int j);

#endif
