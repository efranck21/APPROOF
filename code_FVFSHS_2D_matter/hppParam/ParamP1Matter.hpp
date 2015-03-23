#ifndef PARAMP1MATTER_HPP
#define PARAMP1MATTER_HPP

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
    - sigma: table of sigma scattering coefficient for each cell
    - sigmaA: table of sigma absortion coefficient for each cell
    - eps_value : epsilon value (-100000 if eps is variable)
    - sig_value : sigma value (-100000 if sig is variable)
    - sigA_value : sigmaA value (-100000 if sig is variable)
    - Cv_value : RhoCv value for temperature equation
    - a_value : Stefan constant value
    - TempBB : temperature of the black body use in the boundary condition
**/

class ParamP1Matter {
public:
  int nbvar;
  double *eps;
  double *sigma;
  double * sigmaA;
  int nbcell;
  double eps_value;
  double sig_value;
  double sigA_value;
  double a_value;
  double Cv_value;
  double TempBB;
  
  ParamP1Matter(){
    nbvar=4;
    eps_value=0;
    sig_value=0;
    sigA_value=0;
    nbcell=0;
    eps=NULL;
    sigma=NULL;
    sigmaA=NULL;
    a_value=0;
    Cv_value=1;
    TempBB=0;
  }
  
  ParamP1Matter(Data & d,Mesh & Mh){ 
    nbvar=3;
    nbcell=Mh.nc;
    eps_value=0;
    sig_value=0;
    sigA_value=0;
    eps=new double[nbcell];
    sigma=new double[nbcell];
    sigmaA=new double[nbcell];
    a_value=0;
    Cv_value=1;
    TempBB=0; 
  }

  ParamP1Matter(const ParamP1Matter & ppm){
    nbvar=ppm.nbvar;
    nbcell=ppm.nbcell;
    eps_value=ppm.eps_value;
    sig_value=ppm.sig_value;
    sigA_value=ppm.sigA_value;
    a_value=ppm.a_value;
    Cv_value=ppm.Cv_value;
    TempBB=ppm.TempBB;
    eps = new double[nbcell];
    for(int i=0;i<nbcell;i++){
      eps[i]= ppm.eps[i];
    }    
    sigma = new double[nbcell];
    for(int i=0;i<nbcell;i++){
      sigma[i]= ppm.sigma[i];
    }
     sigmaA = new double[nbcell];
    for(int i=0;i<nbcell;i++){
      sigmaA[i]= ppm.sigma[i];
    }  
  }
  
  ~ParamP1Matter(){
      if(eps!=NULL){
	delete [] eps;
      }     
       if(sigma!=NULL){
	 delete [] sigma;
      }
        if(sigmaA!=NULL){
	 delete [] sigmaA;
      }  
  }
	
 ParamP1Matter & operator=( ParamP1Matter ppm){
    nbvar=ppm.nbvar;
    nbcell=ppm.nbcell;
    eps_value=ppm.eps_value;
    sig_value=ppm.sig_value;
    sigA_value=ppm.sigA_value;
    a_value=ppm.a_value;
    Cv_value=ppm.Cv_value;
    TempBB=ppm.TempBB;
    if(eps==NULL){
      eps = new double[nbcell];
    }
    for(int i=0;i<nbcell;i++){
      eps[i]= ppm.eps[i];
    }    
    if(sigma==NULL){
      sigma = new double[nbcell];
    }
    for(int i=0;i<nbcell;i++){
      sigma[i]= ppm.sigma[i];
    }
     if(sigmaA==NULL){
      sigmaA = new double[nbcell];
    }
    for(int i=0;i<nbcell;i++){
      sigmaA[i]= ppm.sigmaA[i];
    } 
    return *this;
 }    

}; 

void ParamP1Matter_Init(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1Matter & P1M, const char * filename);

void ParamP1Matter_Initwithfile(ParamP1Matter & P1M, const char * filename);

void ParamP1Matter_InitTab(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1Matter & P1M);

void ParamP1Matter_InitPhysicalTestcase(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1Matter & P1M);

double AverageSigNodeP1Matter(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1Matter & P1M,int numGr,int j);

double AverageSigANodeP1Matter(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1Matter & P1M,int numGr,int j);

double AverageEpsNodeP1Matter(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1Matter & P1M,int numGr,int j);

#endif
