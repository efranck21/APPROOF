#ifndef PARAMP1COMPTON_HPP
#define PARAMP1COMPTON_HPP

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

class ParamP1Compton {
public:
  int nbvar;
  int nb_group;
  int nb_moment;
  double *sigma;
  double * sigmaA;
  int nbcell;
  double sig_value;
  double sigA_value;
  double a_value;
  double Cv_value;
  double TempBB;
  double mc2;
  
  ParamP1Compton(){
    nb_moment=3;
    nb_group=1;
    nbvar=nb_moment*nb_group+1;
    sig_value=0;
    sigA_value=0;
    nbcell=0;
    sigma=NULL;
    sigmaA=NULL;
    a_value=0;
    Cv_value=1;
    TempBB=0;
    mc2=0;
  }
  
  ParamP1Compton(Data & d,Mesh & Mh){
    nb_moment=3;
    nb_group=d.ngroup;
    nbvar=nb_moment*nb_group+1;
    nbcell=Mh.nc;
    sig_value=0;
    sigA_value=0;
    sigma=new double[nbcell];
    sigmaA=new double[nbcell];
    a_value=1;
    Cv_value=1;
    TempBB=0;
    mc2=0;
  }

  ParamP1Compton(const ParamP1Compton & ppc){
    nb_group=ppc.nb_group;
    nb_moment=ppc.nb_moment;
    nbvar=ppc.nbvar;
    nbcell=ppc.nbcell;
    sig_value=ppc.sig_value;
    sigA_value=ppc.sigA_value;
    a_value=ppc.a_value;
    Cv_value=ppc.Cv_value;
    TempBB=ppc.TempBB;
    mc2=ppc.mc2;
      
    sigma = new double[nbcell];
    for(int i=0;i<nbcell;i++){
      sigma[i]= ppc.sigma[i];
    }
     sigmaA = new double[nbcell];
    for(int i=0;i<nbcell;i++){
      sigmaA[i]= ppc.sigma[i];
    }  
  }
  
  ~ParamP1Compton(){
   
       if(sigma!=NULL){
	 delete [] sigma;
      }
        if(sigmaA!=NULL){
	 delete [] sigmaA;
      }  
  }
	
 ParamP1Compton & operator=( ParamP1Compton ppc){
   nb_group=ppc.nb_group;
   nb_moment=ppc.nb_moment;
    nbvar=ppc.nbvar;
    nbcell=ppc.nbcell;
    sig_value=ppc.sig_value;
    sigA_value=ppc.sigA_value;
    a_value=ppc.a_value;
    Cv_value=ppc.Cv_value;
    TempBB=ppc.TempBB;
    mc2=ppc.mc2;

    if(sigma==NULL){
      sigma = new double[nbcell];
    }
    for(int i=0;i<nbcell;i++){
      sigma[i]= ppc.sigma[i];
    }
     if(sigmaA==NULL){
      sigmaA = new double[nbcell];
    }
    for(int i=0;i<nbcell;i++){
      sigmaA[i]= ppc.sigmaA[i];
    } 
    return *this;
 }    

}; 

void ParamP1Compton_Init(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1Compton & P1C, const char * filename);

void ParamP1Compton_Initwithfile(ParamP1Compton & P1C, const char * filename);

void ParamP1Compton_InitTab(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1Compton & P1C);

void ParamP1Compton_InitPhysicalTestcase(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1Compton & P1C);

double AverageSigNodeP1Compton(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1Compton & P1C,int numGr,int j);

double AverageSigANodeP1Compton(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1Compton & P1C,int numGr,int j);


#endif
