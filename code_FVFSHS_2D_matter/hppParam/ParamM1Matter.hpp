#ifndef PARAMM1M_HPP
#define PARAMM1M_HPP

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
 - sigma_A: table of sigma absorption coefficient for each cell
 - eps_value : epsilon value (-100000 if eps is variable)
 - sig_value : sigma value (-100000 if sig is variable)
 - mc2value: value of 1/mc2 for compton scattering
 - OrderAdv : order of remap scheme (MUSCL for the second order)
 **/


class ParamM1Matter {
public:
    int nbvar;
    double *eps;
    double *sigma;
    double *sigmaA ;
    int nbcell;
    double eps_value;
    double sig_value;
    double sigA_value;
    double a_value;
    double mc2value ;
    int OrderAdv;
    double TempBB;
    
    ParamM1Matter(){
        nbvar=4;
        eps_value=0;
        sig_value=0;
        sigA_value=0;
        mc2value=0;
        a_value=0.;
        nbcell=0;
        OrderAdv=0;
        eps=NULL;
        sigma=NULL;
        sigmaA=NULL;
        TempBB=0.;
    }
    
    ParamM1Matter(Data & d,Mesh & Mh){
        nbvar=4;
        nbcell=Mh.nc;
        OrderAdv=1;
        eps_value=0;
        sig_value=0;
        sigA_value=0;
        mc2value=0;
        a_value=0.;
        TempBB=0.;
        eps=new double[nbcell];
        sigma=new double[nbcell];
        sigmaA=new double[nbcell];
    }
    
    ParamM1Matter(const ParamM1Matter & pm){
        nbvar=pm.nbvar;
        nbcell=pm.nbcell;
        OrderAdv=pm.OrderAdv;
        eps_value=pm.eps_value;
        sig_value=pm.sig_value;
        sigA_value=pm.sigA_value;
        mc2value=pm.mc2value;
        a_value=pm.a_value;
        TempBB=pm.TempBB;
        eps = new double[nbcell];
        for(int i=0;i<nbcell;i++){
            eps[i]= pm.eps[i];
        }
        sigma = new double[nbcell];
        sigmaA = new double[nbcell];
        for(int i=0;i<nbcell;i++){
            sigma[i]= pm.sigma[i];
            sigmaA[i]= pm.sigmaA[i];
        }
    }
    
    ~ParamM1Matter(){
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
    
    ParamM1Matter & operator=( ParamM1Matter pm){
        nbvar=pm.nbvar;
        nbcell=pm.nbcell;
        OrderAdv=pm.OrderAdv;
        eps_value=pm.eps_value;
        sig_value=pm.sig_value;
        sigA_value=pm.sigA_value;
        mc2value=pm.mc2value;
        a_value=pm.a_value;
        TempBB=pm.TempBB;
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
        if(sigmaA==NULL){
            sigmaA = new double[nbcell];
        }
        for(int i=0;i<nbcell;i++){
            sigmaA[i]= pm.sigmaA[i];
        }
        return *this;
    }
    
};

void ParamM1M_Init(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamM1Matter & M1M, const char * filename);

void ParamM1M_Initwithfile(ParamM1Matter & M1M, const char * filename);

void ParamM1M_InitTab(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamM1Matter & M1M);

void ParamM1M_InitPhysicalTestcase(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamM1Matter & M1M);

double AverageSigNodeM1M(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamM1Matter & M1M,int numGr,int j);

double AverageEpsNodeM1M(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamM1Matter & M1M,int numGr,int j);

double AverageEnergyNodeM1M(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamM1Matter & M1M,int numGr,int j);

R2 AverageFluxNodeM1M(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamM1Matter & M1M,int numGr,int j);

#endif

