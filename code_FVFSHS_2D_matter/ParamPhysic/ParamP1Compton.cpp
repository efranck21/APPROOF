#include <cmath>
#include <cstdlib>
#include <iostream> 
#include "ParamP1Compton.hpp"
#include <cassert>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;


/** File which contained some functions to intitialize the paramater of P1 model with matter **/ 

void ParamP1Matter_Init(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1Compton & P1C, const char * filename){
/** Initialization of P1 with matter parameter **/
  ParamP1Compton_Initwithfile(P1C,filename);
  ParamP1Compton_InitTab(d,Mh,v,tab,P1C);
  ParamP1Compton_InitPhysicalTestcase(d,Mh,v,tab,P1C); 

}

void ParamP1Compton_Initwithfile(ParamP1Compton & P1C, const char * filename){
 /** Use the file ParamP1M.dat to initialize the P1 with matter parameter **/

  
    FILE* f = fopen(filename,"rb");
    if(!f){cout << "Error when we open file"<<filename<<endl;exit(1);}
     cout <<"Reading of file\"" << filename << "\""<<endl;
    char string[255];

  /* sigma*/
  fgets(string,255,f);  fscanf(f,"%lf\n",&P1C.sig_value);

  /* sigma absorbtion*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&P1C.sigA_value);

   /* a (stefan constant)*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&P1C.a_value);

  /* Cv (matter constant)*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&P1C.Cv_value); 

   /* Temp_BB (Black Body temperature)*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&P1C.TempBB);

   /* mc2*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&P1C.mc2);
  
  cout << "End of reading Paramp1M.dat "<<endl;
 
  }

void ParamP1Compton_InitTab(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1Compton & P1C){
/** Function which initialize the table of P1 matter quantities using the given values (paramP1M.dat)
or the test case if the eps, sigma are variables **/
  
   if (P1C.sig_value != -1){
     for(int i=0;i<P1C.nbcell;i++){
      P1C.sigma[i]=P1C.sig_value;
    }
  }
    if (P1C.sig_value != -1){
     for(int i=0;i<P1C.nbcell;i++){
      P1C.sigmaA[i]=P1C.sigA_value;
    }
  }

     
}


void ParamP1Compton_InitPhysicalTestcase(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1Compton & P1C){
  /** This function modify the parameter to have parameter consistant with the test case **/
 
  
 }

double AverageSigNodeP1Compton(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamP1Compton & P1C,int numGr,int j){
  /** Function which compute an average of sigma at the node **/
  double sigma=0;
  int jG1;
  int Typemean=1;

  double m1=0;

  if(Typemean ==1 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
       jG1=tab.TabInv[numGr].TabCell[i];
       m1=m1+P1C.sigma[jG1]; 
     }
    
    sigma=m1/tab.TabInv[numGr].taille; 
  }
  
  m1=100000000000000;
  if(Typemean ==2 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
      jG1=tab.TabInv[numGr].TabCell[i];
      m1=min(m1,P1C.sigma[jG1]); 
    }
    sigma=m1; 
  }
  
  return sigma;
}

double AverageSigANodeP1Compton(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamP1Compton & P1C,int numGr,int j){
  /** Function which compute an average of sigma absortion at the node **/
  double sigma=0;
  int jG1;
  int Typemean=1;

  double m1=0;

  if(Typemean ==1 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
       jG1=tab.TabInv[numGr].TabCell[i];
       m1=m1+P1C.sigmaA[jG1]; 
     }
    
    sigma=m1/tab.TabInv[numGr].taille; 
  }
  
  m1=100000000000000;
  if(Typemean ==2 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
      jG1=tab.TabInv[numGr].TabCell[i];
      m1=min(m1,P1C.sigmaA[jG1]); 
    }
    sigma=m1; 
  }
  
  return sigma;
}

