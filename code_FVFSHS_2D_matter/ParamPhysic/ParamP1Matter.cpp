#include <cmath>
#include <cstdlib>
#include <iostream> 
#include "ParamP1Matter.hpp"
#include <cassert>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;


/** File which contained some functions to intitialize the paramater of P1 model with matter **/ 

void ParamP1Matter_Init(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1Matter & P1M, const char * filename){
/** Initialization of P1 with matter parameter **/
  ParamP1Matter_Initwithfile(P1M,filename);
  ParamP1Matter_InitTab(d,Mh,v,tab,P1M);
  ParamP1Matter_InitPhysicalTestcase(d,Mh,v,tab,P1M); 

}

void ParamP1Matter_Initwithfile(ParamP1Matter & P1M, const char * filename){
 /** Use the file ParamP1M.dat to initialize the P1 with matter parameter **/

  
    FILE* f = fopen(filename,"rb");
    if(!f){cout << "Error when we open file"<<filename<<endl;exit(1);}
     cout <<"Reading of file\"" << filename << "\""<<endl;
    char string[255];

  /* sigma*/
  fgets(string,255,f);  fscanf(f,"%lf\n",&P1M.sig_value);

   /* eps*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&P1M.eps_value);

  /* sigma absorbtion*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&P1M.sigA_value);

   /* a (stefan constant)*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&P1M.a_value);

  /* Cv (matter constant)*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&P1M.Cv_value); 

   /* Temp_BB (Black Body temperature)*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&P1M.TempBB);
  
  cout << "End of reading Paramp1M.dat "<<endl;
 
  }

void ParamP1Matter_InitTab(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1Matter & P1M){
/** Function which initialize the table of P1 matter quantities using the given values (paramP1M.dat)
or the test case if the eps, sigma are variables **/
  
  if (P1M.eps_value != -1){
     for(int i=0;i<P1M.nbcell;i++){
      P1M.eps[i]=P1M.eps_value;
    }
  }
   if (P1M.sig_value != -1){
     for(int i=0;i<P1M.nbcell;i++){
      P1M.sigma[i]=P1M.sig_value;
    }
  }
    if (P1M.sig_value != -1){
     for(int i=0;i<P1M.nbcell;i++){
      P1M.sigmaA[i]=P1M.sigA_value;
    }
  }

    if (d.nTest == 1){
      double c=29979245800;  
      for(int i=0;i<P1M.nbcell;i++){
	P1M.eps[i]=sqrt(3.)/c;
	P1M.sigmaA[i]=3./c;
	P1M.sigma[i]=0;
	   }   
    }
     if (d.nTest == 2){
      double c=299.79;  
      for(int i=0;i<P1M.nbcell;i++){
	P1M.eps[i]=sqrt(3.)/c;
	P1M.sigmaA[i]=(3./c)*(1./pow(v.var[3][i],3.));
	P1M.sigma[i]=0;
      }
     }
     
     
}


void ParamP1Matter_InitPhysicalTestcase(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1Matter & P1M){
  /** This function modify the parameter to have parameter consistant with the test case **/
 
  if (d.nTest == 1){
    P1M.Cv_value=0.001;
    P1M.TempBB=1000;
    P1M.a_value=0.00000000000000755555;
    cout<<"For this test case rhoCv = "<<P1M.Cv_value<<" TempBB = "<<P1M.TempBB<<" and a="<<P1M.a_value<<endl;
    cout<<"The test case is compute with this prepared values"<<endl;
  }
  if (d.nTest == 2){
    P1M.Cv_value=0.38214*0.14361;
    P1M.TempBB=1.;
    P1M.a_value=0.01372;
    cout<<"For this test case rhoCv = "<<P1M.Cv_value<<" TempBB = "<<P1M.TempBB<<" and a="<<P1M.a_value<<endl;
    cout<<"The test case is compute with this prepared values"<<endl;
  }
  
  
 }

double AverageSigNodeP1Matter(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamP1Matter & P1M,int numGr,int j){
  /** Function which compute an average of sigma at the node **/
  double sigma=0;
  int jG1;
  int Typemean=1;

  double m1=0;

  if(Typemean ==1 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
       jG1=tab.TabInv[numGr].TabCell[i];
       m1=m1+P1M.sigma[jG1]; 
     }
    
    sigma=m1/tab.TabInv[numGr].taille; 
  }
  
  m1=100000000000000;
  if(Typemean ==2 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
      jG1=tab.TabInv[numGr].TabCell[i];
      m1=min(m1,P1M.sigma[jG1]); 
    }
    sigma=m1; 
  }
  
  return sigma;
}

double AverageSigANodeP1Matter(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamP1Matter & P1M,int numGr,int j){
  /** Function which compute an average of sigma absortion at the node **/
  double sigma=0;
  int jG1;
  int Typemean=1;

  double m1=0;

  if(Typemean ==1 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
       jG1=tab.TabInv[numGr].TabCell[i];
       m1=m1+P1M.sigmaA[jG1]; 
     }
    
    sigma=m1/tab.TabInv[numGr].taille; 
  }
  
  m1=100000000000000;
  if(Typemean ==2 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
      jG1=tab.TabInv[numGr].TabCell[i];
      m1=min(m1,P1M.sigmaA[jG1]); 
    }
    sigma=m1; 
  }
  
  return sigma;
}

double AverageEpsNodeP1Matter(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamP1Matter & P1M,int numGr,int j){
  /** Function which compute an average of eps at the node **/
  double eps=0;
  int jG1;

  double m1=0;
  int Typemean=1;


  if(Typemean ==1 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
      jG1=tab.TabInv[numGr].TabCell[i];
      m1=m1+P1M.eps[jG1]; 
    }
      
    eps=m1/tab.TabInv[numGr].taille; 
  }


  if(Typemean ==2 ){
    m1=100000000000000;
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
      jG1=tab.TabInv[numGr].TabCell[i];
      m1=min(m1,P1M.eps[jG1]); 
    }
    eps=m1; 
  }
  
  return eps;
}
