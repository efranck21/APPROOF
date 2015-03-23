#include <cmath>
#include <cstdlib>
#include <iostream> 
#include "ParamM1Matter.hpp"
#include <cassert>
#include <fstream>
#include <string.h>
#include <stdlib.h>

using namespace std;

/** File which contained some functions to intitialize the paramater of M1 model **/ 

void ParamM1M_Init(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamM1Matter & M1M, const char * filename){
/** Initialization of M1M parameter **/
  ParamM1M_Initwithfile(M1M,filename);
  ParamM1M_InitTab(d,Mh,v,tab,M1M);
  ParamM1M_InitPhysicalTestcase(d,Mh,v,tab,M1M);

}

void ParamM1M_Initwithfile(ParamM1Matter & M1M, const char * filename){
 /** Use the file ParamM1M.dat to initialize the M1 with matter parameters **/
  
    FILE* f = fopen(filename,"rb");
    if(!f){cout << "Error when we open file"<<filename<<endl;exit(1);}
     cout <<"Reading of file\"" << filename << "\""<<endl;
    char string[255];

  /* sigma*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&M1M.sig_value);

  /* sigma_A*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&M1M.sigA_value);

  /* eps*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&M1M.eps_value);

  /* order advection*/
  fgets(string,255,f);
  fscanf(f,"%d\n",&M1M.OrderAdv);
  
  /* valeur de 1/mc2, propre au Compton */
  fgets(string,255,f);
  fscanf(f,"%lf\n",&M1M.mc2value);
    
  /* Temp_BB (Black Body temperature)*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&M1M.TempBB);
    
  /* a (stefan constant)*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&M1M.a_value);


  cout << "End of reading ParamM1.dat "<<endl;
 
  }

void ParamM1M_InitTab(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamM1Matter & M1M){
/** Function which initialize the table of M1 quantities using the given value (paramM1M.dat)
or the test case if the eps, sigma are variables **/
  
  if (M1M.eps_value != -1){
     for(int i=0;i<M1M.nbcell;i++){
      M1M.eps[i]=M1M.eps_value;
    }
  }
   if (M1M.sig_value != -1){
     for(int i=0;i<M1M.nbcell;i++){
      M1M.sigma[i]=M1M.sig_value;
      M1M.sigmaA[i]=M1M.sigA_value;
    }
  }
 }

void ParamM1M_InitPhysicalTestcase(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamM1Matter & M1M){
 /** This function modify the parameter to have parameter consistant with the test case **/
 
  if (d.nTest == 1){
    M1M.sig_value=1.;
    M1M.sigA_value=1.;
    M1M.eps_value=1.;
    cout<<"For this test case sig = "<<M1M.sig_value<<" and eps = "<<M1M.eps_value<<endl;
    cout<<"The test case is compute with this prepared values"<<endl; 
    
  }
  if (d.nTest == 2){
    M1M.sig_value=0.;
    M1M.eps_value=1.;
    cout<<"For this test case sig = "<<M1M.sig_value<<" and eps = "<<M1M.eps_value<<endl;
    cout<<"The test case is compute with this prepared values"<<endl; 
    
  }
  if (d.nTest == 3 && M1M.sig_value ==0 ){
    M1M.sig_value=1.;
    cout<<"For this test case sigma >0, we take sig = "<<M1M.sig_value<<endl;
    cout<<"The test case is compute with this prepared values"<<endl; 
  
  }
  if (d.nTest == 4 && M1M.sig_value ==0 ){
    M1M.sig_value=1.;
    cout<<"For this test case sigma >0, we take sig = "<<M1M.sig_value<<endl;
    cout<<"The test case is compute with this prepared values"<<endl;
  }
  if(d.nTest == 4 && d.Tx != d.Ty){
    cout<<"The test case is defined on a square with periodic condition (2*2 or 4*4 etc)"<<endl;
    d.Ty=2;
    d.Tx=2;
    cout<<"We run the test with = "<<d.Tx<<" and "<<d.Ty<<endl;
  } 
  
  
 }

double AverageSigNodeM1M(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamM1Matter & M1M,int numGr,int j){
     /** Function which compute an average of sigma at the node **/

  double sigma=0;
  int jG1;
  int Typemean=1;

  double m1=0;

  if(Typemean ==1 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
       jG1=tab.TabInv[numGr].TabCell[i];
       m1=m1+M1M.sigma[jG1];
     }
    
    sigma=m1/tab.TabInv[numGr].taille; //moyenne arithmetique
  }
  
  m1 = 100000000000000 ;
  if(Typemean ==2 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
      jG1=tab.TabInv[numGr].TabCell[i];
      m1=min(m1,M1M.sigma[jG1]);
    }
    sigma=m1; // moyenne harmonique
  }
  
  return sigma;
}

double AverageEpsNodeM1M(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamM1Matter & M1M,int numGr,int j){
     /** Function which compute an average of eps at the node **/

  double eps=0;
  int jG1;
  int Typemean=1;
  
  double m1=0;
  if(Typemean ==1 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
      jG1=tab.TabInv[numGr].TabCell[i];
      m1=m1+M1M.eps[jG1];
    }
    
    eps=m1/tab.TabInv[numGr].taille; 
  }
  if(Typemean ==2 ){
    m1=100000000000000;
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
      jG1=tab.TabInv[numGr].TabCell[i];
      m1=min(m1,M1M.eps[jG1]);
    }
    eps=m1; 
   }
  
  return eps;
}

double AverageEnergyNodeM1M(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamM1Matter & M1M,int numGr,int j){
     /** Function which compute an average of Energy at the node **/

  double Er=0;
  int jG1;
  double m1=0;
  int Typemean=1;
  
  if(Typemean ==1 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
      jG1=tab.TabInv[numGr].TabCell[i];
      m1=m1+v.var[0][jG1]; 
    }
    
    Er=m1/tab.TabInv[numGr].taille; 
  }
  if(Typemean ==2 ){
    m1=100000000000000;
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
      jG1=tab.TabInv[numGr].TabCell[i];
      m1=min(m1,v.var[0][jG1]); 
    }
    Er=m1;
  }
  
  return Er;
}

R2 AverageFluxNodeM1M(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamM1Matter & M1M,int numGr,int j){
     /** Function which compute an average of Flux at the node **/

  R2 Fr(0,0);
  int jG1;
  double m1=0,m2=0;
  int Typemean=1;
  
  if(Typemean ==1 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
      jG1=tab.TabInv[numGr].TabCell[i];
      m1=m1+v.var[1][jG1];
      m2=m2+v.var[2][jG1];  
    }
    
    Fr.x=m1/tab.TabInv[numGr].taille; 
    Fr.y=m2/tab.TabInv[numGr].taille;
   }
  
  if(Typemean ==2 ){
    m1=100000000000000;
    m2=100000000000000;
    
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
      jG1=tab.TabInv[numGr].TabCell[i];
      m1=min(m1,v.var[1][jG1]);  
      m2=min(m2,v.var[2][jG1]); 
    }
    Fr.x=m1; 
    Fr.y=m2;
  }
  
  return Fr;
}
