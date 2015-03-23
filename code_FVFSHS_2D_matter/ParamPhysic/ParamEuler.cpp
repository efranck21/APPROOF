#include <cmath>
#include <cstdlib>
#include <iostream> 
#include "ParamEuler.hpp"
#include <cassert>
#include <fstream>
#include <string.h>
#include <stdlib.h>

using namespace std;

/** File which contained some functions to intitialize the paramater of Euler model **/ 

void ParamEuler_Init(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab, ParamEuler & Euler, const char * filename){
/** Initialization of Euler parameter **/
   ParamEuler_Initwithfile(Euler,filename);
   ParamEuler_InitPhysicalTestcase(d,Mh,v,tab,Euler);
   ParamEuler_InitTab(d,Mh,v,tab,Euler);
}

void ParamEuler_Initwithfile(ParamEuler & Euler, const char * filename){
 /** Use the file ParamEuler.dat to initialize the Euler parameter **/
  
    FILE* f = fopen(filename,"rb");
    if(!f){cout << "Error when we open file"<<filename<<endl;exit(1);}
     cout <<"Reading of file\"" << filename << "\""<<endl;
    char string[255];

  /* Pressure law*/
  fgets(string,255,f);
  fscanf(f,"%d\n",&Euler.lp);

   /* gamma*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&Euler.gamma);

   /* Local high order*/
  fgets(string,255,f);
  fscanf(f,"%d\n",&Euler.LHO);

  /* sigma*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&Euler.sig_value);

   /* eps*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&Euler.eps_value);

   /* gx */
  fgets(string,255,f);
  fscanf(f,"%lf\n",&Euler.gx);

  /* gy */
  fgets(string,255,f);
  fscanf(f,"%lf\n",&Euler.gy);

  /* Order advection */
  fgets(string,255,f);
  fscanf(f,"%d\n",&Euler.OrderAdv); 
  
  cout << "End of reading ParamEuler.dat "<<endl;

  fclose(f); 
  
  }

void ParamEuler_InitTab(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamEuler & Euler){
/** Function which initialize the table of Euler quantities using the given value (paramEuler.dat)
or the test case if the eps, sigma or phi are variables **/
 
  if (Euler.eps_value != -1){
     for(int i=0;i<Mh.nc;i++){
      Euler.eps[i]=Euler.eps_value;
    }
  }
   if (Euler.sig_value != -1){
     for(int i=0;i<Mh.nc;i++){
      Euler.sigma[i]=Euler.sig_value;
    }
  }

   for(int i=0;i<Mh.nc;i++){
       Euler.phi[i]=0;
    }
   
  if (Euler.gx != -100000 || Euler.gy != -100000){
     for(int i=0;i<Mh.nc;i++){
       Euler.phi[i]=1;
    }
  }


 
  
}  

void ParamEuler_InitPhysicalTestcase(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamEuler & Euler){
  /** This function modify the parameter to have parameter consistant with the test case **/
  
  if (d.nTest < 7){
    Euler.gamma=1.4;
    Euler.sig_value=0.;
    Euler.lp =1;
    
   
    cout<<"For this test case sig = "<<Euler.sig_value<<" lp = "<<Euler.lp<<" and gamma="<<Euler.gamma<<endl;
    cout<<"The test case is compute with this prepared values"<<endl;
  }
   if (d.nTest == 7){
    Euler.gamma=1.4;
    Euler.eps_value=1.;
    Euler.lp =1;
   
    cout<<"For this test case eps = "<<Euler.eps_value<<" lp = "<<Euler.lp<<" and gamma="<<Euler.gamma<<endl;
    cout<<"The test case is compute with this prepared values"<<endl;
  }
   if (d.nTest == 8){
    Euler.gamma=1.4;
    Euler.lp =1;
    cout<<"For this test case lp = "<<Euler.lp<<" and gamma="<<Euler.gamma<<endl;
    cout<<"The test case is compute with this prepared values"<<endl;
  }
     
  if (d.nTest == 9){
    Euler.gamma=1.4;
    Euler.eps_value =1;
    cout<<"For this test case lp = "<<Euler.lp<<", gamma="<<Euler.gamma<<" and eps"<<Euler.eps_value<<endl;
    cout<<"The test case is compute with this prepared values"<<endl;
  }
}



double AverageSigNodeEuler(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamEuler & Euler,int numGr,int j){
   /** Function which compute an average of sigma at the node **/
  double sigma=0;
  int jG1; 
  int Typemean =1;

  double m1=0;
  if(Typemean ==1 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
       jG1=tab.TabInv[numGr].TabCell[i];
       m1=m1+Euler.sigma[jG1]; 
     }
    
    
    sigma=m1/tab.TabInv[numGr].taille; 
  }
  if(Typemean ==2 ){  
     m1=100000000000000;
     for(int i=0;i<tab.TabInv[numGr].taille;i++){
       jG1=tab.TabInv[numGr].TabCell[i];
       m1=min(m1,Euler.sigma[jG1]); 
     }
     sigma=m1; 
  }
  
  return sigma;
}

double AverageEpsNodeEuler(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamEuler & Euler,int numGr,int j){
  /** Function which compute an average of eps at the node **/
  double eps=0;
  int jG1;
  int Typemean =1;

  double m1=0;
  if(Typemean ==1 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
       jG1=tab.TabInv[numGr].TabCell[i];
       m1=m1+Euler.eps[jG1]; 
    }
  
  eps=m1/tab.TabInv[numGr].taille;
   }
  
  if(Typemean ==2 ){   
    m1=100000000000000;
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
      jG1=tab.TabInv[numGr].TabCell[i];
      m1=min(m1,Euler.eps[jG1]); 
    }
    eps=m1; 
  }
     
  return eps;
}



