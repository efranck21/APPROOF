#include <cmath>
#include <cstdlib>
#include <iostream> 
#include "ParamP1.hpp"
#include <cassert>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;

/** File which contained some functions to intitialize the paramater of P1 model **/ 

void ParamP1_Init(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1 & P1, const char * filename){
/** Initialization of P1 parameter **/
  ParamP1_Initwithfile(P1,filename);
  ParamP1_InitPhysicalTestcase(d,Mh,v,tab,P1);
  ParamP1_InitTab(d,Mh,v,tab,P1);

}

void ParamP1_Initwithfile(ParamP1 & P1, const char * filename){
 /** Use the file ParamP1.dat to initialize the P1 parameter **/
  
    FILE* f = fopen(filename,"rb");
    if(!f){cout << "Error when we open file"<<filename<<endl;exit(1);}
     cout <<"Reading of file\"" << filename << "\""<<endl;
    char string[255];

  /* sigma*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&P1.sig_value);

   /* eps*/
  fgets(string,255,f);
  fscanf(f,"%lf\n",&P1.eps_value);

  
  
  cout << "End of reading ParamP1.dat  "<<endl;
 
  }

void ParamP1_InitTab(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1 & P1){
/** Function which initialize the table of P1 quantities using the given values (paramP1.dat)
or the test case if the eps, sigma or phi are variables **/
  
  if (P1.eps_value != -1){
     for(int i=0;i<P1.nbcell;i++){
      P1.eps[i]=P1.eps_value;
    }
  }
   if (P1.sig_value != -1){
     for(int i=0;i<P1.nbcell;i++){
      P1.sigma[i]=P1.sig_value;
    }
  }
   if (P1.sig_value == -1 && d.nTest == 5){
     for(int i=0;i<P1.nbcell;i++){
       if(Mh.cells[i].lab==0){
	 P1.sigma[i]=1;
       }
       if(Mh.cells[i].lab > 0){
	 P1.sigma[i]=10000;
       }
     }
   }
   if (P1.sig_value == -1 && d.nTest == 6){
     for(int i=0;i<P1.nbcell;i++){
       if(Mh.cells[i].lab!=1 ){
	 P1.sigma[i]=3;
       }
       if(Mh.cells[i].lab ==1){
	 P1.sigma[i]=0;
       }
     }
   }
   
   
 }

void ParamP1_InitPhysicalTestcase(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamP1 & P1){
 /** This function modify the parameter to have parameter consistant with the test case **/
  
  if (d.nTest == 1){
    P1.sig_value=0.;
    P1.eps_value=1.;
    cout<<"For this test case sig = "<<P1.sig_value<<" and eps = "<<P1.eps_value<<endl;
    cout<<"The test case is compute with this prepared values"<<endl;
  }
  if(d.nTest == 1 && d.Tx != d.Ty){
      cout<<"The test case is defined on a square with periodic condition (2*2 or 4*4 etc)"<<endl;
      d.Ty=2;
      d.Tx=2;
      cout<<"We run the test with = "<<d.Tx<<" and "<<d.Ty<<endl;
    }
  
  if (d.nTest == 2){
    P1.sig_value=1.;
    P1.eps_value=1.; 
    cout<<"For this test case sig = "<<P1.sig_value<<" and eps = "<<P1.eps_value<<endl;
    cout<<"The test case is compute with this prepared values"<<endl;    
  }
  if(d.nTest == 2 && d.Tx != d.Ty){
      cout<<"The test case is defined on a square with periodic condition (2*2 or 4*4 etc)"<<endl;
      d.Ty=2;
      d.Tx=2;
      cout<<"We run the test with = "<<d.Tx<<" and "<<d.Ty<<endl;
    }

  
  if (d.nTest == 3 && P1.sig_value ==0 ){
    P1.sig_value=1.;    
    cout<<"For this test case sigma >0 we take sig = "<<P1.sig_value<<endl;
    cout<<"The test case is compute with this prepared values"<<endl; 
    
  }

  if (d.nTest == 4 && P1.sig_value ==0 ){
    P1.sig_value=1.;    
    cout<<"For this test case sigma >0 we take sig = "<<P1.sig_value<<endl;
    cout<<"The test case is compute with this prepared values"<<endl; 
    
  }
  if(d.nTest == 4 && d.Tx != d.Ty){
    cout<<"The test case is defined on a square with periodic condition (2*2 or 4*4 etc)"<<endl;
    d.Ty=2;
    d.Tx=2;
    cout<<"We run the test with = "<<d.Tx<<" and "<<d.Ty<<endl;
  }
  if(d.nTest == 5){
    cout<<"The test case is defined for the DoubleSquareMeshQ 1*1, eps=1 and sigma variable "<<endl;
    d.Ty=1;
    d.Tx=1;
    P1.eps_value=1;
    P1.sig_value=-1;
    strcpy(d.NameMesh,"DoubleSquareMeshQ");
    cout<<"We run the test with = "<<d.Tx<<" and "<<d.Ty<<endl;
  }
  if(d.nTest == 6){
    cout<<"The test case is defined for the Check Mesh 7*7, eps=sqrt(3.) and sigma variable "<<endl;
    d.Ty=7;
    d.Tx=7;
    P1.eps_value=sqrt(3.);
    P1.sig_value=-1;
    strcpy(d.NameMesh,"CheckMesh");
    cout<<"We run the test with = "<<d.Tx<<" and "<<d.Ty<<endl;
  }
  
}

double AverageSigNodeP1(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamP1 & P1,int numGr,int j){
   /** Function which compute an average of sigma at the node **/
  double sigma=0;
  int jG1;
  int Typemean=1;

  double m1=0;

  if(Typemean ==1 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
       jG1=tab.TabInv[numGr].TabCell[i];
       m1=m1+P1.sigma[jG1]; 
     }
    
    sigma=m1/tab.TabInv[numGr].taille; //moyenne arithmetique
  }
  
  m1=100000000000000;
  if(Typemean ==2 ){
    for(int i=0;i<tab.TabInv[numGr].taille;i++){
      jG1=tab.TabInv[numGr].TabCell[i];
      m1=min(m1,P1.sigma[jG1]); 
    }
    sigma=m1; // moyenne harmonique
  }
  
  return sigma;
}

double AverageEpsNodeP1(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamP1 & P1,int numGr,int j){
   /** Function which compute an average of eps at the node **/
  double eps=0;
  int jG1;

  double m1=0;
  for(int i=0;i<tab.TabInv[numGr].taille;i++){
    jG1=tab.TabInv[numGr].TabCell[i];
    m1=m1+P1.eps[jG1]; 
  }
      
  eps=m1/tab.TabInv[numGr].taille; //moyenne arithmetique

  
  m1=100000000000000;
  for(int i=0;i<tab.TabInv[numGr].taille;i++){
    jG1=tab.TabInv[numGr].TabCell[i];
    m1=min(m1,P1.eps[jG1]); 
  }
  eps=m1; // moyenne harmonique
  
  
  return eps;
}
