#include <cmath>
#include <cstdlib>
#include <iostream> 
#include "ParamDiff.hpp"
#include <cassert>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;

/** File which contained some functions to intitialize the paramater of diffusion model **/ 

void ParamDiff_Init(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab, ParamDiff & Diff, const char * filename){
/** Initialization of diffusion parameter **/
  ParamDiff_Initwithfile(Diff,filename);
  ParamDiff_InitTab(d,Mh,v,tab,Diff);
  ParamDiff_InitPhysicalTestcase(d,Mh,v,tab,Diff);

}

void ParamDiff_Initwithfile(ParamDiff & Diff, const char * filename){
  /** Use the file ParamDiff.dat to initialize the diffusion parameter **/
  
    FILE* f = fopen(filename,"rb");
    if(!f){cout << "Error when we open file"<<filename<<endl;exit(1);}
     cout <<"Reading of file\"" << filename << "\""<<endl;
    char string[255];

  /* coef diff */
  fgets(string,255,f);
  fscanf(f,"%d\n",&Diff.Coefdiff_value);
  
  cout << "End of reading ParamDiff.dat "<<endl;
 
  }

void ParamDiff_InitTab(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamDiff & Diff){
/** Function which initialize the table of diffusion coefficient using the given value (paramDiff.dat)
or the test case if the diffusion coefficient is variable **/
  
  if (Diff.Coefdiff_value != -1){
     for(int i=0;i<Mh.nc;i++){
      Diff.D[i]=Diff.Coefdiff_value;
    }
  }
 }


void ParamDiff_InitPhysicalTestcase(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamDiff & Diff){
  // This function modify the parameter to have parameter consistant with the test case
 
  if (d.nTest == 1 && Diff.Coefdiff_value==0.){
    Diff.Coefdiff_value=1.; 
    cout<<"For this test the diffusion coeffcient is positive. We take = "<<Diff.Coefdiff_value<<endl;
    cout<<"The test case is compute with this prepared values"<<endl; 
    
  }
  if (d.nTest == 2 && Diff.Coefdiff_value==0.){
    Diff.Coefdiff_value=1.; 
    cout<<"For this test the diffusion coeffcient is positive. We take = "<<Diff.Coefdiff_value<<endl;
    cout<<"The test case is compute with this prepared values"<<endl; 
    
  }
  if (d.nTest == 3 && Diff.Coefdiff_value==0.){
    Diff.Coefdiff_value=1.; 
    cout<<"For this test the diffusion coeffcient is positive. We take = "<<Diff.Coefdiff_value<<endl;
    cout<<"The test case is compute with this prepared values"<<endl;
  }
  if(d.nTest == 3 && d.Tx != d.Ty){
    cout<<"The test case is defined on a square with periodic condition (2*2 or 4*4 etc)"<<endl;
    d.Ty=2;
    d.Tx=2;
    cout<<"We run the test with = "<<d.Tx<<" and "<<d.Ty<<endl;
  } 
  
  
 }


double AverageCoefDiffNode(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamDiff & Diff,int numGr,int j){
  /** Function which compute an average of the diffusion coefficient at the node **/
  double Average_Diff=0;
  int jG1;
  int Typemean =1;
  double m1=0;

    if(Typemean ==1 ){
     for(int i=0;i<tab.TabInv[numGr].taille;i++){
       jG1=tab.TabInv[numGr].TabCell[i];
       m1=m1+Diff.D[jG1]; 
     }
      
     Average_Diff=m1/tab.TabInv[numGr].taille;
     }

    if(Typemean ==2 ){
      m1=100000000000000;
      for(int i=0;i<tab.TabInv[numGr].taille;i++){
	jG1=tab.TabInv[numGr].TabCell[i];
	m1=min(m1,Diff.D[jG1]); 
      }
    Average_Diff=m1; 
    }
  
  return Average_Diff;
}
