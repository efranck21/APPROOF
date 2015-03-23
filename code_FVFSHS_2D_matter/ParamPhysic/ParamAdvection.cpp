#include <cmath>
#include <cstdlib>
#include <iostream> 
#include "ParamAdvection.hpp"
#include <cassert>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;

/** File which contained some functions to intitialize the paramater of advection model **/ 

void ParamAdv_Init(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamAdvection & Adv, const char * filename){
  /** Initialization of advection parameter **/
  ParamAdv_Initwithfile(Adv,filename);
  ParamAdv_InitTab(d,Mh,v,tab,Adv);

}

void ParamAdv_Initwithfile(ParamAdvection & Adv, const char * filename){
  /** Use the file ParamAdv.dat to initialize the advection parameter **/
  
    FILE* f = fopen(filename,"rb");
    if(!f){cout << "Error when we open file"<<filename<<endl;exit(1);}
     cout <<"Reading of file\"" << filename << "\""<<endl;
    char string[255];

  /* velocity_y */
  fgets(string,255,f);
  fscanf(f,"%lf\n",&Adv.velocity.x);

  /* velocity_y */
  fgets(string,255,f);
  fscanf(f,"%lf\n",&Adv.velocity.y);

  /* order */
  fgets(string,255,f);
  fscanf(f,"%d\n",&Adv.OrderAdv);

  /*limitor */
  fgets(string,255,f);
  fscanf(f,"%d\n",&Adv.Nlim);
  
  cout << "End of reading ParamAdv.dat "<<endl;
 
  }

void ParamAdv_InitTab(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamAdvection & Adv){
/** Function which initialize the table of advection velocity using the given value (paramAdv.dat)
or the test case if the velocity is variable **/
  
  if (Adv.velocity.x != -100000 || Adv.velocity.y != -100000){
     for(int i=0;i<Mh.nc;i++){
      Adv.a[i].x=Adv.velocity.x;
      Adv.a[i].y=Adv.velocity.y;
    }
  }
 }




R2 AverageVelocityNode(Data & d,Mesh & Mh,variable  & v,TabConnecInv & tab,ParamAdvection & Adv,int numGr,int j){
  /** Function which compute an average of the velocity at the node **/
  R2 Average_velocity(0,0);
  int jG1;
  int Typemean =1;
  double m1=0,m2=0;

    if(Typemean ==1 ){
     for(int i=0;i<tab.TabInv[numGr].taille;i++){
       jG1=tab.TabInv[numGr].TabCell[i];
       m1=m1+Adv.a[jG1].x; 
       m2=m2+Adv.a[jG1].y; 
     }
      
     Average_velocity.x=m1/tab.TabInv[numGr].taille;
     Average_velocity.y=m2/tab.TabInv[numGr].taille;
     }

    if(Typemean ==2 ){
      m1=100000000000000;
      m2=100000000000000;
      for(int i=0;i<tab.TabInv[numGr].taille;i++){
	jG1=tab.TabInv[numGr].TabCell[i];
	m1=min(m1,Adv.a[jG1].x); 
	m2=min(m1,Adv.a[jG1].y); 
      }
    Average_velocity.x=m1;
    Average_velocity.y=m2; 
    }
  
  return Average_velocity;
}
