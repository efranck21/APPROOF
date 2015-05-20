#ifndef DATA_HPP
#define DATA_HPP

#include <cmath>
#include <cstdlib>
#include <iostream> 
#include "R2.hpp"
#include <cassert>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;

/**
 @brief Class: this class contains the parameter of the code.
**/

class Data {
  public :
  /**Number of cells in the x direction**/
  int Nx;
  /**Number of cells in the y direction **/
  int Ny;
  /** Lenght of the domain in the x direction **/
  float Tx;
  /** Lenght of the domain in the y direction **/
  float Ty;
  
  /** Scheme number **/
  int scheme;
   /** Type of scheme (N for nodal scheme, E for edge scheme) **/
  char Typescheme;
  /** Test case number**/
  int nTest; 
  /** Final time**/
  float Tf;
  /** CFL **/
  float CFL;
  /** total order **/
  int Total_order;
  
  /** Mesh name **/
  char *NameMesh;  
  /** the number of variable that we plot **/
  char *NumVar;
  /** Suffix for the output file **/
  char *suffixe;
  
  /** Model :
   - Diffusion : diffusion equation
   - P1 : P1 model
   - M1 : M1 odel
   - Advection : advection
   - Euler : Euler
   - P1Matter : P1 model with matter**/
  char *Typemodel;
  
  /** Time discretization
      - E : explicit
      - S : Semi-implicit.**/
  char Typetime;


  /** Mesh type
   - T : triangular
   - Q : quandrangular **/
  char Typemesh;

  
  /** Time step for the plot **/
  float dtAnim;
  /** Plot solutions: yes (y) or not (n)**/
  char Anim;
  
  /** plot using the cell (N) or the center of the cell (C)**/
  char Typewrite;
  /** plot in 1D or 2D **/
  int dimsave;
  int restart;

  /** nbgroup **/
  int ngroup;

  Data(const char * filename){
    FILE* f = fopen(filename,"rb");
    if(!f){cout << "Erreur d'ouverture du fichier"<<filename<<endl;exit(1);}
     cout <<"Lecture du fichier\"" << filename << "\""<<endl;
     char string[255];
     Typemodel=(char*)malloc(100*sizeof(char));
     NameMesh=(char*)malloc(20*sizeof(char));
     NumVar=(char*)malloc(300*sizeof(char));
     suffixe=(char*)malloc(10*sizeof(char));
     

  /* Nx */
  fgets(string,255,f);
  fscanf(f,"%d\n",&Nx);
  
  /* Ny */
  fgets(string,255,f);
  fscanf(f,"%d\n",&Ny);
  
  /* Tx */
  fgets(string,255,f);
  fscanf(f,"%f\n",&Tx);

  /* Ty */
  fgets(string,255,f);
  fscanf(f,"%f\n",&Ty);
  
  /* Test case number */
  fgets(string,255,f);
  fscanf(f,"%d\n",&nTest);

/* Scheme number*/
  fgets(string,255,f);
  fscanf(f,"%d\n",&scheme);

 /* Geometrical discretization scheme*/
  fgets(string,255,f);
  fscanf(f,"%c\n",&Typescheme);
  
 /* CFL */
  fgets(string,255,f);
  fscanf(f,"%f\n",&CFL);
  
  /* Final time */
  fgets(string,255,f);
  fscanf(f,"%f\n",&Tf);

  /* Total order */
  fgets(string,255,f);
  fscanf(f,"%d\n",&Total_order);
  
  /* Mesh name */
  fgets(string,255,f);
  fscanf(f,"%s\n",NameMesh);
 
  /* Model */
  fgets(string,255,f);
  fscanf(f,"%s\n",Typemodel);

     /* nb group */
  fgets(string,255,f);
  fscanf(f,"%d\n",&ngroup);
  
 /* Time scheme */
  fgets(string,255,f);
  fscanf(f,"%c\n",&Typetime);

  /* Type mesh */
  fgets(string,255,f);
  fscanf(f,"%c\n",&Typemesh);

  /* variable that we plot */
  fgets(string,255,f);
  fscanf(f,"%s\n",NumVar);	

  /* Time step for plot */
  fgets(string,255,f);
  fscanf(f,"%f\n",&dtAnim);
  
  /* Plot yes or no */
  fgets(string,255,f);
  fscanf(f,"%c\n",&Anim);
  
   /* Suffix for plot */
  fgets(string,255,f);
  fscanf(f,"%s\n",suffixe);

  /* Type writting */
  fgets(string,255,f);
  fscanf(f,"%c\n",&Typewrite);

  /* Dimension for plot */
  fgets(string,255,f);
  fscanf(f,"%d\n",&dimsave);


  /* restart (1 for yes and 0 for no) */
  fgets(string,255,f);
  fscanf(f,"%d\n",&restart);
  cout << "end of constrction for Data "<<endl;

  fclose(f); 
  
  }
  Data(const Data & d){
    Nx=d.Nx;
    Ny=d.Ny;
    Tx=d.Tx;
    Ty=d.Ty;
    nTest=d.nTest;
    scheme=d.scheme;
    Tf=d.Tf; 
    CFL=d.CFL;
    NameMesh=strdup(d.NameMesh);
    NumVar=strdup(d.NumVar);
    suffixe=strdup(d.suffixe);
    Typemodel=strdup(d.Typemodel);
    Typetime=d.Typetime;
    Typemesh=d.Typemesh;
    dtAnim=d.dtAnim;
    Anim=d.Anim;
    Typewrite=d.Typewrite;
    dimsave=d.dimsave;
    Typescheme=d.Typescheme;
    restart=d.restart;
    Total_order=d.Total_order;
    ngroup=d.ngroup;

	
  }
  Data & operator =(const Data & d){
     Nx=d.Nx;
    Ny=d.Ny;
    Tx=d.Tx;
    Ty=d.Ty;
    nTest=d.nTest;
    scheme=d.scheme;
    Tf=d.Tf; 
    CFL=d.CFL;
    NameMesh=strdup(d.NameMesh);
    NumVar=strdup(d.NumVar);
    suffixe=strdup(d.suffixe);
    Typemodel=strdup(d.Typemodel);
    Typetime=d.Typetime;
    Typemesh=d.Typemesh;
    dtAnim=d.dtAnim;
    Anim=d.Anim;
    Typewrite=d.Typewrite;
    dimsave=d.dimsave;
    Typescheme=d.Typescheme;
    restart=d.restart;
    Total_order=d.Total_order;
    ngroup=d.ngroup;
    

    return *this;
  }
  ~Data(){
    free(NameMesh);
    free(NumVar);
    free(suffixe);
    free(Typemodel);
  }
  
};
    
#endif
