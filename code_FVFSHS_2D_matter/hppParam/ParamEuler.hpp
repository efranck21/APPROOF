#ifndef PARAMEULER_HPP
#define PARAMEULER_HPP

#include <cmath>
#include <cstdlib>
#include <iostream> 
#include "R2.hpp"
#include "Data.hpp"
#include "ClassMesh.hpp"
#include "FunctionsConnect.hpp"
///#include "HOfunctions.hpp"
#include <cassert>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;

/** Class: physical parameter for Euler model 
    - nbvar: number of variable for the model
    - nbcell: number of cells
    - eps: table of epsilon coefficient for each cell
    - sigma: table of sigma coefficient for each cell
    - eps_value : epsilon value (-100000 if eps is variable)
    - sig_value : sigma value (-100000 if sig is variable)
    - gx, gy : gravity vector (-100000 if the gravity is variable)
    - gamma : gamma value
    - OrderAdv : order of remap scheme (MUSCL for the second order)
    - lp : number of pressure law (1 for perfect gas law)
    - LHO : order of hydrostatic equilibrium reconstruction  
**/

class LSProblem {
public:
   int NColum;
  int NLign;
  double * A;
  double * b_rho;
  double * b_pressure;
  double * Sol_rho;
  double * Sol_pressure;
  int *Stencil;

  
  LSProblem(){
    NColum=0;
    NLign=0;
    A=NULL;
    b_rho=NULL;
    b_pressure=NULL;
    Sol_rho=NULL;
    Sol_pressure=NULL;
    Stencil=NULL;
  }
  

  LSProblem(const LSProblem & ls){
    NColum=ls.NColum;
    NLign=ls.NLign;
    
    A = new double[ls.NLign*ls.NColum];
    for(int i=0;i<ls.NLign*ls.NColum;i++){
      A[i]= ls.A[i];
    }
    
    b_rho = new double[ls.NLign];
    for(int i=0;i<ls.NLign;i++){
      b_rho[i]= ls.b_rho[i];
    }
    
    Sol_rho = new double[ls.NColum];
    for(int i=0;i<ls.NColum;i++){
      Sol_rho[i]= ls.Sol_rho[i];
    }
    
     b_pressure = new double[ls.NLign];
    for(int i=0;i<ls.NLign;i++){
      b_pressure[i]= ls.b_pressure[i];
    }
    
    Sol_pressure = new double[ls.NColum];
    for(int i=0;i<ls.NColum;i++){
      Sol_pressure[i]= ls.Sol_pressure[i];
    }
    
    Stencil= new int[ls.NLign];
    for(int i=0;i<ls.NLign;i++){
      Stencil[i]= ls.Stencil[i];
    }
  }
  
  ~LSProblem(){
  
      if(A!=NULL){
	delete [] A;
      }     
       if(b_rho!=NULL){
	 delete [] b_rho;
      }  
        if(Sol_rho!=NULL){
	 delete [] Sol_rho;
      }
	if(b_pressure!=NULL){
	 delete [] b_pressure;
      }  
        if(Sol_pressure!=NULL){
	 delete [] Sol_pressure;
      }
	if(Stencil!=NULL){
	 delete [] Stencil;
      }  
  }
	
 LSProblem & operator=(const LSProblem & ls){
   NColum=ls.NColum;
   NLign=ls.NLign;
   /**
   if(A==NULL){
     A = new double[ls.NLign*ls.NColum];
   }

   for(int i=0;i<ls.NLign*ls.NColum;i++){
     A[i]= ls.A[i];
     }
   
   if(b_rho==NULL){
     b_rho = new double[ls.NLign];
   }
   for(int i=0;i<ls.NLign;i++){
     b_rho[i]= ls.b_rho[i];
   }
   
    if(b_pressure==NULL){
     b_pressure = new double[ls.NLign];
   }
   for(int i=0;i<ls.NLign;i++){
     b_pressure[i]= ls.b_pressure[i];
     }
   
   if(Sol_rho==NULL){
     Sol_rho = new double[ls.NColum];
   }
   for(int i=0;i<ls.NColum;i++){
     Sol_rho[i]= ls.Sol_rho[i];
   }
   
   if(Sol_pressure==NULL){
     Sol_pressure = new double[ls.NColum];
   }
   for(int i=0;i<ls.NColum;i++){
     Sol_pressure[i]= ls.Sol_pressure[i];
   }
   
   if(Stencil==NULL){
     Stencil= new int[ls.NLign];
   }
   for(int i=0;i<ls.NLign;i++){
     Stencil[i]= ls.Stencil[i];
   }
   **/
   
   if(A!=NULL){
     delete [] A;
    }
     
   if(ls.A !=NULL){
     A = new double[ls.NLign*ls.NColum];
     for(int i=0;i<ls.NLign*ls.NColum;i++){
       A[i]= ls.A[i];
     }
   }

   if(b_rho!=NULL){
     delete [] b_rho;
    }
     
   if(ls.b_rho !=NULL){
     b_rho = new double[ls.NLign];
     for(int i=0;i<ls.NLign;i++){
       b_rho[i]= ls.b_rho[i];
     }
   }

   if(b_pressure!=NULL){
     delete [] b_pressure;
    }
     
   if(ls.b_pressure !=NULL){
     b_pressure = new double[ls.NLign];
     for(int i=0;i<ls.NLign;i++){
       b_pressure[i]= ls.b_pressure[i];
     }
   }

   if(Sol_rho!=NULL){
     delete [] Sol_rho;
    }
     
   if(ls.Sol_rho !=NULL){
     Sol_rho = new double[ls.NColum];
     for(int i=0;i<ls.NColum;i++){
       Sol_rho[i]= ls.Sol_rho[i];
     }
   }

   if(Sol_pressure!=NULL){
     delete [] Sol_pressure;
    }
     
   if(ls.Sol_pressure !=NULL){
     Sol_pressure = new double[ls.NColum];
     for(int i=0;i<ls.NColum;i++){
       Sol_pressure[i]= ls.Sol_pressure[i];
     }
   }

   if(Stencil!=NULL){
     delete [] Stencil;
    }
     
   if(ls.Stencil !=NULL){
     Stencil = new int[ls.NLign];
     for(int i=0;i<ls.NLign;i++){
       Stencil[i]= ls.Stencil[i];
     }
   }
    
    return *this;}
   
}; 

class ParamEuler {
public:
  int nbvar;
  int nbcell;
  int nbnode;
  int nbnodeinternal;
  double *eps;
  double *sigma;
  double eps_value;
  double sig_value;
  double gx;
  double gy;
  int lp;
  double gamma;
  double * phi;
  int LHO;
  int OrderAdv;
  LSProblem * TabLS;
  R2 * TabRhoG;
  
  ParamEuler(){
    nbvar=4;
    nbcell=0;
    nbnode=0;
    nbnodeinternal=0;
    lp = 0 ;
    gamma =0;
    LHO=1;
    eps_value=0;
    sig_value=0;
    gx=0;
    gy=0;
    OrderAdv=0;
    eps=NULL;
    sigma=NULL;
    phi=NULL;
    TabLS=NULL;
    TabRhoG=NULL;
  }
  
  ParamEuler(Data & d,Mesh & Mh){ 
    nbvar=3;
    nbcell=Mh.nc;
    nbnode=Mh.nv;
    nbnodeinternal=(d.Nx+1)*(d.Ny+1);
    lp = 0 ;
    gamma =0;
    LHO=1;
    eps_value=0;
    sig_value=0;
    gx=0;
    gy=0;
    OrderAdv=1;
    eps=new double[Mh.nc];
    sigma=new double[Mh.nc];
    phi=new double[Mh.nc];
    TabLS=new LSProblem[nbnode];
    TabRhoG=new R2[Mh.nv];
  }

  ParamEuler(const ParamEuler & pe){
    nbvar=pe.nbvar;
    nbcell=pe.nbcell;
    nbnode=pe.nbnode;
    nbnodeinternal=pe.nbnodeinternal;
    lp=pe.lp;
    gamma=pe.gamma; 
    LHO=pe.LHO;
    eps_value=pe.eps_value;
    sig_value=pe.sig_value;
    gx=0;
    gy=0; 
    OrderAdv=pe.OrderAdv;
    eps = new double[nbcell];
    for(int i=0;i<nbcell;i++){
      eps[i]= pe.eps[i];
    }    
    sigma = new double[nbcell];
    for(int i=0;i<nbcell;i++){
      sigma[i]= pe.sigma[i];
    }  
    phi = new double[nbcell];
    for(int i=0;i<nbcell;i++){
      phi[i]= pe.phi[i];
    }
    TabLS =new LSProblem[nbnode];
    for(int i=0;i<nbnode;i++){
      TabLS[i]= pe.TabLS[i];
    }
    TabRhoG =new R2[nbnode];
    for(int i=0;i<nbnode;i++){
      TabRhoG[i]= pe.TabRhoG[i];
    }
  }
  
  ~ParamEuler(){
      if(eps!=NULL){
	delete [] eps;
      }     
      if(sigma!=NULL){
	delete [] sigma;
      }  
      if(phi!=NULL){
	delete [] phi;
      }
       if(TabLS!=NULL){
      	delete [] TabLS;
      }
      if(TabRhoG!=NULL){
	delete [] TabRhoG;
      }
  }
	
 ParamEuler & operator=(const ParamEuler & pe){
    nbvar=pe.nbvar;
    lp=pe.lp;
    nbcell=pe.nbcell;
    nbnode=pe.nbnode;
    nbnodeinternal=pe.nbnodeinternal;
    gamma=pe.gamma;
    LHO=pe.LHO;
    eps_value=pe.eps_value;
    sig_value=pe.sig_value;
    gx=pe.gx;
    gy=pe.gy;
    OrderAdv=pe.OrderAdv;
    
    /**if(eps==NULL){
      eps = new double[nbcell];
    }
    for(int i=0;i<nbcell;i++){
      eps[i]= pe.eps[i];
      }

    
    if(sigma==NULL){
      sigma = new double[nbcell];
    }
    for(int i=0;i<nbcell;i++){
      sigma[i]= pe.sigma[i];
    }
    
     if(phi==NULL){
      phi = new double[nbcell];
    }
    for(int i=0;i<nbcell;i++){
      phi[i]= pe.phi[i];
      }
    
     if(TabLS==NULL){
      TabLS = new LSProblem[nbnode];
    }
    for(int i=0;i<nbnode;i++){
      TabLS[i]= pe.TabLS[i];
    }
    
    if(TabRhoG==NULL){
      TabRhoG = new R2[nbnode];
    }
    for(int i=0;i<nbnode;i++){
      TabRhoG[i]= pe.TabRhoG[i];
    }
    **/

    if(eps!=NULL){
     delete [] eps;
    }
     
   if(pe.eps !=NULL){
     eps = new double[nbcell];
     for(int i=0;i<nbcell;i++){
       eps[i]= pe.eps[i];
     }
   }

   if(sigma!=NULL){
     delete [] sigma;
    }
     
   if(pe.sigma !=NULL){
     sigma = new double[nbcell];
     for(int i=0;i<nbcell;i++){
       sigma[i]= pe.sigma[i];
     }
   }

   if(phi!=NULL){
     delete [] phi;
    }
     
   if(pe.phi !=NULL){
     phi = new double[nbcell];
     for(int i=0;i<nbcell;i++){
       phi[i]= pe.phi[i];
     }
   }

    if(TabRhoG !=NULL){
     delete [] TabRhoG;
    }
     
   if(pe.TabRhoG !=NULL){
     TabRhoG = new R2[nbnode];
     for(int i=0;i<nbnode;i++){
       TabRhoG[i].x= pe.TabRhoG[i].x;
       TabRhoG[i].y= pe.TabRhoG[i].y;
     }
   }

    if(TabLS !=NULL){
     delete [] TabLS;
    }
     
   if(pe.TabLS !=NULL){
     TabLS= new LSProblem[nbnode];
     for(int i=0;i<nbnode;i++){
       TabLS[i]=pe.TabLS[i];
     }
   }
    
    return *this;
 }    
}; 


void ParamEuler_Init(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab, ParamEuler & Euler, const char * filename);

void ParamEuler_Initwithfile(ParamEuler & Euler, const char * filename);

void ParamEuler_InitTab(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamEuler & Euler);

void ParamEuler_InitPhysicalTestcase(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamEuler & Euler);

double AverageSigNodeEuler(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamEuler & Euler,int numGr,int j);

double AverageEpsNodeEuler(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamEuler & Euler,int numGr,int j);



#endif
