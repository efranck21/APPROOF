#ifndef VARIABLES_HPP
#define VARIABLES_HPP
#include <cassert>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Data.hpp"

using namespace std;



class variable {
public:
  int nbvar;
  int nbcell;
  double **var; 

  variable(){
    nbvar=0;
    nbcell=0;
    var=NULL;
  }
  
  variable(Data d,int nc,int nv){
    nbvar=0;
    if(!strcmp(d.Typemodel,"Diffusion")){nbvar=1;}
    if(!strcmp(d.Typemodel,"Advection")){nbvar=1;}
    if(!strcmp(d.Typemodel,"P1")) {nbvar=3;}
    if(!strcmp(d.Typemodel,"Euler")) {nbvar=4;}
    if(!strcmp(d.Typemodel,"M1")) {nbvar=3;}
    if(!strcmp(d.Typemodel,"P1Matter")) {nbvar=4;}
    if(!strcmp(d.Typemodel,"M1Matter")) {nbvar=4;}
    if(!strcmp(d.Typemodel,"P1Compton")) {nbvar=d.ngroup*3+1;}

    nbcell=nc;
     
    var = new double*[nbvar];
    for(int i=0;i<nbvar;i++){
      var[i]= new double[nbcell];
    }
  
  }

  variable(const variable & v){
    nbvar=v.nbvar;
    nbcell=v.nbcell;
  
    var = new double*[nbvar];
    for(int i=0;i<nbvar;i++){
      var[i]= new double[nbcell];
    }  
    
    for(int j=0;j<nbcell;j++){
      for(int i=0;i<nbvar;i++){
	var[i][j]=v.var[i][j];
      }
    }
  }
  
  ~variable(){
    if(var!=NULL){
      for(int i=0;i<nbvar;i++){
	if (var[i]!=NULL){
	  delete [] var[i];}
	
      }
    delete [] var;
    }
  }
    
	
 variable & operator=(variable v){
    nbvar=v.nbvar;
    nbcell=v.nbcell;
    
    nbvar=v.nbvar;
    nbcell=v.nbcell;
    if(var==NULL){
    var = new double*[nbvar];
    }
   
    for(int i=0;i<nbvar;i++){
      if(var[i]==NULL){
      var[i]= new double[nbcell];
      }
    }
    
    for(int j=0;j<nbcell;j++){
      for(int i=0;i<nbvar;i++){
	var[i][j]=v.var[i][j];
      }
    }
    
    return *this;
 }    
}; 
    

class vectorflux {
public:
  int nbvar;
  double *vflux;
 

  vectorflux(){
    nbvar=0;   
    vflux=NULL;
  }
  
  vectorflux(Data d){
    nbvar=0;
    if(!strcmp(d.Typemodel,"Diffusion")){nbvar=1;}
    if(!strcmp(d.Typemodel,"Advection")){nbvar=1;}
    if(!strcmp(d.Typemodel,"P1")) {nbvar=3;}
    if(!strcmp(d.Typemodel,"Euler")) {nbvar=4;}
    if(!strcmp(d.Typemodel,"M1")) {nbvar=3;}
    if(!strcmp(d.Typemodel,"P1Matter")) {nbvar=4;}
    if(!strcmp(d.Typemodel,"M1Matter")) {nbvar=4;}
    if(!strcmp(d.Typemodel,"P1Compton")) {nbvar=d.ngroup*3+1;}
    vflux = new double[nbvar];
    
  }

  vectorflux(int nbv){
    nbvar=nbv;
    vflux = new double[nbvar];
  }

  vectorflux(const vectorflux & v){
    nbvar=v.nbvar;
   
     vflux = new double[nbvar];
    
     for(int i=0;i<nbvar;i++){
     
       vflux[i]=v.vflux[i];
     }
  }
  
  ~vectorflux(){
    if(vflux!=NULL){
      delete [] vflux;
    }
  }
	
 vectorflux & operator=( vectorflux v){
    nbvar=v.nbvar;
    if(vflux==NULL){
      vflux = new double[nbvar];
    }
    
    for(int i=0;i<nbvar;i++){
     
      vflux[i]=v.vflux[i];
    }  
    return *this;
 }    

}; 
#endif
