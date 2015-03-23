

#include <iostream>
#include <cmath>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include "Functions.hpp"
#include "Data.hpp"
#include "Functionsgeo.hpp"
#include "Tensor.hpp"

/** file which contains some functions uselful **/



void WriteDataToMain(Data & d,ParamPhysic & Param){
/** function which plot the parameters at the beginning of the run **/
  
  if(!strcmp(d.Typemodel,"Diffusion")){
    cout<<"//////////Solve diffusion equation//////////"<<endl;
    if(d.scheme==1 && d.Typescheme =='N'){
      cout<<"//////////////////// Linear nodal scheme ////////////////////"<<endl;
    }
    if(d.scheme==1 && d.Typescheme =='E'){
      cout<<"////////////////// Linear VF4 Edge scheme ////////////////////"<<endl;
    }
    if(d.scheme==2 && d.Typescheme =='N'){
      cout<<"///////// Nonlinear nodal diffusion scheme ///////////"<<endl;}
  }
  
  if(!strcmp(d.Typemodel,"P1")){
      cout<<"////////////// Solve P1 model ///////////////"<<endl;
      if(d.scheme==1){
      cout<<"////////////////// Nodal JL-(a) scheme /////////////////"<<endl;}
      if(d.scheme==2){
	cout<<"//////////////// Nodal JL-(b) scheme /////////////////"<<endl;}
      if(d.scheme==3){
	cout<<"///////// Nodal JL-(b) with local source term //////////"<<endl;}    
        
      cout<<"//////////////// Eps : "<<Param.P1.eps_value<<"///////////////"<<endl;
      cout<<"////////////// sigma : "<<Param.P1.sig_value<<"///////////////"<<endl;
  }
    
  if(!strcmp(d.Typemodel,"Euler")){
    cout<<"////////////// Solve Euler model ///////////////"<<endl;
    if(d.scheme==1){
      cout<<"////////////// Nodal Lagrange+remap JL-(a) scheme ///////////////"<<endl;}
    if(d.scheme==2){
      cout<<"////////////// Nodal Lagrange+remap JL-(b) scheme /////////////////"<<endl;}
    if(d.scheme==3){
      cout<<"///// Nodal Lagrange+remap JL-(b) scheme with local source term ////////"<<endl;}    
    
    cout<<"//////////////// Eps : "<<Param.Euler.eps_value<<"///////////////"<<endl;
    cout<<"////////////// sigma : "<<Param.Euler.sig_value<<"///////////////"<<endl;
    cout<<"////////////// gamma : "<<Param.Euler.gamma<<"///////////////"<<endl;
    cout<<"/////// Pressure law : "<<Param.Euler.lp<<"///////////////"<<endl;
    cout<<"////// gravity vector: "<<Param.Euler.gx<<" and "<<Param.Euler.gy<<" ///////"<<endl;
     cout<<"/////// Local WB order  : "<<Param.Euler.LHO<<"///////////////"<<endl;
  }
  
   if(!strcmp(d.Typemodel,"Advection")){
     cout<<"////////// Solve advection equation //////////"<<endl;
     if(d.scheme==1 && d.Typescheme =='N'){
       cout<<"//////////////////// Linear nodal scheme ////////////////////"<<endl;
     }
     if(d.scheme==1 && d.Typescheme =='E'){
       cout<<"////////////////// Linear Edge scheme ////////////////////"<<endl;
     }
     cout<<"////////// Order :"<<Param.Adv.OrderAdv <<"///////////////"<<endl;
    
  }
  if(!strcmp(d.Typemodel,"M1")){
        cout<<"////////////// Solve M1 model ///////////////"<<endl;
      if(d.scheme==1){
      cout<<"//////////// Nodal Lagrange+remap scheme /////////////////"<<endl;}
      if(d.scheme==2){
	cout<<"////////// Nodal Lagrange+remap JL-(b) scheme ///////////////"<<endl;}
      if(d.scheme==3){
	cout<<"///// Noda Lagrange+remap JL-(b) scheme with local source term ///////"<<endl;}    
        
      cout<<"//////////////// Eps : "<<Param.M1.eps_value<<"///////////////"<<endl;
      cout<<"////////////// sigma : "<<Param.M1.sig_value<<"///////////////"<<endl;
  }

  if(!strcmp(d.Typemodel,"M1Matter")){
    cout<<"////////////// Solve M1 model with matter ///////////////"<<endl;
    if(d.scheme==1){
      cout<<"//////////// Nodal Lagrange+remap scheme /////////////////"<<endl;}
    if(d.scheme==2){
      cout<<"////////// Nodal Lagrange+remap JL-(b) scheme ///////////////"<<endl;}
    if(d.scheme==3){
      cout<<"///// Noda Lagrange+remap JL-(b) scheme with local source term ///////"<<endl;}
      cout<<"//////////////// Eps : "<<Param.M1.eps_value<<"///////////////"<<endl;
      cout<<"////////////// sigma : "<<Param.M1.sig_value<<"///////////////"<<endl;
  }

  if(!strcmp(d.Typemodel,"P1Matter")){
      cout<<"////////////// Solve P1 model ///////////////"<<endl;
      if(d.scheme==1){
      cout<<"////////////////// Nodal JL-(a) scheme /////////////////"<<endl;}
      if(d.scheme==2){
	cout<<"//////////////// Nodal JL-(b) scheme /////////////////"<<endl;}
      if(d.scheme==3){
	cout<<"///////// Nodal JL-(b) with local source term //////////"<<endl;}    
        
      cout<<"//////////////// Eps : "<<Param.P1M.eps_value<<"///////////////"<<endl;
      cout<<"////////////// sigma : "<<Param.P1M.sig_value<<"///////////////"<<endl;
      cout<<"//////////////// Cv : "<<Param.P1M.Cv_value<<"///////////////"<<endl;
      cout<<"///////////////// a : "<<Param.P1M.a_value<<"///////////////"<<endl;
  }
  
  cout<<"/////////////////////////////////////////////////////////"<<endl;
  cout<<"///////////////Mesh :"<<d.NameMesh<<"////////////////"<<endl;
  cout<<"///////////// Cell Mesh : "<<d.Nx<<" * "<<d.Ny<<"/////////////"<<endl;
  cout<<"///////////// Domain : "<<d.Tx<<" * "<<d.Ty<<"/////////////"<<endl;
  cout<<"/////////////// Final time : "<<d.Tf<<"//////////////"<<endl;
}


	 


void echanger(int * tableau, int i, int j){
  /** function which interchange two numbers (integer) in a table **/
   int temp;
   temp = tableau[i];
   tableau[i] = tableau[j];
   tableau[j] = temp;
}


void echangerd(double * tableau, int i, int j){
    /** function which interchange two numbes (double) in a table **/
   double temp;
   temp = tableau[i];
   tableau[i] = tableau[j];
   tableau[j] = temp;
}


void Tri(int * tableau, int longueur){
  /** tri function for integer table **/
   int i, maximum;
   while(longueur>0)
   {
      maximum = 0;
      for(i=1; i<longueur; i++)
         if(tableau[i]>tableau[maximum])
            maximum = i;
      echanger(tableau, maximum, longueur-1);
      longueur--;
   }
}


void Tridouble(double * tableau, int longueur){
  /** tri function for double table **/
  int i, maximum;
   while(longueur>0)
   {
      maximum = 0;
      for(i=1; i<longueur; i++)
         if(tableau[i]>tableau[maximum])
            maximum = i;
      echangerd(tableau, maximum, longueur-1);
      longueur--;
   }
}


double Valabs(double x){
  /** absolute value **/
  if (x<0.){
    return -x;
  }
  else return x;
}





    

