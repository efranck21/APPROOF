
#include <iostream>
#include <cmath>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include "Functions.hpp"
#include "Data.hpp"
#include "Functionsgeo.hpp"
#include "Tensor.hpp"
#include "Functionsinit.hpp"
#include "Flux.hpp"
#include "ParamPhys.hpp"

/** File with fucntions which allows to initialize some quantities**/



double CFL(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamPhysic & Param){
  /** Function which initiliaze the time step for your schemeand test case **/
double dt=0;
 double hx=0,hy=0,h=0;
  hx=1/double(d.Nx);
  hy=1/double(d.Ny);
  h=min(hx,hy);

 if(d.Typetime=='S'){
   if(Param.Model == 3 &&  d.scheme!=3){ //P1 model
     dt=d.CFL*Param.P1.eps[0]*(pow(h,d.Total_order));
   }
   if(Param.Model == 3  &&  d.scheme==3){ // P1 Model
     dt=d.CFL*Param.P1.sigma[0]*h*h;
   }

   if(Param.Model == 4){ // P1 matter Model
     dt=d.CFL*Param.P1M.eps[0]*h;
   }
   
   if(Param.Model == 5  &&  d.scheme==3){ // Euler Model
     dt=d.CFL*h*h;
   }
  if(Param.Model == 5  &&  d.scheme==4){  // Euler Model
     dt=d.CFL*Param.Euler.eps[0]*h;
   }
  if(Param.Model == 6 && d.nTest <3 ){  // M1 Model
     dt=d.CFL*h*h;
   }
  if (Param.Model == 7) // M1 Model with matter
  {
    cout << "M1 Model with matter not define in semi implicite configuration"<< endl ;
    exit(1) ;
  }
 }

 if(d.Typetime=='E'){
 
  if(Param.Model ==1){
    dt=(d.CFL*(Param.Diff.D[0])*h*h); // Diffusion model
  }

  if(Param.Model == 2){
    dt=d.CFL*pow(h,Param.Adv.OrderAdv); //Advection model-
  }
  
  if(Param.Model == 3){
    dt=d.CFL*Param.P1.eps[0]*(pow(h,d.Total_order)); // P1 model
  }

  if(Param.Model == 4){
    dt=d.CFL*Param.P1M.eps[0]*hx; // P1 Matter model
 }
 
  if(Param.Model == 5){ 
    dt=d.CFL*Param.Euler.eps[0]*h; // Euler model
  }
  
  if(Param.Model == 6 && d.nTest >2 ){
    dt=d.CFL*pow(h,Param.M1.OrderAdv); // M1 model
  }
  if(Param.Model == 6 && d.nTest <3 ){
    dt=d.CFL*Param.M1.eps_value*h; // M1 model
  }
  if(Param.Model == 7 && d.nTest >2 ){
    dt=d.CFL*pow(h,Param.M1M.OrderAdv); // M1 model with matter
  }
  if(Param.Model == 7 && d.nTest <3 ){
    dt=d.CFL*pow(h,Param.M1M.OrderAdv);
  }

 }
 return dt;
}





