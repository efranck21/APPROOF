
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "TimeDiscret.hpp"
#include <fstream>
#include "Flux.hpp"
#include "BoundaryCondition.hpp"
#include "AnalyticSolutions.hpp"
#include "Diagnostics.hpp"
#include "FunctionsAdvection.hpp"
#include "HOfunctions.hpp"
#include <assert.h>

using namespace std;


int SIschemeImplemented(Data & d,Mesh & Mh, variable & v,TabConnecInv & TConnectInv, ParamPhysic & Param){
  /** function which allows to known if the semi implicit scheme exist or not **/
  int res=0;
  
  /** if res=0 the semi implicit scheme does not exist, if res=1 the semi implicit scheme exist **/

  
  if(Param.Model == 1) {
    res=0;  
  }

  if(Param.Model == 2) {
    res=0;
  }

  if(Param.Model == 3) {
     switch(d.scheme){
    case 1:
      res=1;
      break;
    case 2:
      res=0;
      break;  
    case 3:
      res=1;
      break;  
    default: cout <<" the scheme associated with the number"<<d.scheme<<" does not exist "<<endl;
    }
  }

  if(Param.Model == 4) {
       switch(d.scheme){
    case 1:
      res=1;
      break;
    case 2:
      res=0;
      break;  
    case 3:
      res=1;
      break;  
    default: cout <<" the scheme associated with the number"<<d.scheme<<" does not exist "<<endl;
    }

  }

  if(Param.Model == 5) {
       switch(d.scheme){
    case 1:
      res=1;
      break;
    case 2:
      res=0;
      break;  
    case 3:
      res=1;
      break;
    case 4:
      res=1;
      break;  
    default: cout <<" the scheme associated with the number"<<d.scheme<<" does not exist "<<endl;
    }

  }

  if(Param.Model == 6) {
       switch(d.scheme){
    case 1:
      res=1;
      break;
    case 2:
      res=0;
      break;  
    case 3:
      res=1;
      break;  
    default: cout <<" the scheme associated with the number"<<d.scheme<<" does not exist "<<endl;
    }

  }
  return res;

}

void SemiImplicit(Data & d,Mesh & Mh, variable & v,TabConnecInv & TConnectInv, ParamPhysic & Param,double dt,int & ntAnim, double & time){

  int nt=0;
  variable temp(d,Mh.nc,Mh.nv);
  vectorflux flux(d);
  vectorflux source(d);
  R2 res(0,0);
  R2 x(0,0);
  temp=v;
  R2 *ur =new R2[Mh.nv];
  double c=0;
 
  /** computaion of the table for variables parameters  **/ 
  ParamPhysic_InitTab(d,Mh,v,TConnectInv,Param);

  if(Param.Model == 5){
    for(int r=0;r<Mh.nv;r++){
      if(Mh.xr(r).lab>-1){
	LSInit(d,Mh,v,TConnectInv,Param.Euler,Param.Euler.TabLS[r],Param.Euler.LHO,r);
      }    
    }
  }

   while(time<d.Tf)
     {
       if(time+dt>d.Tf) {
	 dt=d.Tf-time;
       }

       if(d.Anim=='y' && time>=d.dtAnim*ntAnim){
         SaveRestart(d,Mh,v,Param,nt,time);
	 SaveData(d,Mh,v,Param,ntAnim);
 
	 ntAnim++;
       } 

       /** apply boundary condition (phantom cell) **/ 
       BoundaryCondition(d,Mh,v,Param,time);              

       /** Diagnostic : masse, velocity etc **/
       Diagnostics_Quantities(d,Mh,v,Param,nt,time);   
  
       int j=0;     
       /** Computation of the nodal flux solving nodal system **/
        if(strcmp(d.Typemodel,"Advection") && d.Typescheme == 'N' ){	 
	 int r=0;
         #pragma omp parallel for
	 for(r=0;r<Mh.nv;r++){
	   if(Param.Model == 5){
	     InitTab_rhoGravity(d,Mh,v,TConnectInv,Param.Euler,r);
	   }
	   ur[r]=SolveurNodal(d,r,Mh,v,TConnectInv,Param);
	 } 
       }
  
       /** Computation of the flux and source term and update of the time solution **/
#pragma omp parallel for private(flux,source,x)  
       
       for( j=0;j<Mh.nc;j++){
	 if(Mh.cells[j].lab!=-1){
	   
	   if(d.Typescheme == 'N' ){
	     flux=ChoiceFluxN(d,j,Mh,v,TConnectInv,Param,ur);
	   }             
	   if(d.Typescheme == 'E' ){
	     flux=ChoiceFluxE(d,j,Mh,v,TConnectInv,Param);
	   }
	   source=ChoiceSource(d,j,Mh,v,TConnectInv,Param,ur);

	   for(int i=0;i<v.nbvar;i++){	
	     temp.var[i][j]=v.var[i][j]+(dt/Mh.area(j))*(flux.vflux[i]+source.vflux[i]);
	   }
	 	  
	   x=SemiImplicitPart(d,j,Mh,v,temp,TConnectInv,Param,dt);	   
	   temp.var[1][j]=x.x;
	   temp.var[2][j]=x.y;
		    
	   if (temp.var[0][j]<0){
	     c++; 
	   }
	 }	 
       }

        /** Computation of the source terms for the matter relaxation and update of the time solution **/
       if(!strcmp(d.Typemodel,"P1Matter")){
	 for(j=0;j<Mh.nc;j++){
	   if(Mh.cells[j].lab!=1){
	     res=SplittingCalculMatter(d,Param,j,Mh,temp,v,TConnectInv,dt,0.00000001);
	     temp.var[0][j]=res.x;
	     temp.var[v.nbvar-1][j]=res.y;
	   }
	 }
       }
  		
	v=temp;
	/** computaion of the table for variables parameters  **/ 
	ParamPhysic_InitTab(d,Mh,v,TConnectInv,Param);
	time+=dt;
	nt++;
	  
     }
   delete [] ur;
    cout<<"negative coefficients " <<c<<" "<<" dt "<<dt<<endl;
}




