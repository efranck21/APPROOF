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

/** Explicit scheme 
 **/
void ExplicitDiscret(Data & d,Mesh & Mh, variable & v,TabConnecInv & TConnectInv, ParamPhysic & Param,double dt,int & ntAnim,double & time){
  
  int nt=0;
  variable temp(d,Mh.nc,Mh.nv);
  vectorflux flux(d);
  vectorflux source(d);
  temp=v;
  R2 **ur;
  double c=0;

  ur = new R2*[d.ngroup];
  for(int i=0;i<d.ngroup;i++){
    ur[i]= new R2[Mh.nv];
  }

  /** computaion of the table for variables parameters  **/ 
  ParamPhysic_InitTab(d,Mh,v,TConnectInv,Param);


  if(Param.Model == 5 && Param.Euler.LHO >1 ){
    for(int r=0;r<Mh.nv;r++){
      if(Mh.xr(r).lab>-1){
	/** We constructt the mean square problem for high order reconstruction 
        only for the node of the physical domain (not for the ghost node)**/
	LSInit(d,Mh,v,TConnectInv,Param.Euler,Param.Euler.TabLS[r],Param.Euler.LHO,r);
      }    
    }
  } 
   
   while(time<d.Tf)
     {
       if(time+dt>d.Tf) {
	 dt=d.Tf-time;
       }
      

       /** apply boundary condition (phantom cell) **/ 
       BoundaryCondition(d,Mh,v,Param,time);
  
       if(d.Anim=='y' && time>=d.dtAnim*ntAnim){
   
	 SaveRestart(d,Mh,v,Param,nt,time);   
	 SaveData(d,Mh,v,TConnectInv,Param,ntAnim);
	    
	 ntAnim++;
       }  
       /** Diagnostic : masse, velocity etc **/
        Diagnostics_Quantities(d,Mh,v,Param,nt,time);   

       /** Computation of the nodal flux solving nodal system **/
       if(Param.Model !=1 && d.Typescheme == 'N' ){	 
	  int r=0;
	 int g=0;
	 for(g=0;g<d.ngroup;g++){
	   for(r=0;r<Mh.nv;r++){
	     if(Param.Model == 5){
	       InitTab_rhoGravity(d,Mh,v,TConnectInv,Param.Euler,r);
	     }
	     ur[g][r]=SolveurNodal(d,r,Mh,v,TConnectInv,Param,g);
	   } 
	 }
	}

    
       int j=0;	
       /** Computation of the flux and source term and update of the time solution **/
       #pragma omp parallel for private(flux,source) 
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
	   //cout<<" ppp "<<flux.vflux[i]<<" "<<source.vflux[i]<<endl;
	   if (v.var[0][j]<0){
	     c++; 
	   }   
	    }
	 }
       }
       
       v=temp;
       
       /** computaion of the table for variables parameters  **/ 
       ParamPhysic_InitTab(d,Mh,v,TConnectInv,Param);
   
       time+=dt;
       nt++;
	
     }
  for(int i=0;i<d.ngroup;i++){
    delete [] ur[i];
  }
  delete [] ur;
   cout<<"negative coefficients " <<c<<" "<<" dt "<<dt<<endl;

}
