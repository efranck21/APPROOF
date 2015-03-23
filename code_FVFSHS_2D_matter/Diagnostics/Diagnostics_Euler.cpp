
#include <iostream>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include "Functions.hpp"
#include "Diagnostics.hpp"


/** File which contains the function associated to the diagnostics (norm Lp, error Lp, masse etc) specific to the Euler equations **/

double NormLP_Velocity(Data & d, Mesh & Mh,variable & v,ParamPhysic & Param,double p,double t){
  /** Function which compute the Lp norm associated with the discrete velocity for the Euler equations **/
  double res=0.;
  int in=0;
  double hx;
  double a1=0,a2=0;
  double rho=0;
  
  hx=1/double(d.Nx);
       if(Mh.nbnodelocal==4){
	   in=(Mh.nc+d.Nx)/2;}
	 if(Mh.nbnodelocal==3){
	   in=((Mh.nc+d.Nx*2)/2);}
	 if(!strcmp(d.NameMesh,"KershawMesh")){
	   for(int j=0;j<Mh.nc;j++){
	     if(Mh.xj(j).x<=1+d.Tx*(1./d.Nx) && Mh.xj(j).x>=1-d.Tx*(1./d.Nx)){
	       if(Mh.xj(j).y<=1+d.Ty*(1./d.Ny) && Mh.xj(j).y>=1-d.Ty*(1./d.Ny)){
		 in=j;
	       }
	     }
	   }
	 }

	 if(d.dimsave==2){

	     for(int j=0;j<Mh.nc;j++){
	        if(Mh.cells[j].lab!=-1){
		  rho=v.var[0][j];
		  a1=Mh.area(j)*(pow(Valabs(v.var[1][j]/rho),p));
		  a2=Mh.area(j)*(pow(Valabs(v.var[2][j]/rho),p));
	       res=res+a1+a2;

		}

	   }
	 }
  
  res=pow(res,1/p);
     return res;

}


double NormLP_Error_All_variables_Euler(Data & d, Mesh & Mh,variable & v,ParamPhysic & Param,double p,double t){
    /** Function which compute the Lp error associated with the discrete quantities and the analytical solutions
 (rho, u, E) for the Euler equations **/
  double res=0.;
  int in=0;
  double hx;
  double a1=0,a2=0,a3=0,a4=0;
  double rhoNum=0,rhoExact=0;
  
  hx=1/double(d.Nx);
       if(Mh.nbnodelocal==4){
	   in=(Mh.nc+d.Nx)/2;}
	 if(Mh.nbnodelocal==3){
	   in=((Mh.nc+d.Nx*2)/2);}
	 if(!strcmp(d.NameMesh,"KershawMesh")){
	   for(int j=0;j<Mh.nc;j++){
	     if(Mh.xj(j).x<=1+d.Tx*(1./d.Nx) && Mh.xj(j).x>=1-d.Tx*(1./d.Nx)){
	       if(Mh.xj(j).y<=1+d.Ty*(1./d.Ny) && Mh.xj(j).y>=1-d.Ty*(1./d.Ny)){
		 in=j;
	       }
	     }
	   }
	 }

	 if(d.dimsave==2){
	     for(int j=0;j<Mh.nc;j++){
	        if(Mh.cells[j].lab!=-1){
		  if(Param.Euler.LHO < 2) {
		    rhoNum=v.var[0][j];
		    rhoExact=SolFon(d,Param,Mh.xj(j),t,Mh.xj(in),0);
		    a1=Mh.area(j)*(pow(Valabs(SolFon(d,Param,Mh.xj(j),t,Mh.xj(in),0)-v.var[0][j]),p));
		    a2=Mh.area(j)*(pow(Valabs(SolFon(d,Param,Mh.xj(j),t,Mh.xj(in),1)-v.var[1][j]),p));
		    a3=Mh.area(j)*(pow(Valabs(SolFon(d,Param,Mh.xj(j),t,Mh.xj(in),2)-v.var[2][j]),p));
		    a4=Mh.area(j)*(pow(Valabs(SolFon(d,Param,Mh.xj(j),t,Mh.xj(in),3)-v.var[3][j]),p));
		    res=res+a1+a2+a3+a4;
		  }
		  else  {
		    rhoNum=v.var[0][j];
		    rhoExact=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),t,Mh.xj(in),j,0);
		    a1=Mh.area(j)*(pow(Valabs(SolFonEulerAverage(d,Mh,Param,Mh.xj(j),t,Mh.xj(in),j,0)-v.var[0][j]),p));
		    a2=Mh.area(j)*(pow(Valabs(SolFonEulerAverage(d,Mh,Param,Mh.xj(j),t,Mh.xj(in),j,1)-v.var[1][j]),p));
		    a3=Mh.area(j)*(pow(Valabs(SolFonEulerAverage(d,Mh,Param,Mh.xj(j),t,Mh.xj(in),j,2)-v.var[2][j]),p));
		    a4=Mh.area(j)*(pow(Valabs(SolFonEulerAverage(d,Mh,Param,Mh.xj(j),t,Mh.xj(in),j,3)-v.var[3][j]),p));
		    res=res+a1+a2+a3+a4;
		  }
		}	
	   }
	 }

  res=pow(res,1/p);
     return res;	
}


