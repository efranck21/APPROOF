
#include <iostream>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include "Functions.hpp"
#include "Diagnostics.hpp"
#include "ParamPhys.hpp"

/** File which contains the function associated to the diagnostics (norm Lp, error Lp, masse etc) for all the models **/



double NormLP_Error_All_variables(Data & d, Mesh & Mh,variable & v,ParamPhysic & Param,double p,double t){
       /** Function which compute the Lp error associated with all discrete quantities and 
all analytical solutions for the different models **/
  double res=0.;
  int in=0;
  double hx;
  double a=0;
  
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
	  for(int k=0;k<v.nbvar;k++){
	    a=Mh.area(j)*(pow(Valabs(SolFon(d,Param,Mh.xj(j),t,Mh.xj(in),k)-v.var[k][j]),p)); 
	    res=res+a;
	  }
	}
      }
      
  }
  
  res=pow(res,1/p);
  return res;	
}

double NormLP_Error_One_variable(Data & d, Mesh & Mh,variable & v,ParamPhysic & Param,double p,double t,int var){
      /** Function which compute the Lp error associated with one discrete quantity and 
one analytical solution for the different models **/
  double res=0.;
  int in=0;
  double hx;
  
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
		res=res+Mh.area(j)*pow(Valabs(SolFon(d,Param,Mh.xj(j),t,Mh.xj(in),var)-v.var[var][j]),p);
	      }
	   }
	 }
	
	 if(d.dimsave==1){
	   for(int j=0;j<Mh.nc;j++){
	      if(Mh.cells[j].lab!=-1){	     
		res=res+hx*pow(Valabs(SolFon(d,Param,Mh.xj(j),t,Mh.xj(in),var)-v.var[var][j]),p);
	      }
	   }
	 }
res=pow(res,1/p);
return res;
	 
}

double NormLP_Stability(Data & d, Mesh & Mh,variable & v,ParamPhysic & Param,double p){
   /** Function which compute the Lp norm associated with one discrete quantity for the different models **/
  double res=0.;
  double sum=0;
  

  for(int j=0;j<Mh.nc;j++){
     if(Mh.cells[j].lab!=-1){
    for(int i=0;i<v.nbvar;i++){
      sum=sum+pow(Valabs(v.var[i][j]),p);
    }
    res=res+Mh.area(j)*sum;
    sum=0;
     }
  }		  
  res=pow(res,1/p);
     return res;	 
}


double Mass(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,double time){
  /** Function which compute the mass associated with the first quantity for the different models **/
  double mass=0;
  double p=1.;
  
      for(int j=0;j<Mh.nc;j++){
	 if(Mh.cells[j].lab!=-1){
	   mass=mass+Mh.area(j)*pow(v.var[0][j],p);
	 }   
    }
      return mass;
}


void Diagnostics_Quantities(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,int nt,double time){
  /** Function which gives the diagnostics for the numerical stability **/
  double Mass_0=0;
  double p=1.;
  
    if((nt%20)==0) {
      Mass_0 = Mass(d,Mh,v,Param,time);
      cout<<"Iteration n "<<nt<<" correspondant au temps :"<<time<<endl;
      cout<<"masse "<<pow(Mass_0,1./p)<<endl;

      if(!strcmp(d.Typemodel,"Euler")){
	cout<<"velocity :"<<NormLP_Velocity(d,Mh,v,Param,2,time)<<endl;
      }
       
    }
 
}

double Diagnostics_Error(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,double time, double p){
  double res=0;
  /** Function which gives the errors for the test case and the model choosen in Data **/
  if(!strcmp(d.Typemodel,"Euler")){
    res=NormLP_Error_All_variables_Euler(d,Mh,v,Param,p,time);
  }
  if(Param.Model == 1){
    res=NormLP_Error_One_variable(d,Mh,v,Param,p,time,0);
  }
  if(Param.Model == 2){
    res=NormLP_Error_One_variable(d,Mh,v,Param,p,time,0);
  }
   if(Param.Model == 3){
     if(d.nTest == 3 ){
       res=NormLP_Error_One_variable(d,Mh,v,Param,p,time,0);
     } 
     else  
       {
	 res=NormLP_Error_All_variables(d,Mh,v,Param,p,time);
       }
   }

   if(Param.Model == 6){
       if(d.nTest != 3 ){
       res=NormLP_Error_One_variable(d,Mh,v,Param,p,time,0);
     } 
     else  
       {
	 res=NormLP_Error_All_variables(d,Mh,v,Param,p,time);
       }
   }
   return res;
}
