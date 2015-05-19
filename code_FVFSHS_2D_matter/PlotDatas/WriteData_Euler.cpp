#include <iostream>
#include <cassert>
#include "Data.hpp"
#include "WriteData.hpp"
#include "Functions.hpp"
#include "AnalyticSolutions.hpp"
#include "Diagnostics.hpp"
#include "Flux.hpp"

/** File which contains the functions for the plot of Euler quantities **/


void SaveInterEner(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamPhysic & Param,int ntAnim){
/** Function which plot the internal energy for Eler model **/
  char string[255];
  FILE *fileresult=NULL;
  int j,n;
  
  if(d.Anim=='y'){ sprintf(string,"%s%s%04d.%s","./DATA/","InterEner",ntAnim,d.suffixe);}

  else { sprintf(string,"%s%s%s","./DATA/","var",d.suffixe);}
  fileresult=fopen(string,"w");
  for(j=0;j<Mh.nc;j++){
     if(Mh.cells[j].lab!=-1){
      
      if(Mh.nbnodelocal==4)
	{n=Mh(j,0);
	  fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,EnergyIn(d,Mh,v,Param.Euler,j));
	n=Mh(j,1);
	fprintf(fileresult,"%e %e %e\n\n",Mh.xr(n).x,Mh.xr(n).y,EnergyIn(d,Mh,v,Param.Euler,j));
	n=Mh(j,3);
	fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,EnergyIn(d,Mh,v,Param.Euler,j));
	n=Mh(j,2);
	fprintf(fileresult,"%e %e %e\n\n\n",Mh.xr(n).x,Mh.xr(n).y,EnergyIn(d,Mh,v,Param.Euler,j));
       
	
	}
      if(Mh.nbnodelocal==3)
	{n=Mh(j,0);
	  fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,EnergyIn(d,Mh,v,Param.Euler,j));
	n=Mh(j,1);
	fprintf(fileresult,"%e %e %e\n\n",Mh.xr(n).x,Mh.xr(n).y,EnergyIn(d,Mh,v,Param.Euler,j));
	n=Mh(j,2);
	fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,EnergyIn(d,Mh,v,Param.Euler,j));
	n=Mh(j,2);
	fprintf(fileresult,"%e %e %e\n\n\n",Mh.xr(n).x,Mh.xr(n).y,EnergyIn(d,Mh,v,Param.Euler,j));       
	
	}

      }

  }

  fclose(fileresult);
}


void SaveFonction2DEuler(Data & d,Mesh & Mh,TabConnecInv & tab,ParamPhysic & Param,float t,Vertex c,int Num){
/** Function which plot the non conservative quantities for Euler**/
  char string[255];
  FILE *fileresult=NULL;
  int j,n;
  double rho;
  int in=0;

   sprintf(string,"%s%s%d","./DATA/","Solfond",Num);
  fileresult=fopen(string,"w");

  for(j=0;j<Mh.nc;j++){
     if(Mh.cells[j].lab!=-1){
       if(d.nTest !=1 && d.nTest !=3 && d.nTest !=5){
       if(Num ==0 ) {
	 rho=1;
       }
       else {
	 rho=SolFonEuler(d,Param,Mh.xj(j),t,c,0);
       }
      if(Mh.nbnodelocal==4)
	{n=Mh(j,0);
	  fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,SolFon(d,Param,Mh.xj(j),t,c,Num)/rho);
	n=Mh(j,1);
	fprintf(fileresult,"%e %e %e\n\n",Mh.xr(n).x,Mh.xr(n).y,SolFon(d,Param,Mh.xj(j),t,c,Num)/rho);
	n=Mh(j,3);
	fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,SolFon(d,Param,Mh.xj(j),t,c,Num)/rho);
	n=Mh(j,2);
	fprintf(fileresult,"%e %e %e\n\n\n",Mh.xr(n).x,Mh.xr(n).y,SolFon(d,Param,Mh.xj(j),t,c,Num)/rho);
       
	
	}
      if(Mh.nbnodelocal==3)
	{n=Mh(j,0);
	  fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,SolFon(d,Param,Mh.xj(j),t,c,Num)/rho);
	n=Mh(j,1);
	fprintf(fileresult,"%e %e %e\n\n",Mh.xr(n).x,Mh.xr(n).y,SolFon(d,Param,Mh.xj(j),t,c,Num)/rho);
	n=Mh(j,2);
	fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,SolFon(d,Param,Mh.xj(j),t,c,Num)/rho);
	n=Mh(j,2);
	fprintf(fileresult,"%e %e %e\n\n\n",Mh.xr(n).x,Mh.xr(n).y,SolFon(d,Param,Mh.xj(j),t,c,Num)/rho);
       
	
	}  
       }
       else{
       if(Num ==0 ) {
	 rho=1;
       }
       else {
	 rho=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),t,Mh.xj(in),j,0);
       }
      if(Mh.nbnodelocal==4)
	{n=Mh(j,0);
	  fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,SolFonEulerAverage(d,Mh,Param,Mh.xj(j),t,Mh.xj(in),j,Num)/rho);
	n=Mh(j,1);
	fprintf(fileresult,"%e %e %e\n\n",Mh.xr(n).x,Mh.xr(n).y,SolFonEulerAverage(d,Mh,Param,Mh.xj(j),t,Mh.xj(in),j,Num)/rho);
	n=Mh(j,3);
	fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,SolFonEulerAverage(d,Mh,Param,Mh.xj(j),t,Mh.xj(in),j,Num)/rho);
	n=Mh(j,2);
	fprintf(fileresult,"%e %e %e\n\n\n",Mh.xr(n).x,Mh.xr(n).y,SolFonEulerAverage(d,Mh,Param,Mh.xj(j),t,Mh.xj(in),j,Num)/rho);
       
	
	}
      if(Mh.nbnodelocal==3)
	{n=Mh(j,0);
	  fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,SolFonEulerAverage(d,Mh,Param,Mh.xj(j),t,Mh.xj(in),j,Num)/rho);
	n=Mh(j,1);
	fprintf(fileresult,"%e %e %e\n\n",Mh.xr(n).x,Mh.xr(n).y,SolFonEulerAverage(d,Mh,Param,Mh.xj(j),t,Mh.xj(in),j,Num)/rho);
	n=Mh(j,2);
	fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,SolFonEulerAverage(d,Mh,Param,Mh.xj(j),t,Mh.xj(in),j,Num)/rho);
	n=Mh(j,2);
	fprintf(fileresult,"%e %e %e\n\n\n",Mh.xr(n).x,Mh.xr(n).y,SolFonEulerAverage(d,Mh,Param,Mh.xj(j),t,Mh.xj(in),j,Num)/rho);
       
	
	}  
       }
     }
       
 }
  fclose(fileresult);

}



