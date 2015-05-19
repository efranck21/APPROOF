#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include "Initialisation.hpp"
#include <math.h>
#include "AnalyticSolutions.hpp"
#include "Diagnostics.hpp"
#include "Functions.hpp"
#include "Flux.hpp"
#include "ParamPhys.hpp"

  void InitP1(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param){
     /** function which construct the initial datas 
      for the P1 model
  **/
  
    int in=0;
    double a=0;
    in=centerdomain(d,Mh);
     switch(d.nTest)
    {
   
       case 1 : 	
      for(int j=0;j<Mh.nc;j++){
	v.var[0][j]=CasTestOnde(d,Param,Mh.xj(j),0,0);
	v.var[1][j]=CasTestOnde(d,Param,Mh.xj(j),0,1);
	v.var[2][j]=CasTestOnde(d,Param,Mh.xj(j),0,2);
      }
      break;

       case 2 : 
		
      for(int j=0;j<Mh.nc;j++){
	v.var[0][j]=CasTestTelegraph(d,Param,Mh.xj(j),0,0);
	v.var[1][j]=CasTestTelegraph(d,Param,Mh.xj(j),0,1);
	v.var[2][j]=CasTestTelegraph(d,Param,Mh.xj(j),0,2);
      }
      break;
      
      
    case 3 : 
	
	for(int j=0;j<Mh.nc;j++){
	  v.var[0][j]=SolFonChaleur(d,Param,Mh.xj(j),0,Mh.xj(in));
	  v.var[1][j]=0.;
	  v.var[2][j]=0.; 
	}
      break;
      
     case 4 :
      for(int j=0;j<Mh.nc;j++){	
	v.var[0][j]=UniformSolution(d,Param,Mh.xj(j),0,0);
	v.var[1][j]=UniformSolution(d,Param,Mh.xj(j),0,1);
	v.var[2][j]=UniformSolution(d,Param,Mh.xj(j),0,2);
      }
      break;

      case 5 :
      double norm;
      for(int j=0;j<Mh.nc;j++){
	norm=sqrt(pow(Mh.xj(j).x,2.)+pow(Mh.xj(j).y,2.));
	v.var[0][j]=0.001+100*exp(-pow((norm/0.1),2.));
	v.var[1][j]=0;
	v.var[2][j]=0; /// A verifier
      }     
      break;
      
      case 6 :
      for(int j=0;j<Mh.nc;j++){
	v.var[0][j]=0;
	v.var[1][j]=0;
	v.var[2][j]=0; /// A verifier
      }     
      break;

       case 7 :
      for(int j=0;j<Mh.nc;j++){
	if( ((Mh.xj(j).x -1.0) * (Mh.xj(j).y -1.0))>0.0){
	  a=1;
	}  
	else  {
	  a=-1;
	} 
	v.var[0][j]=0;
	v.var[1][j]=0;
	v.var[2][j]=a; /// A verifier
      }     
      break;

       case 8 :
      for(int j=0;j<Mh.nc;j++){
	if( ((Mh.xj(j).x -1.0) * (Mh.xj(j).y -1.0))>0.0){
		  a=1;
	}  
	else  {
	  a=-1;
	} 
	v.var[0][j]=0;
	v.var[1][j]=a;
	v.var[2][j]=a; /// A verifier
      }     
      break;

       case 9 :
      for(int j=0;j<Mh.nc;j++){
	if( ((Mh.xj(j).x -1.0) * (Mh.xj(j).y -1.0))>0.0){
	  a=1;
	}  
	else  {
	  a=-1;
	}     
	v.var[0][j]=a;
	v.var[1][j]=a;
	v.var[2][j]=a; /// A verifier
      }     
      break;

      case 10 :
      for(int j=0;j<Mh.nc;j++){
	if( ((Mh.xj(j).x -1.0) * (Mh.xj(j).y -1.0))>0.0 && Mh.xj(j).x>1.0){
	  a=1;
	}  
	else  {
	  a=0;
	}     
	v.var[0][j]=a;
	v.var[1][j]=0;
	v.var[2][j]=0; /// A verifier
      }     
      break;

  }
}







