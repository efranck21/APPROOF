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

void InitM1Matter(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param){
    /** function which construct the initial datas
     for the p1 model with matter
     **/
    double T=0;
    
    switch(d.nTest)
    {
            
        case 1 :
            for(int j=0;j<Mh.nc;j++){
	      if(Mh.xj(j).x < 0.5*d.Tx){
		v.var[0][j]=2;
		v.var[1][j]=0;
		v.var[2][j]=0;
		v.var[3][j]=1.e-14;
	      }
	      else {
		v.var[0][j]=1.e-14;
		v.var[1][j]=0;
		v.var[2][j]=0;
		v.var[3][j]=2;
	      }
            }
            break;
            
        case 2 :
            // Oslon
            T=0.056234;
            for(int j=0;j<Mh.nc;j++){
                v.var[0][j]=Param.M1M.a_value*pow(T,4.);
                v.var[1][j]=0;
                v.var[2][j]=0;
                v.var[3][j]=T;
            }
            break;
            
    }
}








