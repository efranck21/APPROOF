#include <iostream>
#include <cmath>
#include <cassert>
#include "Variables.hpp"
#include "Flux.hpp"
#include "Source.hpp"
#include "BoundaryCondition.hpp"
#include "Tensor.hpp"
#include "FunctionsAdvection.hpp"



vectorflux SourceM1Matter(Data & d,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab, ParamPhysic & Param){
    /** Function wich compute the explicit source terms for M1 **/
    
    
    vectorflux res(d);
    double kj=0;
    double Mr[2][2];
    R2 uj(0,0);
    
    int compton_scattering = 1 ; //(1 if taken into account, 0 if not)
    
 
    res.vflux[0] = 0. ;
    res.vflux[1] = -(4./3.)*((Param.M1M.sigma[numCell])/(Param.M1M.eps[numCell]*Param.M1M.eps[numCell]))*Mh.area(numCell)*v.var[1][numCell];
    res.vflux[2] = -(4./3.)*((Param.M1M.sigma[numCell])/(Param.M1M.eps[numCell]*Param.M1M.eps[numCell]))*Mh.area(numCell)*v.var[2][numCell];
    res.vflux[3] = 0. ;
            
    if (compton_scattering == 1)
      {
                
	double T = 1. ; // temperature electron
	double sigma_a = Param.M1M.sigma[numCell] ;
        
	double temp = 0. ;
	double source_E = 0. ;
	double source_F = 0. ;
	R2 g ;
        
	g.x = v.var[1][numCell] / v.var[0][numCell] ;
	g.y = v.var[2][numCell] / v.var[0][numCell] ;
        
	temp = 2. + sqrt(4.-3.*(g,g)) ;
        
	source_E = sqrt(sqrt( v.var[0][numCell] )) * temp*sqrt(sqrt(temp)) * ( (4.+(g,g))/(2.*temp-(g,g)) ) * sqrt(sqrt( (temp-(g,g))/(1.-(g,g)) )) ;
	
	source_F = sqrt(sqrt( v.var[0][numCell] )) * sqrt(sqrt(temp)) * ( (20.+(g,g))/(2.*temp+(g,g)) ) * sqrt(sqrt( (temp-(g,g))/(1.-(g,g)) )) ;
        
	res.vflux[0] += Mh.area(numCell)*sigma_a*(pow(v.var[3][numCell],4.) - v.var[0][numCell]) + Mh.area(numCell)*Param.M1M.sigma[numCell]*v.var[0][numCell]*( 4.*v.var[3][numCell] - source_E ) ;
	res.vflux[1] += Mh.area(numCell)*( -sigma_a -Param.M1M.sigma[numCell]/3. +Param.M1M.sigma[numCell]*( 4.*v.var[3][numCell] - source_F ) ) * v.var[1][numCell] ;
	res.vflux[2] += Mh.area(numCell)*( -sigma_a -Param.M1M.sigma[numCell]/3. +Param.M1M.sigma[numCell]*( 4.*v.var[3][numCell] - source_F ) ) * v.var[2][numCell] ;
	
	res.vflux[3] += - Mh.area(numCell)*sigma_a*(pow(v.var[3][numCell],4.) - v.var[0][numCell]) + Mh.area(numCell)*Param.M1M.sigma[numCell]*v.var[0][numCell]*( 4.*v.var[3][numCell] - source_E ) ;
      }
    


    return res;
}







