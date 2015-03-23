#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <iostream>
#include <cmath>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include "Variables.hpp"
#include "Functionsgeo.hpp"
#include "FunctionsConnect.hpp" 

#include "ParamPhys.hpp"



void WriteDataToMain(Data & d,ParamPhysic & Param);

void Tri(int * tableau, int longueur);

void Tridouble(double * tableau, int longueur);

void echanger(int * tableau, int i, int j);

void echangerd(double * tableau, int i, int j);

double positivepart(double x);


double Valabs(double x);



#endif
      
	  
	 
