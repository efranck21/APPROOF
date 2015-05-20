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
#include "WriteData.hpp"

using namespace std;

void ChoiceInit(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,double & InitTime){
  /** functions which call the functions which construct the initial datas 
      for the dfferent model
  **/
  if(d.restart==0) {
    if(Param.Model == 1) {InitDiff(d,Mh,v,Param);}
    
    if(Param.Model == 2) {Initadv(d,Mh,v,Param);}
    
    if(Param.Model == 3) {InitP1(d,Mh,v,Param);}

    if(Param.Model == 4) {InitP1Matter(d,Mh,v,Param);}
    
    if(Param.Model == 5) {InitE(d,Mh,v,Param);}
   
    if(Param.Model == 6) {InitM1(d,Mh,v,Param);}
      
    if(Param.Model == 7) {InitM1Matter(d,Mh,v,Param);}

    if(Param.Model == 8) {InitP1Compton(d,Mh,v,Param);}


  }
  else{
    LoadRestart(d,Mh,v,Param,InitTime);
  }
}











