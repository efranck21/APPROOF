#include <cmath>
#include <cstdlib>
#include <iostream> 
#include "ParamPhys.hpp"
#include <cassert>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;


/** File which contained some generic functions to intitialize the paramater of the model **/ 

void ParamPhysic_Init(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamPhysic & Param){
 /** Initialization of the model parameters **/
  if (Param.Model == 1){
    ParamDiff_Init(d,Mh,v,tab,Param.Diff,"ParamPhysic/ParamDiff.dat");
  }
  if (Param.Model == 2){
    ParamAdv_Init(d,Mh,v,tab,Param.Adv,"ParamPhysic/ParamAdv.dat");
  }  
  if (Param.Model == 3){
    ParamP1_Init(d,Mh,v,tab,Param.P1,"ParamPhysic/ParamP1.dat");
  } 
  if (Param.Model == 4){
    ParamP1Matter_Init(d,Mh,v,tab,Param.P1M,"ParamPhysic/ParamP1Matter.dat");
  } 
  if (Param.Model == 5){
    ParamEuler_Init(d,Mh,v,tab,Param.Euler,"ParamPhysic/ParamEuler.dat");
  } 
  if (Param.Model == 6){
    ParamM1_Init(d,Mh,v,tab,Param.M1,"ParamPhysic/ParamM1.dat");
  }
  if (Param.Model == 7){
    ParamM1M_Init(d,Mh,v,tab,Param.M1M,"ParamPhysic/ParamM1Matter.dat");
  }
  
}

void ParamPhysic_InitTab(Data & d,Mesh & Mh,variable &v,TabConnecInv & tab, ParamPhysic & Param){
  /** Function which initialize the table of quantities using the given values **/
  if (Param.Model == 1){
    ParamDiff_InitTab(d,Mh,v,tab,Param.Diff);
  }
  if (Param.Model == 2){
    ParamAdv_InitTab(d,Mh,v,tab,Param.Adv);
  }  
  if (Param.Model == 3){
    ParamP1_InitTab(d,Mh,v,tab,Param.P1);
  } 
  if (Param.Model == 4){
    ParamP1Matter_InitTab(d,Mh,v,tab,Param.P1M);
  } 
  if (Param.Model == 5){
    ParamEuler_InitTab(d,Mh,v,tab,Param.Euler);
  } 
  if (Param.Model == 6){
    ParamM1_InitTab(d,Mh,v,tab,Param.M1);
  }
  if (Param.Model == 7){
    ParamM1M_InitTab(d,Mh,v,tab,Param.M1M);
  }
  
}
