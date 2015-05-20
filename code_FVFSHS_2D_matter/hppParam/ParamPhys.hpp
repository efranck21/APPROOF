#ifndef PARAMPHYS_HPP
#define PARAMPHYS_HPP

#include <cmath>
#include <cstdlib>
#include <iostream> 
#include "ParamP1.hpp"
#include "ParamP1Matter.hpp"
#include "ParamDiff.hpp"
#include "ParamAdvection.hpp"
#include "ParamEuler.hpp"
#include "ParamM1.hpp"
#include "ParamM1Matter.hpp"
#include "ParamP1Compton.hpp"
#include "ClassMesh.hpp"
#include <cassert>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;


/** Class : generic class which contain the parameter associated with the model. 
All the subclass are construct but only the sub class associated with the model use 
is initialize.

Model is an integer associated with the model. This integer is use to avoid the test 
with string.

Model :
-Diffusion : Model = 1
-Advection : Model = 2
-P1 : Model = 3
-P1 with matter : Model = 4
-Euler : Model = 5
-M1 : Model = 6
-M1 with matter : Model = 7
-P1 with compton : Model =8
 **/

class ParamPhysic {
public:
  ParamAdvection Adv;
  ParamDiff Diff;
  ParamP1 P1;
  ParamEuler Euler; 
  ParamM1 M1;
  ParamM1Matter M1M ;
  ParamP1Matter P1M;
  ParamP1Compton P1C;
  int Model;

   ParamPhysic(){
     ParamP1();
     ParamDiff();
     ParamAdvection();
     ParamEuler();
     ParamM1();
     ParamM1Matter();
     ParamP1Matter();
     ParamP1Compton();
     Model=0; 
  }

  ParamPhysic(Data & d, Mesh & Mh){
    if (!strcmp(d.Typemodel,"Diffusion")){
      ParamDiff Diff_temp(d,Mh);
      Diff = Diff_temp;
      Model=1;
    }
    if (!strcmp(d.Typemodel,"Advection")){
      ParamAdvection Adv_temp(d,Mh);
      Adv = Adv_temp;
      Model=2;
    }  
    if (!strcmp(d.Typemodel,"P1")){
      ParamP1 P1_temp(d,Mh);
      P1 = P1_temp;
      Model=3;
    }
    if (!strcmp(d.Typemodel,"P1Matter")){
      ParamP1Matter P1M_temp(d,Mh);
      P1M = P1M_temp;
      Model=4;
    } 
    
    if (!strcmp(d.Typemodel,"Euler")){
      ParamEuler Euler_temp(d,Mh);
      Euler = Euler_temp;
      Model=5;
    } 
    if (!strcmp(d.Typemodel,"M1")){
      ParamM1 M1_temp(d,Mh);
      M1 = M1_temp;
      Model=6;
    } 
    if (!strcmp(d.Typemodel,"M1Matter")){
      ParamM1Matter M1M_temp(d,Mh);
      M1M = M1M_temp;
      Model=7;
    }

    if (!strcmp(d.Typemodel,"P1Compton")){
      ParamP1Compton P1C_temp(d,Mh);
      P1C = P1C_temp;
      Model=8;
    }

      
  }

   ParamPhysic(const ParamPhysic & ppy){
     Diff=ppy.Diff;
     P1=ppy.P1;
     Adv=ppy.Adv;
     Euler=ppy.Euler;
     M1=ppy.M1;
     M1M=ppy.M1M;
     P1M=ppy.P1M;
     P1C=ppy.P1C;
     Model=ppy.Model;
  }

  ~ParamPhysic(){
  }


  ParamPhysic & operator=( ParamPhysic ppy){
    Diff=ppy.Diff;
    P1=ppy.P1;
    Adv=ppy.Adv;
    Euler=ppy.Euler; 
    M1=ppy.M1;
    M1M=ppy.M1M;
    P1M=ppy.P1M;
    P1C=ppy.P1C;
    Model=ppy.Model;
    return *this;
 }    

}; 

void ParamPhysic_Init(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamPhysic  & Param);

void ParamPhysic_InitTab(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamPhysic  & Param);

#endif
