
#include <iostream>
#include <cmath>
#include <cassert>
#include "Flux.hpp"
#include "Variables.hpp"
#include "Functions.hpp"
#include <stdio.h>
#include "Tensor.hpp"
#include "BoundaryCondition.hpp"
#include "Functionsinit.hpp"
#include "Source.hpp"

vectorflux SourceDiff(Data & d,int numCell,Mesh & Mh, variable & v,ParamPhysic & Param){
  /** Function wich compute the explicit source term for diffusion **/
  vectorflux res(1);

  res.vflux[0]=0;
  return res;
}
