#include <iostream>
#include <cmath>
#include <cassert>
#include "Flux.hpp"
#include "Source.hpp"
#include "Variables.hpp"
#include "Functions.hpp"
#include <stdio.h>
#include "Tensor.hpp"
#include "BoundaryCondition.hpp"
#include "Functionsinit.hpp"


vectorflux Sourceadv(Data & d,int numCell, Mesh & Mh, variable & v,ParamPhysic & Param) {
  /** Function wich compute the explicit source term for advection **/
  vectorflux res(d);
  res.vflux[0]=0;
 
    return res;
 
}
