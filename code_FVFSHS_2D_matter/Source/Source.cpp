
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Source.hpp"


/** File which contains the functions for the source terms **/


vectorflux ChoiceSource(Data & d ,int numCell,Mesh & Mh, variable & v, TabConnecInv & tab,ParamPhysic & Param,R2* ur){
 /**  Function which gives the explicit source term associated with your model and your scheme (in Data) **/
vectorflux res(d);

  if(Param.Model ==1){      
    res=SourceDiff(d,numCell,Mh,v,Param);    /** Diffusion source term **/
  }

  if(Param.Model ==2){      
    res=Sourceadv(d,numCell,Mh,v,Param);    /** Advection source term **/
  }
    
  if(Param.Model ==3){  
    res=SourceP1(d,numCell,Mh,v,tab,Param); /** P1 source term **/
  }

 if(Param.Model ==4){
    res=SourceP1Matter(d,numCell,Mh,v,tab,Param); /** P1 with matter source term **/
  }
  
  if(Param.Model ==5){
    res=SourceE(d,numCell,Mh,v,tab,Param,ur); /** Euler source term **/
  }

  if(Param.Model ==6){
    res=SourceM1(d,numCell,Mh,v,tab,Param); /** M1 source term **/
  }

  if(Param.Model ==7){
    res=SourceM1Matter(d,numCell,Mh,v,tab,Param); /** M1 with Matter source term **/
  }

  
  return res;
}


R2 SemiImplicitPart(Data & d ,int numCell,Mesh & Mh, variable & v,variable & temp, TabConnecInv & tab,ParamPhysic & Param,double dt){
/**  Function which gives the implicit source term (scattering or friction source term)
 associated with your model and your scheme (in Data) **/
  R2 res(0,0);

  if(Param.Model ==3){
    res=SolveurImpliciteP1(d,numCell,Mh,v,tab,Param,dt,temp.var[1][numCell],temp.var[2][numCell]);
  }
  
  if(Param.Model ==4){
      res=SolveurImpliciteP1Matter(d,numCell,Mh,v,tab,Param,dt,temp.var[1][numCell],temp.var[2][numCell]);
  }
  
  if(Param.Model ==6){
    res=SolveurImpliciteM1(d,numCell,Mh,v,tab,Param,dt,temp.var[1][numCell],temp.var[2][numCell]);	   
  }
  if(Param.Model ==5){
    res=SolveurImpliciteE(d,numCell,Mh,v,tab,Param,dt,temp.var[1][numCell],temp.var[2][numCell],temp.var[0][numCell]);
  }

  return res;
}
