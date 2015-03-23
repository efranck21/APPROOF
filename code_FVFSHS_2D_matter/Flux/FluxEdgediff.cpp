
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

vectorflux FluxVF4Diff(Data & d,int numCell, Mesh & Mh,  variable & v,  TabConnecInv & tab,ParamPhysic & Param){
   /**  Edge fluxes for linear diffusion scheme VF4   **/
    double dis=0.;
  double s=0.;
  int k,A,B,C,D;
  vectorflux res(1);
  int **tabCellEdge=NULL;// table with the neighbour cell for each edge//  
  Vertex p1,p2;


  //We find the neighbour cell for each edge //
  if(Mh.nbnodelocal==3){

    A=Mh(numCell,0);
    B=Mh(numCell,1);
    C=Mh(numCell,2);

    tabCellEdge = new int*[3];
    for(int i=0;i<3;i++){
      tabCellEdge[i]=new int[3];
    }
    
    tabCellEdge[0][0]=InverseEdge(Mh,A,B,numCell,tab);
    tabCellEdge[0][1]=A;
    tabCellEdge[0][2]=B;
    tabCellEdge[1][0]=InverseEdge(Mh,B,C,numCell,tab);
    tabCellEdge[1][1]=B;
    tabCellEdge[1][2]=C;
    tabCellEdge[2][0]=InverseEdge(Mh,C,A,numCell,tab);
    tabCellEdge[2][1]=C;
    tabCellEdge[2][2]=A;
  }
  else if(Mh.nbnodelocal==4){
    A=Mh(numCell,0);
    B=Mh(numCell,1);
    C=Mh(numCell,2);
    D=Mh(numCell,3);

    tabCellEdge = new int*[4];
    for(int i=0;i<4;i++){
      tabCellEdge[i]=new int[3];
    }
    
    tabCellEdge[0][0]=InverseEdge(Mh,A,B,numCell,tab);
    tabCellEdge[0][1]=A;
    tabCellEdge[0][2]=B;
    tabCellEdge[1][0]=InverseEdge(Mh,B,C,numCell,tab);
    tabCellEdge[1][1]=B;
    tabCellEdge[1][2]=C;
    tabCellEdge[2][0]=InverseEdge(Mh,C,D,numCell,tab);
    tabCellEdge[2][1]=C;
    tabCellEdge[2][2]=D;
    tabCellEdge[3][0]=InverseEdge(Mh,D,A,numCell,tab);
    tabCellEdge[3][1]=D;
    tabCellEdge[3][2]=A;
  }  
  else
    { cout << "problem in the flux edge diff"<<endl;}

  for(int i=0;i<Mh.nbnodelocal;i++){
    k=tabCellEdge[i][0];
        
      dis=sqrt(pow(Mh.xj(k).x-Mh.xj(numCell).x,2)+pow(Mh.xj(k).y-Mh.xj(numCell).y,2));
      s=s+Mh.ljk(numCell,i)*(v.var[0][k]-v.var[0][numCell])/dis; 
  }
 
  s=s*(Param.Diff.D[numCell]);    
  for(int i=0;i<Mh.nbnodelocal;i++){
    delete [] tabCellEdge[i];
  }
  delete [] tabCellEdge;
  res.vflux[0]=s;
      
  
  
  return res;      
  
}
