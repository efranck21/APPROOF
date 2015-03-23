
#include <cassert>
#include <iostream>
#include <cmath>
#include<cstdlib>
#include "Data.hpp"
#include "ClassMesh.hpp"
#include "Tensor.hpp"
#include "Functions.hpp"
#include "Functionsinit.hpp"
#include "Flux.hpp"

tensor inittensor(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,char c,int numGr,int numCell){
  tensor t(d,c);
   R2 xj1(0,0);
  R2 xj2(0,0);
  R2 w(0,0);
  R2 vv(0,0);
  double T;
  double coP1=1;
  int xg,xd,e1,e2,r;
  int ts,th;
  double ld=geomCorrection(d,Mh,numCell);
  ts = 3;
  th = 2; 

  if (c=='s'){ t.num = ts; }
  if (c=='h'){ t.num = th; }
   
  r=NodeGtoL(Mh,numCell,numGr);


  if(c=='s'){

    if(ts==1){
      t.ten[0][0]=coP1*(Mh.xr(numCell,r).x-Mh.xj(numCell).x)*Mh.njr(numCell,r).x;
      t.ten[0][1]=coP1*(Mh.xr(numCell,r).y-Mh.xj(numCell).y)*Mh.njr(numCell,r).x;
      t.ten[1][0]=coP1*(Mh.xr(numCell,r).x-Mh.xj(numCell).x)*Mh.njr(numCell,r).y;
      t.ten[1][1]=coP1*(Mh.xr(numCell,r).y-Mh.xj(numCell).y)*Mh.njr(numCell,r).y;
      
    }

      if(ts==2){   

	/////////////// Cas quadrangle //////////////
	T=Vjr(Mh,numCell,numGr);
	
	t.ten[0][0]=T/Mh.ljr(numCell,r);
	t.ten[0][1]=0;
	t.ten[1][0]=0;
	t.ten[1][1]=T/Mh.ljr(numCell,r);      
      }
      if(ts==3){
	double a11=0,a12=0,a21=0,a22=0;
	int r,jG;
	for(int j=0;j<tab.TabInv[numGr].taille;j++){
	  jG=tab.TabInv[numGr].TabCell[j];   
	  r=NodeGtoL(Mh,jG,numGr);  

	  a11=a11+Mh.ljr(jG,r)*(Mh.xr(jG,r).x-Mh.xj(jG).x)*Mh.njr(jG,r).x;
	  a12=a12+Mh.ljr(jG,r)*(Mh.xr(jG,r).y-Mh.xj(jG).y)*Mh.njr(jG,r).x;
	  a21=a21+Mh.ljr(jG,r)*(Mh.xr(jG,r).x-Mh.xj(jG).x)*Mh.njr(jG,r).y;
	  a22=a22+Mh.ljr(jG,r)*(Mh.xr(jG,r).y-Mh.xj(jG).y)*Mh.njr(jG,r).y; 
	}
	
	 r=NodeGtoL(Mh,numCell,numGr);
	t.ten[0][0]=(Vjr(Mh,numCell,numGr)*a11)/(Vr(Mh,numGr,tab)*Mh.ljr(numCell,r));
	t.ten[0][1]=(Vjr(Mh,numCell,numGr)*a12)/(Vr(Mh,numGr,tab)*Mh.ljr(numCell,r));
	t.ten[1][0]=(Vjr(Mh,numCell,numGr)*a21)/(Vr(Mh,numGr,tab)*Mh.ljr(numCell,r));
	t.ten[1][1]=(Vjr(Mh,numCell,numGr)*a22)/(Vr(Mh,numGr,tab)*Mh.ljr(numCell,r));
	
      }//////// fin tenseur 4

  }

  if(c=='h'){
    if(th==1){
      t.ten[0][0]=Mh.njr(numCell,r).x*Mh.njr(numCell,r).x;
      t.ten[0][1]=Mh.njr(numCell,r).x*Mh.njr(numCell,r).y;
      t.ten[1][0]=Mh.njr(numCell,r).y*Mh.njr(numCell,r).x;
      t.ten[1][1]=Mh.njr(numCell,r).y*Mh.njr(numCell,r).y;
      
      
     }

     if(th==3){
      t.ten[0][0]=ld*Mh.njr(numCell,r).x*Mh.njr(numCell,r).x;
      t.ten[0][1]=ld*Mh.njr(numCell,r).x*Mh.njr(numCell,r).y;
      t.ten[1][0]=ld*Mh.njr(numCell,r).y*Mh.njr(numCell,r).x;
      t.ten[1][1]=ld*Mh.njr(numCell,r).y*Mh.njr(numCell,r).y;
      
      
     }

      if(th==2){
	xd=NextNodeLocal(Mh,r);
	xg=PreviousNodeLocal(Mh,r);
	// recuperation du numéro des arete	
	e1=xg;
	e2=r;

	t.ten[0][0]=(Mh.ljk(numCell,e1)*Mh.njk(numCell,e1).x*Mh.njk(numCell,e1).x+Mh.ljk(numCell,e2)*Mh.njk(numCell,e2).x*Mh.njk(numCell,e2).x)/(2*Mh.ljr(numCell,r));
	t.ten[0][1]=(Mh.ljk(numCell,e1)*Mh.njk(numCell,e1).x*Mh.njk(numCell,e1).y+Mh.ljk(numCell,e2)*Mh.njk(numCell,e2).x*Mh.njk(numCell,e2).y)/(2*Mh.ljr(numCell,r));
	t.ten[1][0]=(Mh.ljk(numCell,e1)*Mh.njk(numCell,e1).y*Mh.njk(numCell,e1).x+Mh.ljk(numCell,e2)*Mh.njk(numCell,e2).y*Mh.njk(numCell,e2).x)/(2*Mh.ljr(numCell,r));
	t.ten[1][1]=(Mh.ljk(numCell,e1)*Mh.njk(numCell,e1).y*Mh.njk(numCell,e1).y+Mh.ljk(numCell,e2)*Mh.njk(numCell,e2).y*Mh.njk(numCell,e2).y)/(2*Mh.ljr(numCell,r));

      }
      t.ten[0][0]=t.ten[0][0];
      t.ten[0][1]=t.ten[0][1];
      t.ten[1][0]=t.ten[1][0];
      t.ten[1][1]=t.ten[1][1];
  }
  return t;
}
