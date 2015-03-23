#include <iostream>
#include <cmath>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include "Functionsgeo.hpp"
#include "Data.hpp"

/** File with Geometric functions **/


R2 GravityCenter(cell cc,int nbnode){
  /** function which compute the gravity center of the cell "cc" **/
  R2 res;
  Vertex A,B,C,D;
  if(nbnode==4){
    A=cc.vertices[0];
    B=cc.vertices[1];
    C=cc.vertices[2];
    D=cc.vertices[3];

    res.x=(1/(6*cc.area))*((A.x+B.x)*(A^B)+(B.x+C.x)*(B^C)+(C.x+D.x)*(C^D)+(D.x+A.x)*(D^A));
    res.y=(1/(6*cc.area))*((A.y+B.y)*(A^B)+(B.y+C.y)*(B^C)+(C.y+D.y)*(C^D)+(D.y+A.y)*(D^A));
  }
  if(nbnode==3){
    A=cc.vertices[0];
    B=cc.vertices[1];
    C=cc.vertices[2];

    res.x=(1/(6*cc.area))*((A.x+B.x)*(A^B)+(B.x+C.x)*(B^C)+(C.x+A.x)*(C^A));                        
    res.y=(1/(6*cc.area))*((A.y+B.y)*(A^B)+(B.y+C.y)*(B^C)+(C.y+A.y)*(C^A));
  }
  return res;
}






R2 ChoiceCenterCell(Data & d,cell cc, int nbnode){
  R2 res;

     res=GravityCenter(cc,nbnode);
  return res;
}


int centerdomain(Data & d,Mesh & Mh){
  /** Function which compute the center of mesh for the gaussian function **/
  int in=0;
 if(Mh.nbnodelocal==4){
       in=(Mh.nc+d.Nx)/2;}
  if(Mh.nbnodelocal==3){
    in=(Mh.nc+2*d.Nx)/2;}
   if(!strcmp(d.NameMesh,"KershawMesh")){
	   for(int j=0;j<Mh.nc;j++){
	     if(Mh.cells[j].centercell.x<=1+d.Tx*(1./d.Nx) && Mh.cells[j].centercell.x>=1-d.Tx*(1./d.Nx)){
	       if(Mh.cells[j].centercell.y<=1+d.Ty*(1./d.Ny) && Mh.cells[j].centercell.y>=1-d.Ty*(1./d.Ny)){
		 in=j;
	       }
	     }
	   }
	 }

   return in;
}

double Vjr(Mesh & Mh,int j,int numGr){
  /** Fonction which compute the area at the sub control volume for a node r and the cell j**/
    double T;
    int r;
    R2 xj1(0,0);
    R2 xj2(0,0);
    
    r=NodeGtoL(Mh,j,numGr);
    
    if(r==Mh.nbnodelocal-1){
      xj2.x=(Mh.xr(j,r).x+Mh.xr(j,0).x)/2;
      xj2.y=(Mh.xr(j,r).y+Mh.xr(j,0).y)/2;
	   
    }
    else{
      xj2.x=(Mh.xr(j,r).x+Mh.xr(j,r+1).x)/2;
      xj2.y=(Mh.xr(j,r).y+Mh.xr(j,r+1).y)/2;
    }

    if(r==0){
      xj1.x=(Mh.xr(j,r).x+Mh.xr(j,Mh.nbnodelocal-1).x)/2;
      xj1.y=(Mh.xr(j,r).y+Mh.xr(j,Mh.nbnodelocal-1).y)/2;	    
    }
    else{
      xj1.x=(Mh.xr(j,r).x+Mh.xr(j,r-1).x)/2;
      xj1.y=(Mh.xr(j,r).y+Mh.xr(j,r-1).y)/2;}

    /////////////// Cas quadrangle ///////////////
    T=0.5*((Mh.xr(numGr)^xj2)+(xj2^Mh.xj(j))+(Mh.xj(j)^xj1)+(xj1^Mh.xr(numGr)));
    return T;
 }


double Vr(Mesh & Mh,int numGr,TabConnecInv & tab){
  /** Function which compute the area at the control volume for a node r**/
  int jG;
  double V=0;

for(int j=0;j<tab.TabInv[numGr].taille;j++){
      jG=tab.TabInv[numGr].TabCell[j];
      V=V+Vjr(Mh,jG,numGr);
 }
 return V;
}


cell ControlVolume(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,int numGr){
  /** Function which compute the control volum for node numGr**/
  int jG=0;
  cell CVol=0;
  int k=0;
  int r=0;
  R2 xj1(0,0);
  R2 xj2(0,0);
  R2 *verticesTemp;

  
  CVol.nbnodelocal=2*tab.TabInv[numGr].taille;
  verticesTemp = new R2[CVol.nbnodelocal];
  CVol.vertices =new Vertex[CVol.nbnodelocal];
  
  for(int j=0;j<tab.TabInv[numGr].taille;j++){
      jG=tab.TabInv[numGr].TabCell[j];

       r=NodeGtoL(Mh,jG,numGr);
      if(r==Mh.nbnodelocal-1){
      xj2.x=(Mh.xr(jG,r).x+Mh.xr(jG,0).x)/2;
      xj2.y=(Mh.xr(jG,r).y+Mh.xr(jG,0).y)/2;
	   
    }
    else{
      xj2.x=(Mh.xr(jG,r).x+Mh.xr(jG,r+1).x)/2;
      xj2.y=(Mh.xr(jG,r).y+Mh.xr(jG,r+1).y)/2;
    }

    if(r==0){
      xj1.x=(Mh.xr(jG,r).x+Mh.xr(jG,Mh.nbnodelocal-1).x)/2;
      xj1.y=(Mh.xr(jG,r).y+Mh.xr(j,Mh.nbnodelocal-1).y)/2;	    
    }
    else{
      xj1.x=(Mh.xr(jG,r).x+Mh.xr(jG,r-1).x)/2;
      xj1.y=(Mh.xr(jG,r).y+Mh.xr(jG,r-1).y)/2;
    }

  CVol.vertices[k].x=xj2.x;
  CVol.vertices[k].y=xj2.y;
  CVol.vertices[k+1].x=Mh.xj(jG).x;
  CVol.vertices[k+1].y=Mh.xj(jG).y;
  k=k+2;
}

   

  CVol.area=Vr(Mh,numGr,tab);
  CVol.lab=0;
  return CVol;
}
  


cell SubControlVolume(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,int numGr,int j){
  /** Function which compute the control volum for node numGr**/
  int jG=0;
  cell CVol=0;
  int r=0;
  R2 xj1(0,0);
  R2 xj2(0,0);
  R2 *verticesTemp;

  
  CVol.nbnodelocal=tab.TabInv[numGr].taille;
  verticesTemp = new R2[CVol.nbnodelocal];
  CVol.vertices =new Vertex[CVol.nbnodelocal];

  jG=tab.TabInv[numGr].TabCell[j];
  r=NodeGtoL(Mh,jG,numGr);

  if(j==0){
  
  CVol.vertices[0].x=Mh.xr(numGr).x;
  CVol.vertices[0].y=Mh.xr(numGr).y;

  r=NodeGtoL(Mh,jG,numGr);
  if(r==Mh.nbnodelocal-1){
    xj2.x=(Mh.xr(jG,r).x+Mh.xr(jG,0).x)/2;
    xj2.y=(Mh.xr(jG,r).y+Mh.xr(jG,0).y)/2;
    
    }
    else{
      xj2.x=(Mh.xr(jG,r).x+Mh.xr(jG,r+1).x)/2;
      xj2.y=(Mh.xr(jG,r).y+Mh.xr(jG,r+1).y)/2;
    }

  if(r==0){
      xj1.x=(Mh.xr(jG,r).x+Mh.xr(jG,Mh.nbnodelocal-1).x)/2;
      xj1.y=(Mh.xr(jG,r).y+Mh.xr(jG,Mh.nbnodelocal-1).y)/2;	    
  }
  else{
    xj1.x=(Mh.xr(jG,r).x+Mh.xr(jG,r-1).x)/2;
    xj1.y=(Mh.xr(jG,r).y+Mh.xr(jG,r-1).y)/2;
  }

  CVol.vertices[1].x=xj2.x;
  CVol.vertices[1].y=xj2.y;
  
  CVol.vertices[2].x=Mh.xj(jG).x;
  CVol.vertices[2].y=Mh.xj(jG).y;
  
  CVol.vertices[3].x=xj1.x;
  CVol.vertices[3].y=xj1.y;
 
  CVol.area=Vjr(Mh,jG,numGr);
   }

  if(j==1){
  
  CVol.vertices[1].x=Mh.xr(numGr).x;
  CVol.vertices[1].y=Mh.xr(numGr).y;

  r=NodeGtoL(Mh,jG,numGr);
  if(r==Mh.nbnodelocal-1){
    xj2.x=(Mh.xr(jG,r).x+Mh.xr(jG,0).x)/2;
    xj2.y=(Mh.xr(jG,r).y+Mh.xr(jG,0).y)/2;
    
    }
    else{
      xj2.x=(Mh.xr(jG,r).x+Mh.xr(jG,r+1).x)/2;
      xj2.y=(Mh.xr(jG,r).y+Mh.xr(jG,r+1).y)/2;
    }

  if(r==0){
      xj1.x=(Mh.xr(jG,r).x+Mh.xr(jG,Mh.nbnodelocal-1).x)/2;
      xj1.y=(Mh.xr(jG,r).y+Mh.xr(jG,Mh.nbnodelocal-1).y)/2;	    
  }
  else{
    xj1.x=(Mh.xr(jG,r).x+Mh.xr(jG,r-1).x)/2;
    xj1.y=(Mh.xr(jG,r).y+Mh.xr(jG,r-1).y)/2;
  }

  CVol.vertices[2].x=xj2.x;
  CVol.vertices[2].y=xj2.y;
  
  CVol.vertices[3].x=Mh.xj(jG).x;
  CVol.vertices[3].y=Mh.xj(jG).y;
  
  CVol.vertices[0].x=xj1.x;
  CVol.vertices[0].y=xj1.y;
 
  CVol.area=Vjr(Mh,jG,numGr);
   }

  if(j==2){
  
  CVol.vertices[2].x=Mh.xr(numGr).x;
  CVol.vertices[2].y=Mh.xr(numGr).y;

  r=NodeGtoL(Mh,jG,numGr);
  if(r==Mh.nbnodelocal-1){
    xj2.x=(Mh.xr(jG,r).x+Mh.xr(jG,0).x)/2;
    xj2.y=(Mh.xr(jG,r).y+Mh.xr(jG,0).y)/2;
    
    }
    else{
      xj2.x=(Mh.xr(jG,r).x+Mh.xr(jG,r+1).x)/2;
      xj2.y=(Mh.xr(jG,r).y+Mh.xr(jG,r+1).y)/2;
    }

  if(r==0){
      xj1.x=(Mh.xr(jG,r).x+Mh.xr(jG,Mh.nbnodelocal-1).x)/2;
      xj1.y=(Mh.xr(jG,r).y+Mh.xr(jG,Mh.nbnodelocal-1).y)/2;	    
  }
  else{
    xj1.x=(Mh.xr(jG,r).x+Mh.xr(jG,r-1).x)/2;
    xj1.y=(Mh.xr(jG,r).y+Mh.xr(jG,r-1).y)/2;
  }

  CVol.vertices[3].x=xj2.x;
  CVol.vertices[3].y=xj2.y;
  
  CVol.vertices[0].x=Mh.xj(jG).x;
  CVol.vertices[0].y=Mh.xj(jG).y;
  
  CVol.vertices[1].x=xj1.x;
  CVol.vertices[1].y=xj1.y;
 
  CVol.area=Vjr(Mh,jG,numGr);
   }

  if(j==3){
  
  CVol.vertices[3].x=Mh.xr(numGr).x;
  CVol.vertices[3].y=Mh.xr(numGr).y;

  r=NodeGtoL(Mh,jG,numGr);
  if(r==Mh.nbnodelocal-1){
    xj2.x=(Mh.xr(jG,r).x+Mh.xr(jG,0).x)/2;
    xj2.y=(Mh.xr(jG,r).y+Mh.xr(jG,0).y)/2;
    
    }
    else{
      xj2.x=(Mh.xr(jG,r).x+Mh.xr(jG,r+1).x)/2;
      xj2.y=(Mh.xr(jG,r).y+Mh.xr(jG,r+1).y)/2;
    }

  if(r==0){
      xj1.x=(Mh.xr(jG,r).x+Mh.xr(jG,Mh.nbnodelocal-1).x)/2;
      xj1.y=(Mh.xr(jG,r).y+Mh.xr(jG,Mh.nbnodelocal-1).y)/2;	    
  }
  else{
    xj1.x=(Mh.xr(jG,r).x+Mh.xr(jG,r-1).x)/2;
    xj1.y=(Mh.xr(jG,r).y+Mh.xr(jG,r-1).y)/2;
  }

  CVol.vertices[0].x=xj2.x;
  CVol.vertices[0].y=xj2.y;
  
  CVol.vertices[1].x=Mh.xj(jG).x;
  CVol.vertices[1].y=Mh.xj(jG).y;
  
  CVol.vertices[2].x=xj1.x;
  CVol.vertices[2].y=xj1.y;
 
  CVol.area=Vjr(Mh,jG,numGr);
   }
 
  CVol.lab=0;
  return CVol;
}


double AverageQuantity(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,ParamPhysic & Param,int numGr,int var){
  /** Function which compute an average of the quantity var(var) around the node numGr**/
  int nbmeans;
  int jG=0;
  double m=0;

  nbmeans=2;
  if(nbmeans==1){
    for(int j=0;j<tab.TabInv[numGr].taille;j++){
      jG=tab.TabInv[numGr].TabCell[j];
      m=m+Vjr(Mh,jG,numGr)*v.var[var][jG];
    }

    m=m/Vr(Mh,numGr,tab);
  }
  if(nbmeans==2){
    for(int j=0;j<tab.TabInv[numGr].taille;j++){
      jG=tab.TabInv[numGr].TabCell[j];
      m=m+v.var[var][jG];
    }
    m=m/tab.TabInv[numGr].taille;
  }

 
  
  return m;
}

double geomCorrection(Data & d,Mesh & Mh,int j){
  /** Geometric correction for the GLACE nodal tensor**/
  double res;
  double num=0;
  double a11=0,a12=0,a21=0,a22=0;
  double tr=0,Det=0,mu=0;
  double a=0;
  
  for(int r=0;r<Mh.nbnodelocal;r++){
    num=num+sqrt(pow(Mh.ljr(j,r)*Mh.njr(j,r).x,2.)+pow(Mh.ljr(j,r)*Mh.njr(j,r).y,2.));
    a11=a11+Mh.ljr(j,r)*Mh.njr(j,r).x*Mh.njr(j,r).x;
    a12=a12+Mh.ljr(j,r)*Mh.njr(j,r).x*Mh.njr(j,r).y;
    a21=a21+Mh.ljr(j,r)*Mh.njr(j,r).y*Mh.njr(j,r).x;
    a22=a22+Mh.ljr(j,r)*Mh.njr(j,r).y*Mh.njr(j,r).y;
 
  }
  tr=a11+a22;
 
   Det=a11*a22-a12*a21;
  if(tr*tr-4*Det<0){
    a=0;
  }
  else{
    a=tr*tr-4*Det;
  }
  mu=0.5*(tr+sqrt(a));
  
  res=sqrt(num/mu);
  return res;
}


void OneLayerStencilNodeQ(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,int *& Stencil,int NStencil, int numGr){
  /** function which compute the list of the cell directly around a node for a quandragular mesh**/
  
  if(NStencil != tab.TabInv[numGr].taille){cout<<"problem in the one layer node stencil for quadrangular mesh"<<endl; exit(0); }
 
  for(int j=0;j<tab.TabInv[numGr].taille;j++){
    Stencil[j]=tab.TabInv[numGr].TabCell[j];
  }


}  
void TwoLayerStencilNodeQ(Data & d,Mesh & Mh,variable & v,TabConnecInv & tab,int *& Stencil,int NStencil, int numGr){
  /** function which compute the list of the cell directly around a node for a quandragular mesh**/
  int * tabStencilCell=NULL;
  int NtabStencilCell=0;
  int k=0,jG=0,present=0;

  
  for(int i=0;i<NStencil;i++){
	Stencil[i]=-1;
  }
  
  
  if(NStencil != 4*tab.TabInv[numGr].taille){cout<<"problem in the two layer node stencil for quadrangular mesh"<<endl; exit(0); }
  
  for(int j=0;j<tab.TabInv[numGr].taille;j++){
    jG=tab.TabInv[numGr].TabCell[j];
    
    NtabStencilCell=tabVF9(Mh,tab,tabStencilCell,jG);
    
    
    for(int i=0;i<NtabStencilCell;i++){
      present=PresenceNodetoTab(Stencil,tabStencilCell[i],NStencil);
	if( present == 0 ){
	  Stencil[k]=tabStencilCell[i];
	  k++;
	}
    }
    
  }
  
  delete [] tabStencilCell;  
}
