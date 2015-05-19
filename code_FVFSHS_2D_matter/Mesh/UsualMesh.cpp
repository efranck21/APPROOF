
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "UsualMesh.hpp"
#include "time.h"
#include "math.h"
using namespace std;

/***  function which construct the mesh (using the mesh choose in DATA )***/

Mesh ChoiceMesh(Data & d){
  int nb=0;
  if(d.Typemesh=='T')
    nb=3;
     else if(d.Typemesh=='Q'){
       nb=4;
     }
 
  Mesh Mh(d);
  
  if(!strcmp(d.NameMesh,"CartesianMesh"))
      Mh=CartesianMesh(d);
  
  if(!strcmp(d.NameMesh,"RandomCartesianMesh"))
    Mh=RandomCartesianMesh(d);
  
  if(!strcmp(d.NameMesh,"SmoothCartesianMesh"))
    Mh=SmoothCartesianMesh(d);
  
  if(!strcmp(d.NameMesh,"KershawMesh"))
    Mh=KershawMesh(d);
  
  if(!strcmp(d.NameMesh,"CartesianMeshT"))
    Mh=CartesianMeshT(d);
  
  if(!strcmp(d.NameMesh,"RandomCartesianMeshT"))
    Mh=RandomCartesianMeshT(d);

   if(!strcmp(d.NameMesh,"DoubleSquareMeshQ"))
    Mh=DoubleSquareMeshQ(d);

  if(!strcmp(d.NameMesh,"CheckMesh"))
    Mh=CheckMesh(d);

  /** lab = 0: nodes in the domain **/
  /** lab = -1 : nodes out of the domain **/
  /** lab = 2 : nodes at the boundary left **/
  /** lab = 3 : nodes at the boundary top **/
  /** lab = 4 : nodes at the boundary right **/
  /** lab = 5 : nodes at the boundary bottom **/
  /** lab = 1 : nodes at the corner of the domain **/

  /** lab = 0: cells in the domain **/
  /** lab = -1 : ghost cells **/
  /** lab >0 : cells in a subdomain **/
 
  return Mh;
	
}


/***  Construction of Cartesian Mesh ***/
Mesh CartesianMesh(Data & d){

double Dx=0.,Dy=0.;
 Mesh Mh(d);

 
Dx=d.Tx/(double)d.Nx; // uniform step mesh in the x direction
Dy=d.Ty/(double)d.Ny;


 
Mh.nv=(d.Nx+3)*(d.Ny+3);
Mh.nc=(d.Nx+2)*(d.Ny+2);


 Vertex iniver(d);
Mh.vertices = new Vertex[Mh.nv];
 for(int i=0;i<Mh.nv;i++){
   Mh.vertices[i]=iniver;
 }

 cell ini(Mh.nbnodelocal);
 Mh.cells = new cell[Mh.nc];
 for(int i=0;i<Mh.nc;i++){
   Mh.cells[i]=ini;
 }
 
 
int k=0;
 for(int iy=0;iy<d.Ny+3;iy++){
   for(int ix=0;ix<d.Nx+3;ix++){
     Mh.vertices[k].x=(ix-1)*Dx;
     Mh.vertices[k].y=(iy-1)*Dy;

    
    if(ix==0) {Mh.vertices[k].lab=-1;}
    if(iy==d.Ny+2) {Mh.vertices[k].lab=-1;}
    if(ix==d.Nx+2) {Mh.vertices[k].lab=-1;}
    if(iy==0) {Mh.vertices[k].lab=-1;}
    if((ix==0 && iy==0) || (ix==d.Nx+2 && iy==0) || (iy==d.Ny+2 && ix==0) || (iy==d.Ny+2 && ix==d.Nx+2)){
      Mh.vertices[k].lab=-1;}

    if(ix==1 && iy>0 && iy<d.Ny+2) {Mh.vertices[k].lab=2;}
    if(iy==d.Ny+1 && ix>0 && ix<d.Nx+2) {Mh.vertices[k].lab=3;}
    if(ix==d.Nx+1 && iy>0 && iy<d.Ny+2) {Mh.vertices[k].lab=4;}
    if(iy==1 && ix>0 && ix<d.Nx+2) {Mh.vertices[k].lab=5;}
    if((ix==1 && iy==1) || (ix==d.Nx+1 && iy==1) || (iy==d.Ny+1 && ix==1) || (iy==d.Ny+1 && ix==d.Nx+1)){
      Mh.vertices[k].lab=1;}
    
    if(ix!=0 && iy!=0 && ix!=d.Nx+2 && iy!=d.Ny+2 && ix!=1 && iy!=1 && ix!=d.Nx+1 && iy!=d.Ny+1) {
      Mh.vertices[k].lab=0;

    }
     k++;
  }
}

 
k=0;
 for(int iy=0;iy<d.Ny+2;iy++){
  for(int ix=0;ix<d.Nx+2;ix++){
    if(Mh.vertices[(d.Nx+3)*iy+ix].lab>-1 && Mh.vertices[(d.Nx+3)*iy+ix+1].lab>-1 && Mh.vertices[(d.Nx+3)*(iy+1)+ix+1].lab>-1 && Mh.vertices[(d.Nx+3)*(iy+1)+ix].lab>-1){
     
     Mh.cells[k].init(Mh.vertices,(d.Nx+3)*iy+ix,(d.Nx+3)*iy+ix+1,(d.Nx+3)*(iy+1)+ix+1,(d.Nx+3)*(iy+1)+ix,0);
      }
      else{
      Mh.cells[k].init(Mh.vertices,(d.Nx+3)*iy+ix,(d.Nx+3)*iy+ix+1,(d.Nx+3)*(iy+1)+ix+1,(d.Nx+3)*(iy+1)+ix,-1);
      }	
    k++;	
  }
 }

return Mh;
}

/***  Construction of the Random Cartesian Mesh ***/
Mesh RandomCartesianMesh(Data & d){

double Dx=0.,Dy=0.;
 Mesh Mh(d);
  double theta=0;
Dx=d.Tx/(double)d.Nx; 
Dy=d.Ty/(double)d.Ny;

Mh.nv=(d.Nx+3)*(d.Ny+3);
Mh.nc=(d.Nx+2)*(d.Ny+2);


 Vertex iniver(d);
Mh.vertices = new Vertex[Mh.nv];
 for(int i=0;i<Mh.nv;i++){
   Mh.vertices[i]=iniver;
 }

 cell ini(Mh.nbnodelocal);
 Mh.cells = new cell[Mh.nc];
 for(int i=0;i<Mh.nc;i++){
   Mh.cells[i]=ini;
 }

 
srand(152);

int k=0;
 for(int iy=0;iy<d.Ny+3;iy++){
   for(int ix=0;ix<d.Nx+3;ix++){
     theta=(rand()/(float)RAND_MAX)*2.0-1.0;
     Mh.vertices[k].x=(ix-1)*Dx;
     Mh.vertices[k].y=(iy-1)*Dy;  

      if(ix==0) {Mh.vertices[k].lab=-1;}
    if(iy==d.Ny+2) {Mh.vertices[k].lab=-1;}
    if(ix==d.Nx+2) {Mh.vertices[k].lab=-1;}
    if(iy==0) {Mh.vertices[k].lab=-1;}
    if((ix==0 && iy==0) || (ix==d.Nx+2 && iy==0) || (iy==d.Ny+2 && ix==0) || (iy==d.Ny+2 && ix==d.Nx+2)){
      Mh.vertices[k].lab=-1;}

      if(ix==1 && iy>0 && iy<d.Ny+2) {Mh.vertices[k].lab=2;}
    if(iy==d.Ny+1 && ix>0 && ix<d.Nx+2) {Mh.vertices[k].lab=3;}
    if(ix==d.Nx+1 && iy>0 && iy<d.Ny+2) {Mh.vertices[k].lab=4;}
    if(iy==1 && ix>0 && ix<d.Nx+2) {Mh.vertices[k].lab=5;}
    if((ix==1 && iy==1) || (ix==d.Nx+1 && iy==1) || (iy==d.Ny+1 && ix==1) || (iy==d.Ny+1 && ix==d.Nx+1)){
      Mh.vertices[k].lab=1;}
  
    if(ix!=0 && iy!=0 && ix!=d.Nx+2 && iy!=d.Ny+2 && ix!=1 && iy!=1 && ix!=d.Nx+1 && iy!=d.Ny+1) {
      //  if(ix!=0 && iy!=0 && ix!=d.Nx+2 && iy!=d.Ny+2) {
      Mh.vertices[k].x=(ix-1)*Dx+theta*(Dx/5);
      Mh.vertices[k].y=(iy-1)*Dy+theta*(Dy/5);
    }
    if(ix!=0 && iy!=0 && ix!=d.Nx+2 && iy!=d.Ny+2 && ix!=1 && iy!=1 && ix!=d.Nx+1 && iy!=d.Ny+1){
          Mh.vertices[k].lab=0;
    }
     k++;
  }
}

 
k=0;
 for(int iy=0;iy<d.Ny+2;iy++){
  for(int ix=0;ix<d.Nx+2;ix++){
    if(Mh.vertices[(d.Nx+3)*iy+ix].lab>-1 && Mh.vertices[(d.Nx+3)*iy+ix+1].lab>-1 && Mh.vertices[(d.Nx+3)*(iy+1)+ix+1].lab>-1 && Mh.vertices[(d.Nx+3)*(iy+1)+ix].lab>-1){
     Mh.cells[k].init(Mh.vertices,(d.Nx+3)*iy+ix,(d.Nx+3)*iy+ix+1,(d.Nx+3)*(iy+1)+ix+1,(d.Nx+3)*(iy+1)+ix,0);
    } 
    else{
      Mh.cells[k].init(Mh.vertices,(d.Nx+3)*iy+ix,(d.Nx+3)*iy+ix+1,(d.Nx+3)*(iy+1)+ix+1,(d.Nx+3)*(iy+1)+ix,-1);
    }	
    k++;	
  }
 }

return Mh;
}

/***  Construction of Collela Mesh (smooth deformation of a Cartesian mesh) ***/
Mesh SmoothCartesianMesh(Data & d){

double Dx=0.,Dy=0.;
 Mesh Mh(d);
 double a=0.1;
Dx=d.Tx/(double)d.Nx; 
Dy=d.Ty/(double)d.Ny;

Mh.nv=(d.Nx+3)*(d.Ny+3);
Mh.nc=(d.Nx+2)*(d.Ny+2);

 Vertex iniver(d);
Mh.vertices = new Vertex[Mh.nv];
 for(int i=0;i<Mh.nv;i++){
   Mh.vertices[i]=iniver;
 }

 cell ini(Mh.nbnodelocal);
 Mh.cells = new cell[Mh.nc];
 for(int i=0;i<Mh.nc;i++){
   Mh.cells[i]=ini;
 }

 
srand(52);

int k=0;
 for(int iy=0;iy<d.Ny+3;iy++){
   for(int ix=0;ix<d.Nx+3;ix++){
     Mh.vertices[k].x=(ix-1)*Dx;
     Mh.vertices[k].y=(iy-1)*Dy;  

     if(ix==0) {Mh.vertices[k].lab=-1;}
     if(iy==d.Ny+2) {Mh.vertices[k].lab=-1;}
     if(ix==d.Nx+2) {Mh.vertices[k].lab=-1;}
     if(iy==0) {Mh.vertices[k].lab=-1;}
     if((ix==0 && iy==0) || (ix==d.Nx+2 && iy==0) || (iy==d.Ny+2 && ix==0) || (iy==d.Ny+2 && ix==d.Nx+2)){
       Mh.vertices[k].lab=-1;}
     
       if(ix==1 && iy>0 && iy<d.Ny+2) {Mh.vertices[k].lab=2;}
    if(iy==d.Ny+1 && ix>0 && ix<d.Nx+2) {Mh.vertices[k].lab=3;}
    if(ix==d.Nx+1 && iy>0 && iy<d.Ny+2) {Mh.vertices[k].lab=4;}
    if(iy==1 && ix>0 && ix<d.Nx+2) {Mh.vertices[k].lab=5;}
    if((ix==1 && iy==1) || (ix==d.Nx+1 && iy==1) || (iy==d.Ny+1 && ix==1) || (iy==d.Ny+1 && ix==d.Nx+1)){
      Mh.vertices[k].lab=1;}
    
    if(ix!=0 && iy!=0 && ix!=d.Nx+2 && iy!=d.Ny+2 && ix!=1 && iy!=1 && ix!=d.Nx+1 && iy!=d.Ny+1) {
      Mh.vertices[k].x=(ix-1)*Dx+a*sin(2*M_PI*((ix-1)*Dx))*sin(2*M_PI*((iy-1)*Dy));
      Mh.vertices[k].y=(iy-1)*Dy+a*sin(2*M_PI*((ix-1)*Dx))*sin(2*M_PI*((iy-1)*Dy));;
    }
    if(ix!=0 && iy!=0 && ix!=d.Nx+2 && iy!=d.Ny+2 && ix!=1 && iy!=1 && ix!=d.Nx+1 && iy!=d.Ny+1){
          Mh.vertices[k].lab=0;
    }
     k++;
  }
}

 
k=0;
 for(int iy=0;iy<d.Ny+2;iy++){
  for(int ix=0;ix<d.Nx+2;ix++){
    if(Mh.vertices[(d.Nx+3)*iy+ix].lab>-1 && Mh.vertices[(d.Nx+3)*iy+ix+1].lab>-1 && Mh.vertices[(d.Nx+3)*(iy+1)+ix+1].lab>-1 && Mh.vertices[(d.Nx+3)*(iy+1)+ix].lab>-1){
     Mh.cells[k].init(Mh.vertices,(d.Nx+3)*iy+ix,(d.Nx+3)*iy+ix+1,(d.Nx+3)*(iy+1)+ix+1,(d.Nx+3)*(iy+1)+ix,0);
    } 
    else{
      Mh.cells[k].init(Mh.vertices,(d.Nx+3)*iy+ix,(d.Nx+3)*iy+ix+1,(d.Nx+3)*(iy+1)+ix+1,(d.Nx+3)*(iy+1)+ix,-1);
    }	
    k++;	
  }
 }

return Mh;
}



/***  Construction of the Kershaw Mesh (strong and nonsmooth deformation of the Cartesian mesh) ***/
Mesh KershawMesh(Data & d){
double Dx=0.,Dy=0.;
 Mesh Mh(d);
 double Ks=1.0;
 
Dx=d.Tx/(double)d.Nx; 
Dy=d.Ty/(double)d.Ny;

Mh.nv=(d.Nx+3)*(d.Ny+3);
Mh.nc=(d.Nx+2)*(d.Ny+2);


 Vertex iniver(d);
Mh.vertices = new Vertex[Mh.nv];
 for(int i=0;i<Mh.nv;i++){
   Mh.vertices[i]=iniver;
 }


 cell ini(Mh.nbnodelocal);
 Mh.cells = new cell[Mh.nc];
 for(int i=0;i<Mh.nc;i++){
   Mh.cells[i]=ini;
 }

 
int k=0;
 for(int iy=0;iy<d.Ny+3;iy++){
   for(int ix=0;ix<d.Nx+3;ix++){
     Mh.vertices[k].x=(ix-1)*Dx;
     Mh.vertices[k].y=(iy-1)*Dy;  

      if(ix==0) {Mh.vertices[k].lab=-1;}
     if(iy==d.Ny+2) {Mh.vertices[k].lab=-1;}
     if(ix==d.Nx+2) {Mh.vertices[k].lab=-1;}
     if(iy==0) {Mh.vertices[k].lab=-1;}
     if((ix==0 && iy==0) || (ix==d.Nx+2 && iy==0) || (iy==d.Ny+2 && ix==0) || (iy==d.Ny+2 && ix==d.Nx+2)){
       Mh.vertices[k].lab=-1;}
     
       if(ix==1 && iy>0 && iy<d.Ny+2) {Mh.vertices[k].lab=2;}
    if(iy==d.Ny+1 && ix>0 && ix<d.Nx+2) {Mh.vertices[k].lab=3;}
    if(ix==d.Nx+1 && iy>0 && iy<d.Ny+2) {Mh.vertices[k].lab=4;}
    if(iy==1 && ix>0 && ix<d.Nx+2) {Mh.vertices[k].lab=5;}
    if((ix==1 && iy==1) || (ix==d.Nx+1 && iy==1) || (iy==d.Ny+1 && ix==1) || (iy==d.Ny+1 && ix==d.Nx+1)){
      Mh.vertices[k].lab=1;}

     
    if(ix!=0 && iy!=0 && ix!=d.Nx+2 && iy!=d.Ny+2 && ix!=1 && iy!=1 && ix!=d.Nx+1 && iy!=d.Ny+1)  {Mh.vertices[k].lab=0;}
    

    if(Mh.vertices[k].x>=0 && Mh.vertices[k].x<0.5*d.Tx)
       {
         if (Mh.vertices[k].y>-1 && Mh.vertices[k].y<=0.25*d.Ty)
           {
             Mh.vertices[k].x+=(Mh.vertices[k].y-0.125*d.Ty)*Ks*Mh.vertices[k].x/(0.25*d.Tx);
           }

         if (Mh.vertices[k].y>0.25*d.Ty && Mh.vertices[k].y<=0.5*d.Ty)
           {
             Mh.vertices[k].x-=(Mh.vertices[k].y-0.125*d.Ty-0.25*d.Ty)*Ks*Mh.vertices[k].x/(0.25*d.Tx);
           }

         if (Mh.vertices[k].y>0.5*d.Ty && Mh.vertices[k].y<=0.75*d.Ty)
           {
             Mh.vertices[k].x+=(Mh.vertices[k].y-0.125*d.Ty-0.50*d.Ty)*Ks*Mh.vertices[k].x/(0.25*d.Tx);
           }

         if (Mh.vertices[k].y>0.75*d.Ty) //&& Mh.vertices[k].y<1*d.Ty)
           {
	     Mh.vertices[k].x-=(Mh.vertices[k].y-0.125*d.Ty-0.75*d.Ty)*Ks*Mh.vertices[k].x/(0.25*d.Tx);
           }
       }
     else
       {
         if (Mh.vertices[k].y>-1 && Mh.vertices[k].x>=0 && Mh.vertices[k].x<=d.Tx && Mh.vertices[k].y<=0.25*d.Ty)
           {
             Mh.vertices[k].x+=(Mh.vertices[k].y-0.125*d.Ty)*Ks*(d.Tx*1.0-Mh.vertices[k].x)/(0.25*d.Tx);
           }

         if (Mh.vertices[k].x>=0 && Mh.vertices[k].x<=d.Tx && Mh.vertices[k].y>0.25*d.Ty && Mh.vertices[k].y<=0.5*d.Ty)
           {
             Mh.vertices[k].x-=(Mh.vertices[k].y-0.125*d.Ty-0.25*d.Ty)*Ks*(d.Tx*1.0-Mh.vertices[k].x)/(0.25*d.Tx);
           }

         if (Mh.vertices[k].x>=0 && Mh.vertices[k].x<=d.Tx && Mh.vertices[k].y>d.Ty*0.5 && Mh.vertices[k].y<=d.Ty*0.75)
           {
             Mh.vertices[k].x+=(Mh.vertices[k].y-d.Ty*0.125-d.Ty*0.50)*Ks*(d.Tx*1.0-Mh.vertices[k].x)/(0.25*d.Tx);
           }

         if (Mh.vertices[k].x>=0 && Mh.vertices[k].x<=d.Tx && Mh.vertices[k].y>d.Ty*0.75)// && Mh.vertices[k].y<d.Ty*1)
           {
	     Mh.vertices[k].x-=(Mh.vertices[k].y-d.Ty*0.125-d.Ty*0.75)*Ks*(d.Tx*1.0-Mh.vertices[k].x)/(0.25*d.Tx);
           }
       }
    k++;

   } 
 }


k=0;
 for(int iy=0;iy<d.Ny+2;iy++){
  for(int ix=0;ix<d.Nx+2;ix++){
    if(Mh.vertices[(d.Nx+3)*iy+ix].lab>-1 && Mh.vertices[(d.Nx+3)*iy+ix+1].lab>-1 && Mh.vertices[(d.Nx+3)*(iy+1)+ix+1].lab>-1 && Mh.vertices[(d.Nx+3)*(iy+1)+ix].lab>-1){
     Mh.cells[k].init(Mh.vertices,(d.Nx+3)*iy+ix,(d.Nx+3)*iy+ix+1,(d.Nx+3)*(iy+1)+ix+1,(d.Nx+3)*(iy+1)+ix,0);
    } 
    else{
      Mh.cells[k].init(Mh.vertices,(d.Nx+3)*iy+ix,(d.Nx+3)*iy+ix+1,(d.Nx+3)*(iy+1)+ix+1,(d.Nx+3)*(iy+1)+ix,-1);
    }	
    k++;	
  }
 }

return Mh;
}




/***  Construction of the regular triangular Mesh ***/
Mesh CartesianMeshT(Data & d){
double Dx=0.,Dy=0.;
 Mesh Mh(d);
 int lab=0,n1=0,n2=0,n3=0;
 
Dx=d.Tx/(double)d.Nx; 
Dy=d.Ty/(double)d.Ny;

Mh.nv=(d.Nx+3)*(d.Ny+3);
 Mh.nc=2*((d.Nx+2)*(d.Ny+2));

 Vertex iniver(d);
Mh.vertices = new Vertex[Mh.nv];
 for(int i=0;i<Mh.nv;i++){
   Mh.vertices[i]=iniver;
 }

 cell ini(Mh.nbnodelocal);
 Mh.cells = new cell[Mh.nc];
 for(int i=0;i<Mh.nc;i++){
   Mh.cells[i]=ini;
 }
 

int k=0;
 for(int iy=0;iy<d.Ny+3;iy++){
  for(int ix=0;ix<d.Nx+3;ix++){
      Mh.vertices[k].x=(ix-1)*Dx;
     Mh.vertices[k].y=(iy-1)*Dy;  

      if(ix==0) {Mh.vertices[k].lab=-1;}
     if(iy==d.Ny+2) {Mh.vertices[k].lab=-1;}
     if(ix==d.Nx+2) {Mh.vertices[k].lab=-1;}
     if(iy==0) {Mh.vertices[k].lab=-1;}
     if((ix==0 && iy==0) || (ix==d.Nx+2 && iy==0) || (iy==d.Ny+2 && ix==0) || (iy==d.Ny+2 && ix==d.Nx+2)){
       Mh.vertices[k].lab=-1;}
     
       if(ix==1 && iy>0 && iy<d.Ny+2) {Mh.vertices[k].lab=2;}
    if(iy==d.Ny+1 && ix>0 && ix<d.Nx+2) {Mh.vertices[k].lab=3;}
    if(ix==d.Nx+1 && iy>0 && iy<d.Ny+2) {Mh.vertices[k].lab=4;}
    if(iy==1 && ix>0 && ix<d.Nx+2) {Mh.vertices[k].lab=5;}
    if((ix==1 && iy==1) || (ix==d.Nx+1 && iy==1) || (iy==d.Ny+1 && ix==1) || (iy==d.Ny+1 && ix==d.Nx+1)){
      Mh.vertices[k].lab=1;}
     
     if(ix!=0 && iy!=0 && ix!=d.Nx+2 && iy!=d.Ny+2 && ix!=1 && iy!=1 && ix!=d.Nx+1 && iy!=d.Ny+1)  {Mh.vertices[k].lab=0;}
     
     k++;
  }
 }


 
k=0;
 for(int iy=0;iy<d.Ny+2;iy++){
  for(int ix=0;ix<d.Nx+2;ix++){
    
    if((iy % 2)==0){
      if((ix % 2)==0){
	n1=(d.Nx+3)*iy+ix;
	n2=(d.Nx+3)*iy+ix+1;
	n3=(d.Nx+3)*(iy+1)+ix;
	if(Mh.vertices[n1].lab>-1 && Mh.vertices[n2].lab>-1 && Mh.vertices[n3].lab>-1){
	  lab=0;
	}
	else{
	  lab=-1;
	}
	Mh.cells[k].init(Mh.vertices,n1,n2,n3,0,lab);
	k++;

	n1=(d.Nx+3)*iy+ix+1;
	n2=(d.Nx+3)*(iy+1)+ix+1;
	n3=(d.Nx+3)*(iy+1)+ix;
	if(Mh.vertices[n1].lab>-1 && Mh.vertices[n2].lab>-1 && Mh.vertices[n3].lab>-1){
	  lab=0;
	}
	else{
	  lab=-1;
	}
	Mh.cells[k].init(Mh.vertices,n1,n2,n3,0,lab);
	k++;}
	 
      else {

	n1=(d.Nx+3)*iy+ix;
	n2=(d.Nx+3)*(iy+1)+ix+1;
	n3=(d.Nx+3)*(iy+1)+ix;
	if(Mh.vertices[n1].lab>-1 && Mh.vertices[n2].lab>-1 && Mh.vertices[n3].lab>-1){
	  lab=0;
	}
	else{
	  lab=-1;
	}
	Mh.cells[k].init(Mh.vertices,n1,n2,n3,0,lab);
	k++;

	n1=(d.Nx+3)*iy+ix;
	n2=(d.Nx+3)*iy+ix+1;
	n3=(d.Nx+3)*(iy+1)+ix+1;
	if(Mh.vertices[n1].lab>-1 && Mh.vertices[n2].lab>-1 && Mh.vertices[n3].lab>-1){
	  lab=0;
	}
	else{
	  lab=-1;
	}
	Mh.cells[k].init(Mh.vertices,n1,n2,n3,0,lab);
	k++;
      }
    }
    else {
      if((ix % 2)==0){

	n1=(d.Nx+3)*iy+ix;
	n2=(d.Nx+3)*(iy+1)+ix+1;
	n3=(d.Nx+3)*(iy+1)+ix;
	if(Mh.vertices[n1].lab>-1 && Mh.vertices[n2].lab>-1 && Mh.vertices[n3].lab>-1){
	  lab=0;
	}
	else{
	  lab=-1;
	}
	Mh.cells[k].init(Mh.vertices,n1,n2,n3,0,lab);
	k++;

	n1=(d.Nx+3)*iy+ix;
	n2=(d.Nx+3)*iy+ix+1;
	n3=(d.Nx+3)*(iy+1)+ix+1;
	if(Mh.vertices[n1].lab>-1 && Mh.vertices[n2].lab>-1 && Mh.vertices[n3].lab>-1){
	  lab=0;
	}
	else{
	  lab=-1;
	}
	Mh.cells[k].init(Mh.vertices,n1,n2,n3,0,lab);
	k++;
	
      }
      else {
	n1=(d.Nx+3)*iy+ix;
	n2=(d.Nx+3)*iy+ix+1;
	n3=(d.Nx+3)*(iy+1)+ix;
	if(Mh.vertices[n1].lab>-1 && Mh.vertices[n2].lab>-1 && Mh.vertices[n3].lab>-1){
	  lab=0;
	}
	else{
	  lab=-1;
	}
	Mh.cells[k].init(Mh.vertices,n1,n2,n3,0,lab);
	k++;

	n1=(d.Nx+3)*iy+ix+1;
	n2=(d.Nx+3)*(iy+1)+ix+1;
	n3=(d.Nx+3)*(iy+1)+ix;
	if(Mh.vertices[n1].lab>-1 && Mh.vertices[n2].lab>-1 && Mh.vertices[n3].lab>-1){
	  lab=0;
	}
	else{
	  lab=-1;
	}
	Mh.cells[k].init(Mh.vertices,n1,n2,n3,0,lab);
	k++;
	
      }
    }  
    
  }
 }

return Mh;
}

/***  Construction of random triangular Mesh ***/
Mesh RandomCartesianMeshT(Data & d){
  double Dx=0.,Dy=0.,theta=0.;
 Mesh Mh(d);
 int lab=0,n1=0,n2=0,n3=0;
 
Dx=d.Tx/(double)d.Nx; 
Dy=d.Ty/(double)d.Ny;
 
Mh.nv=(d.Nx+3)*(d.Ny+3);
 Mh.nc=2*((d.Nx+2)*(d.Ny+2));

 Vertex iniver(d);
Mh.vertices = new Vertex[Mh.nv];
 for(int i=0;i<Mh.nv;i++){
   Mh.vertices[i]=iniver;
 }

 cell ini(Mh.nbnodelocal);
 Mh.cells = new cell[Mh.nc];
 for(int i=0;i<Mh.nc;i++){
   Mh.cells[i]=ini;
 }
 

int k=0;
 for(int iy=0;iy<d.Ny+3;iy++){
  for(int ix=0;ix<d.Nx+3;ix++){
    theta=(rand()/(float)RAND_MAX)*2.0-1.0;
    Mh.vertices[k].x=(ix-1)*Dx;
    Mh.vertices[k].y=(iy-1)*Dy;  


    if(ix==0) {Mh.vertices[k].lab=-1;}
     if(iy==d.Ny+2) {Mh.vertices[k].lab=-1;}
     if(ix==d.Nx+2) {Mh.vertices[k].lab=-1;}
     if(iy==0) {Mh.vertices[k].lab=-1;}
     if((ix==0 && iy==0) || (ix==d.Nx+2 && iy==0) || (iy==d.Ny+2 && ix==0) || (iy==d.Ny+2 && ix==d.Nx+2)){
       Mh.vertices[k].lab=-1;}
     
       if(ix==1 && iy>0 && iy<d.Ny+2) {Mh.vertices[k].lab=2;}
    if(iy==d.Ny+1 && ix>0 && ix<d.Nx+2) {Mh.vertices[k].lab=3;}
    if(ix==d.Nx+1 && iy>0 && iy<d.Ny+2) {Mh.vertices[k].lab=4;}
    if(iy==1 && ix>0 && ix<d.Nx+2) {Mh.vertices[k].lab=5;}
    if((ix==1 && iy==1) || (ix==d.Nx+1 && iy==1) || (iy==d.Ny+1 && ix==1) || (iy==d.Ny+1 && ix==d.Nx+1)){
      Mh.vertices[k].lab=1;}

    if(ix!=0 && iy!=0 && ix!=d.Nx+2 && iy!=d.Ny+2 && ix!=1 && iy!=1 && ix!=d.Nx+1 && iy!=d.Ny+1) {
      Mh.vertices[k].x=(ix-1)*Dx+theta*(Dx/5);
      Mh.vertices[k].y=(iy-1)*Dy+theta*(Dy/5);
    }
     
     if(ix!=0 && iy!=0 && ix!=d.Nx+2 && iy!=d.Ny+2 && ix!=1 && iy!=1 && ix!=d.Nx+1 && iy!=d.Ny+1)  {Mh.vertices[k].lab=0;}
     
     k++;
  }
 }


 
k=0;
 for(int iy=0;iy<d.Ny+2;iy++){
  for(int ix=0;ix<d.Nx+2;ix++){
    
    if((iy % 2)==0){
      if((ix % 2)==0){
	n1=(d.Nx+3)*iy+ix;
	n2=(d.Nx+3)*iy+ix+1;
	n3=(d.Nx+3)*(iy+1)+ix;
	if(Mh.vertices[n1].lab>-1 && Mh.vertices[n2].lab>-1 && Mh.vertices[n3].lab>-1){
	  lab=0;
	}
	else{
	  lab=-1;
	}
	Mh.cells[k].init(Mh.vertices,n1,n2,n3,0,lab);
	k++;

	n1=(d.Nx+3)*iy+ix+1;
	n2=(d.Nx+3)*(iy+1)+ix+1;
	n3=(d.Nx+3)*(iy+1)+ix;
	if(Mh.vertices[n1].lab>-1 && Mh.vertices[n2].lab>-1 && Mh.vertices[n3].lab>-1){
	  lab=0;
	}
	else{
	  lab=-1;
	}
	Mh.cells[k].init(Mh.vertices,n1,n2,n3,0,lab);
	k++;}
	 
      else {

	n1=(d.Nx+3)*iy+ix;
	n2=(d.Nx+3)*(iy+1)+ix+1;
	n3=(d.Nx+3)*(iy+1)+ix;
	if(Mh.vertices[n1].lab>-1 && Mh.vertices[n2].lab>-1 && Mh.vertices[n3].lab>-1){
	  lab=0;
	}
	else{
	  lab=-1;
	}
	Mh.cells[k].init(Mh.vertices,n1,n2,n3,0,lab);
	k++;

	n1=(d.Nx+3)*iy+ix;
	n2=(d.Nx+3)*iy+ix+1;
	n3=(d.Nx+3)*(iy+1)+ix+1;
	if(Mh.vertices[n1].lab>-1 && Mh.vertices[n2].lab>-1 && Mh.vertices[n3].lab>-1){
	  lab=0;
	}
	else{
	  lab=-1;
	}
	Mh.cells[k].init(Mh.vertices,n1,n2,n3,0,lab);
	k++;
      }
    }
    else {
      if((ix % 2)==0){

	n1=(d.Nx+3)*iy+ix;
	n2=(d.Nx+3)*(iy+1)+ix+1;
	n3=(d.Nx+3)*(iy+1)+ix;
	if(Mh.vertices[n1].lab>-1 && Mh.vertices[n2].lab>-1 && Mh.vertices[n3].lab>-1){
	  lab=0;
	}
	else{
	  lab=-1;
	}
	Mh.cells[k].init(Mh.vertices,n1,n2,n3,0,lab);
	k++;

	n1=(d.Nx+3)*iy+ix;
	n2=(d.Nx+3)*iy+ix+1;
	n3=(d.Nx+3)*(iy+1)+ix+1;
	if(Mh.vertices[n1].lab>-1 && Mh.vertices[n2].lab>-1 && Mh.vertices[n3].lab>-1){
	  lab=0;
	}
	else{
	  lab=-1;
	}
	Mh.cells[k].init(Mh.vertices,n1,n2,n3,0,lab);
	k++;
	
      }
      else {
	n1=(d.Nx+3)*iy+ix;
	n2=(d.Nx+3)*iy+ix+1;
	n3=(d.Nx+3)*(iy+1)+ix;
	if(Mh.vertices[n1].lab>-1 && Mh.vertices[n2].lab>-1 && Mh.vertices[n3].lab>-1){
	  lab=0;
	}
	else{
	  lab=-1;
	}
	Mh.cells[k].init(Mh.vertices,n1,n2,n3,0,lab);
	k++;

	n1=(d.Nx+3)*iy+ix+1;
	n2=(d.Nx+3)*(iy+1)+ix+1;
	n3=(d.Nx+3)*(iy+1)+ix;
	if(Mh.vertices[n1].lab>-1 && Mh.vertices[n2].lab>-1 && Mh.vertices[n3].lab>-1){
	  lab=0;
	}
	else{
	  lab=-1;
	}
	Mh.cells[k].init(Mh.vertices,n1,n2,n3,0,lab);
	k++;
	
      }
    }  
    
  }
 }

return Mh;
}



Mesh DoubleSquareMeshQ(Data & d){
// Nx et Ny le nombre de maille dans les directions x et y
 double theta=0;
 int label=0;
 int indic=1;
 if( d.Tx!=1 || d.Ty!=1){
   cout <<" la taille du carré doit etre 1"<<endl;exit(1);
 }


 double Dx=0.,Dy=0.;
 Mesh Mh(d);

 
Dx=d.Tx/(double)d.Nx; // uniform step mesh in the x direction
Dy=d.Ty/(double)d.Ny;

Mh.nv=(d.Nx+3)*(d.Ny+3);
Mh.nc=(d.Nx+2)*(d.Ny+2);


 Vertex iniver(d);
Mh.vertices = new Vertex[Mh.nv];
 for(int i=0;i<Mh.nv;i++){
   Mh.vertices[i]=iniver;
 }

 cell ini(Mh.nbnodelocal);
 Mh.cells = new cell[Mh.nc];
 for(int i=0;i<Mh.nc;i++){
   Mh.cells[i]=ini;
 }
 
 srand(52);
int k=0;
 for(int iy=0;iy<d.Ny+3;iy++){
   for(int ix=0;ix<d.Nx+3;ix++){
      theta=(rand()/(float)RAND_MAX)*2.0-1.0;
     Mh.vertices[k].x=(ix-1)*Dx;
     Mh.vertices[k].y=(iy-1)*Dy;
     indic=1;
    
    if(ix==0) {Mh.vertices[k].lab=-1;}
    if(iy==d.Ny+2) {Mh.vertices[k].lab=-1;}
    if(ix==d.Nx+2) {Mh.vertices[k].lab=-1;}
    if(iy==0) {Mh.vertices[k].lab=-1;}
    if((ix==0 && iy==0) || (ix==d.Nx+2 && iy==0) || (iy==d.Ny+2 && ix==0) || (iy==d.Ny+2 && ix==d.Nx+2)){
      Mh.vertices[k].lab=-1;}

    if(ix==1 && iy>0 && iy<d.Ny+2) {Mh.vertices[k].lab=2;}
    if(iy==d.Ny+1 && ix>0 && ix<d.Nx+2) {Mh.vertices[k].lab=3;}
    if(ix==d.Nx+1 && iy>0 && iy<d.Ny+2) {Mh.vertices[k].lab=4;}
    if(iy==1 && ix>0 && ix<d.Nx+2) {Mh.vertices[k].lab=5;}
    if((ix==1 && iy==1) || (ix==d.Nx+1 && iy==1) || (iy==d.Ny+1 && ix==1) || (iy==d.Ny+1 && ix==d.Nx+1)){
      Mh.vertices[k].lab=1;}
     
    if(ix!=0 && iy!=0 && ix!=d.Nx+2 && iy!=d.Ny+2 && ix!=1 && iy!=1 && ix!=d.Nx+1 && iy!=d.Ny+1) {
      Mh.vertices[k].lab=0;
      ///////////////: construction carré  1//////////////
      if(Mh.vertices[k].x<=(3./16.) && Mh.vertices[k].x+Dx > (3./16.) && Mh.vertices[k].y>=(9./16.) && Mh.vertices[k].y<=(13./16.)){
	Mh.vertices[k].x=(3./16.);
	Mh.vertices[k].lab=6;
	indic=0;
      }
      if(Mh.vertices[k].y<=(9./16.) && Mh.vertices[k].y+Dy > (9./16.) && Mh.vertices[k].x>=(3./16.) && Mh.vertices[k].x<=(7./16.)){
	Mh.vertices[k].y=(9./16.);
	Mh.vertices[k].lab=6;
	indic=0;
      }
      if(Mh.vertices[k].x>=(7./16.) && Mh.vertices[k].x-Dx < (7./16.) && Mh.vertices[k].y>=(9./16.) && Mh.vertices[k].y<=(13./16.)){
	Mh.vertices[k].x=(7./16.);
	Mh.vertices[k].lab=6;
	indic=0;
      }
      if(Mh.vertices[k].y>=(13./16.) && Mh.vertices[k].y-Dy < (13./16.) && Mh.vertices[k].x>=(3./16.) && Mh.vertices[k].x<=(7./16.)){
	Mh.vertices[k].y=(13./16.);
	Mh.vertices[k].lab=6;
	indic=0;
      }
      ////////////////////////////////////////////// coin ///////////////////////////////////:
      if(Mh.vertices[k].x<=(3./16.) && Mh.vertices[k].x+Dx > (3./16.) && Mh.vertices[k].y<=(9./16.) && Mh.vertices[k].y+Dy>(9./16.)){
	Mh.vertices[k].x=(3./16.);
	Mh.vertices[k].y=(9./16.);
	Mh.vertices[k].lab=6;
	indic=0;
      }
      if(Mh.vertices[k].x<=(3./16.) && Mh.vertices[k].x+Dx > (3./16.) && Mh.vertices[k].y>=(13./16.) && Mh.vertices[k].y-Dy<(13./16.)){
	Mh.vertices[k].x=(3./16.);
	Mh.vertices[k].y=(13./16.);
	Mh.vertices[k].lab=6;
	indic=0;
      }
       if(Mh.vertices[k].x>=(7./16.) && Mh.vertices[k].x-Dx < (7./16.) && Mh.vertices[k].y<=(9./16.) && Mh.vertices[k].y+Dy>(9./16.)){
	Mh.vertices[k].x=(7./16.);
	Mh.vertices[k].y=(9./16.);
	Mh.vertices[k].lab=6;
	indic=0;
      }
      if(Mh.vertices[k].x>=(7./16.) && Mh.vertices[k].x-Dx < (7./16.) && Mh.vertices[k].y>=(13./16.) && Mh.vertices[k].y-Dy<(13./16.)){
	Mh.vertices[k].x=(7./16.);
	Mh.vertices[k].y=(13./16.);
	Mh.vertices[k].lab=6;
	indic=0;
      }

      ///////////// Constrution carré 2////////////
       if(Mh.vertices[k].x<=(9./16.) &&  Mh.vertices[k].x+Dx >(9./16.) && Mh.vertices[k].y>=(3./16.) && Mh.vertices[k].y<=(7./16.)){
	Mh.vertices[k].x=(9./16.);
	Mh.vertices[k].lab=7;
	indic=0;
      }
      if(Mh.vertices[k].y<=(3./16.) && Mh.vertices[k].y+Dy>(3./16.) && Mh.vertices[k].x>=(9./16.) && Mh.vertices[k].x<=(13./16.)){
	Mh.vertices[k].y=(3./16.);
	Mh.vertices[k].lab=7;
	indic=0;
      }
      if(Mh.vertices[k].x>=(13./16.) && Mh.vertices[k].x-Dx <(13./16.) && Mh.vertices[k].y>=(3./16.) && Mh.vertices[k].y<=(7./16.)){
	Mh.vertices[k].x=(13./16.);
	Mh.vertices[k].lab=7;
	indic=0;
      }
      if(Mh.vertices[k].y>=(7./16.) && Mh.vertices[k].y-Dy<(7./16.) && Mh.vertices[k].x>=(9./16.) && Mh.vertices[k].x<=(13./16.)){
	Mh.vertices[k].y=(7./16.);
	Mh.vertices[k].lab=7;
	indic=0;
      }
      ////////////////////////////////////////////// coin ///////////////////////////////////:
      if(Mh.vertices[k].x<=(9./16.) && Mh.vertices[k].x+Dx > (9./16.) && Mh.vertices[k].y<=(3./16.) && Mh.vertices[k].y+Dy>(3./16.)){
	Mh.vertices[k].x=(9./16.);
	Mh.vertices[k].y=(3./16.);
	Mh.vertices[k].lab=7;
	indic=0;
      }
      if(Mh.vertices[k].x<=(9./16.) && Mh.vertices[k].x+Dx > (9./16.) && Mh.vertices[k].y>=(7./16.) && Mh.vertices[k].y-Dy<(7./16.)){
	Mh.vertices[k].x=(9./16.);
	Mh.vertices[k].y=(7./16.);
	Mh.vertices[k].lab=7;
	indic=0;
      }
       if(Mh.vertices[k].x>=(13./16.) && Mh.vertices[k].x-Dx < (13./16.) && Mh.vertices[k].y<=(3./16.) && Mh.vertices[k].y+Dy>(3./16.)){
	Mh.vertices[k].x=(13./16.);
	Mh.vertices[k].y=(3./16.);
	Mh.vertices[k].lab=7;
	indic=0;
      }
      if(Mh.vertices[k].x>=(13./16.) && Mh.vertices[k].x-Dx < (13./16.) && Mh.vertices[k].y>=(7./16.) && Mh.vertices[k].y-Dy<(7./16.)){
	Mh.vertices[k].x=(13./16.);
	Mh.vertices[k].y=(7./16.);
	Mh.vertices[k].lab=7;
	indic=0;
      }
      ////////////////////
      
      Mh.vertices[k].x=Mh.vertices[k].x+indic*theta*(Dx/6);
      Mh.vertices[k].y=Mh.vertices[k].y+indic*theta*(Dy/6);
       if(((3./16.)<=Mh.vertices[k].x) && (Mh.vertices[k].x<=(7./16.)) && ((9./16.)<=Mh.vertices[k].y) && (Mh.vertices[k].y<=(13./16.))){
      	Mh.vertices[k].lab=6;
       }
      if(((9./16.)<=Mh.vertices[k].x) && (Mh.vertices[k].x<=(13./16.)) && ((3./16.)<=Mh.vertices[k].y) && (Mh.vertices[k].y<=(7./16.))){
      	Mh.vertices[k].lab=7;
	
	}
    }
    indic=1;
     k++;
  }
}

 
 k=0;
 for(int iy=0;iy<d.Ny+2;iy++){
   for(int ix=0;ix<d.Nx+2;ix++){
     if(Mh.vertices[(d.Nx+3)*iy+ix].lab>-1 && Mh.vertices[(d.Nx+3)*iy+ix+1].lab>-1 && Mh.vertices[(d.Nx+3)*(iy+1)+ix+1].lab>-1 && Mh.vertices[(d.Nx+3)*(iy+1)+ix].lab>-1){
       
       label=0;
       
       if(Mh.vertices[(d.Nx+3)*iy+ix].lab==6 && Mh.vertices[(d.Nx+3)*iy+ix+1].lab==6 && Mh.vertices[(d.Nx+3)*(iy+1)+ix+1].lab==6 && Mh.vertices[(d.Nx+3)*(iy+1)+ix].lab==6){
	 label=1;
       }
       if(Mh.vertices[(d.Nx+3)*iy+ix].lab==7 && Mh.vertices[(d.Nx+3)*iy+ix+1].lab==7 && Mh.vertices[(d.Nx+3)*(iy+1)+ix+1].lab==7 && Mh.vertices[(d.Nx+3)*(iy+1)+ix].lab==7){
	 label=2;
       }        
       Mh.cells[k].init(Mh.vertices,(d.Nx+3)*iy+ix,(d.Nx+3)*iy+ix+1,(d.Nx+3)*(iy+1)+ix+1,(d.Nx+3)*(iy+1)+ix,label);
     }
     else{
       Mh.cells[k].init(Mh.vertices,(d.Nx+3)*iy+ix,(d.Nx+3)*iy+ix+1,(d.Nx+3)*(iy+1)+ix+1,(d.Nx+3)*(iy+1)+ix,-1);
     }	
     k++;	
   }
 }
 
 return Mh;
}


Mesh CheckMesh(Data & d){
  // Nx et Ny le nombre de maille dans les directions x et y
  double Dx=0.,Dy=0.;
  double theta=0;
  int label=0;
  int indic=1;
  int mod=0;
  Mesh M(d);
  
  if(d.Tx!=7 || d.Ty!=7){
    cout<<" Le domaine fait 7*7"<<endl;exit(1);
  }
  if((d.Nx%70==0 && d.Ny%70==0)){
   
    if (d.Nx==70){ mod=10;}
    if (d.Nx==140){ mod=20;}
    if (d.Nx==210){ mod=30;}
    if (d.Nx==280){ mod=40;}
    if (d.Nx==350){ mod=50;}
    if (d.Nx==420){ mod=60;}
    if (d.Nx==490){ mod=70;}
    if (d.Nx==560){ mod=80;}
    if (d.Nx==630){ mod=90;}
    if (d.Nx==700){ mod=100;}
    
     Mesh Mh(d);
    Dx=d.Tx/(double)d.Nx; // uniform step mesh in the x direction
    Dy=d.Ty/(double)d.Ny;
    
    Mh.nv=(d.Nx+3)*(d.Ny+3);
    Mh.nc=(d.Nx+2)*(d.Ny+2);
    

    Vertex iniver(d);
    Mh.vertices = new Vertex[Mh.nv];
    for(int i=0;i<Mh.nv;i++){
      Mh.vertices[i]=iniver;
    }

 cell ini(Mh.nbnodelocal);
 Mh.cells = new cell[Mh.nc];
 for(int i=0;i<Mh.nc;i++){
   Mh.cells[i]=ini;
 }
 
 srand(52);
 int k=0;
 for(int iy=0;iy<d.Ny+3;iy++){
   for(int ix=0;ix<d.Nx+3;ix++){
     // si on veut changer srand(time(NULL));
     theta=(rand()/(float)RAND_MAX)*2.0-1.0;
     Mh.vertices[k].x=ix*Dx;
     Mh.vertices[k].y=iy*Dy;
       
     ///////////////////////////////////// construction des carré////////////////////
           
     Mh.vertices[k].lab=0;
     ///////////////: construction carré  1//////////////
     if((ix%mod==0 ) || (iy%mod==0)){
	indic=0;
     }
     if((ix>=3*mod && ix<=4*mod) && (iy>=5*mod && iy<=6*mod)){
       indic=1;//Carré manquant en haut//
     }
     if((iy>6*mod) || (iy<mod)){
       indic=1;
     }
     
     if(ix==0) {Mh.vertices[k].lab=-1;}
     if(iy==d.Ny+2) {Mh.vertices[k].lab=-1;}
     if(ix==d.Nx+2) {Mh.vertices[k].lab=-1;}
     if(iy==0) {Mh.vertices[k].lab=-1;}
     if((ix==0 && iy==0) || (ix==d.Nx+2 && iy==0) || (iy==d.Ny+2 && ix==0) || (iy==d.Ny+2 && ix==d.Nx+2)){
       Mh.vertices[k].lab=-1;}

     if(ix==1 && iy>0 && iy<d.Ny+2) {Mh.vertices[k].lab=2;}
     if(iy==d.Ny+1 && ix>0 && ix<d.Nx+2) {Mh.vertices[k].lab=3;}
     if(ix==d.Nx+1 && iy>0 && iy<d.Ny+2) {Mh.vertices[k].lab=4;}
     if(iy==1 && ix>0 && ix<d.Nx+2) {Mh.vertices[k].lab=5;}
     if((ix==1 && iy==1) || (ix==d.Nx+1 && iy==1) || (iy==d.Ny+1 && ix==1) || (iy==d.Ny+1 && ix==d.Nx+1)){
       Mh.vertices[k].lab=1;}
     
    
    if(iy>=mod && iy<=2*mod){
      if((ix>=mod && ix<=2*mod) || (ix>=3*mod && ix<=4*mod) || (ix>=5*mod && ix<=6*mod)){
	Mh.vertices[k].lab=6;
      }
    }
    if((iy>=2*mod && iy<=3*mod) ||(iy>=4*mod && iy<=5*mod)){
      if((ix>=2*mod && ix<=3*mod) || (ix>=4*mod && ix<=5*mod)){
	Mh.vertices[k].lab=6;
      }
    }
    if(iy>=3*mod && iy<=4*mod){
      if((ix>=mod && ix<=2*mod) || (ix>=5*mod && ix<=6*mod)){
	Mh.vertices[k].lab=6;
      }
      if(ix>=3*mod && ix<=4*mod){
	Mh.vertices[k].lab=7;
      }
    }
    if(iy>=5*mod && iy<=6*mod){
      if((ix>=mod && ix<=2*mod) || (ix>=5*mod && ix<=6*mod)){
	Mh.vertices[k].lab=6;
      }
    }

     if(ix!=0 && iy!=0 && ix!=d.Nx+2 && iy!=d.Ny+2 && ix!=1 && iy!=1 && ix!=d.Nx+1 && iy!=d.Ny+1) {
       Mh.vertices[k].x=Mh.vertices[k].x+indic*theta*(Dx/4);
       Mh.vertices[k].y=Mh.vertices[k].y+indic*theta*(Dy/4);
     }
      
    indic=1;
     k++;
  }
 }

 
k=0;
 for(int iy=0;iy<d.Ny+2;iy++){
  for(int ix=0;ix<d.Nx+2;ix++){ 
    label=0;
      if(Mh.vertices[(d.Nx+3)*iy+ix].lab>-1 && Mh.vertices[(d.Nx+3)*iy+ix+1].lab>-1 && Mh.vertices[(d.Nx+3)*(iy+1)+ix+1].lab>-1 && Mh.vertices[(d.Nx+3)*(iy+1)+ix].lab>-1){
 
	if(Mh.vertices[(d.Nx+3)*iy+ix].lab==6 && Mh.vertices[(d.Nx+3)*iy+ix+1].lab==6 && Mh.vertices[(d.Nx+3)*(iy+1)+ix+1].lab==6 && Mh.vertices[(d.Nx+3)*(iy+1)+ix].lab==6){
	  label=1;
	}
	if(Mh.vertices[(d.Nx+3)*iy+ix].lab==7 && Mh.vertices[(d.Nx+3)*iy+ix+1].lab==6 && Mh.vertices[(d.Nx+3)*(iy+1)+ix+1].lab==6 && Mh.vertices[(d.Nx+3)*(iy+1)+ix].lab==6){
	  label=1;
	}
	if(Mh.vertices[(d.Nx+3)*iy+ix].lab==6 && Mh.vertices[(d.Nx+3)*iy+ix+1].lab==7 && Mh.vertices[(d.Nx+3)*(iy+1)+ix+1].lab==6 && Mh.vertices[(d.Nx+3)*(iy+1)+ix].lab==6){
	  label=1;
	}
	if(Mh.vertices[(d.Nx+3)*iy+ix].lab==6 && Mh.vertices[(d.Nx+3)*iy+ix+1].lab==6 && Mh.vertices[(d.Nx+3)*(iy+1)+ix+1].lab==7 && Mh.vertices[(d.Nx+3)*(iy+1)+ix].lab==6){
	  label=1;
	}
	if(Mh.vertices[(d.Nx+3)*iy+ix].lab==6 && Mh.vertices[(d.Nx+3)*iy+ix+1].lab==6 && Mh.vertices[(d.Nx+3)*(iy+1)+ix+1].lab==6 && Mh.vertices[(d.Nx+3)*(iy+1)+ix].lab==7){
	  label=1;
	}
	if(Mh.vertices[(d.Nx+3)*iy+ix].lab==7 && Mh.vertices[(d.Nx+3)*iy+ix+1].lab==7 && Mh.vertices[(d.Nx+3)*(iy+1)+ix+1].lab==7 && Mh.vertices[(d.Nx+3)*(iy+1)+ix].lab==7){
	  label=2;
	}   
	Mh.cells[k].init(Mh.vertices,(d.Nx+3)*iy+ix,(d.Nx+3)*iy+ix+1,(d.Nx+3)*(iy+1)+ix+1,(d.Nx+3)*(iy+1)+ix,label);
      }
      else {
	Mh.cells[k].init(Mh.vertices,(d.Nx+3)*iy+ix,(d.Nx+3)*iy+ix+1,(d.Nx+3)*(iy+1)+ix+1,(d.Nx+3)*(iy+1)+ix,-1);
	}
      k++;	
    }
 }
 
 return Mh;
 }
 else{
   return M;
   cout<< " Le maillage doit faire 70 ou 140 maille par direction "<<endl;
 }
 
}
