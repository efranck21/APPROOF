#ifndef CLASSMESH_HPP
#define CLASSMESH_HPP

#include <cmath>
#include <cstdlib>
#include <iostream> 
#include "R2.hpp"
#include <cassert>
#include <fstream>
#include "Variables.hpp"
using namespace std;
/**
Class : Vertex

contains:
- lab : Vertex label (uselful for different sub domain)
- R2 which correspond to the vertex

**/
class Vertex : public R2{
public:
  int lab;
  Vertex(){ R2(); lab=0; }

  Vertex(R2 P, int r,Data d){
    vectorflux vec(d);
    x=P.x;
    y=P.y; 
    lab=r;
  }
  Vertex(R2 P,Data d){
    vectorflux vec(d);
    x=P.x;
    y=P.y;
    lab=0;
  }
  Vertex(Data d){
    R2(0,0);
    lab=0;
  }
  ~Vertex(){}

  Vertex(const Vertex & v){
    lab=v.lab;
    x=v.x;
    y=v.y;
  }

   Vertex & operator =( Vertex  v){
    lab=v.lab;
    x=v.x;
    y=v.y;
    return *this; 
  }
  
};

/**
Class : cell

- nbnodelocal : number of vertex in the cell
- vertices : table of vertex
- numver : local intex of vertex
- area : area of the cell
- normevertex : table of nodal normal njr for the cell j
- lvertex : table of nodal lenght ljr for the cell j
- normeedge : table of edge normal njk for the cell j
- ledge : table of edge normal ljk for the cell j
- lab : label of the cell (for subdomain)

**/

class cell : public R2{
public:
  int nbnodelocal;
  Vertex *vertices;
  int *numver;
  R area;
  Vertex centercell;
  R2 *normvertex;     
  double *lvertex;          
  R2 *normedge;       
  double *ledge;
  int lab;
  cell(){
    nbnodelocal=0;
    vertices=NULL;
    area=0;
    normvertex=NULL;
    lvertex=NULL;
    normedge=NULL;
    ledge=NULL;
    numver=NULL;
    lab=0;
  };
  cell(int nbnod){
    nbnodelocal=nbnod;
    vertices=NULL;
    area=0;
    normvertex=NULL;
    lvertex=NULL;
    normedge=NULL;
    ledge=NULL;
    lab=0;
    numver=NULL;
  };
  cell( const cell & c){
    nbnodelocal=c.nbnodelocal;
    area=c.area;
    lab=c.lab;
    centercell=c.centercell;
    vertices =new Vertex[nbnodelocal];
    normvertex =new R2[nbnodelocal];
    lvertex =new double[nbnodelocal];
    normedge=new R2[nbnodelocal];
    ledge=new double[nbnodelocal];
    numver = new int[nbnodelocal];
    if(c.vertices!=NULL) {
     for(int i=0;i<nbnodelocal;i++){
       vertices[i]=c.vertices[i];}
     }
     if(c.normvertex!=NULL) {
     for(int i=0;i<nbnodelocal;i++){
       normvertex[i]=c.normvertex[i];}
     }
     if(c.lvertex!=NULL) {
     for(int i=0;i<nbnodelocal;i++){
       lvertex[i]=c.lvertex[i];}
     }
      if(c.normedge!=NULL) {
     for(int i=0;i<nbnodelocal;i++){
       normedge[i]=c.normedge[i];}
     }
     if(c.ledge!=NULL) {
     for(int i=0;i<nbnodelocal;i++){
       ledge[i]=c.ledge[i];}
     }

     if(c.numver!=NULL) {
     for(int i=0;i<nbnodelocal;i++){
       numver[i]=c.numver[i];}
     }
     }
     
     

     cell & operator = (cell c){
    nbnodelocal=c.nbnodelocal;
    area=c.area;
    lab=c.lab;
    centercell=c.centercell;

       if(vertices!=NULL)
       {delete[] vertices;}
     if(normvertex!=NULL)
       {delete[] normvertex;}
     if(lvertex!=NULL){
       delete[] lvertex;}
     if(normedge!=NULL)
       {delete[] normedge;}
     if(ledge!=NULL)
       {delete[] ledge;}
     if(numver!=NULL)
       {delete[] numver;}

    vertices =new Vertex[nbnodelocal];
    normvertex =new R2[nbnodelocal];
    lvertex =new double[nbnodelocal];
    normedge=new R2[nbnodelocal];
    ledge=new double[nbnodelocal];
    numver=new int[nbnodelocal];
    if(c.vertices!=NULL) {
     for(int i=0;i<nbnodelocal;i++){
       vertices[i]=c.vertices[i];}
     }
     if(c.normvertex!=NULL) {
     for(int i=0;i<nbnodelocal;i++){
       normvertex[i]=c.normvertex[i];}
     }
     if(c.lvertex!=NULL) {
     for(int i=0;i<nbnodelocal;i++){
       lvertex[i]=c.lvertex[i];}
     }
      if(c.normedge!=NULL) {
     for(int i=0;i<nbnodelocal;i++){
       normedge[i]=c.normedge[i];}
     }
     if(c.ledge!=NULL) {
     for(int i=0;i<nbnodelocal;i++){
       ledge[i]=c.ledge[i];}
     }
      if(c.numver!=NULL) {
     for(int i=0;i<nbnodelocal;i++){
       numver[i]=c.numver[i];}
      }
     return *this;
     }
    
  ~cell(){
    if(vertices!=NULL)
       {delete[] vertices;}
     if(normvertex!=NULL)
       {delete[] normvertex;}
     if(lvertex!=NULL){
       delete[] lvertex;}
     if(normedge!=NULL)
       {delete[] normedge;}
     if(ledge!=NULL)
       {delete[] ledge;}
     if(numver!=NULL)
       {delete[] numver;}
    
     };
  /**
     functions which initialize the cell using the global index of then ode
     @param v : table of all the vertices
     @param i0,i1,i2,i3 global index of the vertex in the tavle v
     @param ir : label of the cell.
  **/

  void init(Vertex * v,int i0,int i1,int i2,int i3,int ir){
    if(vertices!=NULL)
       {delete[] vertices;}
     if(normvertex!=NULL)
       {delete[] normvertex;}
     if(lvertex!=NULL){
       delete[] lvertex;}
     if(normedge!=NULL)
       {delete[] normedge;}
     if(ledge!=NULL)
       {delete[] ledge;}
     if(numver!=NULL)
       {delete[] numver;}

    vertices=new Vertex[nbnodelocal];
    normvertex=new R2[nbnodelocal];
    lvertex=new double[nbnodelocal];
    normedge=new R2[nbnodelocal];
    ledge=new double[nbnodelocal];
    numver=new int[nbnodelocal];
    
    if(nbnodelocal==3){
  vertices[0]=v[i0]; vertices[1]=v[i1]; vertices[2]=v[i2];
  numver[0]=i0; numver[1]=i1; numver[2]=i2;

    lvertex[0]=0.5*sqrt(pow((vertices[1].x-vertices[2].x),2)+pow((vertices[1].y-vertices[2].y),2));
    lvertex[1]=0.5*sqrt(pow((vertices[2].x-vertices[0].x),2)+pow((vertices[2].y-vertices[0].y),2));
    lvertex[2]=0.5*sqrt(pow((vertices[0].x-vertices[1].x),2)+pow((vertices[0].y-vertices[1].y),2));

    normvertex[0].x=(1/(2*lvertex[0]))*(-vertices[2].y+vertices[1].y);
    normvertex[0].y=(1/(2*lvertex[0]))*(vertices[2].x-vertices[1].x);
    normvertex[1].x=(1/(2*lvertex[1]))*(-vertices[0].y+vertices[2].y);
    normvertex[1].y=(1/(2*lvertex[1]))*(vertices[0].x-vertices[2].x);
    normvertex[2].x=(1/(2*lvertex[2]))*(-vertices[1].y+vertices[0].y);
    normvertex[2].y=(1/(2*lvertex[2]))*(vertices[1].x-vertices[0].x);
    
    ledge[0]=2*lvertex[2]; 
    ledge[1]=2*lvertex[0];
    ledge[2]=2*lvertex[1];
 
    normedge[0].x=-normvertex[2].x;
    normedge[0].y=-normvertex[2].y;
    normedge[1].x=-normvertex[0].x;
    normedge[1].y=-normvertex[0].y;
    normedge[2].x=-normvertex[1].x; 
    normedge[2].y=-normvertex[1].y;
   
    area=0.5*((vertices[0]^vertices[1])+(vertices[1]^vertices[2])+(vertices[2]^vertices[0]));
    } 
    else if(nbnodelocal==4){
    
      
 vertices[0]=v[i0]; vertices[1]=v[i1]; vertices[2]=v[i2]; vertices[3]=v[i3];
 numver[0]=i0; numver[1]=i1; numver[2]=i2; numver[3]=i3;
    lvertex[0]=0.5*sqrt(pow((vertices[1].x-vertices[3].x),2)+pow((vertices[1].y-vertices[3].y),2));
    lvertex[1]=0.5*sqrt(pow((vertices[2].x-vertices[0].x),2)+pow((vertices[2].y-vertices[0].y),2));
    lvertex[2]=0.5*sqrt(pow((vertices[3].x-vertices[1].x),2)+pow((vertices[3].y-vertices[1].y),2));
    lvertex[3]=0.5*sqrt(pow((vertices[0].x-vertices[2].x),2)+pow((vertices[0].y-vertices[2].y),2));
  

    normvertex[0].x=(1./(2*lvertex[0]))*(-vertices[3].y+vertices[1].y);
    normvertex[0].y=(1./(2*lvertex[0]))*(vertices[3].x-vertices[1].x);
    normvertex[1].x=(1./(2*lvertex[1]))*(-vertices[0].y+vertices[2].y);
    normvertex[1].y=(1./(2*lvertex[1]))*(vertices[0].x-vertices[2].x);
    normvertex[2].x=(1./(2*lvertex[2]))*(-vertices[1].y+vertices[3].y);
    normvertex[2].y=(1./(2*lvertex[2]))*(vertices[1].x-vertices[3].x);
    normvertex[3].x=(1./(2*lvertex[3]))*(-vertices[2].y+vertices[0].y);
    normvertex[3].y=(1./(2*lvertex[3]))*(vertices[2].x-vertices[0].x); 
 
    
    ledge[0]=sqrt(pow((vertices[1].x-vertices[0].x),2)+pow((vertices[1].y-vertices[0].y),2));   
    ledge[1]=sqrt(pow((vertices[2].x-vertices[1].x),2)+pow((vertices[2].y-vertices[1].y),2));
    ledge[2]=sqrt(pow((vertices[3].x-vertices[2].x),2)+pow((vertices[3].y-vertices[2].y),2));
    ledge[3]=sqrt(pow((vertices[0].x-vertices[3].x),2)+pow((vertices[0].y-vertices[3].y),2));
 
    normedge[0].x=(1/ledge[0])*(vertices[1].y-vertices[0].y);    
    normedge[0].y=(1/ledge[0])*(vertices[0].x-vertices[1].x);
    normedge[1].x=(1/ledge[1])*(vertices[2].y-vertices[1].y);
    normedge[1].y=(1/ledge[1])*(vertices[1].x-vertices[2].x);
    normedge[2].x=(1/ledge[2])*(vertices[3].y-vertices[2].y);
    normedge[2].y=(1/ledge[2])*(vertices[2].x-vertices[3].x);
    normedge[3].x=(1/ledge[3])*(vertices[0].y-vertices[3].y);
    normedge[3].y=(1/ledge[3])*(vertices[3].x-vertices[0].x);
  
    area=0.5*((vertices[0]^vertices[1])+(vertices[1]^vertices[2])+(vertices[2]^vertices[3])+(vertices[3]^vertices[0]));
    }
    else {cout << "problem on the number of vertex"<<endl;exit(1);}
    lab=ir;
	
	
      assert(area>=0);
  }
  /**function which compute the center of the cell **/
   void initCC(R2 a,int r){
      centercell.x=a.x;
      centercell.y=a.y;
      centercell.lab=r;
    }

  };

   

/** Class: Mesh
    - nv : Vertices number
    - nc : Cell number
    - nbnodelocal : number of nodes by cell
    - vertices : table of vertices
    - cells : tables of cells
     
**/
class Mesh { public :
  int nv,nc,nbnodelocal;
   Vertex *vertices;
   cell *cells;
   
    
  Mesh(Data d){
    if(d.Typemesh=='T') {
      nbnodelocal=3;
    }
    if(d.Typemesh=='Q'){
      nbnodelocal=4;
    }
    nv=0; nc=0;
    vertices=NULL;
    cells=NULL;
    
     
   }
  Mesh(const Mesh &  M){
    
     
     nv=M.nv; nc=M.nc; nbnodelocal=M.nbnodelocal;
     cell ini(nbnodelocal);
     Vertex iniver;

     vertices = new Vertex[nv];
     for(int i=0;i<nv;i++){
       vertices[i]=iniver;}


     cells= new cell[nc]; 
     for(int i=0;i<nc;i++){
       cells[i]=ini;}

     if(M.vertices!=NULL) {
     for(int i=0;i<nv;i++){
       vertices[i]=M.vertices[i];
     }
     }
      if(M.cells!=NULL) {
     for(int i=0;i<nc;i++){
       cells[i]=M.cells[i];
     }
      }
   }
   ~Mesh(){
    
     if(vertices!=NULL)
       {delete[] vertices;}
     if(cells!=NULL)
       {delete[] cells;}
   }
     
 


  Mesh & operator =(Mesh M){
     nv=M.nv; nc=M.nc; nbnodelocal=M.nbnodelocal;

     Vertex iniver;
      if(vertices!=NULL)
       {delete[] vertices;}
      if(cells!=NULL)
       {delete[] cells;}

     
     vertices = new Vertex[nv];
     for(int i=0;i<nv;i++){
       vertices[i]=iniver;
     }

    
     

     cell ini(nbnodelocal);
     cells= new cell[nc]; 
     for(int i=0;i<nc;i++){
       cells[i]=ini;} 

     if(M.vertices!=NULL) {
     for(int i=0;i<nv;i++){
       vertices[i]=M.vertices[i];
     }
     }
      if(M.cells!=NULL) {
     for(int i=0;i<nc;i++){
       cells[i]=M.cells[i];
     }
      }

 
     return *this;
   }

  /** operator which gives the node j in the cell it **/
  int operator ()(int it,int j){ return cells[it].numver[j];}

  /** function which gives the nodal normal for the node r in the cell j **/
  R2 njr(int j,int r){ 
    if(r >=0 && r < nbnodelocal){
      return cells[j].normvertex[r];}
    else{exit(1);}
  }

  /** function which gives the nodal normal for the edge k in the cell j **/
  R2 njk(int j,int k){ 
    if(k >=0 && k < nbnodelocal){
      return cells[j].normedge[k];}
    else{exit(1);}
  }

  /** function which gives the nodal lenght for the edge r in the cell j **/
  double ljr(int j, int r){
    if(r >=0 && r < nbnodelocal){
      return cells[j].lvertex[r];}
    else{exit(1);}
  }

  /** function which gives the nodal lenght for the edge k in the cell j **/
  double ljk(int j, int k) {
      if(k >=0 && k < nbnodelocal){
    return cells[j].ledge[k];}
      else{exit(1);}
  }
  
  /** function which gives the area of the cell j **/
  double area(int j){
    return cells[j].area;
  }

    /** function which gives the lab of the cell j **/
  int celllab(int j){
    return cells[j].lab;
  }

    /** function which gives the center of the cell j **/
  Vertex xj(int j){
    return cells[j].centercell;
  }

    /** function which gives the global node r **/
  Vertex & xr(int r){
    return vertices[r];
  }

    /** function which gives the global index of the node with the local index r in the cell j **/
  Vertex xr(int j,int r){
    if(r >=0 && r < nbnodelocal){
     return cells[j].vertices[r];}
   else{exit(1);}
      }   

 };


    
#endif	       
  
    
    

    


