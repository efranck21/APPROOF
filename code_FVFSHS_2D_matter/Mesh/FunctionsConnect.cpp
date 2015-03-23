#include <iostream>
#include <cmath>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"
#include "Functions.hpp"
#include "Data.hpp"
#include "Functionsgeo.hpp"
#include "FunctionsConnect.hpp"
#include "Tensor.hpp"


/** File which contains some functions linked tothe connectivity of the mesh**/

void CellLocal( Mesh & Mh,int numnode, TabCellLocal & CellOfNode){
  /** give a table with the number of the cell which have the node "numnode" **/

  int c=0;
  int end_list=-1;
  int *head_s = new int[Mh.nv];
  int *next_p= new int[Mh.nc*Mh.nbnodelocal];
    int i,j,k,p;

    for(i=0;i<Mh.nv;i++){
      head_s[i] = end_list;
    }
    for(k=0;k<Mh.nc;k++){     
      for(j=0;j<Mh.nbnodelocal;j++)
	{p = Mh.nbnodelocal*k+j;
	  i = Mh(k,j);
	  next_p[p]=head_s[i];
	  head_s[i]=p;
	}
    }

    for(int p=head_s[numnode]; p!=end_list;p=next_p[p]){
      k=p/Mh.nbnodelocal;
      j=p % Mh.nbnodelocal;
       assert(numnode==Mh(k,j));
       
	c++;
    }

    CellOfNode.taille=c;
    CellOfNode.TabCell= new int[c];
    int cc=0;
    for(int p=head_s[numnode]; p!=end_list;p=next_p[p]){
      k=p/Mh.nbnodelocal;
      j=p % Mh.nbnodelocal;
	CellOfNode.TabCell[cc]=k;
      cc++;
    }
    delete [] head_s;
    delete [] next_p;
    
}


void CreateTabInv(Mesh & Mh, TabConnecInv & tab){
  /** function which create the connectivy table**/
  tab.nbvertex=Mh.nv;
  tab.TabInv = new TabCellLocal[Mh.nv];
  int i;

#pragma omp parallel for
  for( i=0;i<Mh.nv;i++){
    CellLocal(Mh,i,tab.TabInv[i]);
  }
}
    


int InverseEdge( Mesh & Mh,int node1,int node2,int numcell,TabConnecInv & tab){
  /** function which take an Edge (defined by 2 nodes) and a cell and gives the cell which have the same edge.**/ 

  int c=0;
  int cell=0;

  for(int i=0;i<tab.TabInv[node1].taille;i++){
    for(int j=0;j<tab.TabInv[node2].taille;j++){
      if(tab.TabInv[node1].TabCell[i]==tab.TabInv[node2].TabCell[j] && tab.TabInv[node1].TabCell[i]!=numcell )
	{cell=tab.TabInv[node1].TabCell[i];

	  c++;}

    }
  }
  if( c==0) {return -1;}
  if (c==1)  {return cell;}
  else { cout << "problem they are some neighbour for this edge"<<" "<<node1<<" "<<node2<<endl; exit(1);}

}
  
  

int tabVF9(Mesh & Mh, TabConnecInv & tab,int *& tabVF9,int cell){
  /** Give the stencil of the nodal scheme for the cell "cell" (9 cells) **/
  int taille;
  int k,c=0;
  int p=0;
  int numGr=0;
  int * temp=NULL;

  if(Mh.nbnodelocal==4){
    taille=9;
  }
  
    if(Mh.nbnodelocal==3){
    taille=50;
  }

  temp = new int[taille];
  for(int j=0;j<taille;j++){
    temp[j]=-1;
  }

  for(int r=0;r<Mh.nbnodelocal;r++){
    numGr=Mh(cell,r);
     for(int j=0;j<tab.TabInv[numGr].taille;j++){

       k=PresenceNodetoTab(temp,tab.TabInv[numGr].TabCell[j],taille);
       if(k==0){
	 temp[p]=tab.TabInv[numGr].TabCell[j];
	 p++;
       }
  
     }
  }

  for(int j=0;j<taille;j++){
    if (temp[j]!=-1) {
      c++;
	}
  }
  
  tabVF9 = new int[c];
  for(int j=0;j<c;j++){
    tabVF9[j]=temp[j];
  }

  delete [] temp;
  return c;
}






int PresenceNodetoTab(int * tab,int i,int taille) {
  /** verify if the node i in the table "tab" **/
  int res=0; 
  for(int j=0;j<taille;j++){
    
    if(tab[j]==i){ 
      res=1;
      break;}
  }   
return res;
   
  }


int NodeGtoL(Mesh & Mh,int  numCell,int numNode){
  /** function which give the local index of the vertex with the global index ans the cell **/
    int res=0,c=0;
    for(int i=0;i<Mh.nbnodelocal;i++){
      if(Mh(numCell,i)==numNode){
	c++;
	res=i;
      }
    }
    if(c==0){cout<<"The node is not in this cell"<<endl;exit(1);}
    else if(c==1) {return res;}
    else {cout <<" problem "<<endl;exit(1);}

  }  



double StepMesh(Mesh & Mh){
  /** Function which compute the step mesh**/
  double tabmax[Mh.nc];
  int max=0;
  
  for(int j=0;j<Mh.nc;j++){
    for(int i=0;i<Mh.nbnodelocal;i++){
      if(Mh.ljk(j,i)>Mh.ljk(j,max)){
	max=i;}
    }
    tabmax[j]=Mh.ljk(j,max);
  }
  
  max=0;

  for(int j=0;j<Mh.nc;j++){
    if(tabmax[j]>tabmax[max]){
      max=j;}

  }

  return tabmax[max];

}


int NextNodeLocal(Mesh & Mh,int r){
  /** Function which give the following node in the cell**/
  int rp=0;
   if(r==Mh.nbnodelocal-1){
    rp=0;
  }
  else{
    rp=r+1;
  }
   return rp;
}

int PreviousNodeLocal(Mesh & Mh,int r){
  /** Function which give the previous node in the cell**/
  int rm=0;
   if(r==0){
    rm=Mh.nbnodelocal-1;
  }
  else{
    rm=r-1;
  }
   return rm;
}
