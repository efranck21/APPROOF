#ifndef TENSOR_HPP
#define TENSOR_HPP
#include <cassert>
#include <iostream>
#include <cmath>
#include<cstdlib>
#include "Data.hpp"
#include "ClassMesh.hpp"
#include "Functions.hpp"
#include "Variables.hpp"


class tensor {
public:
  double ten[2][2];
  int num;
  char type;

  tensor(){
  
    num=0;
    type='0';
    ten[0][0]=0.;
    ten[0][1]=0.;
    ten[1][0]=0.;
    ten[1][1]=0.;
  }
  
  tensor(Data d,char c){
    if(c=='h' || c=='s'){
    type=c;
    }
    else{ cout <<"probleme dans le type de tenseur demander"<<endl;exit(1);}
    num =0;
  
    ten[0][0]=0.;
    ten[0][1]=0.;
    ten[1][0]=0.;
    ten[1][1]=0.;
    
  }
  
  
  tensor(const tensor & t){
    num=t.num;
    type=t.type;
    
    for(int i=0;i<2;i++){
      for(int j=0;j<2;j++){
	ten[i][j]=t.ten[i][j];
	
      }
      
    }
  }
  
  ~tensor(){}
     
	
 tensor & operator=( tensor t){
    num=t.num;
    type=t.type;
    
    for(int i=0;i<2;i++){
      for(int j=0;j<2;j++){
	ten[i][j]=t.ten[i][j];
      }
      
    }
    

    return *this;
 } 

};   

tensor inittensor(Data & d,Mesh & Mh,variable & v, TabConnecInv & tab,char c,int numGr,int numCell);


#endif
