#ifndef FUNCTIONSCONNECT_HPP
#define FUNCTIONSCONNECT_HPP

#include <iostream>
#include <cmath>
#include <cassert>
#include "ClassMesh.hpp"
#include "math.h"

/** Structure qui contient les mailles associées a un noeud
    @brief Structure qui contient les mailles associées a un noeud**/
struct TabCellLocal{
  int taille;
  int *TabCell;
};

/** Structure pour la table de connectivité 
@brief Structure pour la table de connectivité **/
 struct TabConnecInv{
  int nbvertex;
  TabCellLocal *TabInv;
};

void CellLocal( Mesh & Mh,int numnode, TabCellLocal & CellOfNode);

int InverseEdge( Mesh & Mh,int node1,int node2,int numcell,TabConnecInv & tab);

int tabVF9(Mesh & Mh, TabConnecInv & tab,int *& tabVF9,int cell);

int PresenceNodetoTab(int * tab,int i,int taille); 

void CreateTabInv(Mesh & Mh, TabConnecInv & tab);

int NodeGtoL(Mesh & Mh,int  numCell,int numNode);

double StepMesh(Mesh & Mh);

int NextNodeLocal(Mesh & Mh,int r);

int PreviousNodeLocal(Mesh & Mh,int r);

#endif
