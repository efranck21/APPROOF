#ifndef WRITEDATA_HPP
#define WRITEDATA_HPP

#include <iostream>
#include <cassert>
#include "Data.hpp"
#include "Variables.hpp"
#include "ClassMesh.hpp"
#include "Flux.hpp"

using namespace std;

void SaveData(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,int ntAnim);

void SaveOneData(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,int Num,int ntAnim);

void  SaveOneData1D(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,int Num,int ntAnim);

void SaveOneDatacentercell(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,int Num,int ntAnim);

void SaveMesh(Data & d,Mesh & Mh);

void SaveFonction2D(Data & d,Mesh & Mh,ParamPhysic & Param,float t,Vertex c,int var);

void SaveFonction2DEuler(Data & d,Mesh & Mh,ParamPhysic & Param,float t,Vertex c,int var);

void SaveFonction1D(Data & d,Mesh & Mh,ParamPhysic & Param,float t,Vertex c,int var);

void SaveFonction(Data & d,Mesh & Mh,ParamPhysic & Param,float t,Vertex c);

void SaveInterEner(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,int ntAnim);

void SaveRestart(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,int nstep, double time);

void LoadRestart(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param, double & InitTime);

#endif
