#ifndef HOFUNCTIONS_HPP
#define HOFUNCTIONS_HPP
#include "UsualMesh.hpp"
#include "Data.hpp"
#include "ClassMesh.hpp"
#include "Variables.hpp"
#include "svd_solver.hpp"
#include "ParamEuler.hpp"

using namespace std;

void ConstructionLSMatrix(Data & d, Mesh & Mh, variable & v,TabConnecInv & tab,ParamEuler & Euler, LSProblem & LS,int PolynomialOrder,int numGr);

void ConstructionLSRHS(Data & d, Mesh & Mh, variable & v,TabConnecInv & tab,ParamEuler & Euler, LSProblem & LS,int PolynomialOrder,int numGr);

void SolverLSProblem(Data & d, Mesh & Mh, variable & v,TabConnecInv & tab,ParamEuler & Euler, LSProblem & LS,int PolynomialOrder,int numGr);

void ConstructionLSStencil(Data & d, Mesh & Mh, variable & v,TabConnecInv & tab,ParamEuler & Euler, LSProblem & LS,int PolynomialOrder,int numGr);

void LSInit(Data & d, Mesh & Mh, variable & v,TabConnecInv & tab,ParamEuler & Euler, LSProblem & LS,int PolynomialOrder,int numGr);

void QuadratureUnitSquare(double & p, double & w, int nb, int Order);

double Mapping(Data & d, Mesh & Mh,R2 & p, double & J,double px, double py,cell & C);

double ComputeMatrixcoefficent(Data & d, Mesh & Mh,int d1, int d2, int numCell,int numGr,int Order);

void LocalDegreePolynomial(int & d1, int & d2, int j,int PolynomialOrder);

#endif
