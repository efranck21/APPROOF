#include <cmath>
#include <cstdlib>
#include <iostream>
#include "Flux.hpp"
#include "Initialisation.hpp"
#include "WriteData.hpp"
#include "AnalyticSolutions.hpp"
#include "Diagnostics.hpp"
#include <fstream>
#include "Functions.hpp"
#include "TimeDiscret.hpp"
#include "Functionsgeo.hpp"
#include "Functionsinit.hpp"
#include "ParamPhys.hpp"

using namespace std;

/** Following works :**/
/** 1) Check test case for 2 and 6 for P1, 1 and 2 for P1 matter **/
/** 2) Droniou nonlinear diffusion scheme (Manu) **/
/** 3) P3 model (Manu) **/
/** 4) Generic Sn models  (Manu) **/
/** 5) Second order extension. For now just write for P1 system without source term and works only without limiting **/
/** 6) Multi-groupes models (Thomas) **/



int main() { 
  /** Parameter of the code initialized **/
  Data d("parametre.dat"); 
 
  /** Mesh construction **/
  Mesh Mh(d);    
  Mh=ChoiceMesh(d);

  
  R2 p;
  for(int i=0;i<Mh.nc;i++){
    p=GravityCenter(Mh.cells[i],Mh.nbnodelocal);
    Mh.cells[i].initCC(p,1);
  }
  int in=centerdomain(d,Mh);
  
      
  /** Connectivity table construction **/
  struct TabConnecInv TConnectInv;
  CreateTabInv(Mh,TConnectInv);  

  /** Physical parameters construction **/
  ParamPhysic Param(d,Mh);
  
  /** Variable contruction and initialization **/
  variable v(d,Mh.nc,Mh.nv);
  
  /** Physical parameter initialization **/
  double time=0.;
  double dt=0.;
  ParamPhysic_Init(d,Mh,v,TConnectInv,Param);
  ChoiceInit(d,Mh,v,Param,time);
  cout <<"End of parameters and variables construction"<<endl;

  /** Plot informations of the run and plot mesh **/
  WriteDataToMain(d,Param);
  SaveMesh(d,Mh);
 
  /** definition of the time quantities **/
  int SI=0;
  dt=CFL(d,Mh,v,TConnectInv,Param);
  int ntAnim=0;

 
 


  /** Time loop **/    
  if(d.Typetime=='E') {
    ExplicitDiscret(d,Mh,v,TConnectInv,Param,dt,ntAnim,time);
  }
		

  if(d.Typetime=='S') {
    SI=SIschemeImplemented(d,Mh,v,TConnectInv,Param);
    if(SI==0) {
      cout<<"The semi implicit version of this scheme does not exist" <<endl;
      exit(0);
    }
    else {
      SemiImplicit(d,Mh,v,TConnectInv,Param,dt,ntAnim,time);
    }
  }
      
  cout << "Final time: "<<time<<endl;
  
  

  /** Save final data **/
  SaveRestart(d,Mh,v,Param,0,time);
  SaveData(d,Mh,v,Param,ntAnim);   
 
  /** free memory **/
  for(int i=0;i<TConnectInv.nbvertex;i++){
    delete [] TConnectInv.TabInv[i].TabCell;
  }
  delete [] TConnectInv.TabInv;
 
  
  /** Computation and plot of exact solutions. Computation of errors **/
  if(TraceSolfond(d,Param)==1){
    cout <<"The L"<<1<<" norm is : "<<Diagnostics_Error(d,Mh,v,Param,d.Tf,1.)<<endl;     
    cout <<"The L"<<2<<" norm is : "<<Diagnostics_Error(d,Mh,v,Param,d.Tf,2.)<<endl;
    SaveFonction(d,Mh,Param,d.Tf,Mh.cells[in].centercell); 
  }  
  cout <<" The step mesh is :"<<StepMesh(Mh)<<endl;
 
  
  return 0;
}


/********************** Add Physical parameter ***********************/
/**
Add a physical parameter in the model  "NameModel".
-we modify the file hhpParam/Param"NameModel".hpp adding the parameter in the class Param"NameModel"
-We modify the constructor associated with this class.

-we modify the file ParamPhysic/Param"NameModel".dat adding one line with the name of the paramter 
and one line with the value (if this parameter is integer or double not table)

-We add the reading of this parameter in the function Param"NameModel"_Initwithfile (file ParamPhysic/Param"NameModel".cpp)
(if this parameter is integer or double not table)

-We add the initialization of this parameter in the function Param"NameModel"_InitTab (file ParamPhysic/Param"NameModel".cpp)
(if this parameter is a table)

-In the code we access to the new parameter with Param."NameModel"."NameNewParameter"

**/
/*********************************************************************/





/************************* Add code parameter ************************/
/**
Add a code parameter.
-we modify the file hhpfiles/Data.hpp adding the parameter in the class Data
-We add the reading of this parameter in the consturctor Data(const char * filename) in the file hhpfiles/Data.hpp
-We modify the constructors and operators associated with this class in the files hhpfiles/Data.hpp

-we modify the file Parameter.dat adding one line with the name of the paramter 
and one line with the value

-In the code we access to the new parameter with d."NameNewParameter"
**/
/*********************************************************************/





/*************************** Add Test case ***************************/
/**
Add a test case in the model  "NameModel".
-we modify the file init/Init"NameModel".cpp adding a new initiale values in the switch case (switch(d.nTest)).
In this case we give the initial value for all the variables

-If your test case have an exact solution. We write this exact solution in the file AnaSol/AnalyticSolution"NameModel".cpp
and include the call of this function in the switch case in the function SolFon"NameModel".

-If the test case have an initial data the function TraceSolfond in the file AnaSol/AnalyticSolution.cpp must be return "1" 
and "0" if the test case have not analytical solution
**/
/*********************************************************************/





/***************************** Add scheme ****************************/
/**
Add a scheme in the model  "NameModel".
-we create a function which compute your new flux and write this function in a file Flux/flux"NameModel".cpp

-The scheme are write on the form u_j^(n+1)=u_j^n+dt/Omega_j(Flux+source) with Omega_j the cell volume. 

-We add the call of your new fluxes in the functions ChoicefluxN or ChoicefluxE in the switch case of Flux.cpp

-we create or modify a function which compute the source term and write this function in a file Source/Source"NameModel".cpp

-The scheme are write on the form u_j^(n+1)=u_j^n+dt/Omega_j(Flux+source) with Omega_j the cell volume. 

-We add the call of your new fluxes in the functions ChoiceSource in the switch case of Source.cpp
**/
/*********************************************************************/






/****************************** Add Model ****************************/
/**
Add a new model "NewModel":

1) Physical Parameters:
 -Create a class Param"NewModel" whoch contains the physical parameters and the head of the functions to initialize these parameters
 (Param"NewModel".hpp in the repertory hhpParam)
- Create the functions which initialize the physicam parameters (Param"NewModel".cpp in the repertory ParamPhysic)
- Create the file Param"NewModel".dat which contains the input prameters in the repertory ParamPhysic

2) Variables
-Add in the class variables and vectorflux (hppfiles/Variables.hpp) the number of variables associated with your model

3) Initialization
- Create the file Init"NewModel".cpp in the repertory init and add the different initialization for the different test cases
- Add the in the function "ChoiceInit" (Init/initialisation.cpp) the function which init the variables of your model

4) Analytic Solution
- Create the file AnalyticSolutions"NewModel".cpp in the repertory AnaSol and add the different solutions for the different test cases
- Add the in the function "Solfond" (AnaSol/AnalyticSolutions.cpp) the function which gives your analytic solutions of your model
- (see add test case)

5) Boundary conditions
- Create the file BoundaryCondition"NewModel".cpp in the repertory BC and add the different boundary conditions for the different test cases
- Add the in the function "BoundaryCondition" (BC/BoundaryConditions.cpp) the function which gives your analytic solutions of your model

6) Diagnostic
 - Define in the functions "Diagnostic_Errors" and "Diagnostics quantities" the functions called by the model

7) Flux
- Create the files Flux"NewModel".cpp in the repertory Flux and add the different flux for the different schemes
- Add the in the functions "ChoiceFluxN" or "ChoiceFluxE"  (Flux/flux.cpp) the functions which compute the fluxes of your model

8) Source
- Create the files Source"NewModel".cpp in the repertory Source and add the different source for the different schemes
- Add the in the functions "ChoiceSource"  (Source/Source.cpp) the functions which compute the sources of your model

9)
-Add the plot of parameter for your model in the function "WriteDataToMain" (Functions/functions.cpp)
 **/
/*********************************************************************/






/************************ Current problem ****************************/
/**
-Test case 2 for P1 model does not converge
-Test case 1 and 2 for P1 model with matter does not gives the good results (wrong velocity of propragation)
s **/
/*********************************************************************/
