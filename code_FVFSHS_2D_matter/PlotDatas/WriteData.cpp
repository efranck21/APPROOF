#include <iostream>
#include <cassert>
#include "Data.hpp"
#include "WriteData.hpp"
#include "Functions.hpp"
#include "AnalyticSolutions.hpp"
#include "Diagnostics.hpp"
#include "Flux.hpp"

using namespace std;

/** File which contains the functions for the plot **/


void SaveData(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,int ntAnim){
  /** Function which find the numbers of the variables plotted **/
  int TabNumVar[300];
  int NbVarWrite;
  int c=0,j=0;
  char *copie = new char[4];
  for(int i=0;i<300;i++){
    if(d.NumVar[i]=='E') break;
    else if(d.NumVar[i]==',') {
      copie[c]='\0';
      TabNumVar[j]=atoi(copie);
      j++;
      c=0;}
    
    else {
      copie[c]=d.NumVar[i];
      c++;}
      
  }
    delete [] copie;
    NbVarWrite=j;
    
    for(j=0;j<NbVarWrite;j++){
      if(TabNumVar[j] > v.nbvar) {
	cout<<"the variable that you want plot does not exist"<<endl;
	exit(1);
      } 
      if(d.dimsave==2){
	if(d.Typewrite=='N'){
	  SaveOneData(d,Mh,v,Param,TabNumVar[j],ntAnim);
	}
	if(d.Typewrite=='C'){
	  SaveOneDatacentercell(d,Mh,v,Param,TabNumVar[j],ntAnim);
	}
      }
      if(d.dimsave==1 && d.Ny == 1){
	SaveOneData1D(d,Mh,v,Param,TabNumVar[j],ntAnim);
	
      }
    }

  
  if(Param.Model == 5) {
    SaveInterEner(d,Mh,v,Param,ntAnim); 
  }

}


void SaveOneData(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,int Num,int ntAnim){
 /** Function which plot in 2D the variable "num" using the average in the cell**/
  char string[255];
  FILE *fileresult=NULL;
  int j,n;
  double rhoNum=0,rhoExact=0;
 int in=0;
  
  if(d.Anim=='y'){ sprintf(string,"%s%s%d-%04d.%s","./DATA/","var",Num,ntAnim,d.suffixe);}

  else { sprintf(string,"%s%s%d.%s","./DATA/","var",Num,d.suffixe);}
  fileresult=fopen(string,"w");
  for(j=0;j<Mh.nc;j++){
    if(Mh.cells[j].lab!=-1){
    if(Num==0 || strcmp(d.Typemodel,"Euler")){
      
      if(Mh.nbnodelocal==4)
	{n=Mh(j,0);
	  fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,v.var[Num][j]);
	n=Mh(j,1);
	fprintf(fileresult,"%e %e %e\n\n",Mh.xr(n).x,Mh.xr(n).y,v.var[Num][j]);
	n=Mh(j,3);
	fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,v.var[Num][j]);
	n=Mh(j,2);
	fprintf(fileresult,"%e %e %e\n\n\n",Mh.xr(n).x,Mh.xr(n).y,v.var[Num][j]);
       
	
	}
      if(Mh.nbnodelocal==3)
	{n=Mh(j,0);
	  fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,v.var[Num][j]);
	n=Mh(j,1);
	fprintf(fileresult,"%e %e %e\n\n",Mh.xr(n).x,Mh.xr(n).y,v.var[Num][j]);
	n=Mh(j,2);
	fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,v.var[Num][j]);
	n=Mh(j,2);
	fprintf(fileresult,"%e %e %e\n\n\n",Mh.xr(n).x,Mh.xr(n).y,v.var[Num][j]);       
	
	}
    }
      
    if(Num!=0 && (!strcmp(d.Typemodel,"Euler"))){
      rhoNum=v.var[0][j];
      rhoExact=SolFonEulerAverage(d,Mh,Param,Mh.xj(j),0,Mh.xj(in),j,0);
      if(Mh.nbnodelocal==4)
	{n=Mh(j,0);
	  fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,v.var[Num][j]/rhoNum);
	n=Mh(j,1);
	fprintf(fileresult,"%e %e %e\n\n",Mh.xr(n).x,Mh.xr(n).y,v.var[Num][j]/rhoNum);
	n=Mh(j,3);
	fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,v.var[Num][j]/rhoNum);
	n=Mh(j,2);
	fprintf(fileresult,"%e %e %e\n\n\n",Mh.xr(n).x,Mh.xr(n).y,v.var[Num][j]/rhoNum);
       
	
	}
      if(Mh.nbnodelocal==3)
	{n=Mh(j,0);
	  fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,v.var[Num][j]/rhoNum);
	n=Mh(j,1);
	fprintf(fileresult,"%e %e %e\n\n",Mh.xr(n).x,Mh.xr(n).y,v.var[Num][j]/rhoNum);
	n=Mh(j,3);
	fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,v.var[Num][j]/rhoNum);
	n=Mh(j,2);
	fprintf(fileresult,"%e %e %e\n\n\n",Mh.xr(n).x,Mh.xr(n).y,v.var[Num][j]/rhoNum);       
	
	}
    }
        }

  }
  
  fclose(fileresult);
}



void SaveOneData1D(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,int Num,int ntAnim){
/** Function which plot in 1D the variable "num" using the average in the cell **/
  char string[255];
  FILE *fileresult=NULL;
  int j;
  double rho=0;

  if(d.Anim=='y'){ sprintf(string,"%s%s%d-%04d.%s","./DATA/","var",Num,ntAnim,d.suffixe);}
  else { sprintf(string,"%s%s%d.%s","./DATA/","var",Num,d.suffixe);}
  fileresult=fopen(string,"w");
  
  for(j=0;j<Mh.nc;j++){
    if(Mh.cells[j].lab!=-1){
   
     if(Num!=0 && Param.Model == 5){
       rho=v.var[0][j];
       fprintf(fileresult,"%e %e\n",Mh.xj(j).x,v.var[Num][j]/rho);
     }
     else{
       fprintf(fileresult,"%e %e\n",Mh.xj(j).x,v.var[Num][j]);
     }
    }
  }

  fclose(fileresult);
}



void SaveOneDatacentercell(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,int Num,int ntAnim){
/** Function which plot in 2D the variable "num" using the center of the cell **/
  char string[255];
  FILE *fileresult=NULL;
  int j;
  double rho=0;

  
  if(d.Anim=='y'){ sprintf(string,"%s%s%d-%04d.%s","./DATA/","var",Num,ntAnim,d.suffixe);}
  else { sprintf(string,"%s%s%d.%s","./DATA/","var",Num,d.suffixe);}
  fileresult=fopen(string,"w");
  for(j=0;j<Mh.nc;j++){
    if(Mh.cells[j].lab!=-1){
    if(Num==0 || Param.Model != 5){
	  fprintf(fileresult,"%e %e %e\n",Mh.xj(j).x,Mh.xj(j).y,v.var[Num][j]);	
    }
    if(Num!=0 && Param.Model == 5){
      rho=v.var[0][j];
      fprintf(fileresult,"%e %e %e\n",Mh.xj(j).x,Mh.xj(j).y,v.var[Num][j]/rho);	
    }
  }
  }
  fclose(fileresult);
}


void SaveRestart(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param,int nstep, double time){
/** Function which plot in 2D the variable "num" using the center of the cell **/
  char string[255];
  FILE *fileresult=NULL;
  int j,i;

  for(i=0;i<v.nbvar;i++){

  if(d.Tf== time){ sprintf(string,"%s%s%d.%s","./DATA/","var",i,"restartfinal");}
  else {sprintf(string,"%s%s%d.%s%s%f","./DATA/","var",i,d.suffixe,"restart",time);}
  fileresult=fopen(string,"w");

  fprintf(fileresult," %d %e %e %e\n",nstep,time,0.0,0.0);
  
  for(j=0;j<Mh.nc;j++){
    if(Mh.cells[j].lab!=-1){
      fprintf(fileresult," %d %e\n",j,v.var[i][j]);	
    }
  }
  fclose(fileresult);
}
}

void LoadRestart(Data & d,Mesh & Mh,variable & v,ParamPhysic & Param, double & InitTime){
/** Function which plot in 2D the variable "num" using the center of the cell **/
  char filename[255];
  char string[2555];
  FILE *f=NULL;
  int j,i,ncell=0;
  double value=0;
  int iteration;

  for(i=0;i<v.nbvar;i++){
      sprintf(filename,"%s%s%d.%s","./DATA/","var",i,"restart");
      f = fopen(filename,"rb");
     if(!f){cout << "Erreur d'ouverture du fichier restart "<<filename<<endl;exit(1);}
      cout <<"Lecture du fichier restart\"" << filename << "\""<<endl;

      fscanf(f,"%d",&iteration);
      fscanf(f,"%lf\n",&InitTime);
       fgets(string,255,f);
   
      for(j=0;j<Mh.nc;j++){
	if(Mh.cells[j].lab!=-1){
	  fscanf(f,"%d",&ncell);
	  fseek(f, 1, SEEK_CUR);	
	  fscanf(f,"%lf",&value);
	  v.var[i][ncell]=value;
	   fgets(string,255,f);
	}
      }
  fclose(f);
   }
}

  

void SaveMesh(Data & d,Mesh & Mh){
/** Function which plot the mesh **/
  char string[255];
  FILE *fileresult=NULL;
  int j,r,n;
  
  sprintf(string,"%s%s.%s","./DATA/","Mesh",d.suffixe);
  fileresult=fopen(string,"w");

  for(j=0;j<Mh.nc;j++){
    //  if(Mh.cells[j].lab!=-1){
    for(r=0;r<Mh.nbnodelocal;r++){
      n=Mh(j,r);
      fprintf(fileresult,"%e %e\n",Mh.xr(n).x,Mh.xr(n).y);
    }
    fprintf(fileresult,"\n");
    // }
  }
  fclose(fileresult);
}


void SaveFonction2D(Data & d,Mesh & Mh,ParamPhysic & Param,float t,Vertex c,int Num){
/** Function which plot in 2D the exact solution "num" **/
  char string[255];
  FILE *fileresult=NULL;
  int j,n;


  
  sprintf(string,"%s%s%d","./DATA/","Solfond",Num);
  fileresult=fopen(string,"w");
  for(j=0;j<Mh.nc;j++){
     if(Mh.cells[j].lab!=-1){
      if(Mh.nbnodelocal==4)
	{n=Mh(j,0);
	  fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,SolFon(d,Param,Mh.xj(j),t,c,Num));
	n=Mh(j,1);
	fprintf(fileresult,"%e %e %e\n\n",Mh.xr(n).x,Mh.xr(n).y,SolFon(d,Param,Mh.xj(j),t,c,Num));
	n=Mh(j,3);
	fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,SolFon(d,Param,Mh.xj(j),t,c,Num));
	n=Mh(j,2);
	fprintf(fileresult,"%e %e %e\n\n\n",Mh.xr(n).x,Mh.xr(n).y,SolFon(d,Param,Mh.xj(j),t,c,Num));
       
	
	}
      if(Mh.nbnodelocal==3)
	{n=Mh(j,0);
	  fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,SolFon(d,Param,Mh.xj(j),t,c,Num));
	n=Mh(j,1);
	fprintf(fileresult,"%e %e %e\n\n",Mh.xr(n).x,Mh.xr(n).y,SolFon(d,Param,Mh.xj(j),t,c,Num));
	n=Mh(j,2);
	fprintf(fileresult,"%e %e %e\n",Mh.xr(n).x,Mh.xr(n).y,SolFon(d,Param,Mh.xj(j),t,c,Num));
	n=Mh(j,2);
	fprintf(fileresult,"%e %e %e\n\n\n",Mh.xr(n).x,Mh.xr(n).y,SolFon(d,Param,Mh.xj(j),t,c,Num));
       
	
	}

      

       }
  }
  fclose(fileresult);
}





void SaveFonction1D(Data & d,Mesh & Mh,ParamPhysic & Param,float t,Vertex c,int Num){
/** Function which plot in 1D the exact solution "num" **/
  char string[255];
  FILE *fileresult=NULL;
  int j;


  
  sprintf(string,"%s%s%d","./DATA/","Solfond",Num);
  fileresult=fopen(string,"w");

  for(j=0;j<d.Nx;j++){
    if(Mh.cells[j].lab!=-1){
      fprintf(fileresult,"%e %e\n",Mh.xj(j).x,SolFon(d,Param,Mh.xj(j),t,c,Num)) ;
    } 
  }

  fclose(fileresult);
}

void SaveFonction(Data & d,Mesh & Mh,ParamPhysic & Param,float t,Vertex c){
/** Function which find the numbers of the exact solutions plotted **/
 int TabNumVar[300];
  int NbVarWrite;
  
  
  int cc=0,j=0;
  char *copie = new char[4];
  for(int i=0;i<300;i++){
    if(d.NumVar[i]=='E') break;
    else if(d.NumVar[i]==',') {
      copie[cc]='\0';
      TabNumVar[j]=atoi(copie);
      j++;
      cc=0;}
    
    else {
      copie[cc]=d.NumVar[i];
      cc++;}
      
  }
    delete [] copie;
    NbVarWrite=j;
    
    for(j=0;j<NbVarWrite;j++){
      if(d.dimsave==2){
	if(d.dimsave==2 && Param.Model != 5){
	  SaveFonction2D(d,Mh,Param,t,c,TabNumVar[j]);
	}
	if(d.dimsave==2 && Param.Model == 5){
	  SaveFonction2DEuler(d,Mh,Param,t,c,TabNumVar[j]);
	}
      }
      if(d.dimsave==1 && Param.Model != 5){
	SaveFonction1D(d,Mh,Param,t,c,TabNumVar[j]);
	
      }
    }

}
