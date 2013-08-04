// -*- coding: utf-8 -*-
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define Nx 50
#define Ny 50
#define DT 10000
#define TD (DT)
#define TTD 10
#include "../../libadv2.h"

void Calc(double **u,double **un);
void InitValue(double **u);
void InitialPara();
void OutPara();
double Du;

double c_ux[Nx+1],d_ux[Nx+1],l_ux[Nx+1];
double c_uy[Ny+1],d_uy[Ny+1],l_uy[Ny+1];

int main(int argc, char **argv){
  double **u,**un;
  int i,j;
  u = mallmatrix();
  un = mallmatrix();
  
  InitialPara();
  InitValue(u);
  
  Initial_LUx(Du,c_ux,d_ux,l_ux);
  Initial_LUy(Du,c_uy,d_uy,l_uy);

  OutPara();
  FOR(i,TTD){
    OutPut1(i,u);
    printf("L1 = %.15lf\n",L1(u));
    FOR(j , TD){
      Calc(u,un);
    }
  }
}

void InitValue(double **u){
  int i,j;
  FOR(j,Ny){
    FOR(i,Nx){
      u[j][i] = 0.0;
      if(j*dy >= 0.4 && j*dy <= 0.6 &&
	 i*dx >= 0.4 && i*dx <= 0.6){
	u[j][i] = 1.0;
      }
    }
  }
}
void InitialPara(){
  Lx= 1.0;
  Ly = 1.0;
  dx = Lx/(double)Nx;
  dy = Ly/(double)Ny;
  dx2 = dx*dx;
  dy2 = dy*dy;
  ddxx = 1/dx;
  ddyy = 1/dy;
  ddxx2 = 1/dx2;
  ddyy2 = 1/dy2;
  dt = 1/(double)DT;

  Du = 0.01;
  sprintf(FOLDER,"./data/1/");
}

void OutPara(){
  char path[256];
  FILE *fp;
  int i,j;
  sprintf(path,"%s/parameter.txt",FOLDER);
  if((fp = fopen(path,"wt")) == NULL ){
    printf("Directory not found \n");
    exit(1);
  }
  FSFPRT(fp,"Du = ",Du);
  FSFPRT(fp,"Lx = ",Lx);
  FSFPRT(fp,"Ly = ",Ly);
  FSFPRT(fp,"dx = ",dx);
  FSFPRT(fp,"dy = ",dy);
  FSFPRT(fp,"Nx = ",(double)Nx);
  FSFPRT(fp,"Ny = ",(double)Ny);
  FSFPRT(fp,"dt = ",dt);
  FSPRT(fp,"This file is \t");
  FSPRT(fp,__FILE__);
  fclose(fp);
}

void Calc(double **u,double **un){
  int i,i1,i2,i3,i4;
  int j,j1,j2,j3,j4;
  double uhx[Nx+1],uhy[Ny+1];

  FOR(j,Ny){
    j1 = J1(j);
    j2 = J2(j);
    j3 = J3(j);
    j4 = J4(j);
    FOR(i,Nx){
      i1 = I1(i);
      i2 = I2(i);
      i3 = I3(i);
      i4 = I4(i);
      uhx[i] = u[j][i]+dt*(0.5*Du*DIFFUX(u,i2,i,i1,j));
    }
    LUx(uhx,c_ux,d_ux,l_ux);
    FOR(i,Nx){
      u[j][i] = uhx[i];
    }
  }
  
  FOR(i,Nx){
    i1 = I1(i);
    i2 = I2(i);
    i3 = I3(i);
    i4 = I4(i);
    FOR(j,Ny){
      j1 = J1(j);
      j2 = J2(j);
      j3 = J3(j);
      j4 = J4(j);
      uhy[j] = u[j][i] + dt*(0.5*Du*DIFFUY(u,i,j2,j,j1));
    }    

    LUy(uhy,c_uy,d_uy,l_uy);
    FOR(j,Ny){
      u[j][i] = uhy[j];
    }
  }
  if(isnan(un[0][0])){
    END();
    exit(1);
  }
}
