// -*- coding: utf-8 -*-
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define Nx 50
#define Ny 50
#define DT 10000
#define TD (DT*10)
#define TTD 10
#include "../../libadv2.h"

void Calc(double **u,double **un);
void InitValue(double **u);
void InitialPara();
void OutPara();
double c;

int main(int argc, char **argv){
  double **u,**un;
  int i,j;
  u = mallmatrix();
  un = mallmatrix();
  
  InitialPara();
  InitValue(u);
  OutPara();
  FOR(i,TTD){
    OutPut1(i,u);
    printf("L1 = %.15lf\n",L1(u));
    FOR(j , TD/2){
      Calc(u,un);
      Calc(un,u);
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

  c = 0.01;
  sprintf(FOLDER,"./data/0/");
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
  FSFPRT(fp,"c = ",c);
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
      un[j][i] = u[j][i]+dt*(c*UPWX(-1*c,u,i4,i2,i,i1,i3,j)
			     +c*UPWY(-1*c,u,i,j4,j2,j,j1,j3));
    }
  }
  if(isnan(un[0][0])){
    END();
    exit(1);
  }
}
