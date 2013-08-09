// -*- coding: utf-8 -*-
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define Nx 50
#define Ny 50
#define DT 10000
#define TD (DT*10)
#define TTD 100
#include "../../libadv2.h"

void Calc(double **u,double **v,double **un,double **vn);
void InitValue(double **u,double **v);
void InitialPara();
void OutPara();
double Du,Dv,b,a;

int main(int argc, char **argv){
  double **u,**un;
  double **v,**vn;
  int i,j;
  u = mallmatrix();
  un = mallmatrix();
  
  v = mallmatrix();
  vn = mallmatrix();

  InitialPara();
  InitValue(u,v);
  OutPara();
  FOR(i,TTD){
    OutPut2(i,u,v);
    printf("L1 = %.15lf\n",L1(u));
    FOR(j , TD/2){
      Calc(u,v,un,vn);
      Calc(un,vn,u,v);
    }
  }
}

void InitValue(double **u,double **v){
  int i,j;
  srand(RSEED);
  FOR(j,Ny){
    FOR(i,Nx){
      u[j][i] = 1.0;
      v[j][i] = 0.0;
      if(j*dy >= 0.4 && j*dy <= 0.6 &&
	 i*dx >= 0.4 && i*dx <= 0.6){
	v[j][i] = 1.0;
      }
    }
  }
}
void InitialPara(){
  Lx = 1.0;
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

  Du = 2.0*pow(10.0,-5.0);
  Dv = 1.0*pow(10.0,-5.0);
  a = 0.025;
  b = 0.0542+a;

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

void Calc(double **u,double **v,double **un,double **vn){
  int i,i1,i2,i3,i4;
  int j,j1,j2,j3,j4;
  double uxx,vxx,uyy,vyy;
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
      uxx = DIFFUX(u,i2,i,i1,j);
      vxx = DIFFUX(v,i2,i,i1,j);
      uyy = DIFFUY(u,i,j2,j,j1);
      vyy = DIFFUY(v,i,j2,j,j1);
      un[j][i] = u[j][i]+dt*(Du*(uxx+uyy)-u[j][i]*v[j][i]*v[j][i]+a*(1-u[j][i]));
      vn[j][i] = v[j][i]+dt*(Dv*(vxx+vyy)+u[j][i]*v[j][i]*v[j][i]-b*v[j][i]);
    }
  }
  if(isnan(un[0][0])){
    END();
    exit(1);
  }
}
