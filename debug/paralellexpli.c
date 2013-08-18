// -*- coding: utf-8 -*-
// Explicit method and one factor(variable u)
// Pararell version
// Compile can gcc -fopenmp ./test.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define N 500
#define DT 70000
#define TD (DT*10)
#define TTD 10
#include "../libadv.h"

double Du;

void Calc(double *u,double *un,int first,int end);
void InitValue(double *u);
void InitialPara();
int main(int argc, char **argv){
  double *u,*un;
  double *temp;
  int i,j;
  u = mallmatrix();
  un = mallmatrix();

  InitialPara();
  InitValue(u);  
  for( i = 0; i <= TTD; i++){
    OutPut1(i,u);
    printf("L1 = %.15lf\n",L1(u));
    for( j = 0; j < TD; j++){
#ifndef _OPENMP
      Calc(u,un,0,N);
#else      
#pragma omp parallel
#pragma omp sections
      {
#pragma omp section
	{
	  Calc(u,un,0,250);
	}
#pragma omp section
	{
	  Calc(u,un,251,500);
	}
      }
#endif
      
      temp = un;
      un = u;
      u = temp;
    }
  }
  freematrix(u);
  freematrix(un);
  return 0;
}
void InitialPara(){
  Lx = 1.0;
  Du = 0.1;
  dx = Lx/(double)N;
  dt = 1/(double)DT;
  sprintf(FOLDER,"./data/10/");
  dx2 = dx*dx;
  ddxx = 1/dx2;
}
void InitValue(double *u){
  int i;
  FOR(i,Nx){
    u[i] = 0.0;
    if(i*dx >= 0.4 && i*dx <= 0.6){
      u[i] = 1.0;
    }
  }
}
void Calc(double *u,double *un,int first,int end){
  int i,i1,i2;
  //for(i = first; i<= end;i++){
  REP(i,first,end){
    i1 = I1(i);
    i2 = I2(i);
    un[i] = u[i] + dt*Du*ddxx*(u[i1]-2*u[i]+u[i2]);
  }
  if(isnan(un[0])){
    exit(1);
    END();
  }
}
