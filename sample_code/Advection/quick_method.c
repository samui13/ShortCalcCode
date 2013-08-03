// -*- coding: utf-8 -*-
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define N 100
#define DT 10000
#define TD (DT*10)
#define TTD 10
#include "../../libadv.h"
double c;
void Calc(double *u,double *un);
void InitValue(double *u);
void InitialPara();
int main(int argc, char **argv){
  double *u,*un;
  int i,j;
  u = mallmatrix();
  un = mallmatrix();

  InitialPara();
  InitValue(u);  
  for( i = 0; i <= TTD; i++){
    OutPut1(i,u);
    printf("L1 = %.15lf\n",L1(u));
    for( j = 0; j < TD/2; j++){
      Calc(u,un);
      Calc(un,u);
    }
  }
  freematrix(u);
  freematrix(un);
  return 0;
}
void InitialPara(){
  Lx = 1.0;
  c = 0.01;
  dx = Lx/(double)N;
  dt = 1/(double)DT;
  sprintf(FOLDER,"./data/1/");
  dx2 = dx*dx;
  ddxx = 1/dx;
  ddxx2 = 1/dx2;
}
void InitValue(double *u){
  int i;
  FOR(i,Nx){
    u[i] = 0.0;
    if(i*dx >= Lx-0.2 && i*dx <= Lx-0.05){
      u[i] = 1.0;
    }
  }
}
void Calc(double *u,double *un){
  int i,i1,i2,i3,i4;
  FOR(i,Nx){
    i1 = I1(i);
    i2 = I2(i);
    i3 = I3(i);
    i4 = I4(i);
    un[i] = u[i]+dt*(c*QUICK(-c,u,i4,i2,i,i1,i3));
  }
  if(isnan(un[0])){
    exit(1);
    END();
  }
}
