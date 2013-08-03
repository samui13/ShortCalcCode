// -*- coding: utf-8 -*-
// Explicit method and one factor(variable u)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 500
#define DT 70000
#define TD (DT*10)
#define TTD 10
#include "../../libadv.h"

double Du,Dv;

void Calc(double *u,double *v,
	  double *un,double *vn);
void InitValue(double *u,double *v);
void InitialPara();
int main(int argc, char **argv){
  double *u,*un;
  double *v,vn;
  int i,j;
  u = mallmatrix();
  un = mallmatrix();
  v = mallmatrix();
  vn = mallmatrix();

  InitialPara();
  InitValue(u,v);  
  for( i = 0; i <= TTD; i++){
    OutPut2(i,u,v);
    printf("L1 = %.15lf\n",L1(u));
    for( j = 0; j < TD/2; j++){
      Calc(u,v,un,vn);
      Calc(un,vn,u,v);
    }
  }
  freematrix(u);
  freematrix(un);
  freematrix(v);
  freematrix(vn);
  return 0;
}
void InitialPara(){
  Lx = 1.0;
  Du = 0.1;
  Dv = 0.01;
  dx = Lx/(double)N;
  dt = 1/(double)DT;
  sprintf(FOLDER,"./data/0/");
  dx2 = dx*dx;
  ddxx = 1/dx2;
}
void InitValue(double *u,double *v){
  int i;
  FOR(i,Nx){
    u[i] = 0.0;
    v[i] = 0.0;
    if(i*dx >= 0.4 && i*dx <= 0.6){
      u[i] = 1.0;
      v[i] = 1.0;
    }
  }
}
void Calc(double *u,double *v,double *un,double *vn){
  int i,i1,i2;
  FOR(i,Nx){
    i1 = I1(i);
    i2 = I2(i);
    un[i] = u[i] + dt*Du*ddxx*(u[i1]-2*u[i]+u[i2]);
    vn[i] = v[i] + dt*Dv*DIFFU(u,i2,i,i1);
  }
  if(isnan(un[0])){
    exit(1);
    END();
  }
}
