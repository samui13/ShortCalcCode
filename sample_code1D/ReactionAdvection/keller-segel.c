// -*- coding: utf-8 -*-
// Explicit method and two factor(variable u,v)
//
// u_t = D_u u'' -b(uv')'
// v_t = D_v v'' +u -v
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 500
#define DT 70000
#define TD (DT*10)
#define TTD 10
#include "../../libadv.h"

double Du,b,Dv;

void Calc(double *u,double *v,
	  double *un,double *vn);
void InitValue(double *u,double *v);
void InitialPara();
int main(int argc, char **argv){
  double *u,*un;
  double *v,*vn;
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
  b = 0.5;
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
  }
  v[N/2] = 1.0;
}

void Calc(double *u,double *v,double *un,double *vn){
  int i,i1,i2,i3,i4;
  double vxx;
  double vx;
  FOR(i,Nx){
    i1 = I1(i);
    i2 = I2(i);
    i3 = I3(i);
    i4 = I4(i);
    vxx = DIFFU(v,i2,i,i1);
    vx = DIFFE(v,i2,i,i1);
    un[i] = u[i] + dt*(Du*DIFFU(u,i2,i,i1)+b*UPW(-b*vx,u,i4,i2,i,i1,i3)-b*u[i]*vxx);
    vn[i] = v[i] + dt*(Dv*vxx+u[i]-v[i]);
  }
  if(isnan(un[0])){
    exit(1);
    END();
  }
}
