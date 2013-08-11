// -*- coding: utf-8 -*-
// Explicit method and one factor(variable u)
//
// u_t = D_u u''-uvv+a(1-u)
// v_t = D_v v''+uvv-bv
// 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 500
#define DT 70000
#define TD (DT)
#define TTD 18000
//#define __BOUM__ 0
#include "../../libadv.h"

double Du,Dv,b,a;

void Calc(double *u,double *v,
	  double *un,double *vn);
void InitValue(double *u,double *v);
void InitialPara();
int main(int argc, char **argv){
  double *u,*un;
  double *v,*vn;
  double *uo;
  int i,j;
  u = mallmatrix();
  un = mallmatrix();
  v = mallmatrix();
  vn = mallmatrix();
  
  uo = mallmatrix();


  InitialPara();
  InitValue(u,v);  
  cp1(u,uo);
  for( i = 0; i <= TTD; i++){
    OutPut2(i,u,v);
    printf("L1 = %.15lf, epsL1 = %.15lf\n",L1(u),epsL1(u,uo));
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
  Du = 2.0*pow(10.0,-5.0);
  Dv = 1.0*pow(10.0,-5.0);
  a = 0.04;
  b = 0.10075;
  printf("%lf,%lf\n",Du,Dv);
  dx = Lx/(double)N;
  dt = 1/(double)DT;
  sprintf(FOLDER,"./data/0/");
  dx2 = dx*dx;
  ddxx = 1/dx;
  ddxx2 = 1/dx2;
}
void InitValue(double *u,double *v){
  int i;
  srand(RSEED);
  
  FOR(i,Nx){
    u[i] = 1.0;
    v[i] = 0.0;
  }
  //InitialValue_GaussDistribution(v,sqrt(0.002),Lx/2);
  InitialValue_Palse(v,0.0,0.35);
  /*
  FOR(i,Nx){
    v[i]/=9.0;
  }
  */

}

void Calc(double *u,double *v,double *un,double *vn){
  int i,i1,i2;
  FOR(i,Nx){
    i1 = I1(i);
    i2 = I2(i);
    un[i] = u[i] + dt*(Du*DIFFU(u,i1,i,i2)-u[i]*v[i]*v[i]+a*(1-u[i]));
    vn[i] = v[i] + dt*(Dv*DIFFU(v,i2,i,i1)+u[i]*v[i]*v[i]-b*v[i]);
  }
  if(isnan(un[0])){
    exit(1);
    END();
  }
}
