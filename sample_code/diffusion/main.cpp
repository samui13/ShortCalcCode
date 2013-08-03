#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define N 500
#define DT 70000
#define TD (DT/10)
#define TTD 10
#include "../../libadv.h"
double Du;
void Calc(double u[N+1],double v[N+1],
	  double un[N+1],double vn[N+1]);
void InitValue(double u[N+1],double v[N+1]);
int main(int argc, char **argv){
  Lx = 1.0;
  Du = 0.1;
  dx = Lx/(double)N;
  dt = 1/(double)DT;
  sprintf(FOLDER,"./data/0/");
  dx2 = dx*dx;
  ddxx = 1/dx2;
  double u[N+1],v[N+1];
  double un[N+1],vn[N+1];
  int i,j;
  InitValue(u,v);
  
  for( i = 0; i <= TTD; i++){
    OutPut2(i,u,v);
    for( j = 0; j < TD/2; j++){
      Calc(u,v,un,vn);
      if(isnan(un[0]) || isnan(vn[0]) ) exit(1);
      Calc(un,vn,u,v);
      if(isnan(v[0]) || isnan(u[0])) exit(1);
    }
  }

  return 0;
}
void InitValue(double u[N+1],double v[N+1]){
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
void Calc(double u[N+1],double v[N+1],
	  double un[N+1],double vn[N+1]){
  int i,i1,i2;
  FOR(i,Nx){
    i1 = I1(i);
    i2 = I2(i);
    un[i] = u[i] + dt*Du*ddxx*(u[i1]-2*u[i]+u[i2]);
  }
}
