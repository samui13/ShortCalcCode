#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define N 500
#define DT 70000
#define TD (DT/10)
#define TTD 1
#include "../../libadv.h"
double Du;
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
  Du = 0.1;
  dx = Lx/(double)N;
  dt = 1/(double)DT;
  sprintf(FOLDER,"./data/0/");
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
void Calc(double *u,double *un){
  int i,i1,i2;
  FOR(i,Nx){
    i1 = I1(i);
    i2 = I2(i);
    un[i] = u[i] + dt*Du*ddxx*(u[i1]-2*u[i]+u[i2]);
  }
  if(isnan(un[0])){
    exit(1);
    END();
  }
}
