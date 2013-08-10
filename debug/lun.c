#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 100
#define DT 10000
#define TD DT*10
#define TTD 10
#include "../libadv.h"

double Du;
double *dMM,*cMM,*lMM;

int main(int argc,char **argv){
  double *u,*un;
  int i;
  u = mallmatrix();
  un = mallmatrix();
  
  cMM = (double *)malloc(sizeof(double *)*(N+1));
  dMM = (double *)malloc(sizeof(double *)*(N+1));
  lMM = (double *)malloc(sizeof(double *)*(N+1));
  Lx = 1.0;
  Du = 0.1;
  dx = Lx/(double)N;
  dt = 1/(double)DT;
  sprintf(FOLDER,"./data/1/");

  InitialLUx(Du,cMM,dMM,lMM);
  FOR(i,N){
    u[i] = 0.0;
    if(dx*i >= 0.4 && i*dx <= 0.6){
      u[i] = 1.0;
    }
    printf("%d %lf %lf %lf\n",i,cMM[i],dMM[i],lMM[i]);
  }
  OutPut1(0,u);
  FOR(i,TD/4){
    LUx(u,cMM,dMM,lMM);
  }
  OutPut1(1,u);

  return 0;
}
