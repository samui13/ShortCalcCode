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
  FOR(i,N){
    u[i] = 0.0;
    if(dx*i >= 0.4 && i*dx <= 0.6){
      u[i] = 1.0;
    }
  }
  int j;
  double AlphaDu = Du*dt/dx2;
  double B[(int)(N/2)+1],F[N],D[(int)N/2+1];
  B[0] = -(1+2*AlphaDu)/(-2*AlphaDu)*1;
  B[1] = -(1+2*AlphaDu)/(-1*AlphaDu)*B[0] + 1*1;
  F[0] = -(1+2*AlphaDu)/(-2*AlphaDu)*0-1*(0)*0+u[0]/(-2*AlphaDu);
  F[1] = -(1+2*AlphaDu)/(-1*AlphaDu)*F[0]+-1*(1+2*AlphaDu)/(-1*AlphaDu)+u[1]/(-1*AlphaDu);
  D[N] = -(1+2*AlphaDu)/(-2*AlphaDu)*1+0*1;
  D[N-1] = -(1+2*AlphaDu)/(-1*AlphaDu)*1+1*0;
  F[N] = -(1+2*AlphaDu)/(-2*AlphaDu)*0+0*0+u[N]/(-2*AlphaDu);
  F[N-1] = u[N-1]/(-AlphaDu);

  for(i = 2,i<=Nx-1; i++){
    B[i] = 
  }

  
  /*

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
  */
  return 0;
}
