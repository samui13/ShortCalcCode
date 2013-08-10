#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 500
#define DT 10000
#define TD DT*10
#define TTD 10
#include "../libadv.h"

double Du;
double *dMM,*cMM,*lMM;

int main(int argc,char **argv){
  double *u,*un;
  double *temp;
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
  sprintf(FOLDER,"./data/2/");
  FOR(i,N){
    u[i] = 0.0;
    if(dx*i >= 0.4 && i*dx <= 0.6){
      u[i] = 1.0;
    }
  }
  OutPut1(0,u);
  int j;
  double AlphaDu = Du*dt/(dx*dx);
  double B[(int)(N/2)+1],F[N+1],D[(int)N+1];

  for(j = 0;j<=0;j++){
    B[0] = -(1+2*AlphaDu)/(-2*AlphaDu)*1;
    B[1] = -(1+2*AlphaDu)/(-1*AlphaDu)*B[0]-1;
    F[0] = u[0]/(-2*AlphaDu);
    F[1] = -(1+2*AlphaDu)/(-1*AlphaDu)*F[0]+u[1]/(-1*AlphaDu);

    D[N] = -(1+2*AlphaDu)/(-2*AlphaDu)*1+0*1;
    D[N-1] = -(1+2*AlphaDu)/(-1*AlphaDu)*D[N]-1;
    F[N] = -(1+2*AlphaDu)/(-2*AlphaDu)*0+0*0+u[N]/(-2*AlphaDu);
    F[N-1] = -(1+2*AlphaDu)/(-1*AlphaDu)*F[N]+u[N-1]/(-1*AlphaDu);

    for(i = 2;i<=N/2-1; i++){
      B[i] = -(1+2*AlphaDu)/(-1*AlphaDu)*B[i-1]-1*B[i-2];
      F[i] = -(1+2*AlphaDu)/(-1*AlphaDu)*F[i-1]-1*F[i-2]+u[i]/(-1*AlphaDu);
    }
    for(i = 2;i <= N/2; i++){
      D[N-i] = -(1+2*AlphaDu)/(-1*AlphaDu)*D[N-i+1]-1*D[N-i+2];
      F[N-i] = -(1+2*AlphaDu)/(-1*AlphaDu)*F[N-i+1]-1*F[N-i+2]+u[N-i]/(-1*AlphaDu);
    }

    un[0] = ((F[N/2]-F[N/2-2])*D[N/2+1]+(F[N/2-1]-F[N/2+1])*D[N/2])/(B[N/2-2]*D[N/2+1]-D[N/2]*B[N/2-1]);
    un[N] = ((F[N/2-2]-F[N/2])*B[N/2-1]+(F[N/2+1]-F[N/2-1])*B[N/2-2])/(B[N/2-1]*D[N/2]-D[N/2+1]*B[N/2-2]);

    for(i = 1; i<= N/2-1;i++){
      un[i] = B[i-1]*un[0]+F[i-1];
    }
    for(i = 1;i<=N/2;i++){
      un[N-i] = D[N-i+1]*un[N]+F[N-i+1];
    }
    temp = u;
    u = un;
    un = temp;
  }
  OutPut1(1,u);
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
