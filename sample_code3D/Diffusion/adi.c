// -*- coding:utf-8 -*-
/*
  Diffusion 3D
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define Nx 50
#define Ny 50
#define Nz 50

#define DT 8000
#define TD (DT)
#define TTD INT_MAX

#include "../../libadv3.h"

double Du;
double u1[Nz+1][Ny+1][Nx+1],u2[Nz+1][Ny+1][Nx+1];

void InitialValue(double ***u);
void InitialPara();
void OutPara();
void Calc(double ***u,double ***un);
void End();

double c_ux[Nx+1],d_ux[Nx+1],l_ux[Nx+1];
double c_uy[Ny+1],d_uy[Ny+1],l_uy[Ny+1];
double c_uz[Nz+1],d_uz[Nz+1],l_uz[Nz+1];
int main(int argc,char **argv){
  int i,j,k;
  double ***u,***un;
  double ***v,***vn;
  double ***uo;
  double epL1;
  u = mallmatrix();
  un = mallmatrix();
  uo = mallmatrix();

  InitialPara();
  InitialValue(u);

  Initial_LUx(Du,c_ux,d_ux,l_ux);
  Initial_LUy(Du,c_uy,d_uy,l_uy);
  Initial_LUz(Du,c_uz,d_uz,l_uz);

  OutPara();
  FOR(i,TTD){
    OutPut1(i,u);
    cp1(u,uo);
    FOR(j,TD/2){
      Calc(u,un);
      if(isnan(un[0][0][0]))exit(1);
      Calc(un,u);
      if(isnan(u[0][0][0]))exit(1);
    }
    
    epL1 = epsL1(u,uo);
    printf("L1 = %.15lf\t",epL1);
    if(epL1<=1/(double)EPS && i > 60){
      OutPut1(i+1,u);
      END();
      exit(1);
    }
  }
  freematrix(u);
  freematrix(un);
  freematrix(uo);
}

void OutPara(){
  char path[256];
  FILE *fp;
  sprintf(path,"%s/parameter.txt",FOLDER);
  if((fp = fopen(path,"wt")) == NULL){
    printf("Directory not found\n");
    exit(1);
  }
  FSFPRT(fp,"Du = ",Du);

  FSFPRT(fp,"Lx = ",Lx);
  FSFPRT(fp,"Ly = ",Ly);
  FSFPRT(fp,"Lz = ",Lz);

  FSFPRT(fp,"dx = ",dx);
  FSFPRT(fp,"dy = ",dy);
  FSFPRT(fp,"dz = ",dz);
  FSFPRT(fp,"dt = ",dt);

  FSPRT(fp,"Start Time:");
  //fprintf(fp,"%d",time(NULL));
  fclose(fp);
}
void InitialValue(double ***u){
  int i,j,k;
  FOR(k,Nz){
    FOR(j,Ny){
      FOR(i,Nx){
	u[k][j][i] = 0.0;
	if(i*dx >= 0.4 && i*dx <= 0.6 &&
	   j*dy >= 0.4 && j*dy <= 0.6 &&
	   k*dz >= 0.4 && k*dz <= 0.6){
	  u[k][j][i] = 1.0;
	}
      }
    }
  }
}
void InitialPara(){
  Du = 0.00001;

  Lx = 1.0;
  Ly = 1.0;
  Lz = 1.0;


  dt = 1/(double)DT;
  dx = Lx/(double)Nx;
  dy = Ly/(double)Ny;
  dz = Lz/(double)Nz;
  
  dx2 = dx*dx;
  dy2 = dy*dy;
  dz2 = dz*dz;
  
  ddxx = 1/dx;
  ddyy = 1/dy;
  ddxx = 1/dz;

  ddxx2 = 1/dx2;
  ddyy2 = 1/dy2;
  ddzz2 = 1/dz2;

  sprintf(FOLDER,"./data/0/");
}

void Calc(double ***u,double ***un){
  int i,i1,i2,i3,i4;
  int j,j1,j2,j3,j4;
  int k,k1,k2,k3,k4;

  double u1t[Nx+1],u2t[Ny+1],u3t[Nz+1];
  
  double uxx,uyy,uzz;
  double ux,uy,uz;
  FOR(k,Nz){
    k1 = K1(k);
    k2 = K2(k);
    k3 = K3(k);
    k4 = K4(k);
    FOR(j,Ny){
      j1 = J1(j);
      j2 = J2(j);
      j3 = J3(j);
      j4 = J4(j);
      FOR(i,Nx){
	i1 = I1(i);
	i2 = I2(i);
	i3 = I3(i);
	i4 = I4(i);
	uxx = DIFFUX(u,i2,i,i1,j,k);
	uyy = DIFFUY(u,i,j2,j,j1,k);
	uzz = DIFFUZ(u,i,j,k2,k,k1);
	u1t[i] = u[k][j][i]+(dt)*(Du*(0.5*uxx+uyy+uzz));
      }

      LUx(u1t,c_ux,d_ux,l_ux);

      FOR(i,Nx){
	u1[k][j][i] = u1t[i];
      }
      
    }
    
  }
  
  
  FOR(k,Nz){
    k1 = K1(k);
    k2 = K2(k);
    k3 = K3(k);
    k4 = K4(k);
    FOR(i,Nx){
      i1 = I1(i);
      i2 = I2(i);
      i3 = I3(i);
      i4 = I4(i);
      FOR(j,Ny){
	j1 = J1(j);
	j2 = J2(j);
	j3 = J3(j);
	j4 = J4(j);
	uxx = 0.5*(DIFFUX(u1,i2,i,i1,j,k)+DIFFUX(u,i2,i,1,j,k));
	uyy = DIFFUY(u,i,j2,j,j1,k);
	uzz = DIFFUZ(u,i,j,k2,k,k1);
	u2t[j] = u[k][j][i]+dt*(Du*(uxx+0.5*uyy+uzz));
	
      }
      LUy(u2t,c_uy,d_uy,l_uy);
      FOR(j,Ny){
	u2[k][j][i] = u2t[j];
      }
    }
  }

  FOR(i,Nx){
    i1 = I1(i);
    i2 = I2(i);
    i3 = I3(i);
    i4 = I4(i);
    FOR(j,Ny){
      j1 = J1(j);
      j2 = J2(j);
      j3 = J3(j);
      j4 = J4(j);
      FOR(k,Nz){
	k1 = K1(k);
	k2 = K2(k);
	k3 = K3(k);
	k4 = K4(k);

	uxx = 0.5*(DIFFUX(u1,i2,i,i1,j,k)+DIFFUX(u,i2,i,i1,j,k));
	uyy = 0.5*(DIFFUY(u2,i,j2,j,j1,k)+DIFFUY(u,i,j2,j,j1,k));
	uzz = DIFFUZ(u,i,j,k2,k,k1);
	u3t[k] = u[k][j][i]+dt*(Du*(uxx+uyy+0.5*uzz));
      }
      LUz(u3t,c_uz,d_uz,l_uz);
      FOR(k,Nz){
	un[k][j][i] = u3t[k];
      }
    }
  }
  
}


