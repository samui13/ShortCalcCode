/*
Mimura-Tujikawa. 3D
u_t= Du*u''-b(uv')'+cu(1-u)
v_t=Dv*v''+fu-gv
Designed at 2013.5.13
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define Nx 72
#define Ny 72
#define Nz 40

#define DT 8192
//#define TD (DT/(DT))
#define TD (DT)
#define TTD INT_MAX

#include "./libadv3.h"
//#include "/home/ymnk/template/libadv3.h"
double delta,alpha,beta;
double Du,b,c,Dv,f,g;

double u1[Nz+1][Ny+1][Nx+1],u2[Nz+1][Ny+1][Nx+1];
double v1[Nz+1][Ny+1][Nx+1],v2[Nz+1][Ny+1][Nx+1];

void InitialValue(double ***u,double ***v);
void InitialPara();
void OutPara();
void Calc(double ***u,double ***un,
	  double ***v,double ***vn);
void End();

double c_ux[Nx+1],d_ux[Nx+1],l_ux[Nx+1];
double c_vx[Nx+1],d_vx[Nx+1],l_vx[Nx+1];
double c_uy[Ny+1],d_uy[Ny+1],l_uy[Ny+1];
double c_vy[Ny+1],d_vy[Ny+1],l_vy[Ny+1];
double c_uz[Nz+1],d_uz[Nz+1],l_uz[Nz+1];
double c_vz[Nz+1],d_vz[Nz+1],l_vz[Nz+1];
int main(int argc,char **argv){
  int i,j,k;
  double ***u,***un;
  double ***v,***vn;
  double ***uo;
  double epL1;
  u = mallmatrix();
  un = mallmatrix();
  v = mallmatrix();
  vn = mallmatrix();
  uo = mallmatrix();

  InitialPara();
  InitialValue(u,v);

  Initial_LUx(Du,c_ux,d_ux,l_ux);
  Initial_LUy(Du,c_uy,d_uy,l_uy);
  Initial_LUz(Du,c_uz,d_uz,l_uz);

  Initial_LUx(Dv,c_vx,d_vx,l_vx);
  Initial_LUy(Dv,c_vy,d_vy,l_vy);
  Initial_LUz(Dv,c_vz,d_vz,l_vz);

  OutPara();
  FOR(i,TTD){
    OutPut2(i,u,v);
    cp1(u,uo);
    FOR(j,TD/2){
      Calc(u,un,v,vn);
      if(isnan(un[0][0][0]) || isnan(vn[0][0][0]))exit(1);
      Calc(un,u,vn,v);
      if(isnan(u[0][0][0]) || isnan(v[0][0][0]))exit(1);
    }
    
    epL1 = epsL1(u,uo);
    printf("L1 = %.15lf\t",epL1);
    if(epL1<=1/(double)EPS && i > 60){
      OutPut2(i+1,u,v);
      END();
      exit(1);
    }
  }
  freematrix(u);
  freematrix(un);
  freematrix(v);
  freematrix(vn);
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
  FSFPRT(fp,"b = ",b);
  FSFPRT(fp,"c = ",c);
  FSFPRT(fp,"Dv = ",Dv);
  FSFPRT(fp,"f = ",f);
  FSFPRT(fp,"g = ",g);

  FSFPRT(fp,"alpha = ",alpha);
  FSFPRT(fp,"beta = ",beta);
  FSFPRT(fp,"delta = ",delta);

  FSFPRT(fp,"Lx = ",Lx);
  FSFPRT(fp,"Ly = ",Ly);
  FSFPRT(fp,"Lz = ",Lz);

  FSFPRT(fp,"dx = ",dx);
  FSFPRT(fp,"dy = ",dy);
  FSFPRT(fp,"dz = ",dz);
  FSFPRT(fp,"dt = ",dt);

  FSPRT(fp,"Start Time:");
  fprintf(fp,"%d",time(NULL));
  fclose(fp);
}
void InitialValue(double ***u,double ***v){
  int i,j,k;
  srand(RSEED);
  FOR(k,Nz){
    FOR(j,Ny){
      FOR(i,Nx){
	u[k][j][i] = 1.0+0.00001*RDOUBLE;
	v[k][j][i] = f/g;
      }
    }
  }
}
void InitialPara(){

  double b_c;
  Du = 1.0/16.0;
  c = 9.0/2.0;
  Dv = 16.0;
  f = 1.0;
  g = 32.0;

  delta = sqrt(sqrt((c*g)/(Du*Dv)));
  alpha = (sqrt(3.0)*sqrt(delta))/2.0;
  beta = sqrt(delta)/2.0;

  Lx = M_PI/(delta/2.0);
  Ly = M_PI/(delta/2.0);
  Lz = M_PI/(delta*sqrt(3.0)/2.0);

  b_c=(sqrt(Du*g)+sqrt(c*Dv))*(sqrt(Du*g)+sqrt(c*Dv))/f;         //p=1の時のb_criticalの最小値(+0.0001)
  b = b_c+10.0;
  dt = 1/(double)DT;
  dx = Lx/(double)Nx;
  dy = Ly/(double)Ny;
  dz = Lz/(double)Nz;
  
  dx2 = dx*dx;
  dy2 = dy*dy;
  dz2 = dz*dz;
  sprintf(FOLDER,"./data/5/");
}

void Calc(double ***u,double ***un,
	  double ***v,double ***vn){
  int i,i1,i2,i3,i4;
  int j,j1,j2,j3,j4;
  int k,k1,k2,k3,k4;

  double u1t[Nx+1],u2t[Ny+1],u3t[Nz+1];
  double v1t[Nx+1],v2t[Ny+1],v3t[Nz+1];

  double uxx,uyy,uzz,vxx,vyy,vzz;
  double v2xx,v2yy,v2zz;
  double v3xx,v3yy,v3zz;
  double ux,uy,uz;
  double vx,vy,vz;
  

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
	uxx = (u[k][j][i1]-2*u[k][j][i]+u[k][j][i2])/dx2;
	uyy = (u[k][j1][i]-2*u[k][j][i]+u[k][j2][i])/dy2;
	uzz = (u[k1][j][i]-2*u[k][j][i]+u[k2][j][i])/dz2;
	vxx = (v[k][j][i1]-2*v[k][j][i]+v[k][j][i2])/dx2;
	vyy = (v[k][j1][i]-2*v[k][j][i]+v[k][j2][i])/dy2;
	vzz = (v[k1][j][i]-2*v[k][j][i]+v[k2][j][i])/dz2;

	vx = (v[k][j][i1]-v[k][j][i2])/(2*dx);
	vy = (v[k][j1][i]-v[k][j2][i])/(2*dy);
	vz = (v[k1][j][i]-v[k2][j][i])/(2*dz);
	
	ux = (vx > 0 ? (2*u[k][j][i1]+3*u[k][j][i]-6*u[k][j][i2]+u[k][j][i4]):(-u[k][j][i3]+6*u[k][j][i1]-3*u[k][j][i]-2*u[k][j][i2]))/(6*dx);
	uy = (vy > 0 ? (2*u[k][j1][i]+3*u[k][j][i]-6*u[k][j2][i]+u[k][j4][i]): (-u[k][j3][i]+6*u[k][j1][i]-3*u[k][j][i]-2*u[k][j2][i]))/(6*dy);
	uz = (vz > 0 ? (2*u[k1][j][i]+3*u[k][j][i]-6*u[k2][j][i]+u[k4][j][i]): (-u[k3][j][i]+6*u[k1][j][i]-3*u[k][j][i]-2*u[k2][j][i]))/(6*dz);
	
	u1t[i] = u[k][j][i]+(dt)*(Du*(0.5*uxx+uyy+uzz)-b*u[k][j][i]*(vxx+vyy+vzz)-b*(vx*ux+vy*uy+vz*uz)+c*u[k][j][i]*(1-u[k][j][i]));
	v1t[i] = v[k][j][i]+(dt)*(Dv*((0.5*vxx)+vyy+vzz)+f*u[k][j][i]-g*v[k][j][i]);
      }
      /*
      LUx(u1t,c_ux,d_ux,l_ux);
      LUx(v1t,c_vx,d_vx,l_vx);
      */
      DualLUx(u1t,c_ux,d_ux,l_ux,
	      v1t,c_vx,d_vx,l_vx);
      
      FOR(i,Nx){
	u1[k][j][i] = u1t[i];
	v1[k][j][i] = v1t[i];
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
	uxx = 0.5*((u1[k][j][i1]-2*u1[k][j][i]+u1[k][j][i2])/(dx2)+(u[k][j][i1]-2*u[k][j][i]+u[k][j][i2])/(dx2));
	uyy = (u[k][j1][i]-2*u[k][j][i]+u[k][j2][i])/(dy2);
	uzz = (u[k1][j][i]-2*u[k][j][i]+u[k2][j][i])/(dz2);
	
	v2xx = 0.5*((v1[k][j][i1]-2*v1[k][j][i]+v1[k][j][i2])/(dx2)+(v[k][j][i1]-2*v[k][j][i]+v[k][j][i2])/(dx2));
	vxx = (v[k][j][i1]-2*v[k][j][i]+v[k][j][i2])/(dx2);
	vyy = v2yy = (v[k][j1][i]-2*v[k][j][i]+v[k][j2][i])/(dy2);
	vzz = v2zz = (v[k1][j][i]-2*v[k][j][i]+v[k2][j][i])/(dz2);

	vx = (v[k][j][i1]-v[k][j][i2])/(2*dx);
	vy = (v[k][j1][i]-v[k][j2][i])/(2*dy);
	vz = (v[k1][j][i]-v[k2][j][i])/(2*dz);
	
	ux = (vx > 0 ? (2*u[k][j][i1]+3*u[k][j][i]-6*u[k][j][i2]+u[k][j][i4]):(-u[k][j][i3]+6*u[k][j][i1]-3*u[k][j][i]-2*u[k][j][i2]))/(6*dx);
	uy = (vy > 0 ? (2*u[k][j1][i]+3*u[k][j][i]-6*u[k][j2][i]+u[k][j4][i]): (-u[k][j3][i]+6*u[k][j1][i]-3*u[k][j][i]-2*u[k][j2][i]))/(6*dy);
	uz = (vz > 0 ? (2*u[k1][j][i]+3*u[k][j][i]-6*u[k2][j][i]+u[k4][j][i]): (-u[k3][j][i]+6*u[k1][j][i]-3*u[k][j][i]-2*u[k2][j][i]))/(6*dz);

	u2t[j] = u[k][j][i]+dt*(Du*(uxx+0.5*uyy+uzz)-b*u[k][j][i]*(vxx+vyy+vzz)-b*(vx*ux+vy*uy+vz*uz)+c*u[k][j][i]*(1-u[k][j][i]));
	v2t[j] = v[k][j][i]+dt*(Dv*(v2xx+0.5*v2yy+vzz)+f*u[k][j][i]-g*v[k][j][i]);
	
      }
      /*
      LUy(u2t,c_uy,d_uy,l_uy);
      LUy(v2t,c_vy,d_vy,l_vy);
*/
      DualLUy(u2t,c_uy,d_uy,l_uy,
	      v2t,c_vy,d_vy,l_vy);

      FOR(j,Ny){
	u2[k][j][i] = u2t[j];
	v2[k][j][i] = v2t[j];
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

	uxx = 0.5*((u1[k][j][i1]-2*u1[k][j][i]+u1[k][j][i2])/(dx2)+(u[k][j][i1]-2*u[k][j][i]+u[k][j][i2])/(dx2));
	uyy = 0.5*((u2[k][j1][i]-2*u2[k][j][i]+u2[k][j2][i])/dy2+(u[k][j1][i]-2*u[k][j][i]+u[k][j2][i])/dy2);
	uzz = (u[k1][j][i]-2*u[k][j][i]+u[k2][j][i])/dz2;
	
	v3xx = 0.5*((v1[k][j][i1]-2*v1[k][j][i]+v1[k][j][i2])/(dx2)+(v[k][j][i1]-2*v[k][j][i]+v[k][j][i2])/(dx2));
	vxx =(v[k][j][i1]-2*v[k][j][i]+v[k][j][i2])/(dx2);
	v3yy = 0.5*((v2[k][j1][i]-2*v2[k][j][i]+v2[k][j2][i])/dy2+(v[k][j1][i]-2*v[k][j][i]+v[k][j2][i])/dy2);
	vyy = (v[k][j1][i]-2*v[k][j][i]+v[k][j2][i])/dy2;
	vzz = v3zz = (v[k1][j][i]-2*v[k][j][i]+v[k2][j][i])/dz2;

	vx = (v[k][j][i1]-v[k][j][i2])/(2*dx);
	vy = (v[k][j1][i]-v[k][j2][i])/(2*dy);
	vz = (v[k1][j][i]-v[k2][j][i])/(2*dz);

	ux = (vx > 0 ? (2*u[k][j][i1]+3*u[k][j][i]-6*u[k][j][i2]+u[k][j][i4]):(-u[k][j][i3]+6*u[k][j][i1]-3*u[k][j][i]-2*u[k][j][i2]))/(6*dx);
	uy = (vy > 0 ? (2*u[k][j1][i]+3*u[k][j][i]-6*u[k][j2][i]+u[k][j4][i]): (-u[k][j3][i]+6*u[k][j1][i]-3*u[k][j][i]-2*u[k][j2][i]))/(6*dy);
	uz = (vz > 0 ? (2*u[k1][j][i]+3*u[k][j][i]-6*u[k2][j][i]+u[k4][j][i]): (-u[k3][j][i]+6*u[k1][j][i]-3*u[k][j][i]-2*u[k2][j][i]))/(6*dz);

	u3t[k] = u[k][j][i]+dt*(Du*(uxx+uyy+0.5*uzz)-b*u[k][j][i]*(vxx+vyy+vzz)-b*(vx*ux+vy*uy+vz*uz)+c*u[k][j][i]*(1-u[k][j][i]));
	v3t[k] = v[k][j][i]+dt*(Dv*(v3xx+v3yy+0.5*vzz)+f*u[k][j][i]-g*v[k][j][i]);
      }
      /*
      LUz(u3t,c_uz,d_uz,l_uz);
      LUz(v3t,c_vz,d_vz,l_vz);
      */
      DualLUz(u3t,c_uz,d_uz,l_uz,
	      v3t,c_vz,d_vz,l_vz);

      FOR(k,Nz){
	un[k][j][i] = u3t[k];
	vn[k][j][i] = v3t[k];
      }
    }
  }
}


