#include <limits.h>
#include <math.h>
#include <time.h>
double Lx,Ly,Lz,dt,dx,dy,dz;
double dx2,dy2,dz2;
double ddxx,ddyy,ddzz;
double ddxx2,ddyy2,ddzz2;
char FOLDER[256];
//#define FOLDER "./data/0"
//[Nz+1][Ny+1][Nx+1]
#ifndef Nx
#define Nx N
#define Ny N
#define Nz N
#else
#define N Nx
#endif

#define FOR(i,e) for(i = 0; i <= e; i++)
#define REP(i,s,e) for(i = s; i <= e; i++)

#define SPRT(a) printf("Debug Message ==  %s\n",(a))
#define NPRT(a,b) printf("Debug %s Message == %d\n",(b),(a))
#define DPRT(a,b) printf("Debug %s Message == %.15lf\n",(b),(a))

#define FSPRT(p,s) fprintf(p,"%s",s)
#define FSFPRT(p,s,f) fprintf(p,"%s%lf\n",s,f)
#define FDPRT(p,s) fprintf(p,"%.15lf ",s)

#ifndef __BOUNX__
// Boundary 0 Condition
#define I1(a) ((a)+1 == Nx+1 ? Nx-1 : (a)+1)
#define I2(a) ((a)-1 == -1 ? 1 : (a)-1)
#define I3(a) ((a)+2 == Nx+2 ? Nx-2 : I1((a)+1))
#define I4(a) ((a)-2 == -2 ? 2 : I2((a)-1))
#else
// 周期
#define I1(a) ((a)+1 == Nx+1 ? 0 : a+1)
#define I2(a) ((a)-1 == -1  ? Nx : a-1)
#define I3(a) ((a)+2 == Nx+2 ? 1 : I1(a+1))
#define I4(a) ((a)-2 == Nx-2 ? Nx-1: I2(a-1))
#endif

#ifndef __BOUNY__
// Boundary 0 Condition
#define J1(a) ((a)+1 == Ny+1 ? Ny-1 : (a)+1)
#define J2(a) ((a)-1 == -1 ? 1 : (a)-1)
#define J3(a) ((a)+2 == Ny+2 ? Ny-2 : J1((a)+1))
#define J4(a) ((a)-2 == -2 ? 2 : J2((a)-1))
#else
// 周期
#define J1(a) ((a)+1 == Ny+1 ? 0 : a+1)
#define J2(a) ((a)-1 == -1  ? Ny : a-1)
#define J3(a) ((a)+2 == Ny+2 ? 1 : J1(a+1))
#define J4(a) ((a)-2 == Ny-2 ? Ny-1: J2(a-1))
#endif

#ifndef __BOUNZ__
// Boundary 0 Condition
#define K1(a) ((a)+1 == Nz+1 ? Nz-1 : (a)+1)
#define K2(a) ((a)-1 == -1 ? 1 : (a)-1)
#define K3(a) ((a)+2 == Nz+2 ? Nz-2 : K1((a)+1))
#define K4(a) ((a)-2 == -2 ? 2 : K2((a)-1))
#else
// 周期
#define K1(a) ((a)+1 == Nz+1 ? 0 : a+1)
#define K2(a) ((a)-1 == -1  ? Nz : a-1)
#define K3(a) ((a)+2 == Nz+2 ? 1 : K1(a+1))
#define K4(a) ((a)-2 == Nz-2 ? Nz-1: K2(a-1))
#endif

#ifndef TTD
#define TTD INT_MAX
#endif


//Positive advection UPW
#define PUPWX(u,LLx,Lx,Cx,Rx,RRx,Cy,Cz) (-u[Cz][Cy][RRx]+6*u[Cz][Cy][Rx]-3*u[Cz][Cy][Cx]-2*u[Cz][Cy][Lx])*ddxx/(6.0)
#define PUPWY(u,Cx,LLy,Ly,Cy,Ry,RRy,Cz) (-u[Cz][RRy][Cx]+6*u[Cz][Ry][Cx]-3*u[Cz][Cy][Cx]-2*u[Cz][Ly][Cx])*ddyy/(6.0)
#define PUPWZ(u,Cx,Cy,LLz,Lz,Cz,Rz,RRz) (-u[RRz][Cy][Cx]+6*u[Rz][Cy][Cx]-3*u[Cz][Cy][Cx]-2*u[Lz][Cy][Cx])*ddzz/(6.0)
//Negative advection UPW
#define NUPWX(u,LLx,Lx,Cx,Rx,RRx,Cy,Cz) (2*u[Cz][Cy][Rx]+3*u[Cz][Cy][Cx]-6*u[Cz][Cy][Lx]+u[Cz][Cy][LLx])*ddxx/(6.0)
#define NUPWY(u,Cx,LLy,Ly,Cy,Ry,RRy,Cz) (2*u[Cz][Ry][Cx]+3*u[Cz][Cy][Cx]-6*u[Cz][Ly][Cx]+u[Cz][LLy][Cx])*ddyy/(6.0)
#define NUPWZ(u,Cx,Cy,LLz,Lz,Cz,Rz,RRz) (2*u[Rz][Cy][Cx]+3*u[Cz][Cy][Cx]-6*u[Lz][Cy][Cz]+u[LLz][Cy][Cx])*ddzz/(6.0)
#define UPWX(VELOCI,a,LLx,Lx,Cx,Rx,RRx,Cy,Cz) (VELOCI < 0 ? PUPWX(a,LLx,Lx,Cx,Rx,RRx,Cy,Cz) : NUPWX(a,LLx,Lx,Cx,Rx,RRx,Cy,Cz))
#define UPWY(VELOCI,a,Cx,LLy,Ly,Cy,Ry,RRy,Cz) (VELOCI < 0 ? PUPWY(a,Cx,LLy,Ly,Cy,Ry,RRy,Cz) : NUPWY(a,Cx,LLy,Ly,Cy,Ry,RRy,Cz))
#define UPWZ(VELOCI,a,Cx,Cy,LLz,Lz,Cz,Rz,RRz) (VELOCI < 0 ? PUPWZ(a,Cx,Cy,LLz,Lz,Cz,Rz,RRz) : NUPWZ(a,Cx,Cy,LLz,Lz,Cz,Rz,RRz))

//Positive advection Quick
//i,j(x,y)
#define PQUICKX(u,LLx,Lx,Cx,Rx,RRx,Cy,Cz) (-3*u[Cz][Cy][Cx]+7*u[Cz][Cy][Rx]-3*u[Cz][Cy][Lx]-u[Cz][Cy][RRx])*ddxx*0.125
#define PQUICKY(u,Cx,LLy,Ly,Cy,Ry,RRy,Cz) (-3*u[Cz][Cy][Cx]+7*u[Cz][Ry][Cx]-3*u[Cz][Ly][Cx]-u[Cz][RRy][Cx])*ddyy*0.125
#define PQUICKZ(u,Cx,Cy,LLz,Lz,Cz,Rz,RRz) (-3*u[Cz][Cy][Cx]+7*u[Rz][Cy][Cx]-3*u[Lz][Cy][Cx]-u[RRz][Cy][Cx])*ddzz*0.125
//Negative advection Quick
#define NQUICKX(u,LLx,Lx,Cx,Rx,RRx,Cy,Cz) (3*u[Cz][Cy][Cx]+3*u[Cz][Cy][Rx]-7*u[Cz][Cy][Lx]+u[Cz][Cy][LLx])*ddxx*0.125
#define NQUICKY(u,Cx,LLy,Ly,Cy,Ry,RRy,Cz) (3*u[Cz][Cy][Cx]+3*u[Cz][Ry][Cx]-7*u[Cz][Ly][Cx]+u[Cz][LLy][Cx])*ddyy*0.125
#define NQUICKZ(u,Cx,Cy,LLz,Lz,Cz,Rz,RRz) (3*u[Cz][Cy][Cz]+3*u[Rz][Cy][Cx]-7*u[Lz][Cy][Cx]+u[LLz][Cy][Cx])*ddzz*0.125

#define QUICKX(VELOCI,a,LLx,Lx,Cx,Rx,RRx,Cy,Cz) ( VELOCI < 0 ? PQUICKX(a,LLx,Lx,Cx,Rx,RRx,Cy,Cz) : NQUICKX(a,LLx,Lx,Cx,Rx,RRx,Cy,Cz))
#define QUICKY(VELOCI,a,Cx,LLy,Ly,Cy,Ry,RRy,Cz) ( VELOCI < 0 ? PQUICKY(a,Cx,LLy,Ly,Cy,Ry,RRy,Cz) : NQUICKY(a,Cx,LLy,Ly,Cy,Ry,RRy,Cz))
#define QUICKZ(VELOCI,a,Cx,Cy,LLz,Lz,Cz,Rz,RRz) ( VELOCI < 0 ? PQUICKZ(a,Cx,Cy,LLz,Lz,Cz,Rz,RRz) : NQUICKZ(a,Cx,Cy,LLz,Lz,Cz,Rz,RRz))

#define DIFFUX(a,Lx,Cx,Rx,Cy,Cz) (a[Cz][Cy][Lx]-2*a[Cz][Cy][Cx]+a[Cz][Cy][Rx])*(ddxx2)
#define DIFFUY(a,Cx,Ly,Cy,Ry,Cz) (a[Cz][Ly][Cx]-2*a[Cz][Cy][Cx]+a[Cz][Ry][Cx])*(ddyy2)
#define DIFFUZ(a,Cx,Cy,Lz,Cz,Rz) (a[Lz][Cy][Cx]-2*a[Cz][Cy][Cx]+a[Rz][Cy][Cx])*(ddzz2)
#define DIFFEX(a,Lx,Cx,Rx,Cy,Cz) (a[Cz][Cy][Rx]-a[Cz][Cy][Lx])*0.5*ddxx
#define DIFFEY(a,Cx,Ly,Cy,Ry,Cz) (a[Cz][Ry][Cx]-a[Cz][Ly][Cx])*0.5*ddyy
#define DIFFEZ(a,Cx,Cy,Lz,Cz,Rz) (a[Rz][Cy][Cx]-a[Lz][Cy][Cx])*0.5*ddzz
//EULERは使わないほうがいい。
#define EULER(a,b,c,p,q,r) b[r][q][p] = a[r][q][p]+dt*c

#define AREA(p,q,r) (p == 0 || p == Nx ? dx/2:dx)*(q == 0 || q  == Ny ? dy/2:dy)*(r == 0 || r == Nz ? dz/2 : dz)

#define RSEED (unsigned int)1354682563
#define RDOUBLE (2.0*rand()/(RAND_MAX+1.0)-1.0)
#define EPS 1000000

double L1(double ***u){
  int i,j,k;
  double norm = 0,S;
  FOR(k,Nz){
    FOR(j,Ny){
      FOR(i,Nx){
	S = AREA(i,j,k);
	norm+=S*u[k][j][i];
      }
    }
  }
  return norm;
}

double epsL1(double ***u1,double ***u2){
  int i,j,k;
  double eps = 0,S;
  FOR(k,Nz){
    FOR(j,Ny){
      FOR(i,Nx){
	S = AREA(i,j,k);
	eps+=S*(u1[k][j][i]-u2[k][j][i] > 0 ? u1[k][j][i]-u2[k][j][i]:u2[k][j][i]-u1[k][j][i]);
      }
    }
  }
  return eps;
}
void cp1(double ***u,double ***uo){
  int i,j,k;
  FOR(k,Nz){
    FOR(j,Ny){
      FOR(i,Nx){
	uo[k][j][i] = u[k][j][i];
      }
    }
  }
}

void cp2(double ***u,double ***v,double ***uo,double ***vo){
  int i,j,k;
  FOR(k,Nz){
    FOR(j,Ny){
      FOR(i,Nx){
	uo[k][j][i] = u[k][j][i];
	vo[k][j][i] = v[k][j][i];
      }
    }
  }
}

void cp3(double ***u,double ***v,double ***w,double ***uo,double ***vo,double ***wo){
  int i,j,k;
  FOR(k,Nz){
    FOR(j,Ny){
      FOR(i,Nx){
	uo[k][j][i] = u[k][j][i];
	vo[k][j][i] = v[k][j][i];
	wo[k][j][i] = w[k][j][i];
      }
    }
  }
}

void OutPut1(int t,double ***u){
  char path[256];
  FILE *fp;
  int i,j,k;
  sprintf(path,"%s/data%d.txt",FOLDER,t);
  if((fp = fopen(path,"wt")) == NULL){
    printf("Directory not found\n");
    exit(1);
  }
  FOR(k,Nz){
    FOR(j,Ny){
      FOR(i,Nx){
	FDPRT(fp,dx*i);
	FDPRT(fp,dy*j);
	FDPRT(fp,dz*k);
	FDPRT(fp,u[k][j][i]);
	FSPRT(fp,"\n");
      }
      FSPRT(fp,"\n");
    }
    FSPRT(fp,"\n");
  }
    fclose(fp);
    printf("t = %f\t Output %s\n",dt*t*TD,path);
}
void OutPut2(int t,double ***u,double ***v){
  char path[256];
  FILE *fp;
  int i,j,k;
  sprintf(path,"%s/data%d.txt",FOLDER,t);
  if((fp = fopen(path,"wt")) == NULL){
    printf("Directory not found\n");
    exit(1);
  }
  FOR(k,Nz){
    FOR(j,Ny){
      FOR(i,Nx){
	FDPRT(fp,dx*i);
	FDPRT(fp,dy*j);
	FDPRT(fp,dz*k);
	FDPRT(fp,u[k][j][i]);
	FDPRT(fp,v[k][j][i]);
	FSPRT(fp,"\n");
      }
      FSPRT(fp,"\n");
    }
    FSPRT(fp,"\n");
  }
    fclose(fp);
    printf("t = %f\t Output %s\n",dt*t*TD,path);
}

void OutPut3(int t,double ***u,double ***v,double ***w){
  char path[256];
  FILE *fp;
  int i,j,k;
  sprintf(path,"%s/data%d.txt",FOLDER,t);
  if((fp = fopen(path,"wt")) == NULL){
    printf("Directory not found\n");
    exit(1);
  }
  FOR(k,Nz){
    FOR(j,Ny){
      FOR(i,Nx){
	FDPRT(fp,dx*i);
	FDPRT(fp,dy*j);
	FDPRT(fp,dz*j);
	FDPRT(fp,u[k][j][i]);
	FDPRT(fp,v[k][j][i]);
	FDPRT(fp,w[k][j][i]);
	FSPRT(fp,"\n");
      }
      FSPRT(fp,"\n");
    }
    FSPRT(fp,"\n");
  }
  fclose(fp);
  printf("t = %f\t Output %s\n",dt*t*TD,path);
}

void END(){
  char path[256];
  FILE *fp;
  sprintf(path,"%s/end.txt",FOLDER);
  if((fp = fopen(path,"wt")) == NULL){
    printf("Directory not found");
    exit(1);
  }
  FSPRT(fp,"Compile DATE is\t");
  FSPRT(fp,__DATE__);
  FSPRT(fp,"\n");
  FSPRT(fp,"Compile TIME is\t");
  FSPRT(fp,__TIME__);
  FSPRT(fp,"\n");
  FSPRT(fp,"This File is\t");
  FSPRT(fp,__FILE__);
  FSPRT(fp,"\n");
  FSPRT(fp,"End Time:");
  fprintf(fp,"%d",(unsigned int)time(NULL));
  fclose(fp);
}
void LUx(double *uh,double *cM,double *dM,double *lM){
  int i;
  double z[Nx+1];
  z[0] = uh[0];
  for (i = 1; i<=Nx; i++){
    z[i] = uh[i]-lM[i]*z[i-1];
  }
  uh[Nx] = z[Nx]/dM[Nx];
  for(i = 1; i <= Nx;i++){
    uh[Nx-i] = (z[Nx-i]-cM[Nx-i]*uh[Nx-i+1])/dM[Nx-i];
  }
}
void DualLUx(double *uh,double *ucM,double *udM,double *ulM,
	     double *vh,double *vcM,double *vdM,double *vlM){
  int i;
  double zu[Nx+1],zv[Nx+1];
  zu[0] = uh[0];
  zv[0] = vh[0];
  for (i = 1; i<=Nx; i++){
    zu[i] = uh[i]-ulM[i]*zu[i-1];
    zv[i] = vh[i]-vlM[i]*zv[i-1];
  }
  uh[Nx] = zu[Nx]/udM[Nx];
  vh[Nx] = zv[Nx]/vdM[Nx];
  for(i = 1; i <= Nx;i++){
    uh[Nx-i] = (zu[Nx-i]-ucM[Nx-i]*uh[Nx-i+1])/udM[Nx-i];
    vh[Nx-i] = (zv[Nx-i]-vcM[Nx-i]*vh[Nx-i+1])/vdM[Nx-i];
  }
}
void LUy(double *uh,double *cM,double *dM,double *lM){
  int i;
  double z[Ny+1];
  z[0] = uh[0];
  for (i = 1; i<=Ny; i++){
    z[i] = uh[i]-lM[i]*z[i-1];
  }
  uh[Ny] = z[Ny]/dM[Ny];
  for(i = 1; i <= Ny;i++){
    uh[Ny-i] = (z[Ny-i]-cM[Ny-i]*uh[Ny-i+1])/dM[Ny-i];
  }
}
void DualLUy(double *uh,double *ucM,double *udM,double *ulM,
	     double *vh,double *vcM,double *vdM,double *vlM){
  int i;
  double zu[Ny+1],zv[Ny+1];
  zu[0] = uh[0];
  zv[0] = vh[0];
  for (i = 1; i<=Ny; i++){
    zu[i] = uh[i]-ulM[i]*zu[i-1];
    zv[i] = vh[i]-vlM[i]*zv[i-1];
  }
  uh[Ny] = zu[Ny]/udM[Ny];
  vh[Ny] = zv[Ny]/vdM[Ny];
  for(i = 1; i <= Ny;i++){
    uh[Ny-i] = (zu[Ny-i]-ucM[Ny-i]*uh[Ny-i+1])/udM[Ny-i];
    vh[Ny-i] = (zv[Ny-i]-vcM[Ny-i]*vh[Ny-i+1])/vdM[Ny-i];
  }
}
void LUz(double *uh,double *cM,double *dM,double *lM){
  int i;
  double z[Nz+1];
  z[0] = uh[0];
  for (i = 1; i<=Nz; i++){
    z[i] = uh[i]-lM[i]*z[i-1];
  }
  uh[Nz] = z[Nz]/dM[Nz];
  for(i = 1; i <= Nz;i++){
    uh[Nz-i] = (z[Nz-i]-cM[Nz-i]*uh[Nz-i+1])/dM[Nz-i];
  }
}
void DualLUz(double *uh,double *ucM,double *udM,double *ulM,
	     double *vh,double *vcM,double *vdM,double *vlM){
  int i;
  double zu[Nz+1],zv[Nz+1];
  zu[0] = uh[0];
  zv[0] = vh[0];
  for (i = 1; i<=Nz; i++){
    zu[i] = uh[i]-ulM[i]*zu[i-1];
    zv[i] = vh[i]-vlM[i]*zv[i-1];
  }
  uh[Nz] = zu[Nz]/udM[Nz];
  vh[Nz] = zv[Nz]/vdM[Nz];
  for(i = 1; i <= Nz;i++){
    uh[Nz-i] = (zu[Nz-i]-ucM[Nz-i]*uh[Nz-i+1])/udM[Nz-i];
    vh[Nz-i] = (zv[Nz-i]-vcM[Nz-i]*vh[Nz-i+1])/vdM[Nz-i];
  }
}

void Initial_LUx(double D,double *cM,double *dM,double *lM){
  int i;
  double  R = D*dt/(2*dx*dx);
  
  dM[0] = 1+2*R;
  cM[0] = -2*R;

  lM[0] = -R/dM[0];
  for(i = 1;i < Nx;i++){
    cM[i] = -R;
    lM[i] = -R/dM[i-1];
    dM[i] = (1+2*R)-lM[i]*cM[i-1];
  }
  lM[Nx] = -2*R/dM[Nx-1];
  dM[Nx] = 1+2*R-lM[Nx]*cM[Nx-1];
}
void Initial_LUy(double D,double *cM,double *dM,double *lM){
  int i;
  double  R = D*dt/(2*dy*dy);
  dM[0] = 1+2*R;
  cM[0] = -2*R;
  lM[0] = -R/dM[0];
  for(i = 1;i < Ny;i++){
    cM[i] = -R;
    lM[i] = -R/dM[i-1];
    dM[i] = (1+2*R)-lM[i]*cM[i-1];
  }
  lM[Ny] = -2*R/dM[Ny-1];
  dM[Ny] = 1+2*R-lM[Ny]*cM[Ny-1];
}
void Initial_LUz(double D,double *cM,double *dM,double *lM){
  int i;
  double  R = D*dt/(2*dz*dz);
  dM[0] = 1+2*R;
  cM[0] = -2*R;
  lM[0] = -R/dM[0];
  for(i = 1;i < Nz;i++){
    cM[i] = -R;
    lM[i] = -R/dM[i-1];
    dM[i] = (1+2*R)-lM[i]*cM[i-1];
  }
  lM[Nz] = -2*R/dM[Nz-1];
  dM[Nz] = 1+2*R-lM[Nz]*cM[Nz-1];
}

double*** mallmatrix(){
  int i,j,k;
  double ***matrix;
  matrix = (double ***)malloc(sizeof(double*)*(Nz+1));
  FOR(k,Nz){
    matrix[k] = (double **)malloc(sizeof(double*)*(Ny+1));
    FOR(j,Ny){
      matrix[k][j] = (double*)malloc(sizeof(double)*(Nx+1));
    }
  }
  return matrix;
}

void freematrix(double ***matrix){
  int i,j,k;
  FOR(k,Nz){
      FOR(j,Ny){
	free(matrix[k][j]);
      }
      free(matrix[k]);
  }
  free(matrix);
}
