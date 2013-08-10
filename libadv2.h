#include <limits.h>
double Lx,Ly,dt,dx,dy;
double dx2,dy2;
double ddxx,ddxx2;
double ddyy,ddyy2;
char FOLDER[256];

#ifndef Nx
#define Nx N
#define Ny N
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


#ifndef __BOUNI__
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

#ifndef __BOUNJ__
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

#ifndef TTD
#define TTD INT_MAX
#endif


//Positive advection UPW
#define PUPWX(u,LLx,Lx,Cx,Rx,RRx,Cy) (-u[Cy][RRx]+6*u[Cy][Rx]-3*u[Cy][Cx]-2*u[Cy][Lx])*ddxx/(6.0)
#define PUPWY(u,Cx,LLy,Ly,Cy,Ry,RRy) (-u[RRy][Cx]+6*u[Ry][Cx]-3*u[Cy][Cx]-2*u[Ly][Cx])*ddyy/(6.0)
//Negative advection UPW
#define NUPWX(u,LLx,Lx,Cx,Rx,RRx,Cy) (2*u[Cy][Rx]+3*u[Cy][Cx]-6*u[Cy][Lx]+u[Cy][LLx])*ddxx/(6.0)
#define NUPWY(u,Cx,LLy,Ly,Cy,Ry,RRy) (2*u[Ry][Cx]+3*u[Cy][Cx]-6*u[Ly][Cx]+u[LLy][Cx])*ddyy/(6.0)
#define UPWX(VEL,a,LLx,Lx,Cx,Rx,RRx,Cy) (VEL < 0 ? PUPWX(a,LLx,Lx,Cx,Rx,RRx,Cy) : NUPWX(a,LLx,Lx,Cx,Rx,RRx,Cy))
#define UPWY(VEL,a,Cx,LLy,Ly,Cy,Ry,RRy) (VEL < 0 ? PUPWY(a,Cx,LLy,Ly,Cy,Ry,RRy) : NUPWY(a,Cx,LLy,Ly,Cy,Ry,RRy))

//Positive advection Quick
//i,j(x,y)
#define PQUICKX(u,LLx,Lx,Cx,Rx,RRx,Cy) (-3*u[Cy][Cx]+7*u[Cy][Rx]-3*u[Cy][Lx]-u[Cy][RRx])*0.125*ddxx
#define PQUICKY(u,Cx,LLy,Ly,Cy,Ry,RRy) (-3*u[Cy][Cx]+7*u[Ry][Cx]-3*u[Ly][Cx]-u[RRy][Cx])*0.125*ddyy
//Negative advection Quick
#define NQUICKX(u,LLx,Lx,Cx,Rx,RRx,Cy) (3*u[Cy][Cx]+3*u[Cy][Rx]-7*u[Cy][Lx]+u[Cy][LLx])*0.125*ddxx
#define NQUICKY(u,Cx,LLy,Ly,Cy,Ry,RRy) (3*u[Cy][Cx]+3*u[Ry][Cx]-7*u[Ly][Cx]+u[LLy][Cx])*0.125*ddyy
#define QUICKX(VEL,a,LLx,Lx,Cx,Rx,RRx,Cy) ( VEL < 0 ? PQUICKX(a,LLx,Lx,Cx,Rx,RRx,Cy) : NQUICKX(a,LLx,Lx,Cx,Rx,RRx,Cy))
#define QUICKY(VEL,a,Cx,LLy,Ly,Cy,Ry,RRy) ( VEL < 0 ? PQUICKY(a,Cx,LLy,Ly,Cy,Ry,RRy) : NQUICKY(a,Cx,LLy,Ly,Cy,Ry,RRy))


#define DIFFUX(a,Lx,Cx,Rx,Cy) (a[Cy][Rx]-2*a[Cy][Cx]+a[Cy][Lx])*ddxx2
#define DIFFUY(a,Cx,Ly,Cy,Ry) (a[Ry][Cx]-2*a[Cy][Cx]+a[Ly][Cx])*ddyy2
#define DIFFEX(a,Lx,Cx,Rx,Cy) (a[Cy][Rx]-a[Cy][Lx])*ddxx*0.5
#define DIFFEY(a,Cx,Ly,Cy,Ry) (a[Ry][Cx]-a[Ly][Cx])*ddxx*0.5
//EULERは使わないほうがいい。
#define EULER(a,b,c,p,q) b[q][p] = a[q][p]+dt*c

#define AREA(i,j) (i == 0 || i == Nx ? dx/2:dx)*(j == 0 || j  == Ny ? dy/2:dy)

#define RSEED (unsigned int)1354682563
#define RDOUBLE (2.0*rand()/(RAND_MAX+1.0)-1.0)
#define EPS 1000000

double L1(double **u);
double epsL1(double **u1,double **u2);
void cp1(double **u,double **uo);
void cp2(double **u,double **v,double **uo,double **vo);
void cp3(double **u,double **v,double **w,double **uo,double **vo,double **wo);
void OutPut1(int t,double **u);
void OutPut2(int t,double **u,double **v);
void OutPut3(int t,double **u,double **v,double **w);
void END();
void LUx(double *uh,double *cM,double *dM,double *lM);
void LUy(double *uh,double *cM,double *dM,double *lM);
void Initial_LUx(double D,double *cM,double *dM,double *lM);
void Initial_LUy(double D,double *cM,double *dM,double *lM);
double** mallmatrix();
void freematrix(double **matrix);


double L1(double **u){
  int i,j;
  double norm = 0,S;
  //FOR2(i,j,Nx,Ny){
  FOR(j,Ny){
    FOR(i,Nx){
      S = AREA(i,j);
      norm+=S*u[j][i];
    }
  }
  return norm;
}

double epsL1(double **u1,double **u2){
  int i,j;
  double eps = 0,S;
  FOR(j,Ny){
    FOR(i,Nx){
      S = AREA(i,j);
      eps+=S*(u1[j][i]-u2[j][i] > 0 ? u1[j][i]-u2[j][i]:u2[j][i]-u1[j][i]);
    }
  }
  return eps;
}
void cp1(double **u,double **uo){
  int i,j;
  FOR(j,Ny){
    FOR(i,Nx){
      uo[j][i] = u[j][i];
    }
  }
}

void cp2(double **u,double **v,double **uo,double **vo){
  int i,j;
  FOR(j,Ny){
    FOR(i,Nx){
      uo[j][i] = u[j][i];
      vo[j][i] = v[j][i];
    }
  }
}

void cp3(double **u,double **v,double **w,double **uo,double **vo,double **wo){
  int i,j;
  FOR(j,Ny){
    FOR(i,Nx){
      uo[j][i] = u[j][i];
      vo[j][i] = v[j][i];
      wo[j][i] = w[j][i];
    }
  }
}
void OutPut1(int t,double **u){
  char path[256];
  FILE *fp;
  int i,j;
  sprintf(path,"%s/data%d.txt",FOLDER,t);
  if((fp = fopen(path,"wt")) == NULL){
    printf("Directory not found\n");
    exit(1);
  }
  FOR(j,Ny){
    FOR(i,Nx){
      FDPRT(fp,dx*i);
      FDPRT(fp,dy*j);
      FDPRT(fp,u[j][i]);
      FSPRT(fp,"\n");
    }
    FSPRT(fp,"\n");
  }
  fclose(fp);
  printf("t = %f\t Output %s\n",dt*t*TD,path);
}


void OutPut2(int t,double **u,double **v){
  char path[256];
  FILE *fp;
  int i,j;
  sprintf(path,"%s/data%d.txt",FOLDER,t);
  if((fp = fopen(path,"wt")) == NULL){
    printf("Directory not found\n");
    exit(1);
  }
  FOR(j,Ny){
    FOR(i,Nx){
      FDPRT(fp,dx*i);
      FDPRT(fp,dy*j);
      FDPRT(fp,u[j][i]);
      FDPRT(fp,v[j][i]);
      FSPRT(fp,"\n");
    }
    FSPRT(fp,"\n");
  }
  fclose(fp);
  printf("t = %f\t Output %s\n",dt*t*TD,path);
}


void OutPut3(int t,double **u,double **v,double **w){
  char path[256];
  FILE *fp;
  int i,j;
  sprintf(path,"%s/data%d.txt",FOLDER,t);
  if((fp = fopen(path,"wt")) == NULL){
    printf("Directory not found\n");
    exit(1);
  }
  FOR(j,Ny){
    FOR(i,Nx){
      FDPRT(fp,dx*i);
      FDPRT(fp,dy*j);
      FDPRT(fp,u[j][i]);
      FDPRT(fp,v[j][i]);
      FDPRT(fp,w[j][i]);
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
  FSPRT(fp,"START DATE is\t");
  FSPRT(fp,__DATE__);
  FSPRT(fp,"\n");
  FSPRT(fp,"START TIME is\t");
  FSPRT(fp,__TIME__);
  FSPRT(fp,"\n");
  FSPRT(fp,"This File is\t");
  FSPRT(fp,__FILE__);
  FSPRT(fp,"\n");
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
double** mallmatrix(){
  int j;
  double **matrix;
  matrix =  (double **)malloc(sizeof(double*)*(Ny+1));
  FOR(j,Ny){
    matrix[j] = (double*)malloc(sizeof(double)*(Nx+1));
  }
  return matrix;
}

void freematrix(double **matrix){
  int i,j;
  FOR(j,Ny){
    free(matrix[j]);
  }
  free(matrix);
}

void InitialValue_GaussDistribution(double **ma,double rhoX,double muX,
				    double rhoY,double muY){
  int i,j;
  FOR(j,Ny){
    FOR(i,N){
      ma[j][i] = 1/sqrt(2*M_PI*rhoX*rhoX)*1/sqrt(2*M_PI*rhoY*rhoY)*exp(-(dx*i-muX)*(dx*i-muX)/(2*rhoX*rhoX))*exp(-(j*dy-muY)*(j*dy-muY)/(2*rhoY*rhoY));
    }
  }
}

