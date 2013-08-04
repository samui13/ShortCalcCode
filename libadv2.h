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
#define PUPWX(u,p,q) (-u[q][I3(p)]+6*u[q][I1(p)]-3*u[q][p]-2*u[q][I2(p)])/(6*dx)
//#define PUPWX(u,LL,L,C,R,RR,Q) (-u[Q][L]+6*u[Q][R]-3*u[Q][C]-2*u[Q][L])*ddxx/(6.0)
#define PUPWY(u,p,q) (-u[J3(q)][p]+6*u[J1(q)][p]-3*u[q][p]-2*u[J2(q)][p])/(6*dy)
//Negative advection UPW
#define NUPWX(u,p,q) (2*u[q][I1(p)]+3*u[q][p]-6*u[q][I2(p)]+u[q][I4(p)])/(6*dx)
#define NUPWY(u,p,q) (2*u[J1(q)][p]+3*u[q][p]-6*u[J2(q)][p]+u[J4(q)][p])/(6*dy)
#define UPWX(c,a,p,q) c*(c < 0 ? PUPWX(a,p,q) : NUPWX(a,p,q))
#define UPWY(c,a,p,q) c*(c < 0 ? PUPWY(a,p,q) : NUPWY(a,p,q))

//Positive advection Quick
//i,j(x,y)
#define PQUICKX(u,p,q) (-3*u[q][p]+7*u[q][I1(p)]-3*u[q][I2(p)]-u[q][I3(p)])/(8*dx)
#define PQUICKY(u,p,q) (-3*u[q][p]+7*u[J1(q)][p]-3*u[J2(q)][p]-u[J3(q)][p])/(8*dy)
//Negative advection Quick
#define NQUICKX(u,p,q) (3*u[q][p]+3*u[q][I1(p)]-7*u[q][I2(p)]+u[q][I4(p)])/(8*dx)
#define NQUICKY(u,p,q) (3*u[q][p]+3*u[J1(q)][p]-7*u[J2(q)][p]+u[J4(q)][p])/(8*dy)
#define QUICKX(c,a,p,q) c*( c < 0 ? PQUICKX(a,p,q) : NQUICKX(a,p,q))
#define QUICKY(c,a,p,q) c*( c < 0 ? PQUICKY(a,p,q) : NQUICKY(a,p,q))


#define DIFFUX(a,p,q) (a[q][I1(p)]-2*a[q][p]+a[q][I2(p)])/(dx*dx)
#define DIFFUY(a,p,q) (a[J1(q)][p]-2*a[q][p]+a[J2(q)][p])/(dy*dy)
#define DIFFEX(a,p,q) (a[q][I1(p)]-a[q][I2(p)])/(2*dx)
#define DIFFEY(a,p,q) (a[J1(q)][p]-a[J2(q)][p])/(2*dy)
//EULERは使わないほうがいい。
#define EULER(a,b,c,p,q) b[q][p] = a[q][p]+dt*c

#define AREA(i,j) (i == 0 || i == Nx ? dx/2:dx)*(j == 0 || j  == Ny ? dy/2:dy)

#define RSEED (unsigned int)1354682563
#define RDOUBLE (2.0*rand()/(RAND_MAX+1.0)-1.0)
#define EPS 1000000

double L1(double u[Ny+1][Nx+1]){
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

double epsL1(double u1[Ny+1][Nx+1],double u2[Ny+1][Nx+1]){
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
void cp1(double u[Ny+1][Nx+1],double uo[Ny+1][Nx+1]){
  int i,j;
  FOR(j,Ny){
    FOR(i,Nx){
      uo[j][i] = u[j][i];
    }
  }
}

void cp2(double u[Ny+1][Nx+1],double v[Ny+1][Nx+1],double uo[Ny+1][Nx+1],double vo[Ny+1][Nx+1]){
  int i,j;
  FOR(j,Ny){
    FOR(i,Nx){
      uo[j][i] = u[j][i];
      vo[j][i] = v[j][i];
    }
  }
}

void cp3(double u[Ny+1][Nx+1],double v[Ny+1][Nx+1],double w[Ny+1][Nx+1],double uo[Ny+1][Nx+1],double vo[Ny+1][Nx+1],double wo[Ny+1][Nx+1]){
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


void OutPut2(int t,double u[Ny+1][Nx+1],double v[Ny+1][Nx+1]){
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


void OutPut3(int t,double u[Ny+1][Nx+1],double v[Ny+1][Nx+1],double w[Ny+1][Nx+1]){
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
