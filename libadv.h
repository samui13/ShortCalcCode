#include <limits.h>
#include <math.h>
double Lx,dt,dx,dx2;
double ddxx,ddxx2;
char FOLDER[256];
//#define FOLDER "./data/0"
#ifndef Nx
#define Nx N
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
#ifndef __BOUN__
// Boundary 0 Condition
#define I1(a) ((a)+1 == N+1 ? N-1 : (a)+1)
#define I2(a) ((a)-1 == -1 ? 1 : (a)-1)
#define I3(a) ((a)+2 == N+2 ? N-2 : I1((a)+1))
#define I4(a) ((a)-2 == -2 ? 2 : I2((a)-1))
#else
// 周期
#define I1(a) ((a)+1 == N+1 ? 0 : (a)+1)
#define I2(a) ((a)-1 == -1  ? N : (a)-1)
#define I3(a) ((a)+2 == N+2 ? 1 : I1((a)+1))
#define I4(a) ((a)-2 == N-2 ? N-1: I2((a)-1))
#endif

#ifndef TTD
#define TTD INT_MAX
#endif

//Positive advection UPW
#define PUPW(u,LL,L,C,R,RR) (-u[RR]+6*u[R]-3*u[C]-2*u[L])*ddxx/(6.0)
//Negative advection UPW
#define NUPW(u,LL,L,C,R,RR) (2*u[R]+3*u[C]-6*u[L]+u[LL])*ddxx/(6.0)
#define UPW(c,a,LL,L,C,R,RR) (c < 0 ? PUPW(a,LL,L,C,R,RR) : NUPW(a,LL,L,C,R,RR))

//Positive advection Quick
#define PQUICK(u,LL,L,C,R,RR) (-3*u[C]+7*u[R]-3*u[L]-u[RR])*ddxx*0.125
//Negative advection Quick
#define NQUICK(u,LL,L,C,R,RR) (3*u[C]+3*u[R]-7*u[L]+u[LL])*ddxx*0.125
#define QUICK(c,a,LL,L,C,R,RR) ( c < 0 ? PQUICK(a,LL,L,C,R,RR) : NQUICK(a,LL,L,C,R,RR))


#define DIFFU(a,L,C,R) (a[L]-2*a[C]+a[R])*ddxx2
#define DIFFE(a,L,C,R) (a[L]-a[R])*ddxx*0.5
#define EULER(a,b,TERM,p) b[p] = a[p]+dt*TERM

#define ABS(s) ((s) >0 ? (s) : -1*(s))
#define TWICE(s) (s*s)
#define RSEED (unsigned int)1354682563
#define RDOUBLE (2.0*rand()/(RAND_MAX+1.0)-1.0)
#define EPS 1000000
double L1(double *u);
double epsL1(double *u1,double *u2);
double epsL2(double *u1,double *u2);
void cp1(double *u,double *uo);
double Liapun3L1(double *u,double *v,double *w,
		 double *u1,double *v1,double *w1,double z0);
double Liapun1L1(double *u,double *v,double *w,
		 double *u1,double *v1,double *w1,double z0);
double Liapun1L2(double *u,double *v,double *w,
		 double *u1,double *v1,double *w1,double z0);
void cp3(double *u,double *v,double *w,double *uo,double *vo,double *wo);
void OutPut1(int t,double *u);
void OutPut2(int t,double *u,double *v);
void OutPut3(int t,double *u,double *v,double *w);
void END();
double* mallmatrix();
void freematrix(double *m);
void InitialLUx(double D,double *cM,double *dM,double *lM);
void LUx(double *uh,double *cM,double *dM,double *lM);


double L1(double *u){
  int i;
  double norm = 0;
  REP(i,1,N){
    norm+= dx*0.5*(u[i]+u[I2(i)]);
  }
  return norm;
}

double epsL1(double *u1,double *u2){
  double eps = 0;
  int i;
  FOR(i,N){
    eps+=ABS(u1[i]-u2[i]);
  }
  return eps;
}
double epsL2(double *u1,double *u2){
  double eps = 0;
  int i;
  FOR(i,N){
    eps+=(u1[i]-u2[i])*(u1[i]-u2[i]);
  }
  return sqrt(eps);
}

void cp1(double *u,double *uo){
  int i;
  FOR(i,N){
    uo[i] = u[i];
  }
}
double Liapun3L1(double *u,double *v,double *w,
	       double *u1,double *v1,double *w1,double z0){
  double zn;
  int i;
  zn = epsL1(u,u1)+epsL1(v,v1)+epsL1(w,w1);
  FOR(i,Nx){
    u1[i] = u[i]+(z0/zn)*(u1[i]-u[i]);
    v1[i] = v[i]+(z0/zn)*(v1[i]-v[i]);
    w1[i] = w[i]+(z0/zn)*(w1[i]-w[i]);
  }
  return zn;
}
double Liapun1L1(double *u,double *v,double *w,
	       double *u1,double *v1,double *w1,double z0){
  double zn;
  int i;
  zn = epsL1(w,w1);
  FOR(i,Nx){
    u1[i] = u[i]+(z0/zn)*(u1[i]-u[i]);
    v1[i] = v[i]+(z0/zn)*(v1[i]-v[i]);
    w1[i] = w[i]+(z0/zn)*(w1[i]-w[i]);
  }
  return zn;
}

double Liapun1L2(double *u,double *v,double *w,
	       double *u1,double *v1,double *w1,double z0){
  double zn;
  int i;
  zn = fabs(epsL2(u,u1));
  FOR(i,Nx){
    u1[i] = u[i]+(z0/zn)*(u1[i]-u[i]);
    v1[i] = v[i]+(z0/zn)*(v1[i]-v[i]);
    w1[i] = w[i]+(z0/zn)*(w1[i]-w[i]);
  }
  return zn;
}

void cp3(double *u,double *v,double *w,double *uo,double *vo,double *wo){
  int i;
  FOR(i,N){
    uo[i] = u[i];
    vo[i] = v[i];
    wo[i] = w[i];
  }
}

void OutPut1(int t,double *u){
  char path[256];
  FILE *fp;
  int i,j;
  sprintf(path,"%s/data%d.txt",FOLDER,t);
  if((fp = fopen(path,"wt")) == NULL){
    printf("Directory not found\n");
    exit(1);
  }
  FOR(i,N){
    FDPRT(fp,dx*i);
    FDPRT(fp,u[i]);
    FSPRT(fp,"\n");
  }
  fclose(fp);
  printf("t = %f\t Output %s\n",dt*t*TD,path);
}

void OutPut2(int t,double *u,double *v){
  char path[256];
  FILE *fp;
  int i,j;
  sprintf(path,"%s/data%d.txt",FOLDER,t);
  if((fp = fopen(path,"wt")) == NULL){
    printf("Directory not found\n");
    exit(1);
  }
  FOR(i,N){
    FDPRT(fp,dx*i);
    FDPRT(fp,u[i]);
    FDPRT(fp,v[i]);
    FSPRT(fp,"\n");
  }
  fclose(fp);
  printf("t = %f\t Output %s\n",dt*t*TD,path);
}


void OutPut3(int t,double *u,double *v,double *w){
  char path[256];
  FILE *fp;
  int i,j;
  sprintf(path,"%s/data%d.txt",FOLDER,t);
  if((fp = fopen(path,"wt")) == NULL){
    printf("Directory not found\n");
    exit(1);
  }
  FOR(i,N){
    FDPRT(fp,dx*i);
    FDPRT(fp,u[i]);
    FDPRT(fp,v[i]);
    FDPRT(fp,w[i]);
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
  fclose(fp);
}

double* mallmatrix(){
  int i;
  double *matrix;
  matrix = (double *)malloc(sizeof(double *)*(N+1));
  return matrix;
}
void freematrix(double *m){
  int i;
  free(m);
}

void InitialLUx(double D,double *cM,double *dM,double *lM){
  int i;
  double R = D*dt/(2*dx*dx);
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

void InitialValue_GaussDistribution(double *ma,double rho,double mu){
  int i;
  FOR(i,N){
    ma[i] = 1/sqrt(2*M_PI*rho*rho)*exp(-(dx*i-mu)*(dx*i-mu)/(2*rho*rho));
  }
}
void InitialValue_Palse(double *ma,double rho,double mu){
  int i;
  FOR(i , N){
    if(i*dx >= mu && i*dx <= Lx-mu){
      ma[i] = 1.0;
    }
  }
}
void InitialValue_Discreate(double *ma){
  int i;
  
}
