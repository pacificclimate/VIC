#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include "vicNl.h"
#include "newt_raph_func_fast.h"

#define MAXSIZE  20
#define MAXTRIAL 150
#define TOLX     1e-4
#define TOLF     1e-1
#define R_MAX    2.0
#define R_MIN    -5.0
#define RELAX1   0.9
#define RELAX2   0.7
#define RELAX3   0.2

NewtonRaphsonMethod::NewtonRaphsonMethod(double T_2[], double res[], int n,
    double deltat, int FS_ACTIVE, int NOFLUX, int EXP_TRANS, double *T0,
    double *moist, double *ice, double *kappa, double *Cs,
    const double *max_moist, const double *bubble, const double *expt,
    const double *alpha, const double *beta, const double *gamma,
    const double *Zsum, double Dp, const double *bulk_dens_min,
    const double *soil_dens_min, const double *quartz,
    const double *bulk_density, const double *soil_density,
    const double *organic, const double *depth, int Nlayers) :

    deltat(deltat), FS_ACTIVE(FS_ACTIVE), NOFLUX(NOFLUX), EXP_TRANS(EXP_TRANS), T0(
        T0), moist(moist), ice(ice), kappa(kappa), Cs(Cs), max_moist(max_moist), bubble(
        bubble), expt(expt), alpha(alpha), beta(beta), gamma(gamma), Zsum(Zsum), Dp(
        Dp), bulk_dens_min(bulk_dens_min), soil_dens_min(soil_dens_min), quartz(
        quartz), bulk_density(bulk_density), soil_density(soil_density), organic(
        organic), depth(depth), Nlayers(Nlayers) {

  if (EXP_TRANS) {
    if (!NOFLUX)
      Bexp = log(Dp + 1.) / (double) (n + 1);
    else
      Bexp = log(Dp + 1.) / (double) (n);
  }

  Ts = T0[0];
  if (!NOFLUX)
    Tb = T0[n + 1];
  else
    Tb = T0[n];
  for (int i = 0; i < n; i++)
    T_2[i] = T0[i + 1];

}

int NewtonRaphsonMethod::compute(double x[], int n)
{

/******************************************************************
  newt_raph    2006    Ming Pan      mpan@princeton.EDU

 Newton-Raphson method to solve non-linear system
 adapted from "Numerical Recipes"

 A relaxation factor "RELAX#" is added to help the Newton-Raphson trials
 to converge during the initial formation of ice, where the shape 
 of the residual function becomes very difficult 

******************************************************************/

  int k, i, index[MAXSIZE], Error;
  double errx, errf, d, fvec[MAXSIZE], fjac[MAXSIZE*MAXSIZE], p[MAXSIZE];
  double a[MAXSIZE], b[MAXSIZE], c[MAXSIZE];

  Error = 0;

  for (k=0; k<MAXTRIAL; k++) {

    // calculate function value for all nodes, i.e. focus = -1
    fda_heat_eqn(x, fvec, n, -1);

    // stop if TOLF is satisfied
    errf=0.0;
    for (i=0; i<n; i++) errf+=fabs(fvec[i]);
    if (errf<=TOLF) {
      //fprintf(stderr, "Number of Newton-Raphson trials (F criterium with F error = %g): %d\n", errf, k);
      return (Error);
    }
    
    // calculate the Jacobian
    fdjac3(x, fvec, a, b, c, n);

    for (i=0; i<n; i++) p[i]=-fvec[i];

    // only the tri-diagnal part of the Jacobian matrix is significant
    // and most of the off-belt entries are encessially zeros
    tridiag(a, b, c, p, n);

    errx=0.0;
    for (i=0; i<n; i++) {
      errx+=fabs(p[i]);

      if   (k>10&&k<=20&&x[i]<R_MAX&&x[i]>R_MIN) x[i]+=p[i]*RELAX1;
      else if (k>20&&k<=60&&x[i]<R_MAX&&x[i]>R_MIN) x[i]+=p[i]*RELAX2;
      else if (k>60&&x[i]<R_MAX&&x[i]>R_MIN) x[i]+=p[i]*RELAX3;
      else x[i]+=p[i];

    }

    // stop if TOLX is satisfied
    if (errx<=TOLX) {
      //fprintf(stderr, "Number of Newton-Raphson trials (x criterium with F error = %g): %d\n", errf, k);
      return (Error);
    }
  }
  Error = 1;
#if VERBOSE
  //fprintf(stderr, "WARNING: Maximum number of trials %d reached in Newton-Raphson search for solution (with F error = %g).\n", MAXTRIAL, errf);
  //for (i=0; i<n; i++) 
  //fprintf(stderr,"%d %.2f %.2f %.2f\n",i+1,x[i],fvec[i],p[i]);  
  //fprintf(stderr,"Explicit method with root_brent will be used to solve soil thermal fluxes\n");
#endif
  //vicerror("");

  return (Error);
}

#undef MAXTRIAL
#undef TOLX
#undef TOLF
#undef R_MAX
#undef R_MIN
#undef RELAX1
#undef RELAX2
#undef RELAX3



#define EPS2     1e-4

void NewtonRaphsonMethod::fdjac3(double x[], double fvec[], double a[], double b[], double c[], int n)
{

/******************************************************************
  fdjac3    2006    Ming Pan      mpan@princeton.EDU

 forward difference approx to Jacobian,
 adapted from "Numerical Recipes"

******************************************************************/

  int i, j;
  double h, temp, f[MAXSIZE];

  for (j=0; j<n; j++) {
    temp=x[j];
    h=EPS2*fabs(temp);
    if (h==0) h=EPS2;
    x[j]=temp+h;
    h=x[j]-temp;

    // only update column j-1, j and j+1, caused by change in x[j]
    fda_heat_eqn(x, f, n, j);

    x[j]=temp;

    b[j]=(f[j]-fvec[j])/h;
    if (j!=0) c[j-1]=(f[j-1]-fvec[j-1])/h;
    if (j!=n-1) a[j+1]=(f[j+1]-fvec[j+1])/h;
  }

}

#undef EPS2
#undef MAXSIZE


void tridiag(double a[], double b[], double c[], double r[], int n)
{

/******************************************************************
  tridiag    2006    Ming Pan      mpan@princeton.EDU

 function to solve tridiagonal linear system
 adapted from "Numerical Recipes"

******************************************************************/

  int j;
  double factor;

  /* forward substitution */
  factor=b[0];
  b[0]=1.0;
  c[0]=c[0]/factor;
  r[0]=r[0]/factor;

  for (j=1; j<n; j++) {

    factor=a[j];
    a[j]=a[j]-b[j-1]*factor;
    b[j]=b[j]-c[j-1]*factor;
    r[j]=r[j]-r[j-1]*factor;

    factor=b[j];
    b[j]=1.0;
    c[j]=c[j]/factor;
    r[j]=r[j]/factor;

  }

  /* backward substitution */
  for (j=n-2; j>=0; j--) {

    factor=c[j];
    c[j]=c[j]-b[j+1]*factor;
    r[j]=r[j]-r[j+1]*factor;

    factor=b[j];
    //b[j]=1.0;
    r[j]=r[j]/factor;
  }

}

