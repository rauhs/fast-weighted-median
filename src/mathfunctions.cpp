// This file is obsolete for th medianf project due to the new
// formulas. Let's keep it for now though...
#include "mathfunctions.hpp"
#define _USE_MATH_DEFINES
#include <cmath>

static double betacf(double a, double b, double x);
static double gammaln(double xx);

#ifndef M_SQRT_2PI
# define M_SQRT_2PI 2.506628274631000
#endif

double normpdf(double x, double mu, double sigma)
{
  double z = ((x - mu)/sigma);
  return exp(-0.5*z*z)/(sqrt(2*M_PI)*sigma);
}

double poisspdf(double x, double lambda)
{
  return exp(-lambda + x * log(lambda) - gammaln(x+1));
}

double betai(double a, double b, double x)
{
  double bt;

  if (x < 0.0 || x > 1.0)
    return -1.0; // Error
  if (x==0.0) return 0.0 ;
  if (x==1.0) return 1.0 ;
  /* Factors in front of the continued fraction. */
  bt = exp(gammaln(a+b)-gammaln(a)-gammaln(b)+a*log(x)+b*log(1.0-x));
  if (x < (a+1.0)/(a+b+2.0))              /* Use continued fraction directly. */
    return bt*betacf(a,b,x)/a;
  else                    /* Use continued faction after making */
    return 1.0-bt*betacf(b,a,1.0-x)/b;              /* the symmetry transformation. */
}

/*********************************************************************
   Continued fraction evaluation routine needed for the incomplete beta
   function, I_x(a,b).
   C.A. Bertulani        May/16/2000
*********************************************************************/
static double betacf(double a, double b, double x)
/* Used by betai: Evaluates continued fraction for incomplete beta function
   by modified Lentz's method.   */
{
#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

  int m,m2;
  double aa,c,d,del,h,qab,qam,qap;

  qab=a+b;
  qap=a+1.0;
  qam=a-1.0;
  c=1.0;          /* First step of Lentz's method.  */
  d=1.0-qab*x/qap;
  if (fabs(d) < FPMIN) d=FPMIN;
  d=1.0/d;
  h=d;
  for( m=1; m<=MAXIT; m++ )
  {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d;             /* One step (the even one) of the recurrence. */
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d;             /* Next step of the recurence the odd one) */
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;   /* Are we done? */
  }
  if (m > MAXIT)
    return -1.0; // Error
  return h;
}

double gammaln(double xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,
                        -86.50532032941677,
                        24.01409824083091,
                        -1.231739572450155,
                        0.1208650973866179e-2,
                        -0.5395239384953e-5};
  int j;
  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.000000000190015;
  for(j = 0; j <= 5; j++)
    ser += cof[j] / ++y;
  return -tmp + log(2.5066282746310005 * ser / x);
}
