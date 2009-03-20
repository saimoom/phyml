#include "utilities.h"
#include "rates.h"
#include "eigen.h"
#include "numeric.h"

/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/

void RWRAPPER_Log_Dnorm(phydbl *x, phydbl *mean, phydbl *sd,  phydbl *res)
{
  *res = Log_Dnorm(*x,*mean,*sd);
}

/*********************************************************/

void RWRAPPER_Rnorm_Trunc(phydbl *mean, phydbl *sd, phydbl *min, phydbl *max, phydbl *res)
{
  *res = Rnorm_Trunc(*mean,*sd,*min,*max);
}

void  RWRAPPER_Cholesky_Decomp(double *A, int *dim)
{
  Cholesky_Decomp(A,*dim);
}

/*********************************************************/

void RWRAPPER_Bivariate_Normal_Density(phydbl *x, phydbl *y, phydbl *mux, phydbl *muy, phydbl *sdx, phydbl *sdy, phydbl *rho, phydbl *dens)
{
  *dens = Bivariate_Normal_Density(*x,*y,*mux,*muy,*sdx,*sdy,*rho);
}

/* /\*********************************************************\/ */

/* void RWRAPPER_Dmu2_And_Mu1_Given_Min_N(phydbl *mu1, phydbl *mu2, phydbl *dt1, phydbl *dt2, int *n_min, phydbl *a, phydbl *b, phydbl *lexp, phydbl *dens) */
/* { */
/*   *dens = RATES_Dmu2_And_Mu1_Given_Min_N(*mu1, *mu2, *dt1, *dt2, *n_min, *a, *b, *lexp); */
/* } */
/* /\*********************************************************\/ */

/* void RWRAPPER_Dgamma(phydbl *x, phydbl *shape, phydbl *scale, phydbl *dens) */
/* { */
/*   *dens = Dgamma(*x,*shape,*scale); */
/* } */

/* /\*********************************************************\/ */

/* void RWRAPPER_Dnorm(phydbl *x, phydbl *mean, phydbl *var, double *dens) */
/* { */
/*   *dens = Dnorm_Moments(*x,*mean,*var); */
/* } */
/* /\*********************************************************\/ */

/* void RWRAPPER_Dmu_One(phydbl *mu, phydbl *dt, phydbl *a, phydbl *b, phydbl *lexp, double *dens) */
/* { */
/*   *dens = RATES_Dmu_One(*mu,*dt,*a,*b,*lexp); */
/* } */
/* /\*********************************************************\/ */

/* void RWRAPPER_Dr_X_Dx(double *r, double *mu, double *y, double *dt, double *a, double *b, double *lexp, double *dens) */
/* { */
/*   *dens = RATES_Dr_X_Dx(*r,*mu,*y,*dt,*a,*b,*lexp); */
/* } */

/* /\*********************************************************\/ */

/* void RWRAPPER_Dmu_Given_Y(double *mu, double *y, double *dt, double *a, double *b, double *lexp, double *dens) */
/* { */
/*   *dens = RATES_Dmu_Given_Y(*mu,*y,*dt,*a,*b,*lexp); */
/* } */

/* /\*********************************************************\/ */

/* void RWRAPPER_Dmu2_And_Mu1(double *mu1, double *mu2, double *dt1, double *dt2, double *a, double *b, double *lexp, double *dens) */
/* { */
/*   *dens = RATES_Dmu2_And_Mu1(*mu1,*mu2,*dt1,*dt2,*a,*b,*lexp); */
/* } */

/* /\*********************************************************\/ */

/* void RWRAPPER_Dmu2_Given_Mu1(double *mu1, double *mu2, double *dt1, double *dt2, double *a, double *b, double *lexp, double *dens) */
/* { */
/* /\*   *dens = RATES_Dmu2_Given_Mu1(*mu1,*mu2,*dt1,*dt2,*a,*b,*lexp); *\/ */
/*   *dens = RATES_Dmu2_Given_Mu1_Bis(*mu1,*mu2,*dt1,*dt2,*a,*b,*lexp); */
/* } */

/* /\*********************************************************\/ */

