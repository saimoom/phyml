/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#ifndef STATS_H
#define STATS_H

#include "utilities.h"
#include "free.h"
#include "lk.h"
#include "optimiz.h"
#include "models.h"
#include "eigen.h"


phydbl *Covariance_Matrix(arbre *tree);
phydbl *Hessian(arbre *tree);
void   Recurr_Hessian(node *a, node *b, int plus_minus, phydbl *inc, phydbl *res, int *is_ok, arbre *tree);
double stdnormal_inv(double p);
double Uni();
int    Rand_Int(int min, int max);
double Ahrensdietergamma(double alpha);
double Rgamma(double shape, double scale);
double Rexp(double lambda);
phydbl Bico(int n, int k);
phydbl Factln(int n);
phydbl Gammln(phydbl xx);
phydbl Pbinom(int N, int ni, phydbl p);
phydbl LnGamma (phydbl alpha);
phydbl IncompleteGamma(phydbl x, phydbl alpha, phydbl ln_gamma_alpha);
phydbl PointChi2 (phydbl prob, phydbl v);
phydbl Bivariate_Normal_Density(phydbl x, phydbl y, phydbl mux, phydbl muy, phydbl sdx, phydbl sdy, phydbl rho);
phydbl PointNormal (phydbl prob);
int    DiscreteGamma (phydbl freqK[], phydbl rK[],phydbl alfa, phydbl beta, int K, int median);
phydbl CDF_Normal(phydbl x, phydbl mean, phydbl var);
phydbl Dnorm_Moments(phydbl x, phydbl mean, phydbl var);
phydbl Dnorm(phydbl x, phydbl mean, phydbl sd);
phydbl CDF_Gamma(phydbl x, phydbl shape, phydbl scale);
phydbl Dgamma_Moments(phydbl x, phydbl mean, phydbl var);
phydbl Dgamma(phydbl x, phydbl shape, phydbl scale);
phydbl LnFact(int n);
int    Choose(int n, int k);
phydbl CDF_Pois(phydbl x, phydbl param);
phydbl Dexp(phydbl x, phydbl param);
phydbl Dpois(phydbl x, phydbl param);
phydbl Rand_Normal_Deviate(phydbl mean, phydbl sd);
phydbl Rnorm(phydbl mean, phydbl sd);
phydbl *Rnorm_Multid(phydbl *mu, phydbl *cov, int dim);
phydbl Rnorm_Trunc(phydbl mean, phydbl sd, phydbl min, phydbl max, int *err);
phydbl *Rnorm_Multid_Trunc(phydbl *mean, phydbl *cov, phydbl *min, phydbl *max, int dim);
phydbl *Hessian_Log(arbre *tree);
void Recurr_Hessian_Log(node *a, node *d, int plus_minus, phydbl *inc, phydbl *res, int *is_ok, arbre *tree);
phydbl Log_Det(int *is_ok, arbre *tree);
phydbl Dnorm_Trunc(phydbl x, phydbl mean, phydbl sd, phydbl lo, phydbl up);
phydbl Normal_Trunc_Mean(phydbl mu, phydbl sd, phydbl min, phydbl max);
phydbl Constraint_Normal_Trunc_Mean(phydbl wanted_mu, phydbl sd, phydbl min, phydbl max);
phydbl Dnorm_Multi(phydbl *x, phydbl *mu, phydbl *cov, int size, int _log);
phydbl Dnorm_Multi_Given_InvCov_Det(phydbl *x, phydbl *mu, phydbl *invcov, phydbl det, int size, int _log);
phydbl Prop_Log_Dnorm_Multi_Given_InvCov_Det(phydbl *x, phydbl *mu, phydbl *invcov, phydbl det, int size);
phydbl Log_Dnorm(phydbl x, phydbl mean, phydbl sd);
phydbl tt800();
double Pnorm_Ihaka_Derived_From_Cody(double x);
int Matinv(double *x, int n, int m);
phydbl *Matrix_Mult(phydbl *A, phydbl *B, int nra, int nca, int nrb, int ncb);
phydbl *Matrix_Transpose(phydbl *A, int dim);
void Normal_Conditional(phydbl *mu, phydbl *cov, phydbl *a, int n, short int *is_1, int n1, phydbl *cond_mu, phydbl *cond_var);
void Normal_Conditional_Unsorted(phydbl *mu, phydbl *cov, phydbl *a, int n, short int *is_1, int n1, phydbl *cond_mu, phydbl *cond_cov);
phydbl Matrix_Det(phydbl *A, int size);
void Get_Reg_Coeff(phydbl *mu, phydbl *cov, phydbl *a, int n, short int *is_1, int n1, phydbl *reg_coeff);


#endif
