
#include "rwrapper.h"
/* #include <R.h> */
/* #include <Rmath.h> */

/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/

void RWRAPPER_Min_Number_Of_Tip_Permut(char **tree_file_name, char **coord_file_name, phydbl *res)
{
  t_tree *tree;
  option *io;
  FILE *fp_tree_file, *fp_coord_file;
  int i;


  srand(time(NULL)); rand();

/*   Rprintf("%s\n",tree_file_name[0]); */
/*   Rprintf("%s\n",coord_file_name[0]); */

  fp_tree_file  = (FILE *)fopen(tree_file_name[0],"r");
  fp_coord_file = (FILE *)fopen(coord_file_name[0],"r");

  io = (option *)Make_Input();
  io->fp_in_tree = fp_tree_file;
  Read_Tree_File(io);
  fclose(io->fp_in_tree);
  tree = io->treelist->tree[0];
  tree->io = io;

  tree->io->z_scores = (phydbl *)mCalloc(tree->n_otu,sizeof(phydbl));

  For(i,tree->n_otu) tree->io->z_scores[i] = TIPO_Read_One_Taxon_Zscore(fp_coord_file,tree->noeud[i]->name,1,tree);
  TIPO_Normalize_Zscores(tree);
  TIPO_Get_Min_Number_Of_Tip_Permut(tree);
  res[0] = (phydbl)tree->tip_order_score;

/*   Rprintf("\n>> y score: %f",tree->tip_order_score); */

  For(i,tree->n_otu) tree->io->z_scores[i] = TIPO_Read_One_Taxon_Zscore(fp_coord_file,tree->noeud[i]->name,2,tree);
  TIPO_Normalize_Zscores(tree);
  TIPO_Get_Min_Number_Of_Tip_Permut(tree);
  res[1] = (phydbl)tree->tip_order_score;

/*   Rprintf("\n>> x score: %f",tree->tip_order_score); */

  fclose(fp_tree_file);
  fclose(fp_coord_file);

}

/*********************************************************/

/* void RWRAPPER_Log_Dnorm(phydbl *x, phydbl *mean, phydbl *sd,  phydbl *res) */
/* { */
/*   *res = Log_Dnorm(*x,*mean,*sd); */
/* } */

/* /\*********************************************************\/ */

/* void RWRAPPER_Rnorm_Trunc(phydbl *mean, phydbl *sd, phydbl *min, phydbl *max, phydbl *res) */
/* { */
/*   *res = Rnorm_Trunc(*mean,*sd,*min,*max); */
/* } */

/* void  RWRAPPER_Cholesky_Decomp(double *A, int *dim) */
/* { */
/*   Cholesky_Decomp(A,*dim); */
/* } */

/*********************************************************/

/* void RWRAPPER_Bivariate_Normal_Density(phydbl *x, phydbl *y, phydbl *mux, phydbl *muy, phydbl *sdx, phydbl *sdy, phydbl *rho, phydbl *dens) */
/* { */
/*   *dens = Bivariate_Normal_Density(*x,*y,*mux,*muy,*sdx,*sdy,*rho); */
/* } */

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


