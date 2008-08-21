/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/


/* Routines for molecular clock trees and molecular dating */
#if defined (MC) ||  defined (RWRAPPER)

#include "spr.h"
#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "bionj.h"
#include "models.h"
#include "free.h"
#include "options.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"
#include "mc.h"
#include "m4.h"
#include "draw.h"
#include "rates.h"

#ifdef RWRAPPER
#include <R.h>
#endif



/*********************************************************/

phydbl RATES_Lk_Rates(arbre *tree)
{
/*   PhyML_Printf("\n>.< LnL rates = %f",lnL); */
  phydbl dens,rate0,rate1,dt0,dt1,lnL;
  
  RATES_Set_Node_Times(tree);
  RATES_Init_Triplets(tree);
  if(tree->rates->model == EXPONENTIAL)
    RATES_Adjust_Clock_Rate(tree);

  dens = UNLIKELY;
  tree->rates->c_lnL = .0;

  RATES_Lk_Rates_Pre(tree->n_root,tree->n_root->v[0],NULL,tree);
  RATES_Lk_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);

  dt0   = (tree->rates->t[tree->n_root->v[0]->num] - tree->rates->t[tree->n_root->num]);
  dt1   = (tree->rates->t[tree->n_root->v[1]->num] - tree->rates->t[tree->n_root->num]);
  if(dt0 < MIN_DT) dt0 = MIN_DT;
  if(dt1 < MIN_DT) dt1 = MIN_DT;
  rate0 = tree->n_root->l[0] / (tree->rates->clock_r * dt0);
  rate1 = tree->n_root->l[1] / (tree->rates->clock_r * dt1);

  if(dt0 < 0.0 || dt1 < 0.0)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(rate0 < 1.E-20 || rate1 < 1.E-20)
    {
      printf("\n. r0=%G r1=%G dt0=%f dt1=%f l=%f l0=%f l1=%f",rate0,rate1,dt0,dt1,tree->e_root->l,tree->n_root->l[0],tree->n_root->l[1]);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  lnL = tree->rates->c_lnL;

  if((tree->rates->model == COMPOUND_COR) ||
     (tree->rates->model == COMPOUND_NOCOR))
    {
      dens = RATES_Lk_Rates_Core(rate0,rate1,dt0,dt1,tree);
      dens *= RATES_Dmu(rate0,dt0,tree->rates->alpha,1./tree->rates->alpha,tree->rates->lexp,0);
    }
  else if(tree->rates->model == EXPONENTIAL)
    {
      dens = Dexp(rate0,tree->rates->nu) * Dexp(rate1,tree->rates->nu);
    }

  tree->rates->c_lnL += log(dens);
  tree->rates->triplet[tree->n_root->num] = log(dens);

/*   printf("\n\n. LNL = %f",tree->rates->c_lnL); */

  return tree->rates->c_lnL;
}

/*********************************************************/

void RATES_Lk_Rates_Pre(node *a, node *d, edge *b, arbre *tree)
{
  if(d->tax) return;
  else
    {
      int i;
      phydbl dens,mu1,mu2,dt1,dt2;
      
      mu1 = -1.;

      if(!b)
	{
	  dt1 = tree->rates->t[d->num] - tree->rates->t[tree->n_root->num];
	  if(dt1 < MIN_DT) dt1 = MIN_DT;
	  if(d == tree->n_root->v[0]) mu1 = tree->n_root->l[0] / (dt1 * tree->rates->clock_r);
	  else if(d == tree->n_root->v[1]) mu1 = tree->n_root->l[1] / (dt1 * tree->rates->clock_r);
	  else
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	}
      else
	{
	  dt1 = tree->rates->t[d->num] - tree->rates->t[a->num];
	  if(dt1 < MIN_DT) dt1 = MIN_DT;
	  mu1 = b->l / (dt1 * tree->rates->clock_r);
	}

      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      dt2 = tree->rates->t[d->v[i]->num] - tree->rates->t[d->num];
	      if(dt2 < MIN_DT) dt2 = MIN_DT;
	      mu2 = d->b[i]->l / (dt2 * tree->rates->clock_r);

	      dens = RATES_Lk_Rates_Core(mu1,mu2,dt1,dt2,tree);

	      int hit = 0;
	      if(dens < MDBL_MIN) {hit = 1; dens = 1.E-10;}

	      tree->rates->c_lnL += log(dens);

/* 	      printf("\n. LNL = %G",tree->rates->c_lnL); */

	      tree->rates->triplet[d->num]       += log(dens);
/* 	      tree->rates->triplet[d->v[i]->num] += log(dens); */

	      RATES_Lk_Rates_Pre(d,d->v[i],d->b[i],tree);
	    }
	}
    }
}

/*********************************************************/

phydbl RATES_Lk_Change_One_Rate(edge *b, phydbl new_rate, arbre *tree)
{
  phydbl dt;
  
  if(b == tree->e_root)
    {
      printf("\n. Can't change rate on root...");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
      
  dt = fabs(tree->rates->t[b->left->num] - tree->rates->t[b->rght->num]);
  if(dt < MIN_DT) dt = MIN_DT;
  b->l = dt * new_rate * tree->rates->clock_r; 
  

  RATES_Update_Triplet(b->left,tree);
  RATES_Update_Triplet(b->rght,tree);

  return(tree->rates->c_lnL);
}

/*********************************************************/

phydbl RATES_Lk_Change_One_Time(node *n, phydbl new_t, arbre *tree)
{  
  if(n == tree->n_root)
    {
      printf("\n. Moving the time of the root node is not permitted.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  else
    {
      int i;
      
      tree->rates->t[n->num] = new_t;
/*       RATES_Set_Node_Times(tree);       */

      RATES_Update_Triplet(n,tree);
      
      For(i,3)
	{
	  if(n->b[i] != tree->e_root)
	    {
	      RATES_Update_Triplet(n->v[i],tree);
	    }
	  else
	    {
	      RATES_Update_Triplet(tree->n_root,tree);
	    }
	}
    }
  return(tree->rates->c_lnL);
}

/*********************************************************/

void RATES_Update_Triplet(node *n, arbre *tree)
{
  phydbl curr_triplet,new_triplet;
  phydbl dt0,dt1,dt2;
  phydbl mu1_mu0,mu2_mu0;
  phydbl l0,l1,l2;
  phydbl mu0,mu1,mu2;
  int i;
  node *v1,*v2;


  if(n->tax) return;

  curr_triplet = tree->rates->triplet[n->num];

  l0 = l1 = l2 = -1.0;
  dt0 = dt1 = dt2 = -100.0;

  if(n == tree->n_root)
    {
      phydbl dens;

      l0 = tree->n_root->l[0];
      l1 = tree->n_root->l[1];
      
      dt0 = tree->rates->t[tree->n_root->v[0]->num] - tree->rates->t[tree->n_root->num];
      dt1 = tree->rates->t[tree->n_root->v[1]->num] - tree->rates->t[tree->n_root->num];
      
      if(dt0 < MIN_DT) dt0 = MIN_DT;
      if(dt1 < MIN_DT) dt1 = MIN_DT;

      mu0 = l0 / (dt0 * tree->rates->clock_r);
      mu1 = l1 / (dt1 * tree->rates->clock_r);

      
      if((tree->rates->model == COMPOUND_COR) ||
	 (tree->rates->model == COMPOUND_NOCOR))
	{
	  dens = RATES_Lk_Rates_Core(mu0,mu1,dt0,dt1,tree);
	  dens *= RATES_Dmu(mu0,dt0,tree->rates->alpha,1./tree->rates->alpha,tree->rates->lexp,0);
	}
      else
	{
	  dens = Dexp(mu0,tree->rates->nu) * Dexp(mu1,tree->rates->nu);
	}

      new_triplet = log(dens);
    }
  else
    {
      v1 = v2 = NULL;
      For(i,3)
	{
	  if((n->v[i] != n->anc) && (n->b[i] != tree->e_root))
	    {
	      if(!v1) 
		{
		  v1  = n->v[i]; 
		  dt1 = tree->rates->t[v1->num] - tree->rates->t[n->num];
		  l1  = n->b[i]->l;
		}
	      else
		{
		  v2  = n->v[i]; 
		  dt2 = tree->rates->t[v2->num] - tree->rates->t[n->num];
		  l2  = n->b[i]->l;
		}
	    }
	  else 
	    {
	      if(n->b[i] == tree->e_root)
		{
		  dt0 = tree->rates->t[n->num] - tree->rates->t[tree->n_root->num];
		  if(n == tree->n_root->v[0])      l0 = tree->n_root->l[0];
		  else if(n == tree->n_root->v[1]) l0 = tree->n_root->l[1];
		  else 
		    {
		      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		      Warn_And_Exit("");
		    }
		}	      
	      else
		{
		  dt0 = tree->rates->t[n->num] - tree->rates->t[n->v[i]->num];
		  l0 = n->b[i]->l;
		}
	    }
	}

      if(dt0 < MIN_DT) dt0 = MIN_DT;
      if(dt1 < MIN_DT) dt1 = MIN_DT;
      if(dt2 < MIN_DT) dt2 = MIN_DT;

      mu0 = l0 / (dt0 * tree->rates->clock_r);
      mu1 = l1 / (dt1 * tree->rates->clock_r);
      mu2 = l2 / (dt2 * tree->rates->clock_r);
      
      mu1_mu0 = RATES_Lk_Rates_Core(mu0,mu1,dt0,dt1,tree);
      mu2_mu0 = RATES_Lk_Rates_Core(mu0,mu2,dt0,dt2,tree);
      
      new_triplet = log(mu1_mu0) + log(mu2_mu0);
    }
  tree->rates->c_lnL = tree->rates->c_lnL + new_triplet - curr_triplet;
  tree->rates->triplet[n->num] = new_triplet;
}

/*********************************************************/
/* Returns f(mu2;mu1) */
phydbl RATES_Lk_Rates_Core(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, arbre *tree)
{
  phydbl dens;
  int k;
  phydbl alpha, beta, lexp;

  lexp = tree->rates->lexp;
  alpha = tree->rates->alpha;
  beta = 1./alpha;
  dens = UNLIKELY;

  /* Discretize rates */
  if(mu1 > 1.0)
    {
      k = (int)((mu1 - tree->rates->min_rate)/tree->rates->step_rate);
      mu1 = tree->rates->min_rate + (phydbl)k*tree->rates->step_rate;
    }
  else
    {
      mu1 = 1./mu1;
      k = (int)((mu1 - tree->rates->min_rate)/tree->rates->step_rate);
      mu1 = tree->rates->min_rate + (phydbl)k*tree->rates->step_rate;
      mu1 = 1./mu1;
    }

  if(mu2 > 1.0)
    {
      k = (int)((mu2 - tree->rates->min_rate)/tree->rates->step_rate);
      mu2 = tree->rates->min_rate + (phydbl)k*tree->rates->step_rate;
    }
  else
    {
      mu2 = 1./mu2;
      k = (int)((mu2 - tree->rates->min_rate)/tree->rates->step_rate);
      mu2 = tree->rates->min_rate + (phydbl)k*tree->rates->step_rate;
      mu2 = 1./mu2;
    }


  if(mu1 < tree->rates->min_rate) mu1 = tree->rates->min_rate;
  if(mu1 > tree->rates->max_rate) mu1 = tree->rates->max_rate;
  
  if(mu2 < tree->rates->min_rate) mu2 = tree->rates->min_rate;
  if(mu2 > tree->rates->max_rate) mu2 = tree->rates->max_rate;
 
 if(dt1 < MIN_DT || dt2 < MIN_DT)
   {
      printf("\n. dt1=%f dt2=%f",dt1,dt2);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }


 switch(tree->rates->model)
   {
   case COMPOUND_COR:
     {       
       dens = RATES_Compound_Core(mu1,mu2,dt1,dt2,alpha,beta,lexp,tree->rates->step_rate,tree->rates->approx);
       break;
    }

   case COMPOUND_NOCOR :
     {
       dens = RATES_Dmu(mu2,dt2,alpha,beta,lexp,0);
       break;
     }
						
   case EXPONENTIAL :
     {
       dens = Dexp(mu2,tree->rates->nu);
       break;
     }

   default : 
     {
       PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
       Warn_And_Exit("");
     }     
   }
 return dens;
}

/*********************************************************/

phydbl RATES_Compound_Core(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl alpha, phydbl beta, phydbl lexp, phydbl eps, int approx)
{
  phydbl p0,p1,v0,v1,v2;
  phydbl dmu1;
  int    equ;
  phydbl dens;
  
  v0 = v1 = v2 = 0.0;
  
  /* Probability of 0 and 1 jumps */
  p0   = Dpois(0,lexp*(dt2+dt1));       
  p1   = Dpois(1,lexp*(dt2+dt1));
  
  dmu1 = RATES_Dmu(mu1,dt1,alpha,beta,lexp,0);
  
  /* Are the two rates equal ? */
  equ = 0;
  if(mu1 < 1.0 && mu2 < 1.0)
    {
      mu1 = 1./mu1;
      mu2 = 1./mu2;
      
      if(fabs(mu1-mu2) < eps) equ = 1;
      
      mu1 = 1./mu1;
      mu2 = 1./mu2;
    }
  else if(mu1 > 1.0 && mu2 > 1.0)
    {
      if(fabs(mu1-mu2) < eps) equ = 1;
    }
  
  /* No jump */
      if(equ)
    {
      v0 = 1.0*Dgamma(mu1,alpha,beta)/dmu1;
/*       Rprintf("\n. mu1=%f mu2=%f",mu1,mu2); */
    }
  else
    {
      v0 = 0.0;
    }

  
  /* One jump */
  v1 = RATES_Dmu1_And_Mu2_One_Jump_Two_Intervals(mu1,mu2,dt1,dt2,alpha,beta);
  v1 /= dmu1;
  
  
  /* Two jumps and more (approximation) */
  if(approx == 1)
    {
      v2 =
	RATES_Dmu(mu1,dt1,alpha,beta,lexp,0)*RATES_Dmu(mu2,dt2,alpha,beta,lexp,0) -
	Dpois(0,lexp*dt1) * Dpois(0,lexp*dt2) *
	Dgamma(mu1,alpha,beta) * Dgamma(mu2,alpha,beta);
    }
  else
    {
      v2 = 
	RATES_Dmu_One(mu1,dt1,alpha,beta,lexp) * 
	RATES_Dmu_One(mu2,dt2,alpha,beta,lexp);

      v2 += Dpois(0,lexp*dt1)*Dgamma(mu1,alpha,beta)*RATES_Dmu(mu2,dt2,alpha,beta,lexp,1);
      v2 += Dpois(0,lexp*dt2)*Dgamma(mu2,alpha,beta)*RATES_Dmu(mu1,dt1,alpha,beta,lexp,1);

    }
/*   printf("\n. %f %f %f %f %f ",mu1,mu2,dt1,dt2,v2); */
  v2 /= dmu1;

  dens = p0*v0 + p1*v1 + v2;
/*   dens = p1*v1 + v2; */
  /*       dens = p1*v1 + v2; */
  /*   dens = v0; */
/*   dens *= dmu1; */
  
  if(dens < 0.0)
    {
      printf("\n. dens=%12G mu1=%12G mu2=%12G dt1=%12G dt2=%12G lexp=%12G alpha=%f v0=%f v1=%f v2=%f p0=%f p1=%f p2=%f",
	     dens,
	     mu1,mu2,dt1,dt2,
	     lexp,
	     alpha,
	     v0,v1,v2,
	     p0,p1,1.-p0-p1);
    }
  
  return dens;

}

/**********************************************************/

void RATES_Print_Triplets(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-1) printf("\n. Node %3d t=%f",i,tree->rates->triplet[i]);
}


/**********************************************************/

void RATES_Print_Rates(arbre *tree)
{
  RATES_Print_Rates_Pre(tree->n_root,tree->n_root->v[0],NULL,tree);
  RATES_Print_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);
}

/*********************************************************/

void RATES_Print_Rates_Pre(node *a, node *d, edge *b, arbre *tree)
{

  if(b)
    printf("\n. Edge %3d l=%f rate=%f t_left=%f t_rght=%f",
	   b->num,
	   b->l,
	   b->l / (fabs(tree->rates->t[b->left->num] - tree->rates->t[b->rght->num]) * tree->rates->clock_r),
	   tree->rates->t[b->left->num],tree->rates->t[b->rght->num]);

  if(d->tax) return;
  else
    {
      int i;

      For(i,3) 
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Print_Rates_Pre(d,d->v[i],d->b[i],tree);
	    }
	}
    }
}

/*********************************************************/

phydbl RATES_Check_Mean_Rates(arbre *tree)
{
  phydbl sum;

  sum = 0.0;
  RATES_Check_Mean_Rates_Pre(tree->n_root,tree->n_root->v[0],NULL,&sum,tree);
  RATES_Check_Mean_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,&sum,tree);

  sum += tree->n_root->l[0] / (fabs(tree->rates->t[tree->n_root->num] - tree->rates->t[tree->n_root->v[0]->num]) * tree->rates->clock_r);
  sum += tree->n_root->l[1] / (fabs(tree->rates->t[tree->n_root->num] - tree->rates->t[tree->n_root->v[1]->num]) * tree->rates->clock_r);
  
  return(sum/(phydbl)(2*tree->n_otu-2));
}

/*********************************************************/

void RATES_Check_Mean_Rates_Pre(node *a, node *d, edge *b, phydbl *sum, arbre *tree)
{
  if(b) *sum += b->l / (fabs(tree->rates->t[b->left->num] - tree->rates->t[b->rght->num]) * tree->rates->clock_r);

  if(d->tax) return;
  else
    {
      int i;

      For(i,3) 
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Check_Mean_Rates_Pre(d,d->v[i],d->b[i],sum,tree);
	    }
	}
    }
}

/*********************************************************/

trate *RATES_Make_Rate_Struct(int n_otu)
{
  trate *rates;
  
  rates          = (trate *)mCalloc(1,sizeof(trate));
  rates->br_r    = (phydbl *)mCalloc(2*n_otu-3,sizeof(phydbl));
  rates->mc_mr   = (phydbl **)mCalloc(2*n_otu-1,sizeof(phydbl *));
  rates->t       = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->true_t  = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->dens    = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
  rates->triplet = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  
  return rates;
}

/*********************************************************/

void RATES_Init_Rate_Struct(trate *rates, int n_otu)
{
  int i;

  rates->model         = COMPOUND_COR;
  rates->n_mc_runs     = 1000; 
  rates->clock_r       = 1.E-03;
/*   rates->clock_r       = 1.E-2; */
  rates->curr_mc_run   = 0;
  rates->c_lnL         = -INFINITY;
  rates->adjust_rates  = 0;
  rates->use_rates     = 1;
  rates->lexp          = 0.13;
  rates->alpha         = 2.0;
  rates->birth_rate    = 0.001;
  rates->max_rate      = 10.;
  rates->min_rate      = 0.1;
  rates->step_rate     = 1.E-1;
  rates->nu            = 1.0;
  rates->approx        = 1;

  For(i,2*n_otu-2) rates->br_r[i] = 1.0;

  For(i,2*n_otu-1) 
    {
      rates->t[i]      = 0.0;
      rates->true_t[i] = 0.0;
    }


  For(i,2*n_otu-1)
    rates->mc_mr[i] = (phydbl *)mCalloc(rates->n_mc_runs,sizeof(phydbl));

}

/*********************************************************/

void RATES_Fill_Node_Rates(phydbl *node_r, arbre *tree)
{

  node_r[tree->n_root->v[0]->num] = tree->rates->br_r[tree->e_root->num];
  node_r[tree->n_root->v[1]->num] = tree->rates->br_r[tree->e_root->num];

  RATES_Fill_Node_Rates_Pre(tree->n_root,tree->n_root->v[0],NULL,node_r,tree);
  RATES_Fill_Node_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,node_r,tree);
}

/*********************************************************/

void RATES_Fill_Node_Rates_Pre(node *a, node *d, edge *b, phydbl *node_r, arbre *tree)
{
  if(b) node_r[d->num] = tree->rates->br_r[b->num];

  if(d->tax) return;
  else
    {
      int i;
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Fill_Node_Rates_Pre(d,d->v[i],d->b[i],node_r,tree);
	    }
	}
    }
}

/*********************************************************/

phydbl RATES_Exp_Y(phydbl mu1, phydbl mu2, phydbl dt1, phydbl lexp)
{
  phydbl prob0jump;
  prob0jump = Dpois(0,lexp*dt1); /* probability of having 0 jump in dt1 */        
/*   return prob0jump * mu1 + (1.-prob0jump) * mu2; */
  return prob0jump * mu1 + (1.-prob0jump) * 1.0;
}

/*********************************************************/

void RATES_Bracket_N_Jumps(int *up, int *down, phydbl param)
{
  phydbl cdf,eps,a,b,c;
  int step;

  step = 10;
  eps = 1.E-10;
  cdf = 0.0;
  c = 1;
  
  while(cdf < 1.-eps)
    {
      c = (int)floor(c * step);
      cdf = CDF_Pois(c,param);      
    }
  
  a = 0.0;
  b = (c-a)/2.;
  step = 0;
  do
    {
      step++;
      cdf = CDF_Pois(b,param);
      if(cdf < eps) a = b;
      else 
	{
	  break;
	}
      b = (c-a)/2.;
    }
  while(step < 1000);
  
  if(step == 1000)
    {
      PhyML_Printf("\n. a=%f b=%f c=%f param=%f",a,b,c,param);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  *up = c;
  *down = a;
}

/*********************************************************/
/* 
   mu   : average rate of the time period dt
   dt   : time period to be considered
   a    : rate at a given time point is gamma distributed. a is the shape parameter
   b    : rate at a given time point is gamma distributed. b is the scale parameter
   lexp : the number of rate switches is Poisson distributed with parameter lexp * dt
*/ 
/* compute f(mu;dt,a,b,lexp), the probability density of mu. We need to integrate over the
   possible number of jumps (n) during the time interval dt */
phydbl RATES_Dmu(phydbl mu, phydbl dt, phydbl a, phydbl b, phydbl lexp, int min_n)
{
  phydbl var,cumpoissprob,dens,mean,poissprob,ab2,gammadens,lexpdt,*suminv,b2;
  int n,up,down;

  var          = 0.0;
  cumpoissprob = 0.0;
  dens         = 0.0;
  n            = 0;
  mean         = a*b;
  ab2          = a*b*b;
  lexpdt       = lexp*dt;  
  b2           = b*b;
  suminv       = NULL;

  if(dt < 0.0)
    {
      PhyML_Printf("\n. dt=%f",dt);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(lexpdt < MDBL_MIN)
    {
      PhyML_Printf("\n. lexpdt=%G lexp=%G dt=%G",lexpdt,lexp,dt);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  RATES_Bracket_N_Jumps(&up,&down,lexpdt);
  For(n,MAX(down,min_n)-1) cumpoissprob += Dpois(n,lexpdt);

  for(n=MAX(down,min_n);n<up+1;n++)
    {
      poissprob    = Dpois(n,lexpdt); /* probability of having n jumps */      
      var          = (2./(n+2.))*ab2; /* var(mu|n) = var(mu|n=0) * 2 / (n+2) */
      gammadens    = Dgamma_Moments(mu,mean,var);
      dens         += poissprob * gammadens;
      cumpoissprob += poissprob;
      if(cumpoissprob > 1.-1.E-06) break;
    }

  return(dens);
}

/*********************************************************/

phydbl RATES_Dmu_One(phydbl mu, phydbl dt, phydbl a, phydbl b, phydbl lexp)
{
  phydbl var,cumpoissprob,dens,mean,poissprob,ab2,gammadens,lexpdt,*suminv,b2;
  int n,up,down;
  
  var          = 0.0;
  cumpoissprob = 0.0;
  dens         = 0.0;
  n            = 0;
  mean         = a*b;
  ab2          = a*b*b;
  lexpdt       = lexp*dt;  
  b2           = b*b;
  suminv       = NULL;
  
  if(dt < 0.0)
    {
      PhyML_Printf("\n. dt=%f",dt);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  
  if(lexpdt < MDBL_MIN)
    {
      PhyML_Printf("\n. lexpdt=%G lexp=%G dt=%G",lexpdt,lexp,dt);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  
  RATES_Bracket_N_Jumps(&up,&down,lexpdt);

  For(n,MAX(1,down)-1) cumpoissprob += Dpois(n,lexpdt);

  for(n=MAX(1,down);n<up+1;n++) /* WARNING: we are considering that at least one jump occurs in the interval */
    {
      poissprob    = Dpois(n,lexpdt); /* probability of having n jumps */      
      var          = (n/((n+1)*(n+1)*(n+2)))*(pow(1-a*b,2) + 2/(n+1)*ab2) + 2*n*n*ab2/pow(n+1,3);      
      gammadens    = Dgamma_Moments(mu,mean,var);
      dens         += poissprob * gammadens;
      cumpoissprob += poissprob;
      if(cumpoissprob > 1.-1.E-06) break;
    }

  return(dens);
}

/*********************************************************/
/*
  0     y      t1          x            dt
  *------------*---*------*------*------*

  * are jumps. t1 is the time of the first jump.
  r = t1 / dt
  y is the rate in [0,t1]
  x is the average rate over [t1,dt]
  mu is the average rate over [0,dt]
  This function computes f(mu;r,y) x f(r;lexp,dt) 
*/ 

phydbl RATES_Dr_X_Dx(phydbl r, phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp)
{
  phydbl x,dr,dx;

  if(r > 1.-1.E-10) r = 1.-1.E-10;
  
  x = (mu-r*y)/(1-r);
  
  if(x < 0.0)
    {
      PhyML_Printf("\n. x=%f ,mu=%f r=%f y=%f dt=%f",x,mu,r,y,dt);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  
  if(x < 1.E-10) x = 1.E-10;
  
  dr = lexp*dt*exp(-lexp*r*dt); /* This is the probability density for r. K = lexp*exp(-lexp*r*dt) gives
				   the density of t1. As r = t1/T, the pdf for r is K*T. 
				*/
  dx = RATES_Dmu(x,(1-r)*dt,a,b,lexp,0); /* This is f(x;T*(1-r)) */
  
  return(dx*dr/(1-r)); /* The division by (1-r) comes from the transformation
			  of the random variable x to get the probability density
			  function of mu = r.y + (1-r).x.
		       */
}

/*********************************************************/

phydbl RATES_Dmu_Given_Y_Trpzd(phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp, 
			       int nsteps, phydbl beg, phydbl end, phydbl prevs)
{
  phydbl scurr;

  scurr = -1.E+10;
  if(nsteps == 1)
    {
#ifdef DEBUG
      if(beg*y > mu)
	{
	  PhyML_Printf("\n. beg=%f y=%f mu=%f",beg,y,mu);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      if(end*y > mu)
	{
	  PhyML_Printf("\n. end=%f y=%f mu=%f",end,y,mu);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
#endif

      
      scurr = 0.5*(end-beg)*(RATES_Dr_X_Dx(beg,mu,y,dt,a,b,lexp)+
			     RATES_Dr_X_Dx(end,mu,y,dt,a,b,lexp));
    }
  else
    {
      phydbl h,x,sum;
      int i,npoints;

      h = (end - beg)*pow(2,2-nsteps);
      x = beg+0.5*h;
      sum = 0.0;
      npoints = pow(2,nsteps-2); /* actual formula is 2^(n-1) - 2^(n-2) */      

      
      For(i,npoints)
	{
#ifdef DEBUG
	  if(x*y > mu)
	    {
	      PhyML_Printf("\n. x=%f y=%f mu=%f",x,y,mu);
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	  if(x > 1.-1.E-10 || x < 1.E-10)
	    {
	      PhyML_Printf("\n. x=%f",x);
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
#endif
	  
	  sum += RATES_Dr_X_Dx(x,mu,y,dt,a,b,lexp);
	  x = x+h;
	}
      scurr = prevs/2. + sum*(end-beg)/pow(2,nsteps-1);      
    }
  return scurr;
}

/*********************************************************/

phydbl RATES_Dmu_Given_Y(phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp)
{
  phydbl p_no_jump;

  p_no_jump = Dpois(0,lexp*(dt));

  if(p_no_jump > 0.990)
    {
      if(fabs(mu-y) < 0.001 * MIN(mu,y)) return 1.0;
      else return 0.0;
    }
  else
    {
      return RATES_Dmu_Given_Y_Romb(mu,y,dt,a,b,lexp);
    }
}

/*********************************************************/
/* f(mu;y) = \int_{r=0}^{r=1} f(u;y.r) x f(r) dr
   i.e., this function computes f(mu,y) by intyegrating over r in [0,1].
 */
phydbl RATES_Dmu_Given_Y_Std(phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp)
{
  phydbl s,olds,prevs,beg,end,eps;
  int i;

  s     = 0.0;
  olds  = 0.0;
  prevs = 0.0;

  eps = MIN_CLOCK_RATE;
  if(eps > mu)
    {
      printf("\n. eps=%G mu=%G",eps,mu);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  beg = eps;
  end = MIN(1.,(mu-eps)/y);
  
  for(i=1;i<20;i++)
    {
      s = RATES_Dmu_Given_Y_Trpzd(mu,y,dt,a,b,lexp,i,beg,end,prevs);
      prevs = s;
      if(fabs(s-olds) < 1.E-03*fabs(olds))  return(s);       
      olds = s;
    }
  return(-1.0E+10);
}

/*********************************************************/

/*********************************************************/
/* Romberg method of integration over the time of the first jump. 
   Adapted from "Numerical Recipes in C". 
   Press, Flannery, Teukolsky, Vetterling, 1988. 
*/
phydbl RATES_Dmu_Given_Y_Romb(phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp)
{

  phydbl *s,*h;
  phydbl ss,dss;
  phydbl beg,end,prevs;
  int i,K,NMAXITER;
  phydbl eps;

  NMAXITER = 100;

  s = (phydbl *)mCalloc(NMAXITER  ,sizeof(phydbl));
  h = (phydbl *)mCalloc(NMAXITER+1,sizeof(phydbl));

  K = 5;

  eps = MIN_CLOCK_RATE;
  if(eps > mu)
    {
      printf("\n. eps=%G mu=%G",eps,mu);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }


  beg = eps;
  end = MIN(1.,(mu-eps)/y);

  prevs = 0.0;
  h[1] = 1.0;
  for(i=1;i<=NMAXITER;i++)
    {
      s[i] = RATES_Dmu_Given_Y_Trpzd(mu,y,dt,a,b,lexp,i,beg,end,prevs);
      prevs = s[i];
      
      if(i >= K)
	{
	  Polint(&h[i-K],&s[i-K],K,0.0,&ss,&dss);
	  if(fabs(dss) <= 1.0E-03*fabs(ss)) 
	    {
	      Free(s); 
	      Free(h);
	      return(ss);
	    }
	}
      h[i+1] = 0.25*h[i];
    }
  
  Free(s);
  Free(h);

  return(ss);

}

/*********************************************************/

phydbl RATES_Dy_Given_Mu(phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp)
{
  phydbl A,B;
/*   phydbl C; */
  A = RATES_Dmu_Given_Y_Romb(mu,y,dt,a,b,lexp);
  B = Dgamma(y,a,b);
/*   C = RATES_Dmu(mu,dt,a,b,a,lexp); */
/*   return(A*B/C); */
  return(A*B);
}

/*********************************************************/
/*
  0                dt1               dt2+dt1
  *-----------------|-------------------*
         mu1        y        mu2

   P(mu2|y) x P(y|mu1) 
*/
phydbl RATES_Dmu2_Given_Y_X_Dy_Given_Mu1(phydbl mu1, phydbl mu2, phydbl y, phydbl dt1, phydbl dt2, phydbl a, phydbl b, phydbl lexp)
{
  phydbl A,B;
  B = RATES_Dy_Given_Mu(mu1,y,dt1,a,b,lexp);
  A = RATES_Dmu_Given_Y(mu2,y,dt2,a,b,lexp);
  return(A*B);
}

/*********************************************************/

phydbl RATES_Dmu2_Given_Mu1_Trpz(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl a, phydbl b, phydbl lexp,
				 int nsteps, phydbl beg, phydbl end, phydbl prevs)
{
  phydbl scurr;

  scurr = -1.E+10;
  if(nsteps == 1)
    {
      scurr = 0.5*(end-beg)*(RATES_Dmu2_Given_Y_X_Dy_Given_Mu1(mu1,mu2,beg,dt1,dt2,a,b,lexp)+
			     RATES_Dmu2_Given_Y_X_Dy_Given_Mu1(mu1,mu2,end,dt1,dt2,a,b,lexp));
    }
  else
    {
      phydbl h,x,sum;
      int i,npoints;

      h = (end - beg)*pow(2,2-nsteps);
      x = beg+0.5*h;
      sum = 0.0;
      npoints = pow(2,nsteps-2); /* actual formula is 2^(n-1) - 2^(n-2) */      

      For(i,npoints)
	{
	  sum += RATES_Dmu2_Given_Y_X_Dy_Given_Mu1(mu1,mu2,x,dt1,dt2,a,b,lexp);
	  x = x+h;
	}
      scurr = prevs/2. + sum*(end-beg)/pow(2,nsteps-1);      
    }
  return scurr;
}

/*********************************************************/
phydbl RATES_Dmu2_Given_Mu1(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl a, phydbl b, phydbl lexp)
{
  phydbl p_no_jump;
  
  p_no_jump = Dpois(0,lexp*(dt1+dt2));

/*   Rprintf("p_no_jump dt1+dt2 = %f\n",p_no_jump); */
/*   if(p_no_jump > 0.99) */
/*     { */
/*       if(fabs(mu1-mu2) < 0.01 * MIN(mu1,mu2)) return 1.0; */
/*       else return 0.0; */
/*     } */
/*   else */
/*     { */
      return(RATES_Dmu2_Given_Mu1_Romb(mu1,mu2,dt1,dt2,a,b,lexp));
/*     } */
}
/*********************************************************/
/* P(mu2/mu1) = \int_{y=0}^{y=infty} P(mu2|y) P(y|mu1) */
phydbl RATES_Dmu2_Given_Mu1_Romb(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl a, phydbl b, phydbl lexp)
{
  phydbl *s,*h;
  phydbl ss,dss;
  phydbl beg,end,prevs;
  int i,K,NMAXITER;

  NMAXITER = 100;

  s = (phydbl *)mCalloc(NMAXITER  ,sizeof(phydbl));
  h = (phydbl *)mCalloc(NMAXITER+1,sizeof(phydbl));

  K = 5;

  beg = 0.01;
  end = 10.;

  prevs = 0.0;
  h[1] = 1.0;
  for(i=1;i<=NMAXITER;i++)
    {
      s[i] = RATES_Dmu2_Given_Mu1_Trpz(mu1,mu2,dt1,dt2,a,b,lexp,i,beg,end,prevs);
      prevs = s[i];
      
/*       Rprintf("%f\n",s[i]); */

      if(i >= K)
	{
	  Polint(&h[i-K],&s[i-K],K,0.0,&ss,&dss);
	  if(fabs(dss) <= 1.E-03*fabs(ss)) 
	    {
	      Free(s); 
	      Free(h);
	      return(ss/RATES_Dmu(mu1,dt1,a,b,lexp,0));
	    }
	}
      h[i+1] = 0.25*h[i];
    }
  
  Free(s);
  Free(h);

  return(ss/RATES_Dmu(mu1,dt1,a,b,lexp,0));

}

/*********************************************************/

phydbl RATES_Dmu2_And_Mu1(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl a, phydbl b, phydbl lexp)
{
  phydbl A,B;

  A = RATES_Dmu2_Given_Mu1(mu1,mu2,dt1,dt2,a,b,lexp);
  B = RATES_Dmu(mu1,dt1,a,b,lexp,0);

  return(A*B);
}

/*********************************************************/

/* Given the times of nodes a (ta) and d (td), the shape of the gamma distribution of instantaneous
   rates, the parameter of the exponential distribution of waiting times between rate jumps and the 
   instantaneous rate at node a, this function works out an expected number of (amino-acids or 
   nucleotide) substitutions per site.
*/
phydbl RATES_Expect_Number_Subst(phydbl t_beg, phydbl t_end, phydbl *r_beg, phydbl alpha, phydbl lexp)
{
  phydbl curr_r, mean_r, curr_t, next_t;
  int n_jumps;

  curr_r = *r_beg;
  mean_r = *r_beg;
      
  curr_t = t_beg + Rexp(lexp); /* Exponentially distributed waiting times */
  next_t = curr_t;
  
  n_jumps = 0;
  while(curr_t < t_end)
    {
      curr_r = Rgamma(alpha,1./alpha); /* Gamma distributed random instantaneous rate */
      
      n_jumps++;
      
      next_t = curr_t + Rexp(lexp);
	  
      if(next_t < t_end)
	{
	  mean_r = (1./(next_t - t_beg)) * (mean_r * (curr_t - t_beg) + curr_r * (next_t - curr_t));
	}
      else
	{
	  mean_r = (1./(t_end - t_beg)) * (mean_r * (curr_t - t_beg) + curr_r * (t_end - curr_t));
	}
      curr_t = next_t;
    }

  printf("\n. [%3d]",n_jumps);

  *r_beg = curr_r;
  return mean_r;
}

/*********************************************************/


void RATES_Get_Mean_Rates_Pre(node *a, node *d, edge *b, phydbl *r_a, phydbl alpha, phydbl lexp, arbre *tree)
{
  phydbl a_t,d_t;

  if(b)
    {
      phydbl mean_r;

      a_t = tree->rates->t[a->num];
      d_t = tree->rates->t[d->num];
      
      mean_r = RATES_Expect_Number_Subst(a_t,d_t,r_a,alpha,lexp);

      tree->rates->br_r[b->num] = mean_r;
/*       printf("[%f]",mean_r); */
    }

  
  /*******/
/*   printf("\n. No correlation"); */
/*   *r_a = Rgamma(alpha,1./alpha); /\* Gamma distributed random instantaneous rate *\/ */
  /*******/


  /* Move to the next branches */
  if(d->tax) return;
  else
    {
      int i;
      
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Get_Mean_Rates_Pre(d,d->v[i],d->b[i],r_a,alpha,lexp,tree);
	    }
	}
    }
}

/*********************************************************/

void RATES_Random_Branch_Lengths(arbre *tree)
{
  phydbl r0,alpha,lexp,correction,mean_r0,mean_r1,a_t,d_t;
  int i;

  mean_r0 = mean_r1 = -1.;

   
  alpha = tree->rates->alpha;
  lexp  = tree->rates->lexp;
  
  printf("\n. alpha = %f lexp = %f",alpha,lexp);

  r0 = Rgamma(alpha,1./alpha);

  if(!tree->n_root)
    {
      RATES_Get_Mean_Rates_Pre(tree->noeud[0],tree->noeud[0]->v[0],tree->noeud[0]->b[0], 
			       &r0, alpha,lexp,tree);
    }
  else
    {
      phydbl r_root;

      r_root = r0;
      RATES_Get_Mean_Rates_Pre(tree->n_root,tree->n_root->v[0],NULL,&r_root,alpha,lexp,tree);

      r_root = r0;
      RATES_Get_Mean_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,&r_root,alpha,lexp,tree);

      a_t = tree->rates->t[tree->n_root->num];
      d_t = tree->rates->t[tree->n_root->v[0]->num];
      r_root = r0;
      mean_r0 = RATES_Expect_Number_Subst(a_t,d_t,&r_root,alpha,lexp);

      a_t = tree->rates->t[tree->n_root->num];
      d_t = tree->rates->t[tree->n_root->v[1]->num];
      r_root = r0;
      mean_r1 = RATES_Expect_Number_Subst(a_t,d_t,&r_root,alpha,lexp);
   }
  
  /* Correction of relative rates on branches such that they are centered on 1.0 */
  correction = 0;
  For(i,2*tree->n_otu-3)
    if(tree->t_edges[i] != tree->e_root)
      correction += tree->rates->br_r[tree->t_edges[i]->num];
  correction += (mean_r0 + mean_r1);
  correction /= 2*tree->n_otu-2;
  
  For(i,2*tree->n_otu-3) tree->rates->br_r[tree->t_edges[i]->num] /= correction;
  mean_r0 /= correction;
  mean_r1 /= correction;

  For(i,2*tree->n_otu-3)
    {
      if(tree->t_edges[i] != tree->e_root)
	{
	  tree->t_edges[i]->l = 
 	    fabs(tree->rates->t[tree->t_edges[i]->left->num] -
		 tree->rates->t[tree->t_edges[i]->rght->num]) * 
	    tree->rates->clock_r * tree->rates->br_r[tree->t_edges[i]->num];

/* 	  printf("\n. l=%f mean_r=%f dt=%f", */
/* 		 tree->t_edges[i]->l, */
/* 		 tree->rates->br_r[tree->t_edges[i]->num], 	     */
/* 		 fabs(tree->rates->t[tree->t_edges[i]->left->num] - */
/* 		      tree->rates->t[tree->t_edges[i]->rght->num])); */
	}
    }

  a_t = tree->rates->t[tree->n_root->num];
  d_t = tree->rates->t[tree->n_root->v[0]->num];
  tree->e_root->l = fabs(a_t - d_t) * mean_r0 * tree->rates->clock_r;
  tree->n_root->l[0] = fabs(a_t - d_t) * mean_r0 * tree->rates->clock_r;
  
  a_t = tree->rates->t[tree->n_root->num];
  d_t = tree->rates->t[tree->n_root->v[1]->num];
  tree->e_root->l += fabs(a_t - d_t) * mean_r1 * tree->rates->clock_r;
  tree->n_root->l[1] = fabs(a_t - d_t) * mean_r1 * tree->rates->clock_r;

/*   printf("\n. e_root->l=%f mean_r0=%f mean_r1=%f",tree->e_root->l,mean_r0,mean_r1); */
}

/*********************************************************/

void RATES_Randomize_Node_Times(arbre *tree)
{
  RATES_Randomize_Node_Times_Pre(tree->n_root,tree->n_root->v[0],tree);
  RATES_Randomize_Node_Times_Pre(tree->n_root,tree->n_root->v[1],tree);
}

/*********************************************************/

void RATES_Set_Node_Times(arbre *tree)
{
  RATES_Set_Node_Times_Pre(tree->n_root,tree->n_root->v[0],tree);
  RATES_Set_Node_Times_Pre(tree->n_root,tree->n_root->v[1],tree);
}

/*********************************************************/

void RATES_Init_Triplets(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-1) tree->rates->triplet[i] = 0.0;
}
/*********************************************************/

void RATES_Set_Node_Times_Pre(node *a, node *d, arbre *tree)
{
  if(d->tax) return;
  else
    {
      node *v1, *v2; /* the two sons of d */
      phydbl t_sup, t_inf;
      int i;

      v1 = v2 = NULL;
      For(i,3) if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	{
	  if(!v1) v1 = d->v[i]; 
	  else    v2 = d->v[i];
	}
	  
      t_inf = MIN(tree->rates->t[v1->num],tree->rates->t[v2->num]);
      t_sup = tree->rates->t[a->num];

      if(t_sup > t_inf)
	{
	  printf("\n. t_sup = %f t_inf = %f",t_sup,t_inf);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
            
      if(tree->rates->t[d->num] > t_inf)      tree->rates->t[d->num] = t_inf;
      else if(tree->rates->t[d->num] < t_sup) tree->rates->t[d->num] = t_sup;

      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Set_Node_Times_Pre(d,d->v[i],tree);	      
	    }
	}
    }
}

/*********************************************************/

void RATES_Randomize_Node_Times_Pre(node *a, node *d, arbre *tree)
{
  if(d->tax) return;
  else
    {
      node *v1, *v2; /* the two sons of d */
      phydbl t_sup, t_inf, u;
      int i;

      v1 = v2 = NULL;
      For(i,3) if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	{
	  if(!v1) v1 = d->v[i]; 
	  else    v2 = d->v[i];
	}
	  
      t_inf = MIN(tree->rates->t[v1->num],tree->rates->t[v2->num]);
      t_sup = tree->rates->t[a->num];

      if(t_sup > t_inf)
	{
	  printf("\n. t_sup = %f t_inf = %f",t_sup,t_inf);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      
      u = rand();
      u /= RAND_MAX;

      tree->rates->t[d->num] = t_sup + u * fabs(t_sup-t_inf);

      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Randomize_Node_Times_Pre(d,d->v[i],tree);	      
	    }
	}
    }
}

/*********************************************************/

phydbl RATES_Dmu2_Given_Mu1_Bis(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl alpha, phydbl beta, phydbl lexp)
{
  phydbl p0,p1,p2,v1,v2;
  
  p0 = Dpois(0,lexp*(dt1+dt2));
  p1 = Dpois(1,lexp*(dt1+dt2));
  p2 = Dpois(2,lexp*(dt1+dt2));       
  
  v1 = 
    p0 * Dnorm_Moments(mu2,mu1,0.01) + 
    p1 * Dnorm_Moments(mu2,mu1,0.04) + 
    p2 * Dnorm_Moments(mu2,mu1,0.09) ;
  
  v2 = RATES_Dmu(mu2,dt2,alpha,beta,lexp,0);
  
  return(v1 + (1.-p0-p1-p2) * v2);
}

/*********************************************************/

void RATES_Replace_Br_Lengths_By_Rates(arbre *tree)
{
  phydbl rate;

  RATES_Replace_Br_Lengths_By_Rates_Pre(tree->n_root,tree->n_root->v[0],NULL,tree);
  RATES_Replace_Br_Lengths_By_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);
  
  rate = tree->n_root->l[0] / (tree->rates->clock_r * fabs(tree->rates->t[tree->n_root->num] - tree->rates->t[tree->n_root->v[0]->num]));
/*   rate = fabs(tree->rates->t[tree->n_root->num] - tree->rates->t[tree->n_root->v[0]->num]); */
  tree->n_root->l[0] = rate;

  rate = tree->n_root->l[1] / (tree->rates->clock_r * fabs(tree->rates->t[tree->n_root->num] - tree->rates->t[tree->n_root->v[1]->num]));
/*   rate = fabs(tree->rates->t[tree->n_root->num] - tree->rates->t[tree->n_root->v[1]->num]); */
  tree->n_root->l[1] = rate;

  tree->e_root->l = tree->n_root->l[0] + tree->n_root->l[1];

  tree->n_root_pos = tree->n_root->l[0] / tree->e_root->l;

}

/*********************************************************/

void RATES_Replace_Br_Lengths_By_Rates_Pre(node *a, node *d, edge *b, arbre *tree)
{
  if(b)
    {
      phydbl rate;      
      rate = b->l / (tree->rates->clock_r * fabs(tree->rates->t[a->num] - tree->rates->t[d->num]));
/*       rate = fabs(tree->rates->t[a->num] - tree->rates->t[d->num]); */
/*       printf("\n. dt = %f clock_r = %f r=%f", */
/* 	     fabs(tree->rates->t[a->num] - tree->rates->t[d->num]), */
/* 	     tree->rates->clock_r, */
/* 	     rate); */
      b->l = rate;
    }
 
  if(d->tax) return;
  else
    {
      int i;
      
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Replace_Br_Lengths_By_Rates_Pre(d,d->v[i],d->b[i],tree);
	    }
	}
    }
}

/*********************************************************/

phydbl RATES_Dmu1_And_Mu2_One_Jump_Two_Intervals(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl alpha, phydbl beta)
{
  phydbl dens;

/*   dens = */
/*     (dt1/(dt1+dt2)) * RATES_Dmu1_And_Mu2_One_Jump_One_Interval(mu1,mu2,alpha,beta) + */
/*     (dt2/(dt1+dt2)) * RATES_Dmu1_And_Mu2_One_Jump_One_Interval(mu2,mu1,alpha,beta) ; */

  dens =
    ((dt1/(dt1+dt2)) * RATES_Dmu1_Given_V_And_N(mu1,mu2,1,dt1,alpha,beta) * Dgamma(mu2,alpha,beta)) +
    ((dt2/(dt1+dt2)) * RATES_Dmu1_Given_V_And_N(mu2,mu1,1,dt2,alpha,beta) * Dgamma(mu1,alpha,beta));

/*   dens += */
/*     (1./(dt1+dt2)) * */
/*     Dgamma(mu1,alpha,beta) * */
/*     Dgamma(mu2,alpha,beta); /\* Jump occurred at the junction *\/ */

  return dens;
}

/*********************************************************/

phydbl RATES_Dmu1_And_Mu2_One_Jump_One_Interval(phydbl mu1, phydbl mu2, phydbl a, phydbl b)
{
  phydbl lbda,h,p_unif,x;
  int n_inter,i;

  n_inter = (int)pow(2,2);
  p_unif = 1./(phydbl)(n_inter-1);
  
  /* Lambda never reaches 1.0  
  */
  h = 0;
  lbda = 0.0;
  for(i=0;i<n_inter-1;i++) 
    {
      lbda += 1./(phydbl)n_inter;
/*       x = (mu1 - (1.-lbda)*mu2)/lbda; */
      x = (mu1 - lbda*mu2)/(1.-lbda);
      if(x > 0.0) h += Dgamma(x,a,b) / (1.-lbda);
/*       printf("\n. ** lbda=%f dens=%f h=%f x=%f",lbda,Dgamma(x,a,b)/(1.-lbda),h,x); */
    }

  h *= Dgamma(mu2,a,b);

  return(h*p_unif);
}

/*********************************************************/

phydbl RATES_Dmu1_Given_V_And_N(phydbl mu1, phydbl v, int n, phydbl dt1, phydbl a, phydbl b)
{
  phydbl lbda,dens,h,u;
  phydbl mean,var;
  int n_points,i;
  phydbl ndb;
  phydbl end, beg;

  n_points = 100;

  end = MIN(mu1/v-0.01,0.99);
/*   end = 0.99; */
  beg = 0.01;
  

  mean = a*b;
  var = a*b*b*2./(n+1.);
  ndb = (phydbl)n/dt1;

  h = (end - beg) / (phydbl)n_points;

  lbda = beg;
  dens = 0.0;
  For(i,n_points-1) 
    {
      lbda += h;
      u = (mu1 - lbda*v)/(1.-lbda);
      if(u > 0.0) dens += Dgamma_Moments(u,mean,var) / (1.-lbda) * ndb * pow(1.-lbda,n-1);
    }
  dens *= 2.;

  lbda = beg;
  u = (mu1 - lbda*v)/(1.-lbda);
  if(u>0.0) dens += Dgamma_Moments(u,mean,var) / (1.-lbda) * ndb * pow(1.-lbda,n-1);

  lbda = end;
  u = (mu1 - lbda*v)/(1.-lbda);
  if(u>0.0) dens += Dgamma_Moments(u,mean,var) / (1.-lbda) * ndb * pow(1.-lbda,n-1);

  dens *= (h/2.);
  dens *= dt1;

  return(dens);
}

/*********************************************************/

phydbl RATES_Dmu_Given_V_And_MinN(phydbl mu, phydbl dt, phydbl v, int minn, phydbl a, phydbl b, phydbl lexp)
{
  int n,up,down;
  phydbl dens,cumpoiss,poiss,lexpdt;
  
  lexpdt   = lexp*dt;
  cumpoiss = 0.0;
  dens     = 0.0;
  poiss    = 0.0;

  RATES_Bracket_N_Jumps(&up,&down,lexpdt);

  For(n,MAX(down,minn)-1) cumpoiss += Dpois(n,lexpdt);

  for(n=MAX(down,minn);n<up+1;n++) /* Integrate over the possible number of jumps */
    {
      poiss     = Dpois(n,lexpdt); 
      dens     += poiss * RATES_Dmu1_Given_V_And_N(mu, v, n, dt, a, b);
      cumpoiss += poiss;
      if(cumpoiss > 1.-1.E-06) break;
    }

  return(dens);
}

/*********************************************************/
/* Joint density of mu1 and a minimum number of jumps occuring in the interval dt1+dt2 given mu1. 
   1 jump occurs at the junction of the two intervals, which makes mu1 and mu2 independant */
phydbl RATES_Dmu2_And_Mu1_Given_Min_N(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, int n_min, phydbl a, phydbl b, phydbl lexp)
{
  phydbl density, lexpdt,cumpoiss,poiss;
  int i;
  int up,down;

  density = 0.0;
  lexpdt = lexp * (dt1+dt2);
  cumpoiss = 0.0;

  RATES_Bracket_N_Jumps(&up,&down,lexpdt);

  For(i,MAX(up,n_min)-1)
    {
      poiss = Dpois(i,lexpdt);
      cumpoiss = cumpoiss + poiss;
    }

  for(i=MAX(up,n_min);i<up;i++)
    {
/*       poiss = Dpois(i-1,lexpdt); /\* Complies with the no correlation model *\/ */
      poiss = Dpois(i,lexpdt);
      cumpoiss = cumpoiss + poiss;

/*       density = density + poiss * RATES_Dmu2_And_Mu1_Given_N(mu1,mu2,dt1,dt2,i-1,a,b,lexp); */
      density = density + poiss * RATES_Dmu2_And_Mu1_Given_N_Full(mu1,mu2,dt1,dt2,i,a,b,lexp);
      if(cumpoiss > 1.-1.E-6) break;
    }

  if(density < 0.0)
    {
      printf("\n. density=%f cmpoiss = %f i=%d n_min=%d mu1=%f mu2=%f dt1=%f dt2=%f",
	     density,cumpoiss,i,n_min,mu1,mu2,dt1,dt2);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  return(density);
}

/*********************************************************/
/* Joint density of mu1 and mu2 given the number of jumps (n) in the interval dt1+dt2, which are considered
   as independant. Hence, for n jumps occuring in dt1+dt2, the number of jumps occuring in dt1 and dt2 are
   n and 0, or n-1 and 1, or n-2 and 2 ... or 0 and n. This function sums over all these scenarios, with
   weights corresponding to the probability of each partitition.
 */

phydbl RATES_Dmu2_And_Mu1_Given_N(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, int n, phydbl a, phydbl b, phydbl lexp)
  {
    phydbl density,lexpdt1,lexpdt2,texpn,cumpoiss,poiss,abb,ab,lognf,logdt1,logdt2,nlogdt;
    int i;

    density  = 0.0;
    lexpdt1  = lexp*dt1;
    lexpdt2  = lexp*dt2;
    texpn    = pow(dt1+dt2,n);
    abb      = a*b*b;
    ab       = a*b;
    cumpoiss = 0.0;
    poiss    = 0.0;
    lognf    = LnFact(n);
    logdt1   = log(dt1);
    logdt2   = log(dt2);
    nlogdt   = n*log(dt1+dt2);

    For(i,n+1)
      {
        poiss = lognf - LnFact(i) - LnFact(n-i) + i*logdt1 + (n-i)*logdt2 - nlogdt;
	poiss = exp(poiss);
	cumpoiss = cumpoiss + poiss;
	density = density + poiss * Dgamma_Moments(mu2,ab,2./((n-i)+2.)*abb) * Dgamma_Moments(mu1,ab,2./(i+2.)*abb);
        if(cumpoiss > 1.-1.E-6) break;
      }

    return(density);
  }

/*********************************************************/

phydbl RATES_Dmu2_And_Mu1_Given_N_Full(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, int n, phydbl a, phydbl b, phydbl lexp)
{
  phydbl dens,dens_i,dens_i_v,poiss,lognf,logdt1,logdt2,nlogdt;
  phydbl vval[10];
  int i,j,n_vval;

  /* Quantiles for alpha=2 (beta=1/2) */
  vval[0] = 0.2659058;
  vval[1] = 0.4121942;
  vval[2] = 0.5486746;
  vval[3] = 0.6882107;
  vval[4] = 0.8391735;
  vval[5] = 1.011157;
  vval[6] = 1.219608;
  vval[7] = 1.497154;
  vval[8] = 1.94486;
  vval[9] = 2.;

/*   vval[0] = 0.2659058; */
/*   vval[1] = 0.4121942; */
/*   vval[2] = 0.5486746; */
/*   vval[1] = 0.6882107; */
/*   vval[4] = 0.8391735; */
/*   vval[2] = 1.011157; */
/*   vval[6] = 1.219608; */
/*   vval[3] = 1.497154; */
/*   vval[4] = 1.94486; */
/*   vval[9] = 10.; */

  dens     = 0.0;
  poiss    = 0.0;
  lognf    = LnFact(n);
  logdt1   = log(dt1);
  logdt2   = log(dt2);
  nlogdt   = n*log(dt1+dt2);
  n_vval    = 10;


  n_vval = 1.0; vval[0] = 1.0;

  for(i=1;i<n;i++) /* At least one jump in both intervals */
    {
      poiss = lognf - LnFact(i) - LnFact(n-i) + i*logdt1 + (n-i)*logdt2 - nlogdt;
      poiss = exp(poiss);
      
      dens_i = 0.0;
      for(j=1;j<n_vval-1;j++) /* integrate over v */
	{
	  dens_i_v =
	    RATES_Dmu1_Given_V_And_N(mu1,vval[j],i,dt1,a,b)*    /* integrate over lbda1 */
	    RATES_Dmu1_Given_V_And_N(mu2,vval[j],n-i,dt2,a,b);  /* integrate over lbda2 */

	  dens_i +=  dens_i_v * Dgamma(vval[j],a,b) * (vval[j+1]-vval[j-1]);
	}

      dens_i +=
	RATES_Dmu1_Given_V_And_N(mu1,vval[0],i,dt1,a,b)*
	RATES_Dmu1_Given_V_And_N(mu2,vval[0],n-i,dt2,a,b) *
	Dgamma(vval[0],a,b) * (vval[1]-vval[0]);

      dens_i +=
	RATES_Dmu1_Given_V_And_N(mu1,vval[n_vval-1],i,dt1,a,b)*
	RATES_Dmu1_Given_V_And_N(mu2,vval[n_vval-1],n-i,dt2,a,b) *
	Dgamma(vval[n_vval-1],a,b) * (vval[n_vval-1]-vval[n_vval-2]);
      
      dens_i /= 2.;

      dens += dens_i * poiss;
    }

/*   /\* All jumps in the first interval *\/ */
/*   poiss =  n*logdt1 - nlogdt; */
/*   poiss = exp(poiss); */
/*   dens += poiss * RATES_Dmu1_Given_V_And_N(mu1,mu2,n,dt1,a,b) * Dgamma(mu2,a,b); */

/*   /\* All jumps in the second interval *\/ */
/*   poiss =  n*logdt2 - nlogdt; */
/*   poiss = exp(poiss); */
/*   dens += poiss * RATES_Dmu1_Given_V_And_N(mu2,mu1,n,dt2,a,b) * Dgamma(mu1,a,b); */

  return(dens);
}

/*********************************************************/

/* Formula (9) in Rannala and Yang (1996) */
phydbl RATES_Yule(arbre *tree)
{
  phydbl sumti,density,lambda;
  int n,i;

  sumti = 0.0;
  for(i=tree->n_otu;i<2*tree->n_otu-1;i++) sumti += tree->rates->t[i];
  sumti -= tree->rates->t[i];

  lambda = tree->rates->birth_rate;
  n = tree->n_otu;
  
  density = 
    (n-1.)*log(2.) + 
    (n-2.)*log(lambda) - 
    lambda*sumti - 
    factln(n) - 
    log(n-1.) - 
    (n-2.)*log(1.-exp(-lambda));
  
  density = exp(density);

  return density;
}

/*********************************************************/

/* Set the clock rate such that the relative rates are centered on 1.0 */
void RATES_Adjust_Clock_Rate(arbre *tree)
{
  phydbl ratio;
  int i;

  ratio = 0.0;
  For(i,2*tree->n_otu-3)
    if(tree->t_edges[i] != tree->e_root)
      ratio +=
	tree->t_edges[i]->l / 
	(fabs(tree->rates->t[tree->t_edges[i]->left->num] -
	      tree->rates->t[tree->t_edges[i]->rght->num]));

  ratio += tree->n_root->l[0] / (fabs(tree->rates->t[tree->n_root->v[0]->num] -
				      tree->rates->t[tree->n_root->num]));

  ratio += tree->n_root->l[1] / (fabs(tree->rates->t[tree->n_root->v[1]->num] -
				      tree->rates->t[tree->n_root->num]));

  tree->rates->clock_r = ratio/(phydbl)(2.*tree->n_otu-2);
}


#endif
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
