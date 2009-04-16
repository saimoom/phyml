/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/


/* Routines for molecular clock trees and molecular dating */
#if defined (MC) ||  defined (RWRAPPER)

#include "utilities.h"
#include "spr.h"
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
#include "mcmc.h"
#include "numeric.h"

#ifdef RWRAPPER
#include <R.h>
#endif



/*********************************************************/

phydbl RATES_Lk_Rates(arbre *tree)
{
  
  RATES_Set_Node_Times(tree);
  RATES_Init_Triplets(tree);

  tree->rates->c_lnL = .0;
  RATES_Lk_Rates_Pre(tree->n_root,tree->n_root->v[0],NULL,tree);
  RATES_Lk_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);

  if(isnan(tree->rates->c_lnL))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  return tree->rates->c_lnL;
}

/*********************************************************/

void RATES_Lk_Rates_Pre(node *a, node *d, edge *b, arbre *tree)
{
  int i,n1,n2;
  phydbl log_dens,mu1,mu2,dt1,dt2;
  
  log_dens = -1.;

  if(d->anc != a)
    {
      printf("\n. d=%d d->anc=%d a=%d root=%d",d->num,d->anc->num,a->num,tree->n_root->num);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(a != tree->n_root)
    {
      dt1 = fabs(tree->rates->nd_t[a->num] - tree->rates->nd_t[a->anc->num]);
      mu1 = tree->rates->nd_r[a->num];
      n1  = tree->rates->n_jps[a->num];

      dt2 = fabs(tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num]);
      mu2 = tree->rates->nd_r[d->num];
      n2  = tree->rates->n_jps[d->num];

      log_dens = RATES_Lk_Rates_Core(mu1,mu2,n1,n2,dt1,dt2,tree);
    }
  else
    {
      dt2 = fabs(tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num]);
      mu2 = tree->rates->nd_r[d->num];
      n2  = tree->rates->n_jps[d->num];

      switch(tree->rates->model)
	{
	case COMPOUND_COR: case COMPOUND_NOCOR:
	  {       
	    log_dens = RATES_Dmu(mu2,n2,dt2,tree->rates->alpha,1./tree->rates->alpha,tree->rates->lexp,0,1);
	    log_dens = log(log_dens);
	    break;
	  }
	case EXPONENTIAL:
	  {
	    log_dens = Dexp(mu2,tree->rates->lexp);
	    log_dens = log(log_dens);
	    break;
	  }
	case GAMMA:
	  {
	    log_dens = Dgamma(mu2,tree->rates->alpha,1./tree->rates->alpha);
	    log_dens = log(log_dens);
	    break;
	  }
	case THORNE:
	  {
	    log_dens = Log_Dnorm(mu2,1.0,sqrt(tree->rates->nu*dt2));
	    if(log_dens < -10000.) printf("\n. run=%d",tree->mcmc->run);	      
	    break;
	  }

	default:
	  {	    
	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Exit("\n. Model not implemented yet.\n");
	    break;
	  }
	}
    }

  tree->rates->c_lnL += log_dens;

  if(isnan(tree->rates->c_lnL))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      MCMC_Print_Param(stderr,tree);
      Exit("\n");
    }

  tree->rates->triplet[a->num] += log_dens;

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Lk_Rates_Pre(d,d->v[i],d->b[i],tree);
	    }
	}
    }
}

/*********************************************************/

phydbl RATES_Lk_Change_One_Rate(node *d, phydbl new_rate, arbre *tree)
{
  tree->rates->nd_r[d->num] = new_rate;
  RATES_Update_Triplet(d,tree);
  RATES_Update_Triplet(d->anc,tree);
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
      
      tree->rates->nd_t[n->num] = new_t;

      RATES_Update_Triplet(n,tree);
      
      For(i,3)
	{
	  if(n->b[i] != tree->e_root) RATES_Update_Triplet(n->v[i],tree);
	  else RATES_Update_Triplet(tree->n_root,tree);
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
  phydbl mu0,mu1,mu2;
  int n0,n1,n2;
  int i;
  node *v1,*v2;

  if(n->tax) return;

  curr_triplet = tree->rates->triplet[n->num];

  dt0 = dt1 = dt2 = -100.0;

  if(n == tree->n_root)
    {
      phydbl log_dens;
      
      log_dens = 0.0;

      dt0 = tree->rates->nd_t[tree->n_root->v[0]->num] - tree->rates->nd_t[tree->n_root->num];
      dt1 = tree->rates->nd_t[tree->n_root->v[1]->num] - tree->rates->nd_t[tree->n_root->num];
      
      mu0 = tree->rates->nd_r[tree->n_root->v[0]->num];
      mu1 = tree->rates->nd_r[tree->n_root->v[1]->num];
      
      n0  = tree->rates->n_jps[tree->n_root->v[0]->num];
      n1  = tree->rates->n_jps[tree->n_root->v[1]->num];


      switch(tree->rates->model)
	{
	case COMPOUND_COR : case COMPOUND_NOCOR : 
	  {
	    log_dens  = RATES_Dmu(mu0,n0,dt0,tree->rates->alpha,1./tree->rates->alpha,tree->rates->lexp,0,1);
	    log_dens *= RATES_Dmu(mu1,n1,dt1,tree->rates->alpha,1./tree->rates->alpha,tree->rates->lexp,0,1);
	    log_dens  = log(log_dens);
	    break;
	  }
	case EXPONENTIAL : 
	  {
	    log_dens = Dexp(mu0,tree->rates->lexp) * Dexp(mu1,tree->rates->lexp);
	    log_dens = log(log_dens);
	    break;
	  }
	case GAMMA :
	  {
	    log_dens = Dgamma(mu0,tree->rates->alpha,1./tree->rates->alpha) * Dgamma(mu1,tree->rates->alpha,1./tree->rates->alpha);
	    log_dens = log(log_dens);
	    break;
	  }
	case THORNE :
	  {
	    log_dens = Log_Dnorm(mu0,1.0,sqrt(sqrt(tree->rates->nu*dt0))) + Log_Dnorm(mu1,1.0,sqrt(tree->rates->nu*dt1)); 
	    if(log_dens < -10000.) printf("\n. run=%d",tree->mcmc->run);	      
	    break;
	  }
	default :
	  {
	    Exit("\n. Model not implemented yet.\n");
	    break;
	  }
	}
      new_triplet = log_dens;

      if(isnan(log_dens) || isinf(fabs(log_dens)))
	{
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  MCMC_Print_Param(stderr,tree);
	  Exit("\n");
	}
    }
  else
    {
      mu0 = mu1 = mu2 = -1.;
      n0 = n1 = n2 = -1;

      mu0 = tree->rates->nd_r[n->num];
      dt0 = fabs(tree->rates->nd_t[n->num] - tree->rates->nd_t[n->anc->num]);
      n0  = tree->rates->n_jps[n->num];

      v1 = v2 = NULL;
      For(i,3)
	{
	  if((n->v[i] != n->anc) && (n->b[i] != tree->e_root))
	    {
	      if(!v1)
		{
		  v1  = n->v[i]; 
		  mu1 = tree->rates->nd_r[v1->num];
		  dt1 = fabs(tree->rates->nd_t[v1->num] - tree->rates->nd_t[n->num]);
		  n1  = tree->rates->n_jps[v1->num];
		}
	      else
		{
		  v2  = n->v[i]; 
		  mu2 = tree->rates->nd_r[v2->num];
		  dt2 = fabs(tree->rates->nd_t[v2->num] - tree->rates->nd_t[n->num]);
		  n2  = tree->rates->n_jps[v2->num];
		}
	    }
	}
 
      mu1_mu0 = RATES_Lk_Rates_Core(mu0,mu1,n0,n1,dt0,dt1,tree);
      mu2_mu0 = RATES_Lk_Rates_Core(mu0,mu2,n0,n2,dt0,dt2,tree);
      
      new_triplet = mu1_mu0 + mu2_mu0;
    }

  tree->rates->c_lnL = tree->rates->c_lnL + new_triplet - curr_triplet;
  tree->rates->triplet[n->num] = new_triplet;
}

/*********************************************************/
/* Returns log(f(mu2;mu1)) */
phydbl RATES_Lk_Rates_Core(phydbl mu1, phydbl mu2, int n1, int n2, phydbl dt1, phydbl dt2, arbre *tree)
{
  phydbl log_dens;
  phydbl alpha, beta, lexp;

  lexp = tree->rates->lexp;
  alpha = tree->rates->alpha;
  beta = 1./alpha;
  log_dens = UNLIKELY;

  if(mu1 < tree->rates->min_rate) mu1 = tree->rates->min_rate;
  if(mu1 > tree->rates->max_rate) mu1 = tree->rates->max_rate;
  
  if(mu2 < tree->rates->min_rate) mu2 = tree->rates->min_rate;
  if(mu2 > tree->rates->max_rate) mu2 = tree->rates->max_rate;
 
/*   if(dt1 < tree->rates->min_dt) dt1 = tree->rates->min_dt; */
/*   if(dt2 < tree->rates->min_dt) dt2 = tree->rates->min_dt; */
  
  switch(tree->rates->model)
    {
    case COMPOUND_COR:
      {       
	log_dens = RATES_Compound_Core(mu1,mu2,n1,n2,dt1,dt2,alpha,beta,lexp,tree->rates->step_rate,tree->rates->approx);
	log_dens = log(log_dens);
	break;
      }
      
    case COMPOUND_NOCOR :
      {
	log_dens = RATES_Dmu(mu2,n2,dt2,alpha,beta,lexp,0,1);
	log_dens = log(log_dens);
	break;
      }
      
    case EXPONENTIAL :
      {
	log_dens = Dexp(mu2,tree->rates->lexp);
	log_dens = log(log_dens);
	break;
      }
      
    case GAMMA :
      {
	log_dens = Dgamma(mu2,tree->rates->alpha,1./tree->rates->alpha);
	log_dens = log(log_dens);
	break;
      }
      
    case THORNE :
      {
	log_dens = Log_Dnorm(mu2,mu1,sqrt(tree->rates->nu*dt2));
	if(log_dens < -10000.) printf("\n. run=%d",tree->mcmc->run);	      
	break;
      }

    default : 
      {
	PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	Warn_And_Exit("");
      }
    }

  if(isnan(log_dens) || isinf(fabs(log_dens)))
    {
      PhyML_Printf("\n. Run=%4d mu2=%f mu1=%f dt2=%f dt1=%f nu=%f log_dens=%G sd=%f\n",
		   tree->mcmc->run,
		   mu2,mu1,dt2,dt1,tree->rates->nu,log_dens,
		   sqrt(tree->rates->nu*dt2));
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  return log_dens;
}

/*********************************************************/

phydbl RATES_Compound_Core(phydbl mu1, phydbl mu2, int n1, int n2, phydbl dt1, phydbl dt2, phydbl alpha, phydbl beta, phydbl lexp, phydbl eps, int approx)
{
  if((n1 > -1) && (n2 > -1))
    {
      return RATES_Compound_Core_Joint(mu1,mu2,n1,n2,dt1,dt2,alpha,beta,lexp,eps,approx);
    }
  else
    {
      return RATES_Compound_Core_Marginal(mu1,mu2,dt1,dt2,alpha,beta,lexp,eps,approx);
    }
}

/*********************************************************/

phydbl RATES_Compound_Core_Marginal(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl alpha, phydbl beta, phydbl lexp, phydbl eps, int approx)
{
  phydbl p0,p1,v0,v1,v2;
  phydbl dmu1;
  int    equ;
  phydbl dens;
  
  v0 = v1 = v2 = 0.0;
  
  /* Probability of 0 and 1 jumps */
  p0   = Dpois(0,lexp*(dt2+dt1));       
  p1   = Dpois(1,lexp*(dt2+dt1));
  
  dmu1 = RATES_Dmu(mu1,-1,dt1,alpha,beta,lexp,0,0);
  
  /* Are the two rates equal ? */
  equ = 0;
  if(fabs(mu1-mu2) < eps) equ = 1;
  
  /* No jump */
  if(equ)
    {
      v0 = 1.0*Dgamma(mu1,alpha,beta)/dmu1;
      /*       Rprintf("\n. mu1=%f mu2=%f",mu1,mu2); */
    }
  else
    {
      v0 = 1.E-100;
    }

  /* One jump */
  v1 = RATES_Dmu1_And_Mu2_One_Jump_Two_Intervals(mu1,mu2,dt1,dt2,alpha,beta);
  v1 /= dmu1;
    
  /* Two jumps and more (approximation) */
  if(approx == 1)
    {
      v2 =
	RATES_Dmu(mu1,-1,dt1,alpha,beta,lexp,0,0)*RATES_Dmu(mu2,-1,dt2,alpha,beta,lexp,0,0) -
	Dpois(0,lexp*dt1) * Dpois(0,lexp*dt2) *
	Dgamma(mu1,alpha,beta) * Dgamma(mu2,alpha,beta);
    }
  else
    {
      v2 = 
	RATES_Dmu_One(mu1,dt1,alpha,beta,lexp) * 
	RATES_Dmu_One(mu2,dt2,alpha,beta,lexp);

      v2 += Dpois(0,lexp*dt1)*Dgamma(mu1,alpha,beta)*RATES_Dmu(mu2,-1,dt2,alpha,beta,lexp,1,0);
      v2 += Dpois(0,lexp*dt2)*Dgamma(mu2,alpha,beta)*RATES_Dmu(mu1,-1,dt1,alpha,beta,lexp,1,0);

    }
/*   printf("\n. %f %f %f %f %f ",mu1,mu2,dt1,dt2,v2); */
  v2 /= dmu1;

  dens = p0*v0 + p1*v1 + v2;
/*   dens = p1*v1 + v2; */
  /*       dens = p1*v1 + v2; */
  /*   dens = v0; */
/*   dens *= dmu1; */
  
  if(dens < MDBL_MIN)
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

phydbl RATES_Compound_Core_Joint(phydbl mu1, phydbl mu2, int n1, int n2, phydbl dt1, phydbl dt2, 
				 phydbl alpha, phydbl beta, phydbl lexp, phydbl eps, int approx)
{
  phydbl density;
  phydbl dmu1;
  

  if(n1 < 0 || n2 < 0)
    {
      printf("\n. n1=%d n2=%d",n1,n2);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  dmu1 = RATES_Dmu(mu1,n1,dt1,alpha,beta,lexp,0,0);
  
  if((n1 < 0) || (n2 < 0))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if((n1 == 0) && (n2 == 0))
    {
      if(fabs(mu1-mu2) < eps) { density = Dgamma(mu1,alpha,beta); }
      else                    { density = 1.E-70; }
    }
  else if((n1 == 0) && (n2 == 1))
    {
      density = 
	Dgamma(mu1,alpha,beta) *
	RATES_Dmu1_Given_V_And_N(mu2,mu1,1,dt2,alpha,beta);
    }
  else if((n1 == 1) && (n2 == 0))
    {
      density = 
	Dgamma(mu2,alpha,beta) *
	RATES_Dmu1_Given_V_And_N(mu1,mu2,1,dt1,alpha,beta);
    }
  else /* independent */
    {
      density = 
	RATES_Dmu(mu1,n1,dt1,alpha,beta,lexp,0,0) * 
	RATES_Dmu(mu2,n2,dt2,alpha,beta,lexp,0,0);
    }

  density /= dmu1;

  density *= Dpois(n2,dt2*lexp);
  
  if(density < 1.E-70) density = 1.E-70;

/*   printf("\n. density = %15G mu1=%3.4f mu2=%3.4f dt1=%3.4f dt2=%3.4f n1=%2d n2=%2d",density,mu1,mu2,dt1,dt2,n1,n2); */
  return density;
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

  
  if((d == tree->n_root->v[0] && d->tax) || (d == tree->n_root->v[1] && d->tax))
    printf("\n. a=%3d ++d=%3d rate=%12f t_left=%12f t_rght=%12f ml=%12f l=%12f %12f",
	   a->num,d->num,
	   tree->rates->nd_r[d->num],
	   tree->rates->nd_t[a->num],tree->rates->nd_t[d->num],
	   tree->rates->ml_l[d->num],
	   tree->rates->cur_l[d->num],
	   (tree->rates->nd_t[d->num]-tree->rates->nd_t[a->num])*tree->rates->clock_r*tree->rates->nd_r[d->num]);
  
  else if((d == tree->n_root->v[0]) || (d == tree->n_root->v[1]))
    printf("\n. a=%3d __d=%3d rate=%12f t_left=%12f t_rght=%12f ml=%12f l=%12f %12f",
	   a->num,d->num,
	   tree->rates->nd_r[d->num],
	   tree->rates->nd_t[a->num],tree->rates->nd_t[d->num],
	   tree->rates->ml_l[d->num],
	   tree->rates->cur_l[d->num],
	   (tree->rates->nd_t[d->num]-tree->rates->nd_t[a->num])*tree->rates->clock_r*tree->rates->nd_r[d->num]);
  else 
    printf("\n. a=%3d   d=%3d rate=%12f t_left=%12f t_rght=%12f ml=%12f l=%12f %12f",
	   a->num,d->num,
	   tree->rates->nd_r[d->num],
	   tree->rates->nd_t[a->num],tree->rates->nd_t[d->num],
	   tree->rates->ml_l[d->num],
	   tree->rates->cur_l[d->num],
	   (tree->rates->nd_t[d->num]-tree->rates->nd_t[a->num])*tree->rates->clock_r*tree->rates->nd_r[d->num]);
  
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
  int i;
  
  sum = 0.0;
  For(i,2*tree->n_otu-2) sum += tree->rates->nd_r[i];
  return(sum/(phydbl)(2*tree->n_otu-2));
}

/*********************************************************/

void RATES_Check_Node_Times(arbre *tree)
{
  RATES_Check_Node_Times_Pre(tree->n_root,tree->n_root->v[0],tree);
  RATES_Check_Node_Times_Pre(tree->n_root,tree->n_root->v[1],tree);
}

/*********************************************************/

void RATES_Check_Node_Times_Pre(node *a, node *d, arbre *tree)
{
  if(tree->rates->nd_t[d->num] < tree->rates->nd_t[a->num])
    {
      PhyML_Printf("\n. a->t=%f d->t=%f",tree->rates->nd_t[d->num],tree->rates->nd_t[a->num]);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  if(d->tax) return;
  else
    {
      int i;

      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))	  
	  RATES_Check_Node_Times_Pre(d,d->v[i],tree);
    }
}
/*********************************************************/

trate *RATES_Make_Rate_Struct(int n_otu)
{
  trate *rates;

  rates               = (trate  *)mCalloc(1,sizeof(trate));
  rates->br_r         = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
  rates->old_r        = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->nd_r         = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->true_r       = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->old_t        = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->nd_t         = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->true_t       = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->dens         = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
  rates->triplet      = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->n_jps        = (int    *)mCalloc(2*n_otu-2,sizeof(int));
  rates->t_jps        = (int    *)mCalloc(2*n_otu-2,sizeof(int));
  rates->cov          = (phydbl *)mCalloc((2*n_otu-2)*(2*n_otu-2),sizeof(phydbl));
  rates->invcov       = (phydbl *)mCalloc((2*n_otu-2)*(2*n_otu-2),sizeof(phydbl));
  rates->ml_l         = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
  rates->cur_l        = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
  rates->u_ml_l       = (phydbl *)mCalloc(2*n_otu-3,sizeof(phydbl));
  rates->u_cur_l      = (phydbl *)mCalloc(2*n_otu-3,sizeof(phydbl));
  rates->cov_r        = (phydbl *)mCalloc((2*n_otu-2)*(2*n_otu-2),sizeof(phydbl));
  rates->mean_r       = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
  rates->lca          = (node **)mCalloc((2*n_otu-1)*(2*n_otu-1),sizeof(node *));
      
  return rates;
}

/*********************************************************/

void Free_Rates(trate *rates)
{
  Free(rates->br_r);
  Free(rates->old_r);     
  Free(rates->nd_r);    
  Free(rates->true_r);     
  Free(rates->old_t);   
  Free(rates->nd_t);   
  Free(rates->true_t);     
  Free(rates->dens);   
  Free(rates->triplet);    
  Free(rates->n_jps);  
  Free(rates->t_jps);    
  Free(rates->cov);   
  Free(rates->invcov);   
  Free(rates->ml_l);
  Free(rates->cur_l);
  Free(rates->u_ml_l);
  Free(rates->u_cur_l);
  Free(rates->cov_r);
  Free(rates->lca);
}

/*********************************************************/

void RATES_Init_Rate_Struct(trate *rates, int n_otu)
{
  int i;

  rates->met_within_gibbs = NO;
  rates->model         = THORNE;
  rates->clock_r       = 1.E-3;
  rates->c_lnL         = -INFINITY;
  rates->c_lnL_jps     = -INFINITY;
  rates->adjust_rates  = 0;
  rates->use_rates     = 1;
  rates->lexp          = 1.E-3;
  rates->alpha         = 2.;
  rates->birth_rate    = 0.001;
  rates->max_rate      = 1.E+2;
  rates->min_rate      = 1.E-2;
  rates->min_dt        = 1.E-8;
  rates->step_rate     = 1.E-4;
  rates->nu            = 1.E-3;
  rates->approx        = 1;
  rates->bl_from_rt    = 0;
  rates->lk_approx     = NORMAL;

  For(i,2*n_otu-2) 
    {
      rates->old_r[i]  = 1.0;
      rates->nd_r[i]   = 1.0;
      rates->br_r[i]   = 1.0;
      rates->n_jps[i]  =  -1;
      rates->t_jps[i]  =  -1;
      rates->mean_r[i] = 1.0;
    }

  For(i,2*n_otu-1) 
    {
      rates->nd_t[i]   = 0.0;
      rates->true_t[i] = 0.0;
    }
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
phydbl RATES_Dmu(phydbl mu, int n_jumps, phydbl dt, phydbl a, phydbl b, phydbl lexp, int min_n, int jps_dens)
{
  if(n_jumps < 0) /* Marginal, i.e., the number of jumps is not fixed */
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
      
      RATES_Bracket_N_Jumps(&up,&down,lexpdt);
      For(n,MAX(down,min_n)-1) cumpoissprob += Dpois(n,lexpdt);
      
      for(n=MAX(down,min_n);n<up+1;n++)
	{
	  poissprob    = Dpois(n,lexpdt); /* probability of having n jumps */      
	  var          = (2./(n+2.))*ab2; /* var(mu|n) = var(mu|n=0) * 2 / (n+2) */
	  gammadens    = Dgamma_Moments(mu,mean,var);
	  dens         += poissprob * gammadens;
	  cumpoissprob += poissprob;
	  if(cumpoissprob > 1.-1.E-04) break;
	}
      
      if(dens < 1.E-70) dens = 1.E-70;

      return(dens);      
    }
  else /* Joint, i.e., return P(mu | dt, n_jumps) */
    {
      phydbl mean, var, density;


      mean = 1.0;
      var = (2./(n_jumps+2.))*a*b*b;

      if(jps_dens)
	density = Dgamma_Moments(mu,mean,var) * Dpois(n_jumps,dt*lexp);
      else
	density = Dgamma_Moments(mu,mean,var);
      
      if(density < 1.E-70) density = 1.E-70;

      return density;
    }
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

  if(mu < 1.E-10)
    {
      PhyML_Printf("\n. mu=%G",mu);
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

/* Given the times of nodes a (ta) and d (td), the shape of the gamma distribution of instantaneous
   rates, the parameter of the exponential distribution of waiting times between rate jumps and the 
   instantaneous rate at node a, this function works out an expected number of (amino-acids or 
   nucleotide) substitutions per site.
*/
void RATES_Expect_Number_Subst(phydbl t_beg, phydbl t_end, phydbl r_beg,  int *n_jumps, phydbl *mean_r, phydbl *r_end, trate *rates, arbre *tree)
{
  phydbl curr_r, curr_t, next_t;

  switch(rates->model)
    {
    case COMPOUND_COR:case COMPOUND_NOCOR:
      {
	/* Compound Poisson */
	if(rates->model == COMPOUND_COR)
	  {
	    curr_r  = r_beg;
	    *mean_r = r_beg;
	  }
	else
	  {
	    curr_r  = Rgamma(rates->alpha,1./rates->alpha);;
	    *mean_r = curr_r;
	  }

	curr_t = t_beg + Rexp(rates->lexp); /* Exponentially distributed waiting times */
	next_t = curr_t;
	
	*n_jumps = 0;
	while(curr_t < t_end)
	  {
	    curr_r = Rgamma(rates->alpha,1./rates->alpha); /* Gamma distributed random instantaneous rate */
	    
	    (*n_jumps)++;
	    
	    next_t = curr_t + Rexp(rates->lexp);
	    
	    if(next_t < t_end)
	      {
		*mean_r = (1./(next_t - t_beg)) * (*mean_r * (curr_t - t_beg) + curr_r * (next_t - curr_t));
	      }
	    else
	      {
		*mean_r = (1./(t_end - t_beg)) * (*mean_r * (curr_t - t_beg) + curr_r * (t_end - curr_t));
	      }
	    curr_t = next_t;
	  }
	
	/*   printf("\n. [%3d %f %f]",*n_jumps,*mean_r,r_beg); */
	
	if(*mean_r < rates->min_rate) *mean_r = rates->min_rate;
	if(*mean_r > rates->max_rate) *mean_r = rates->max_rate;

	*r_end = curr_r;
	break;
      }


    case EXPONENTIAL:
      {
	*mean_r = Rexp(rates->lexp);

	if(*mean_r < rates->min_rate) *mean_r = rates->min_rate;
	if(*mean_r > rates->max_rate) *mean_r = rates->max_rate;

	*r_end  = *mean_r;
	break;
      }
    case GAMMA:
      {
	*mean_r = Rgamma(rates->alpha,1./rates->alpha);

	if(*mean_r < rates->min_rate) *mean_r = rates->min_rate;
	if(*mean_r > rates->max_rate) *mean_r = rates->max_rate;

	*r_end  = *mean_r;
	break;
      }
    case THORNE:
      {
	phydbl sd;

	sd = sqrt(rates->nu*fabs(t_beg-t_end));

	*mean_r = Rnorm_Trunc(r_beg,sd,rates->min_rate,rates->max_rate);

	*r_end  = *mean_r;
	
	break;
      }
    default:
      {
	PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	Exit("\n. Model not implemented yet.\n");
	break;
      }
    }
}
  
/*********************************************************/

void RATES_Get_Mean_Rates_Pre(node *a, node *d, edge *b, phydbl r_a, arbre *tree)
{
  phydbl a_t,d_t;
  phydbl mean_r;
  int n_jumps;
  phydbl r_d;


  a_t = tree->rates->nd_t[a->num];
  d_t = tree->rates->nd_t[d->num];
      
  RATES_Expect_Number_Subst(a_t,d_t,r_a,&n_jumps,&mean_r,&r_d,tree->rates,tree);

  tree->rates->nd_r[d->num]   = mean_r;
  tree->rates->true_r[d->num] = mean_r;
  tree->rates->t_jps[d->num]  = n_jumps;


  /* Move to the next branches */
  if(d->tax) return;
  else
    {
      int i;
      
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Get_Mean_Rates_Pre(d,d->v[i],d->b[i],r_d,tree);
	    }
	}
    }

}

/*********************************************************/

void RATES_Random_Branch_Lengths(arbre *tree)
{
  phydbl r0,r_d;
   
  PhyML_Printf("\n> WARNING: ROOT RATE SET TO ONE.");
  PhyML_Printf("\n. File %s -- line %d\n",__FILE__,__LINE__);
  r0 = 1.0;
  r_d = 0;

  tree->rates->nd_r[tree->n_root->num] = r0;

  RATES_Get_Mean_Rates_Pre(tree->n_root,tree->n_root->v[0],NULL,r0,tree);
  RATES_Get_Mean_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,r0,tree);

  RATES_Update_Cur_Bl(tree);
  RATES_Initialize_True_Rates(tree);

  tree->n_root_pos = 
    tree->rates->cur_l[tree->n_root->v[0]->num] /
    (tree->rates->cur_l[tree->n_root->v[0]->num] + tree->rates->cur_l[tree->n_root->v[1]->num]);
  printf("\n. Root rate here = %f",tree->rates->nd_r[tree->n_root->num]);
  printf("\n. Root position: %f\n",tree->n_root_pos);
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
	  
      t_inf = MIN(tree->rates->nd_t[v1->num],tree->rates->nd_t[v2->num]);
      t_sup = tree->rates->nd_t[a->num];

      if(t_sup > t_inf)
	{
	  printf("\n. t_sup = %f t_inf = %f",t_sup,t_inf);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
            
      if(tree->rates->nd_t[d->num] > t_inf)      tree->rates->nd_t[d->num] = t_inf;
      else if(tree->rates->nd_t[d->num] < t_sup) tree->rates->nd_t[d->num] = t_sup;

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


phydbl RATES_Dmu1_And_Mu2_One_Jump_Two_Intervals(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl alpha, phydbl beta)
{
  phydbl dens;

  if(mu2 < 1.E-10)
    {
      printf("\n. mu2=%G",mu2);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(mu1 < 1.E-10)
    {
      printf("\n. mu2=%G",mu1);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  dens =
    ((dt1/(dt1+dt2)) * RATES_Dmu1_Given_V_And_N(mu1,mu2,1,dt1,alpha,beta) * Dgamma(mu2,alpha,beta)) +
    ((dt2/(dt1+dt2)) * RATES_Dmu1_Given_V_And_N(mu2,mu1,1,dt2,alpha,beta) * Dgamma(mu1,alpha,beta));

  return dens;
}

/*********************************************************/

phydbl RATES_Dmu1_Given_V_And_N(phydbl mu1, phydbl v, int n, phydbl dt1, phydbl a, phydbl b)
{
  phydbl lbda,dens,h,u;
  phydbl mean,var;
  int n_points,i;
  phydbl ndb;
  phydbl end, beg;

  n_points = 20;

  end = MIN(mu1/v-0.01,0.99);
  beg = 0.01;
  
  dens = 0.0;
  
  if(end > beg)
    {
      mean = a*b;
      var = a*b*b*2./(n+1.);
      ndb = (phydbl)n/dt1;
      
      h = (end - beg) / (phydbl)n_points;
      
      lbda = beg;
      For(i,n_points-1) 
	{
	  lbda += h;
	  u = (mu1 - lbda*v)/(1.-lbda);
	  
	  if(u < 1.E-10)
	    {
	      printf("\n. u = %G",u);
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	  
	  dens += Dgamma_Moments(u,mean,var) / (1.-lbda) * ndb * pow(1.-lbda,n-1);
	}
      dens *= 2.;
      
      lbda = beg;
      u = (mu1 - lbda*v)/(1.-lbda);
      if(u < 1.E-10)
	{
	  printf("\n. mu1 = %f lambda = %f v = %f u = %G beg = %f end = %f",mu1,lbda,v,u,beg,end);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      
      dens += Dgamma_Moments(u,mean,var) / (1.-lbda) * ndb * pow(1.-lbda,n-1);
      
      lbda = end;
      u = (mu1 - lbda*v)/(1.-lbda);
      if(u < 1.E-10)
	{
	  printf("\n. u = %G",u);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      
      dens += Dgamma_Moments(u,mean,var) / (1.-lbda) * ndb * pow(1.-lbda,n-1);
      
      dens *= (h/2.);
      dens *= dt1;
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

      density = density + poiss * RATES_Dmu2_And_Mu1_Given_N(mu1,mu2,dt1,dt2,i-1,a,b,lexp);
/*       density = density + poiss * RATES_Dmu2_And_Mu1_Given_N_Full(mu1,mu2,dt1,dt2,i,a,b,lexp); */
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

	if(mu2 < 1.E-10)
	  {
	    printf("\n. mu2=%f",mu2);
	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	  }

	if(mu1 < 1.E-10)
	  {
	    printf("\n. mu1=%f",mu1);
	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	  }


	density = density + poiss * Dgamma_Moments(mu2,ab,2./((n-i)+2.)*abb) * Dgamma_Moments(mu1,ab,2./(i+2.)*abb);
        if(cumpoiss > 1.-1.E-6) break;
      }

    return(density);
  }

/*********************************************************/

/* Logarithm of Formula (9) in Rannala and Yang (1996) */
phydbl RATES_Yule(arbre *tree)
{
  phydbl sumti,density,lambda;
  int n,i;

  sumti = 0.0;
  for(i=tree->n_otu;i<2*tree->n_otu-1;i++) sumti += tree->rates->nd_t[i];
  sumti -= tree->rates->nd_t[i];

  lambda = tree->rates->birth_rate;
  n = tree->n_otu;
  
  density = 
    (n-1.)*log(2.) + 
    (n-2.)*log(lambda) - 
    lambda*sumti - 
    Factln(n) - 
    log(n-1.) - 
    (n-2.)*log(1.-exp(-lambda));
  
  /* Ad-hoc */
  if(tree->rates->nd_t[tree->n_root->num] > -100.) density = -INFINITY;

  return density;
}

/*********************************************************/
/* Set the clock rate such that the relative rates are centered on 1.0 */
/* Set the clock rate such that the relative rates are centered on the ML tree size */
phydbl RATES_Adjust_Clock_Rate(arbre *tree)
{
  int i;
  phydbl mean;

/*   mean = 0.0; */
/*   For(i,2*tree->n_otu-2) mean += tree->rates->nd_r[i]; */
/*   mean /= (phydbl)(2*tree->n_otu-2); */
/*   For(i,2*tree->n_otu-2) tree->rates->nd_r[i] /= mean; */
/*   tree->rates->clock_r *= mean; */

  mean = 0.0;
  For(i,2*tree->n_otu-2)
    mean += (tree->rates->nd_t[tree->noeud[i]->num] - tree->rates->nd_t[tree->noeud[i]->anc->num]) * tree->rates->nd_r[tree->noeud[i]->num];

  tree->rates->clock_r = tree->size/mean;

  return tree->rates->clock_r;
}

/*********************************************************/



/*********************************************************/

void RATES_Initialize_True_Rates(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-2) tree->rates->true_r[i] = tree->rates->nd_r[i];
}
/*********************************************************/

void RATES_Record_Rates(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-2) tree->rates->old_r[i] = tree->rates->nd_r[i];
}

/*********************************************************/

void RATES_Reset_Rates(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-2) tree->rates->nd_r[i] = tree->rates->old_r[i];
}

/*********************************************************/

void RATES_Record_Times(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-1) tree->rates->old_t[i] = tree->rates->nd_t[i];
}

/*********************************************************/

void RATES_Reset_Times(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-1) tree->rates->nd_t[i] = tree->rates->old_t[i];
}

/*********************************************************/

void RATES_Get_Rates_From_Bl(arbre *tree)
{
  phydbl dt,cr;
  node *left, *rght;
  int i;

  dt = -1.0;
  cr = tree->rates->clock_r;

  if(tree->n_root)
    {
      dt = fabs(tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->n_root->v[0]->num]);
      tree->rates->nd_r[tree->n_root->v[0]->num] = 0.5 * tree->e_root->l / (dt*cr);
      dt = fabs(tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->n_root->v[1]->num]);
      tree->rates->nd_r[tree->n_root->v[1]->num] = 0.5 * tree->e_root->l / (dt*cr);
    }
  

  For(i,2*tree->n_otu-3)
    {
      if(tree->t_edges[i] != tree->e_root)
	{
	  left = tree->t_edges[i]->left;
	  rght = tree->t_edges[i]->rght;
	  dt = fabs(tree->rates->nd_t[left->num] - tree->rates->nd_t[rght->num]);	  
	  
	  if(left->anc == rght) tree->rates->nd_r[left->num] = tree->t_edges[i]->l / (dt*cr);
	  else                  tree->rates->nd_r[rght->num] = tree->t_edges[i]->l / (dt*cr);
	}
    }

  RATES_Adjust_Clock_Rate(tree);
}

/*********************************************************/

phydbl RATES_Lk_Jumps(arbre *tree)
{
  int i,n_jps;
  phydbl dens,dt,lexp;
  node *n;

  n = NULL;
  lexp = tree->rates->lexp;
  n_jps = 0;
  dt = 0.0;
  dens = 0.0;

  For(i,2*tree->n_otu-2)
    {
      n = tree->noeud[i];
      dt = fabs(tree->rates->nd_t[n->num]-tree->rates->nd_t[n->anc->num]);
      n_jps = tree->rates->n_jps[n->num];
      dens += log(Dpois(n_jps,lexp*dt));
    }

  tree->rates->c_lnL_jps = dens;
  
  return dens;
}

/*********************************************************/

void RATES_Posterior_Rates(arbre *tree)
{
/*   RATES_Posterior_Rates_Pre(tree->n_root,tree->n_root->v[0],tree); */
/*   RATES_Posterior_Rates_Pre(tree->n_root,tree->n_root->v[1],tree); */

  int node_num;
  node_num = Rand_Int(0,2*tree->n_otu-3);
  RATES_Posterior_Rates_Pre(tree->noeud[node_num]->anc,tree->noeud[node_num],tree);
}

/*********************************************************/

void RATES_Posterior_Times(arbre *tree)
{
/*   RATES_Posterior_Times_Pre(tree->n_root,tree->n_root->v[0],tree); */
/*   RATES_Posterior_Times_Pre(tree->n_root,tree->n_root->v[1],tree); */

  int node_num;
  node_num = Rand_Int(tree->n_otu,2*tree->n_otu-3);
  
  if(tree->rates->met_within_gibbs)
    MCMC_Times_Pre(tree->noeud[node_num]->anc,tree->noeud[node_num],1,tree);
  else
    RATES_Posterior_Times_Pre(tree->noeud[node_num]->anc,tree->noeud[node_num],tree);
}

/*********************************************************/

void RATES_Posterior_Rates_Pre(node *a, node *d, arbre *tree)
{
  phydbl like_mean, like_var;
  phydbl prior_mean, prior_var;
  phydbl post_mean, post_var, post_sd;
  phydbl dt,el,vl,ra,rd,min_r,max_r,cr,nu,min_dt,cel,cvl,cer,cvr;
  phydbl new_l;
  int nb;
  short int *is_1;
  phydbl *cond_mu, *cond_cov; 
  phydbl l_opp; /* length of the branch connected to the root, opposite to the one connected to d */
  edge *b;
  int i;
  node *v2,*v3;
  phydbl T0,T1,T2,T3;
  phydbl U0,U1,U2,U3;

/*   if(a != tree->n_root) return; */

  is_1     = (short int *)mCalloc(2*tree->n_otu-2,sizeof(short int));
  cond_mu  = (phydbl *)mCalloc(1,sizeof(phydbl));
  cond_cov = (phydbl *)mCalloc(1,sizeof(phydbl));

  b = NULL;
  if(a == tree->n_root) b = NULL;
  else For(i,3) if(d->v[i] == a) { b = d->b[i]; break; }

  
  v2 = v3 = NULL;
  if(!d->tax)
    {
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  {
	    if(!v2) { v2 = d->v[i]; }
	    else    { v3 = d->v[i]; }
	  }
    }

  T3 = T2 = 0.0;
  T0 = tree->rates->nd_t[a->num];
  T1 = tree->rates->nd_t[d->num];
  U0 = tree->rates->nd_r[a->num];
  U1 = tree->rates->nd_r[d->num];
  U3 = U2 = -1.0;

  if(!d->tax)
    {
      T2  = tree->rates->nd_t[v2->num];
      T3  = tree->rates->nd_t[v3->num];
      U2  = tree->rates->nd_r[v2->num];
      U3  = tree->rates->nd_r[v3->num];
    }

  nb     = 2*tree->n_otu-3;
  dt     = tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num];
  el     = (b)?(tree->rates->u_ml_l[b->num]):(tree->rates->u_ml_l[tree->e_root->num]);
  vl     = (b)?(tree->rates->cov[b->num*nb+b->num]):(tree->rates->cov[tree->e_root->num*nb+tree->e_root->num]);
  ra     = tree->rates->nd_r[a->num];
  rd     = tree->rates->nd_r[d->num];
  min_r  = tree->rates->min_rate;
  max_r  = tree->rates->max_rate;
  cr     = tree->rates->clock_r;
  nu     = tree->rates->nu;
  min_dt = tree->rates->min_dt;
  cel    = el;   /* conditional branch length expectation */
  cvl    = vl;   /* conditional branch length variance */
  cer    = -1.0; /* conditional rate expectation */
  cvr    = -1.0; /* conditional rate variance */

  l_opp  = -1.0;
  if(a == tree->n_root)
    {
      if(d == tree->n_root->v[0])
	{
	  l_opp = tree->rates->cur_l[tree->n_root->v[1]->num];
	}
      else
	{
	  l_opp = tree->rates->cur_l[tree->n_root->v[0]->num];
	}
    }
  
  if(a == tree->n_root) is_1[tree->e_root->num] = 1;
  else                  is_1[b->num] = 1;

  /* Likelihood */
  Normal_Conditional(tree->rates->u_ml_l,tree->rates->cov,tree->rates->u_cur_l,2*tree->n_otu-3,is_1,1,cond_mu,cond_cov);
  
  cel = cond_mu[0];
  cvl = cond_cov[0];
  
  like_mean = cel;
  like_var  = cvl;
  
  /* Prior */
  For(i,2*tree->n_otu-2) is_1[i] = 0;
  is_1[d->num] = 1;
  
  RATES_Covariance_Mu(tree); /* TO DO: only need to do that once for all the rates (i.e., the cov matrix
				stays the same as long as the times do not change */

  Normal_Conditional(tree->rates->mean_r,tree->rates->cov_r,tree->rates->nd_r,2*tree->n_otu-2,is_1,1,cond_mu,cond_cov);

  cer = cond_mu[0];
  cvr = cond_cov[0];

  prior_mean = log(cer) + log(dt) + log(cr);
  prior_mean = exp(prior_mean);
  if(a == tree->n_root) prior_mean += l_opp;

  prior_var = log(cvr) + 2.*log(dt) + 2.*log(cr) ;
  prior_var = exp(prior_var);
  

  /* Posterior */
  post_mean = (prior_mean/prior_var + like_mean/like_var)/(1./prior_var + 1./like_var);

  post_var  = 1./(1./prior_var + 1./like_var);
  post_sd   = sqrt(post_var);


  if(a == tree->n_root)
    new_l = Rnorm_Trunc(post_mean,post_sd,l_opp+BL_MIN,BL_MAX);
  else
    new_l = Rnorm_Trunc(post_mean,post_sd,BL_MIN,BL_MAX);
    
  if(a == tree->n_root)
    rd = (new_l-l_opp)/(dt*cr);
  else
    rd = new_l / (dt*cr);
  
  if(isnan(rd))
    {
      printf("\n. rd=%G new_l=%G prior_mean=%G prior_var=%G like_mean=%G like_var=%G post_mean=%G post_var=%G dt=%f ra=%G cr=%G root=%d",
	     rd,new_l,
	     prior_mean,prior_var,
	     like_mean,like_var,
	     post_mean,post_var,
	     dt,ra,cr,(a==tree->n_root)?(1):(0));
      printf("\n. t(root)=%f t0=%f t1=%f",
	     tree->rates->nd_t[tree->n_root->num],
	     tree->rates->nd_t[tree->n_root->v[0]->num],
	     tree->rates->nd_t[tree->n_root->v[1]->num]);
      
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }
  
  tree->rates->nd_r[d->num] = rd;

  RATES_Update_Cur_Bl(tree);

  /* Optional here but useful for monitoring ESS for likelihoods */
  tree->c_lnL        = Dnorm_Multi_Given_InvCov_Det(tree->rates->u_cur_l,tree->rates->u_ml_l,tree->rates->invcov,tree->rates->covdet,2*tree->n_otu-3,YES);
  tree->rates->c_lnL = RATES_Lk_Rates(tree);

  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc->out_fp,tree);

  Free(is_1);
  Free(cond_mu);
  Free(cond_cov);



/*   if(Log_Dnorm(U1,U0,sqrt(tree->rates->nu*(T1-T0))) < -2.84) */
/*     { */
/*       PhyML_Printf("\n. RATES i.\n"); */
/*       PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO")); */
/*       PhyML_Printf("\n. U0 = %f U1 = %f U2 = %f U3 = %f",U0,U1,U2,U3); */
/*       PhyML_Printf("\n. T0 = %f T1 = %f T1 = %f T2 = %f T3 = %f MIN(T2,T3)=%f",T0,T1,T1,T2,T3,MIN(T2,T3)); */
/*       PhyML_Printf("\n. (T1-T0)=%f (T2-T1)=%f (T1-T0)=%f",T1-T0,T2-T1,T1-T0); */
/*       PhyML_Printf("\n. Low variance detected --> %f",Log_Dnorm(U1,U0,sqrt(tree->rates->nu*(T1-T0)))); */
/*       PhyML_Printf("\n. %f %f %f",U1,U0,sqrt(tree->rates->nu*(T1-T0))); */
/*       Exit("\n"); */
/*     } */

/*   if(!d->tax) */
/*     { */
/*       if(Log_Dnorm(U2,U1,sqrt(tree->rates->nu*(T2-T1))) < -2.84) */
/* 	{ */
/* 	  PhyML_Printf("\n. RATES ii.\n"); */
/* 	  PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO")); */
/* 	  PhyML_Printf("\n. U0 = %f U1 = %f U2 = %f U3 = %f",U0,U1,U2,U3); */
/* 	  PhyML_Printf("\n. T0 = %f T1 = %f T1 = %f T2 = %f T3 = %f MIN(T2,T3)=%f",T0,T1,T1,T2,T3,MIN(T2,T3)); */
/* 	  PhyML_Printf("\n. (T1-T0)=%f (T2-T1)=%f (T2-T1)=%f",T1-T0,T2-T1,T2-T1); */
/* 	  PhyML_Printf("\n. Low variance detected --> %f",Log_Dnorm(U2,U1,sqrt(tree->rates->nu*(T2-T1)))); */
/* 	  PhyML_Printf("\n. %f %f %f",U2,U1,sqrt(tree->rates->nu*(T2-T1))); */
/* 	  Exit("\n"); */
/* 	} */

/*       if(Log_Dnorm(U3,U1,sqrt(tree->rates->nu*(T3-T1))) < -2.84) */
/* 	{ */
/* 	  PhyML_Printf("\n. RATES iii.\n"); */
/* 	  PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO")); */
/* 	  PhyML_Printf("\n. U0 = %f U1 = %f U2 = %f U3 = %f",U0,U1,U2,U3); */
/* 	  PhyML_Printf("\n. T0 = %f T1 = %f T1 = %f T2 = %f T3 = %f MIN(T2,T3)=%f",T0,T1,T1,T2,T3,MIN(T2,T3)); */
/* 	  PhyML_Printf("\n. (T1-T0)=%f (T2-T1)=%f (T3-T1)=%f",T1-T0,T2-T1,T3-T1); */
/* 	  PhyML_Printf("\n. Low variance detected --> %f",Log_Dnorm(U3,U1,sqrt(tree->rates->nu*(T3-T1)))); */
/* 	  PhyML_Printf("\n. %f %f %f",U3,U1,sqrt(tree->rates->nu*(T3-T1))); */
/* 	  Exit("\n"); */
/* 	} */
/*     } */



/*   if(d->tax) return; */
/*   else */
/*     { */
/*       int i; */

/*       For(i,3) */
/* 	if((d->v[i] != a) && (d->b[i] != tree->e_root)) */
/* 	  RATES_Posterior_Rates_Pre(d,d->v[i],tree); */
/*     } */
}

/*********************************************************/

void RATES_Posterior_Times_Pre(node *a, node *d, arbre *tree)
{
  /*
           T0 (a)
            |
            | L1,U1,b1
            |
           T1 (d)
           / \
 L2,U2,b2 /   \ L3,U3,b3
         /     \
       T2       T3
      (v2)     (v3)

   */

  phydbl L1,L2,L3;
  phydbl EL1,EL2,EL3;
  phydbl U0,U1,U2,U3;
  phydbl T0,T1,T2,T3;
  phydbl T1_LIM_INF, T1_LIM_SUP;
  phydbl BL_LIM_INF, BL_LIM_SUP;
  phydbl T1_new;
  phydbl X,Y;
  phydbl cr;
  int    i,j;
  phydbl *mu, *cov;
  phydbl *cond_mu, *cond_cov;
  short int *is_1;
  phydbl sig11, sig1X, sig1Y, sig22, sig2X, sig2Y, sigXX, sigXY, sigYY;
  phydbl cov11,cov12,cov13,cov22,cov23,cov33;
  int nb;
  edge *b1, *b2, *b3;
  node *v2,*v3;
  phydbl *L2XY;
  phydbl l_opp;
  edge *b, *buff_b;
  node *buff_n;

  if(d->tax) return;

  b = NULL;
  if(a == tree->n_root) b = NULL;
  else For(i,3) if(d->v[i] == a) { b = d->b[i]; break; }

  L2XY     = (phydbl *)mCalloc(3,sizeof(phydbl));
  mu       = (phydbl *)mCalloc(3,sizeof(phydbl));
  cov      = (phydbl *)mCalloc(9,sizeof(phydbl));
  cond_mu  = (phydbl *)mCalloc(3,sizeof(phydbl));
  cond_cov = (phydbl *)mCalloc(9,sizeof(phydbl));
  is_1     = (short int *)mCalloc(3,sizeof(short int));

  if(b) b1 = b;
  else  b1 = tree->e_root;
  b2 = b3 = NULL;
  v2 = v3 = NULL;
  For(i,3)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      {
	if(!v2) { v2 = d->v[i]; b2 = d->b[i]; }
	else    { v3 = d->v[i]; b3 = d->b[i]; }
      }

  T2  = tree->rates->nd_t[v2->num];
  T3  = tree->rates->nd_t[v3->num];

  buff_n = NULL;
  buff_b = NULL;
  if(T3 > T2)
    {
      buff_n = v2;
      v2     = v3;
      v3     = buff_n;

      buff_b = b2;
      b2     = b3;
      b3     = buff_b;
    }
  
  nb = 2*tree->n_otu-3;
  
  is_1[0] = 1; is_1[1] = 0; is_1[2] = 0;
  
  T0  = tree->rates->nd_t[a->num];
  T1  = tree->rates->nd_t[d->num];
  T2  = tree->rates->nd_t[v2->num];
  T3  = tree->rates->nd_t[v3->num];
  U0  = tree->rates->nd_r[a->num]*tree->rates->clock_r;
  U1  = tree->rates->nd_r[d->num]*tree->rates->clock_r;
  U2  = tree->rates->nd_r[v2->num]*tree->rates->clock_r;
  U3  = tree->rates->nd_r[v3->num]*tree->rates->clock_r;
  EL1 = tree->rates->u_ml_l[b1->num];
  EL2 = tree->rates->u_ml_l[b2->num];
  EL3 = tree->rates->u_ml_l[b3->num];
  L1  = tree->rates->cur_l[d->num];
  L2  = tree->rates->cur_l[v2->num];
  L3  = tree->rates->cur_l[v3->num];
  cr  = tree->rates->clock_r;

  T1_new = +1;

  T1_LIM_INF = T0 + (1./tree->rates->nu)*pow((U1-U0)/1.96,2);
  T1_LIM_SUP = MIN(T3 - (1./tree->rates->nu)*pow((U1-U3)/1.96,2),T2 - (1./tree->rates->nu)*pow((U1-U2)/1.96,2));

  BL_LIM_INF = BL_LIM_SUP = -1.0;

  if(T1_LIM_SUP < T1_LIM_INF)
    {
      PhyML_Printf("\n. T1_LIM_SUP = %f T1_LIM_INF=%f",T1_LIM_SUP,T1_LIM_INF);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }


  l_opp = -1.;
  if(a == tree->n_root)
    {
      if(d == tree->n_root->v[0])
	l_opp = tree->rates->cur_l[tree->n_root->v[1]->num];
      else
	l_opp = tree->rates->cur_l[tree->n_root->v[0]->num];

      L1 =
	tree->rates->cur_l[tree->n_root->v[0]->num] +
	tree->rates->cur_l[tree->n_root->v[1]->num] ;
    }
   
  /* Not used but here for the sake of clarity (we want X=0 and Y=0) */
  X = L1/U1 + T0 + L2/U2 - T2;
  Y = L1/U1 + T0 + L3/U3 - T3;
  
  L2XY[0] = 0.0; /* does not matter */
  L2XY[1] = 0.0; /* constraint 1 (X=0) */
  L2XY[2] = 0.0; /* constraint 2 (Y=0) */
  

  if((T1 - T0) > (T2 - T1))
    mu[0] = EL1;
  else
    mu[0] = EL2;

  mu[1] = EL1/U1 + T0 + EL2/U2 - T2;
  mu[2] = EL1/U1 + T0 + EL3/U3 - T3;

  if(a == tree->n_root)
    {
      mu[1] = (EL1 - l_opp)/U1 + T0 + EL2/U2 - T2;
      mu[2] = (EL1 - l_opp)/U1 + T0 + EL3/U3 - T3;
    }

  cov11 = tree->rates->cov[b1->num*nb+b1->num];
  cov12 = tree->rates->cov[b1->num*nb+b2->num];
  cov13 = tree->rates->cov[b1->num*nb+b3->num];
  cov23 = tree->rates->cov[b2->num*nb+b3->num];
  cov22 = tree->rates->cov[b2->num*nb+b2->num];
  cov33 = tree->rates->cov[b3->num*nb+b3->num];

  /* printf("\n. T%d %G %G %G",d->num,cov11,cov22,cov33); */

  sig11 = cov11;
  sig22 = cov22;
  sig1X = (1./U1) * cov11 + (1./U2) * cov12;
  sig1Y = (1./U1) * cov11 + (1./U3) * cov13;
  sig2X = (1./U1) * cov12 + (1./U2) * cov22;
  sig2Y = (1./U1) * cov12 + (1./U3) * cov23;
  sigXX = 1./(U1*U1) * cov11 + 1./(U2*U2) * cov22 + 2./(U1*U2) * cov12;
  sigYY = 1./(U1*U1) * cov11 + 1./(U3*U3) * cov33 + 2./(U1*U3) * cov13;
  sigXY = 1./(U1*U1) * cov11 + 1./(U1*U2) * cov12 + 1./(U1*U3) * cov13 + 1./(U2*U3) * cov23;

  if((T1 - T0) > (T2 - T1))
    {
      cov[0*3+0] = sig11; cov[0*3+1] = sig1X; cov[0*3+2] = sig1Y;
      cov[1*3+0] = sig1X; cov[1*3+1] = sigXX; cov[1*3+2] = sigXY;
      cov[2*3+0] = sig1Y; cov[2*3+1] = sigXY; cov[2*3+2] = sigYY;
    }
  else
    {
      cov[0*3+0] = sig22; cov[0*3+1] = sig2X; cov[0*3+2] = sig2Y;
      cov[1*3+0] = sig2X; cov[1*3+1] = sigXX; cov[1*3+2] = sigXY;
      cov[2*3+0] = sig2Y; cov[2*3+1] = sigXY; cov[2*3+2] = sigYY;
    }

  Normal_Conditional(mu,cov,L2XY,3,is_1,1,cond_mu,cond_cov);
  
  if(cond_cov[0*3+0] < 0.0)
    {
      cond_cov[0*3+0] = 1.E-10;

      PhyML_Printf("\n. Conditional mean=%G var=%G",cond_mu[0],cond_cov[0*3+0]);
      PhyML_Printf("\n. T0=%G T1=%f T2=%f L1=%G L2=%G EL1=%G EL2=%G Nu=%G",T0,T1,T2,L1,L2,EL1,EL2,tree->rates->nu);
      
      printf("\n");
      For(i,3)
	{
	  printf(".mu%d=%12lf\t",i,mu[i]);
	  For(j,3)
	    {
	      printf("%12lf ",cov[i*3+j]);
	    }
	  printf("\n");
	}
    }

  
  if((T1 - T0) > (T2 - T1))
    {
      BL_LIM_INF = (T1_LIM_INF - T0) * U1;
      BL_LIM_SUP = (T1_LIM_SUP - T0) * U1;

      if(a == tree->n_root)
	{
	  L1 = Rnorm_Trunc(cond_mu[0],sqrt(cond_cov[0*3+0]),l_opp+BL_LIM_INF,l_opp+BL_LIM_SUP);
	  T1_new = (L1-l_opp)/U1 + T0;
	}
      else
	{
	  L1 = Rnorm_Trunc(cond_mu[0],sqrt(cond_cov[0*3+0]),BL_LIM_INF,BL_LIM_SUP);
	  T1_new = L1/U1 + T0;
	}
    }
  else
    {
      BL_LIM_INF = (T2 - T1_LIM_SUP) * U2;
      BL_LIM_SUP = (T2 - T1_LIM_INF) * U2; 

      L2 = Rnorm_Trunc(cond_mu[0],sqrt(cond_cov[0*3+0]),BL_LIM_INF,BL_LIM_SUP);
      T1_new = -L2/U2 + T2;
    }
  

  if(T1_new < T0)
    {
      PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
      PhyML_Printf("\n. T0 = %f T1_new = %f T1 = %f",T0,T1_new,T1);
      PhyML_Printf("\n. T1_LIM_INF=%f T1_LIM_SUP=%f",T1_LIM_INF,T1_LIM_SUP);
      PhyML_Printf("\n. L1 = %f",L1);
      PhyML_Printf("\n. BL_LIM_INF = %f BL_LIM_SUP = %f",BL_LIM_INF,BL_LIM_SUP);
      PhyML_Printf("\n. (T1-T0)=%f (T2-T1)=%f",T1-T0,T2-T1);
      PhyML_Printf("\n. L1 = %f L2 = %f cov11=%f cov22=%f cov33=%f",L1,L2,cov11,cov22,cov33);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }
  if(T1_new > MIN(T2,T3))
    {
      PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
      PhyML_Printf("\n. T0 = %f T1_new = %f T1 = %f T2 = %f T3 = %f MIN(T2,T3)=%f",T0,T1_new,T1,T2,T3,MIN(T2,T3));
      PhyML_Printf("\n. T1_LIM_INF=%f T1_LIM_SUP=%f",T1_LIM_INF,T1_LIM_SUP);
      PhyML_Printf("\n. L2 = %f",L2);
      PhyML_Printf("\n. BL_LIM_INF = %f BL_LIM_SUP = %f",BL_LIM_INF,BL_LIM_SUP);
      PhyML_Printf("\n. (T1-T0)=%f (T2-T1)=%f",T1-T0,T2-T1);
      PhyML_Printf("\n. L1 = %f L2 = %f cov11=%f cov22=%f cov33=%f",L1,L2,cov11,cov22,cov33);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  if(isnan(T1_new))
    {
      PhyML_Printf("\n. run=%d",tree->mcmc->run);
      PhyML_Printf("\n. mean=%G var=%G",cond_mu[0],cond_cov[0*3+0]);
      PhyML_Printf("\n. T1=%f L1=%G U1=%G T0=%G",T1,L1,U1,T0);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }
  
  if(Log_Dnorm(U1,U0,sqrt(tree->rates->nu*(T1_new-T0))) < -2.84)
    {
      PhyML_Printf("\n. TIMES i.\n");
      PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
      PhyML_Printf("\n. U0 = %f U1 = %f U2 = %f U3 = %f",U0,U1,U2,U3);
      PhyML_Printf("\n. T0 = %f T1_new = %f T1 = %f T2 = %f T3 = %f MIN(T2,T3)=%f",T0,T1_new,T1,T2,T3,MIN(T2,T3));
      PhyML_Printf("\n. T1_LIM_INF=%f T1_LIM_SUP=%f",T1_LIM_INF,T1_LIM_SUP);
      PhyML_Printf("\n. (T1-T0)=%f (T2-T1)=%f (T1_new-T0)=%f",T1-T0,T2-T1,T1_new-T0);
      PhyML_Printf("\n. Low variance detected --> %f",Log_Dnorm(U1,U0,sqrt(tree->rates->nu*(T1_new-T0))));
      PhyML_Printf("\n. %f %f %f",U1,U0,sqrt(tree->rates->nu*(T1_new-T0)));
      Exit("\n");
    }

  if(Log_Dnorm(U2,U1,sqrt(tree->rates->nu*(T2-T1_new))) < -2.84)
    {
      PhyML_Printf("\n. TIMES ii.\n");
      PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
      PhyML_Printf("\n. U0 = %f U1 = %f U2 = %f U3 = %f",U0,U1,U2,U3);
      PhyML_Printf("\n. T0 = %f T1_new = %f T1 = %f T2 = %f T3 = %f MIN(T2,T3)=%f",T0,T1_new,T1,T2,T3,MIN(T2,T3));
      PhyML_Printf("\n. T1_LIM_INF=%f T1_LIM_SUP=%f",T1_LIM_INF,T1_LIM_SUP);
      PhyML_Printf("\n. (T1-T0)=%f (T2-T1)=%f (T2-T1_new)=%f",T1-T0,T2-T1,T2-T1_new);
      PhyML_Printf("\n. Low variance detected --> %f",Log_Dnorm(U2,U1,sqrt(tree->rates->nu*(T2-T1_new))));
      PhyML_Printf("\n. %f %f %f",U2,U1,sqrt(tree->rates->nu*(T2-T1_new)));
      Exit("\n");
    }

  if(Log_Dnorm(U3,U1,sqrt(tree->rates->nu*(T3-T1_new))) < -2.84)
    {
      PhyML_Printf("\n. TIMES iii.\n");
      PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
      PhyML_Printf("\n. U0 = %f U1 = %f U2 = %f U3 = %f",U0,U1,U2,U3);
      PhyML_Printf("\n. T0 = %f T1_new = %f T1 = %f T2 = %f T3 = %f MIN(T2,T3)=%f",T0,T1_new,T1,T2,T3,MIN(T2,T3));
      PhyML_Printf("\n. T1_LIM_INF=%f T1_LIM_SUP=%f",T1_LIM_INF,T1_LIM_SUP);
      PhyML_Printf("\n. (T1-T0)=%f (T2-T1)=%f (T3-T1_new)=%f",T1-T0,T2-T1,T3-T1_new);
      PhyML_Printf("\n. Low variance detected --> %f",Log_Dnorm(U3,U1,sqrt(tree->rates->nu*(T3-T1_new))));
      PhyML_Printf("\n. %f %f %f",U3,U1,sqrt(tree->rates->nu*(T3-T1_new)));
      Exit("\n");
    }
    


  tree->rates->nd_t[d->num] = T1_new;
  RATES_Update_Cur_Bl(tree);
  
  tree->c_lnL        = Dnorm_Multi_Given_InvCov_Det(tree->rates->u_cur_l,tree->rates->u_ml_l,tree->rates->invcov,tree->rates->covdet,2*tree->n_otu-3,YES);
  tree->rates->c_lnL = RATES_Lk_Rates(tree);

  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc->out_fp,tree);

  Free(mu);
  Free(cov);
  Free(is_1);
  Free(L2XY);

  For(i,3)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      RATES_Posterior_Times_Pre(d,d->v[i],tree);
}

/*********************************************************/

void RATES_Update_Cur_Bl(arbre *tree)
{
  RATES_Update_Cur_Bl_Pre(tree->n_root,tree->n_root->v[0],NULL,tree);
  RATES_Update_Cur_Bl_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);

  tree->e_root->l = 
    tree->rates->cur_l[tree->n_root->v[0]->num] +
    tree->rates->cur_l[tree->n_root->v[1]->num] ;

  tree->rates->u_cur_l[tree->e_root->num] = tree->e_root->l;

  if(tree->e_root->l < 0.0)
    {
      PhyML_Printf("\n. l=%f",tree->e_root->l);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }
}

/*********************************************************/

void RATES_Update_Cur_Bl_Pre(node *a, node *d, edge *b, arbre *tree)
{
  phydbl down_t,up_t,dt,rr,cr;
  
  down_t = tree->rates->nd_t[d->num];
  up_t   = tree->rates->nd_t[a->num];
  dt     = down_t - up_t;
  rr     = tree->rates->nd_r[d->num];
  cr     = tree->rates->clock_r;

  tree->rates->cur_l[d->num] = dt*rr*cr;

  if(tree->rates->cur_l[d->num] < 0.0)
    {
      printf("\n. l=%G dt=%G rr=%G cr=%G up_t=%G down_t=%G",tree->rates->cur_l[d->num],dt,rr,cr,up_t,down_t);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  if(b) 
    {
      tree->rates->u_cur_l[b->num] = tree->rates->cur_l[d->num];
      b->l = tree->rates->cur_l[d->num];
    }

  if(b && isnan(b->l))
    {
      printf("\n. dt=%G rr=%G cr=%G",dt,rr,cr);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }
  


  if(d->tax) return;
  else
    {
      int i;
      For(i,3) 
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  RATES_Update_Cur_Bl_Pre(d,d->v[i],d->b[i],tree);
    }
}

/*********************************************************/

void RATES_Bl_To_Ml(arbre *tree)
{
  RATES_Bl_To_Ml_Pre(tree->n_root,tree->n_root->v[0],NULL,tree);
  RATES_Bl_To_Ml_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);
  tree->rates->u_ml_l[tree->e_root->num] = tree->t_edges[tree->e_root->num]->l;
  tree->rates->ml_l[tree->n_root->v[0]->num] = tree->rates->u_ml_l[tree->e_root->num] * tree->n_root_pos;
  tree->rates->ml_l[tree->n_root->v[1]->num] = tree->rates->u_ml_l[tree->e_root->num] * (1. - tree->n_root_pos);
}

/*********************************************************/

void RATES_Bl_To_Ml_Pre(node *a, node *d, edge *b, arbre *tree)
{

  if(b) 
    {
      tree->rates->u_ml_l[b->num] = b->l;
      tree->rates->ml_l[d->num] = tree->rates->cur_l[d->num];
    }

  if(d->tax) return;
  else
    {
      int i;

      For(i,3) 
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  RATES_Bl_To_Ml_Pre(d,d->v[i],d->b[i],tree);
    }
}

/*********************************************************/

void RATES_Get_Cov_Matrix_Rooted(phydbl *unroot_cov, arbre *tree)
{
  int i,dim;

  dim = 2*tree->n_otu-3;

  RATES_Get_Cov_Matrix_Rooted_Pre(tree->n_root,tree->n_root->v[0],NULL,unroot_cov,tree);
  RATES_Get_Cov_Matrix_Rooted_Pre(tree->n_root,tree->n_root->v[1],NULL,unroot_cov,tree);

  For(i,dim+1) tree->rates->cov[i*(dim+1)+tree->n_root->v[0]->num] /= 2.;
  For(i,dim+1) tree->rates->cov[i*(dim+1)+tree->n_root->v[1]->num] /= 2.;
  For(i,dim+1) tree->rates->cov[tree->n_root->v[0]->num*(dim+1)+i] /= 2.;
  For(i,dim+1) tree->rates->cov[tree->n_root->v[1]->num*(dim+1)+i] /= 2.;

}

/*********************************************************/

void RATES_Get_Cov_Matrix_Rooted_Pre(node *a, node *d, edge *b, phydbl *cov, arbre *tree)
{
  int i, dim, done_left;
  node *n;

  dim       = 2*tree->n_otu-3;
  done_left = 0;
  n         = NULL;

  For(i,dim) 
    { 
      if(tree->t_edges[i] != tree->e_root)
	{
	  n = 
	    (tree->t_edges[i]->left->anc == tree->t_edges[i]->rght)?
	    (tree->t_edges[i]->left):
	    (tree->t_edges[i]->rght);

	  if(b)
	    {
	      tree->rates->cov[d->num*(dim+1) + n->num] = cov[b->num*dim + i];
	    }
	  else
	    {
	      tree->rates->cov[d->num*(dim+1) + n->num] = cov[tree->e_root->num*dim + i];
	    }
	}
      else
	{
	  n = tree->e_root->left;
	  if(b)
	    tree->rates->cov[d->num*(dim+1) + n->num] = cov[b->num*dim + i];
	  else
	    tree->rates->cov[d->num*(dim+1) + n->num] = cov[tree->e_root->num*dim + i];

	  n = tree->e_root->rght;
	  if(b)
	    tree->rates->cov[d->num*(dim+1) + n->num] = cov[b->num*dim + i];
	  else
	    tree->rates->cov[d->num*(dim+1) + n->num] = cov[tree->e_root->num*dim + i];
	}
    }


  if(d->tax) return;
  else
    {
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  RATES_Get_Cov_Matrix_Rooted_Pre(d,d->v[i],d->b[i],cov,tree);
    }
}

/*********************************************************/

void RATES_Covariance_Mu(arbre *tree)
{
  int i,j;
  phydbl dt,var;
  int dim;
  int lca_num;

  dim = 2*tree->n_otu-2;

  For(i,dim*dim) tree->rates->cov_r[i] = 0.0;

  dt =  tree->rates->nd_t[tree->n_root->v[0]->num] - tree->rates->nd_t[tree->n_root->num];
  var = dt * tree->rates->nu;
  tree->rates->cov_r[tree->n_root->v[0]->num*dim+tree->n_root->v[0]->num] = var;

  dt = tree->rates->nd_t[tree->n_root->v[1]->num] - tree->rates->nd_t[tree->n_root->num];
  var = dt * tree->rates->nu;
  tree->rates->cov_r[tree->n_root->v[1]->num*dim+tree->n_root->v[1]->num] = var;

  RATES_Variance_Mu_Pre(tree->n_root,tree->n_root->v[0],tree);
  RATES_Variance_Mu_Pre(tree->n_root,tree->n_root->v[1],tree);

  For(i,dim)
    {
      for(j=i+1;j<dim;j++)
	{
	  lca_num = tree->rates->lca[i*(dim+1)+j]->num;
	  if(lca_num < dim) 
	    {
	      tree->rates->cov_r[i*dim+j] = tree->rates->cov_r[lca_num*dim+lca_num];	    
	      tree->rates->cov_r[j*dim+i] = tree->rates->cov_r[i*dim+j];
	    }
	  else if(lca_num == dim)
	    {
	      tree->rates->cov_r[i*dim+j] = 0.0;	    
	      tree->rates->cov_r[j*dim+i] = 0.0;
	    }
	  else
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Exit("\n");
	    }
	}
    }
}

/*********************************************************/

void RATES_Variance_Mu_Pre(node *a, node *d, arbre *tree)
{
  int dim;
  phydbl dt0,var0;
  phydbl dt1,var1;
  phydbl dt2,var2;
  int i;
  int dir1, dir2;

  dim = 2*tree->n_otu-2;

  if(d->tax) return;

  dir1 = dir2 = -1;
  For(i,3)
    {
      if((d->v[i] != a) && (d->b[i] != tree->e_root))
	{
	  if(dir1 < 0) dir1 = i;
	  else         dir2 = i;
	}
    }


  dt0  = tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num];
  var0 = tree->rates->cov_r[d->num*dim+d->num];

  dt1  = tree->rates->nd_t[d->v[dir1]->num] - tree->rates->nd_t[d->num];
  var1 = tree->rates->nu*dt1;

  dt2  = tree->rates->nd_t[d->v[dir2]->num] - tree->rates->nd_t[d->num];
  var2 = tree->rates->nu*dt2;

  tree->rates->cov_r[d->v[dir1]->num*dim+d->v[dir1]->num] = var0+var1;
  tree->rates->cov_r[d->v[dir2]->num*dim+d->v[dir2]->num] = var0+var2;


  For(i,3)
    {
      if((d->v[i] != a) && (d->b[i] != tree->e_root))
	{
	  RATES_Variance_Mu_Pre(d,d->v[i],tree);
	}
    }
}

/*********************************************************/

void RATES_Fill_Lca_Table(arbre *tree)
{
  int i,j;
  int dim;
  
  dim = 2*tree->n_otu-1;

  For(i,dim)
    {
      for(j=i;j<dim;j++)
	{
	  tree->rates->lca[i*dim+j] = Find_Lca(tree->noeud[i],tree->noeud[j],tree);
	  tree->rates->lca[j*dim+i] = tree->rates->lca[i*dim+j];
	}
    }
}

/*********************************************************/
/*********************************************************/
/*********************************************************/

#endif
