/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/


/* Routines for molecular clock trees and molecular dating */
#if defined (MC) ||  defined (RWRAPPER)

#include "rates.h"

#ifdef RWRAPPER
#include <R.h>
#endif



/*********************************************************/

phydbl RATES_Lk_Rates(t_tree *tree)
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

void RATES_Lk_Rates_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  int i,n1,n2;
  phydbl log_dens,mu1,mu2,dt1,dt2;
  
  log_dens = -1.;

  if(d->anc != a)
    {
      PhyML_Printf("\n. d=%d d->anc=%d a=%d root=%d",d->num,d->anc->num,a->num,tree->n_root->num);
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
      
      if(mu2 < tree->rates->min_rate) mu2 = tree->rates->min_rate;
      if(mu2 > tree->rates->max_rate) mu2 = tree->rates->max_rate; 

      if(dt2 < tree->rates->min_dt)   dt2 = tree->rates->min_dt;
  
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
	    int err;
	    log_dens = Log_Dnorm(mu2,1.0,sqrt(tree->rates->nu*dt2),&err);

	    if(err)
	      {
		PhyML_Printf("\n. log_dens=%f",log_dens);
		PhyML_Printf("\n. mu2=%f dt2=%f",mu2,dt2);
		PhyML_Printf("\n. a->t=%f %f",tree->rates->nd_t[a->num],tree->rates->nd_t[d->num]);
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
/* 		Exit("\n"); */
	      }

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
      MCMC_Print_Param(tree->mcmc,tree);
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

phydbl RATES_Lk_Change_One_Rate(t_node *d, phydbl new_rate, t_tree *tree)
{
  tree->rates->nd_r[d->num] = new_rate;
  RATES_Update_Triplet(d,tree);
  RATES_Update_Triplet(d->anc,tree);
  return(tree->rates->c_lnL);
}

/*********************************************************/

phydbl RATES_Lk_Change_One_Time(t_node *n, phydbl new_t, t_tree *tree)
{  
  if(n == tree->n_root)
    {
      PhyML_Printf("\n. Moving the time of the root t_node is not permitted.");
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

void RATES_Update_Triplet(t_node *n, t_tree *tree)
{
  phydbl curr_triplet,new_triplet;
  phydbl dt0,dt1,dt2;
  phydbl mu1_mu0,mu2_mu0;
  phydbl mu0,mu1,mu2;
  int n0,n1,n2;
  int i;
  t_node *v1,*v2;

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
	    int err;

	    log_dens = 
	      Log_Dnorm(mu0,1.0,sqrt(tree->rates->nu*dt0),&err) + 
	      Log_Dnorm(mu1,1.0,sqrt(tree->rates->nu*dt1),&err); 

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
	  MCMC_Print_Param(tree->mcmc,tree);
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
phydbl RATES_Lk_Rates_Core(phydbl mu1, phydbl mu2, int n1, int n2, phydbl dt1, phydbl dt2, t_tree *tree)
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
 
  if(dt1 < tree->rates->min_dt) dt1 = tree->rates->min_dt;
  if(dt2 < tree->rates->min_dt) dt2 = tree->rates->min_dt;
  
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
	int err;

	log_dens = Log_Dnorm(mu2,mu1,sqrt(tree->rates->nu*dt2),&err);

	if(err)
	  {
	    PhyML_Printf("\n. mu1=%f mu2=%f dt2=%f nu=%G",mu1,mu2,dt2,tree->rates->nu);
	    PhyML_Printf("\n. Run: %d",tree->mcmc->run);
	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
/* 	    Exit("\n"); */
	  }

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
/*   PhyML_Printf("\n. %f %f %f %f %f ",mu1,mu2,dt1,dt2,v2); */
  v2 /= dmu1;

  dens = p0*v0 + p1*v1 + v2;
/*   dens = p1*v1 + v2; */
  /*       dens = p1*v1 + v2; */
  /*   dens = v0; */
/*   dens *= dmu1; */
  
  if(dens < MDBL_MIN)
    {
      PhyML_Printf("\n. dens=%12G mu1=%12G mu2=%12G dt1=%12G dt2=%12G lexp=%12G alpha=%f v0=%f v1=%f v2=%f p0=%f p1=%f p2=%f",
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
      PhyML_Printf("\n. n1=%d n2=%d",n1,n2);
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

/*   PhyML_Printf("\n. density = %15G mu1=%3.4f mu2=%3.4f dt1=%3.4f dt2=%3.4f n1=%2d n2=%2d",density,mu1,mu2,dt1,dt2,n1,n2); */
  return density;
}

/**********************************************************/

void RATES_Print_Triplets(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-1) PhyML_Printf("\n. Node %3d t=%f",i,tree->rates->triplet[i]);
}


/**********************************************************/

void RATES_Print_Rates(t_tree *tree)
{
  RATES_Print_Rates_Pre(tree->n_root,tree->n_root->v[0],NULL,tree);
  RATES_Print_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);
}

/*********************************************************/

void RATES_Print_Rates_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{

  
  if((d == tree->n_root->v[0] && d->tax) || (d == tree->n_root->v[1] && d->tax))
    PhyML_Printf("\n. a=%3d ++d=%3d rate=%12f t_left=%12f t_rght=%12f ml=%12f l=%12f %12f",
	   a->num,d->num,
	   tree->rates->nd_r[d->num],
	   tree->rates->nd_t[a->num],tree->rates->nd_t[d->num],
	   tree->rates->ml_l[d->num],
	   tree->rates->cur_l[d->num],
	   (tree->rates->nd_t[d->num]-tree->rates->nd_t[a->num])*tree->rates->clock_r*tree->rates->nd_r[d->num]);
  
  else if((d == tree->n_root->v[0]) || (d == tree->n_root->v[1]))
    PhyML_Printf("\n. a=%3d __d=%3d rate=%12f t_left=%12f t_rght=%12f ml=%12f l=%12f %12f",
	   a->num,d->num,
	   tree->rates->nd_r[d->num],
	   tree->rates->nd_t[a->num],tree->rates->nd_t[d->num],
	   tree->rates->ml_l[d->num],
	   tree->rates->cur_l[d->num],
	   (tree->rates->nd_t[d->num]-tree->rates->nd_t[a->num])*tree->rates->clock_r*tree->rates->nd_r[d->num]);
  else 
    PhyML_Printf("\n. a=%3d   d=%3d rate=%12f t_left=%12f t_rght=%12f ml=%12f l=%12f %12f",
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

phydbl RATES_Check_Mean_Rates(t_tree *tree)
{
  phydbl sum_r,sum_dt;
  phydbl u,t,t_anc;
  int i;
  
  sum_r  = 0.0;
  sum_dt = 0.0;
  For(i,2*tree->n_otu-2) 
    {
      u     = tree->rates->nd_r[i];
      t     = tree->rates->nd_t[i];
      t_anc = tree->rates->nd_t[tree->noeud[i]->anc->num];

      sum_r += u * fabs(t-t_anc);
      sum_dt += fabs(t-t_anc);      
    }

  return(sum_r / sum_dt);
}

/*********************************************************/

phydbl RATES_Check_Mean_Rates_True(t_tree *tree)
{
  phydbl sum;
  int i;
  
  sum = 0.0;
  For(i,2*tree->n_otu-2) sum += tree->rates->true_r[i];
  return(sum/(phydbl)(2*tree->n_otu-2));
}

/*********************************************************/

void RATES_Check_Node_Times(t_tree *tree)
{
  RATES_Check_Node_Times_Pre(tree->n_root,tree->n_root->v[0],tree);
  RATES_Check_Node_Times_Pre(tree->n_root,tree->n_root->v[1],tree);
}

/*********************************************************/

void RATES_Check_Node_Times_Pre(t_node *a, t_node *d, t_tree *tree)
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
  
  rates                 = (trate  *)mCalloc(1,sizeof(trate));
  rates->br_r           = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
  rates->old_r          = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->nd_r           = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->true_r         = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->old_t          = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->nd_t           = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->true_t         = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->t_mean         = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->t_prior        = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->t_prior_min    = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->t_prior_max    = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->t_has_prior    = (short int *)mCalloc(2*n_otu-1,sizeof(short int));
  rates->dens           = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
  rates->triplet        = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->n_jps          = (int    *)mCalloc(2*n_otu-2,sizeof(int));
  rates->t_jps          = (int    *)mCalloc(2*n_otu-2,sizeof(int));
  rates->cov            = (phydbl *)mCalloc((2*n_otu-2)*(2*n_otu-2),sizeof(phydbl));
  rates->invcov         = (phydbl *)mCalloc((2*n_otu-2)*(2*n_otu-2),sizeof(phydbl));
  rates->ml_l           = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
  rates->cur_l          = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
  rates->u_ml_l         = (phydbl *)mCalloc(2*n_otu-3,sizeof(phydbl));
  rates->u_cur_l        = (phydbl *)mCalloc(2*n_otu-3,sizeof(phydbl));
  rates->cov_r          = (phydbl *)mCalloc((2*n_otu-2)*(2*n_otu-2),sizeof(phydbl));
  rates->cond_var       = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
  rates->mean_r         = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
  rates->lca            = (t_node **)mCalloc((2*n_otu-1)*(2*n_otu-1),sizeof(t_node *));
  rates->reg_coeff      = (phydbl *)mCalloc((2*n_otu-2)*(2*n_otu-2),sizeof(phydbl));
  rates->trip_reg_coeff = (phydbl *)mCalloc((2*n_otu-2)*(6*n_otu-9),sizeof(phydbl));
  rates->trip_cond_cov  = (phydbl *)mCalloc((2*n_otu-2)*9,sizeof(phydbl));
  rates->_2n_vect1      = (phydbl *)mCalloc(2*n_otu,sizeof(phydbl));
  rates->_2n_vect2      = (phydbl *)mCalloc(2*n_otu,sizeof(phydbl));
  rates->_2n_vect3      = (phydbl *)mCalloc(2*n_otu,sizeof(phydbl));
  rates->_2n_vect4      = (phydbl *)mCalloc(2*n_otu,sizeof(phydbl));
  rates->_2n_vect5      = (short int *)mCalloc(2*n_otu,sizeof(short int));
  rates->_2n2n_vect1    = (phydbl *)mCalloc(4*n_otu*n_otu,sizeof(phydbl));
  rates->_2n2n_vect2    = (phydbl *)mCalloc(4*n_otu*n_otu,sizeof(phydbl));
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
  Free(rates->t_prior);
  Free(rates->t_mean);
  Free(rates->t_prior_min);
  Free(rates->t_prior_max);
  Free(rates->t_has_prior);
  Free(rates->dens);   
  Free(rates->triplet);    
  Free(rates->n_jps);  
  Free(rates->t_jps);    
  Free(rates->cov);
  Free(rates->cond_var);
  Free(rates->invcov);   
  Free(rates->ml_l);
  Free(rates->cur_l);
  Free(rates->u_ml_l);
  Free(rates->u_cur_l);
  Free(rates->cov_r);
  Free(rates->lca);
  Free(rates->trip_cond_cov);
  Free(rates->trip_reg_coeff);
  Free(rates->_2n_vect1);
  Free(rates->_2n_vect2);
  Free(rates->_2n_vect3);
  Free(rates->_2n_vect4);
  Free(rates->_2n_vect5);
  Free(rates->_2n2n_vect1);
  Free(rates->_2n2n_vect2);
}

/*********************************************************/

void RATES_Init_Rate_Struct(trate *rates, int n_otu)
{
  int i;

  rates->met_within_gibbs = NO;
  rates->model         = THORNE;
  rates->c_lnL         = -INFINITY;
  rates->c_lnL_jps     = -INFINITY;
  rates->adjust_rates  = 0;
  rates->use_rates     = 1;
  rates->lexp          = 1.E-3;
  rates->alpha         = 2.;
  rates->birth_rate    = 0.001;

  rates->max_rate      = 5.;
  rates->min_rate      = 0.2;

  rates->clock_r       = 1.E-3;
  rates->max_clock     = 1.E-0;
  rates->min_clock     = 1.E-9;

  rates->nu            = 1.E-4;
  rates->max_nu        = 1.E-1;
  rates->min_nu        = 1.E-10;
  rates->lbda_nu       = 1.E+3;

  rates->min_dt        = 0.0;

  rates->step_rate     = 1.E-4;
  rates->approx        = 1;
  rates->bl_from_rt    = 0;
  rates->lk_approx     = NO;

  rates->p_max         = 0.01;

  rates->true_tree_size = 0.0;

  For(i,2*n_otu-2) 
    {
      rates->br_r[i]   = 1.0;
      rates->n_jps[i]  =  -1;
      rates->t_jps[i]  =  -1;
      rates->mean_r[i] = 1.0;      
    }

  For(i,2*n_otu-1) 
    {
      rates->old_r[i]  = 1.0;
      rates->nd_r[i]   = 1.0;
      rates->nd_t[i]   = 0.0;
      rates->true_t[i] = 0.0;
      if(i < n_otu)
	{
	  rates->t_has_prior[i] = 1;
	  rates->t_prior_max[i] = 0.0;
	  rates->t_prior_min[i] = 0.0;
	}
      else
	{
	  rates->t_prior_max[i] =  MDBL_MAX;
	  rates->t_prior_min[i] = -MDBL_MAX;
	}
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
      cdf = Ppois(c,param);      
    }
  
  a = 0.0;
  b = (c-a)/2.;
  step = 0;
  do
    {
      step++;
      cdf = Ppois(b,param);
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
   instantaneous rate at t_node a, this function works out an expected number of (amino-acids or 
   nucleotide) substitutions per site.
*/
void RATES_Expect_Number_Subst(phydbl t_beg, phydbl t_end, phydbl r_beg,  int *n_jumps, phydbl *mean_r, phydbl *r_end, trate *rates, t_tree *tree)
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
	
	/*   PhyML_Printf("\n. [%3d %f %f]",*n_jumps,*mean_r,r_beg); */
	
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
	int err;

	sd = sqrt(rates->nu*fabs(t_beg-t_end));
	*mean_r = Rnorm_Trunc(r_beg,sd,rates->min_rate,rates->max_rate,&err);
	if(err) PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run);
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

void RATES_Get_Mean_Rates_Pre(t_node *a, t_node *d, t_edge *b, phydbl r_a, t_tree *tree)
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

void RATES_Random_Branch_Lengths(t_tree *tree)
{
  phydbl r0,r_d;
   
  r0 = 1.0;
  r_d = 0;

  tree->rates->nd_r[tree->n_root->num] = r0;

  RATES_Get_Mean_Rates_Pre(tree->n_root,tree->n_root->v[0],NULL,r0,tree);
  RATES_Get_Mean_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,r0,tree);

/*   /\* !!!!!!!!!!! *\/ */
/*   RATES_Normalise_Rates(tree); */

  RATES_Update_Cur_Bl(tree);
  RATES_Initialize_True_Rates(tree);

  tree->n_root_pos = 
    tree->rates->cur_l[tree->n_root->v[0]->num] /
    (tree->rates->cur_l[tree->n_root->v[0]->num] + tree->rates->cur_l[tree->n_root->v[1]->num]);
}

/*********************************************************/

void RATES_Set_Node_Times(t_tree *tree)
{
  RATES_Set_Node_Times_Pre(tree->n_root,tree->n_root->v[0],tree);
  RATES_Set_Node_Times_Pre(tree->n_root,tree->n_root->v[1],tree);
}

/*********************************************************/

void RATES_Init_Triplets(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-1) tree->rates->triplet[i] = 0.0;
}
/*********************************************************/

void RATES_Set_Node_Times_Pre(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      t_node *v1, *v2; /* the two sons of d */
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
	  PhyML_Printf("\n. t_sup = %f t_inf = %f",t_sup,t_inf);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
            
      if(tree->rates->nd_t[d->num] > t_inf)      
	{
	  tree->rates->nd_t[d->num] = t_inf;
	  PhyML_Printf("\n. Time for t_node %d is larger than should be. Set it to %f",d->num,tree->rates->nd_t[d->num]);
	}
      else if(tree->rates->nd_t[d->num] < t_sup) 
	{
	  tree->rates->nd_t[d->num] = t_sup;
	  PhyML_Printf("\n. Time for t_node %d is lower than should be. Set it to %f",d->num,tree->rates->nd_t[d->num]);
	}

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
      PhyML_Printf("\n. mu2=%G",mu2);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(mu1 < 1.E-10)
    {
      PhyML_Printf("\n. mu2=%G",mu1);
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
	      PhyML_Printf("\n. u = %G",u);
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
	  PhyML_Printf("\n. mu1 = %f lambda = %f v = %f u = %G beg = %f end = %f",mu1,lbda,v,u,beg,end);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      
      dens += Dgamma_Moments(u,mean,var) / (1.-lbda) * ndb * pow(1.-lbda,n-1);
      
      lbda = end;
      u = (mu1 - lbda*v)/(1.-lbda);
      if(u < 1.E-10)
	{
	  PhyML_Printf("\n. u = %G",u);
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
      PhyML_Printf("\n. density=%f cmpoiss = %f i=%d n_min=%d mu1=%f mu2=%f dt1=%f dt2=%f",
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
	    PhyML_Printf("\n. mu2=%f",mu2);
	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	  }

	if(mu1 < 1.E-10)
	  {
	    PhyML_Printf("\n. mu1=%f",mu1);
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
phydbl RATES_Yule(t_tree *tree)
{
  phydbl sumti,density,lambda;
  int n,i;

  sumti = 0.0;
  for(i=tree->n_otu;i<2*tree->n_otu-1;i++) sumti += tree->rates->nd_t[i];
  sumti -= tree->rates->nd_t[i-1];

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

void RATES_Initialize_True_Rates(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-2) tree->rates->true_r[i] = tree->rates->nd_r[i];
}
/*********************************************************/

void RATES_Record_Rates(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-2) tree->rates->old_r[i] = tree->rates->nd_r[i];
}

/*********************************************************/

void RATES_Reset_Rates(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-2) tree->rates->nd_r[i] = tree->rates->old_r[i];
}

/*********************************************************/

void RATES_Record_Times(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-1) tree->rates->old_t[i] = tree->rates->nd_t[i];
}

/*********************************************************/

void RATES_Reset_Times(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-1) tree->rates->nd_t[i] = tree->rates->old_t[i];
}

/*********************************************************/

void RATES_Get_Rates_From_Bl(t_tree *tree)
{
  phydbl dt,cr;
  t_node *left, *rght;
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
}

/*********************************************************/

phydbl RATES_Lk_Jumps(t_tree *tree)
{
  int i,n_jps;
  phydbl dens,dt,lexp;
  t_node *n;

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

void RATES_Posterior_Clock_Rate(t_tree *tree)
{
  if(1)
    MCMC_Clock_Rate(tree);
  else
    {
      int i,j,dim,err;
      phydbl *z,*z_mean,*z_cov;
      phydbl *r,*dt;
      phydbl *cond_mu, *cond_cov;
      short int *is_1;
      int ref_br;

      dim = 2*tree->n_otu-3;
      
      r        = (phydbl *)mCalloc(dim,sizeof(phydbl));
      dt       = (phydbl *)mCalloc(dim,sizeof(phydbl));
      
      z        = tree->rates->_2n_vect1;
      z_mean   = tree->rates->_2n_vect2;
      cond_mu  = tree->rates->_2n_vect3;
      is_1     = tree->rates->_2n_vect5;
      z_cov    = tree->rates->_2n2n_vect1;
      cond_cov = tree->rates->_2n2n_vect2;
      
      For(i,dim) if(tree->rates->u_cur_l[i] > BL_MIN) break;
      ref_br = i;

      For(i,dim)
	{
	  if(tree->t_edges[i] != tree->e_root)
	    {
	      r[i]  = tree->rates->nd_r[Edge_Num_To_Node_Num(i,tree)];
	      dt[i] = fabs(tree->rates->nd_t[tree->t_edges[i]->left->num] - tree->rates->nd_t[tree->t_edges[i]->rght->num]);
	    }
	  else
	    {
	      r[i]  =
		tree->rates->nd_r[tree->n_root->v[0]->num] * (tree->rates->nd_t[tree->n_root->v[0]->num] - tree->rates->nd_t[tree->n_root->num]) +
		tree->rates->nd_r[tree->n_root->v[1]->num] * (tree->rates->nd_t[tree->n_root->v[1]->num] - tree->rates->nd_t[tree->n_root->num]);
	      
	      dt[i] =
		(tree->rates->nd_t[tree->n_root->v[0]->num] - tree->rates->nd_t[tree->n_root->num]) +
		(tree->rates->nd_t[tree->n_root->v[1]->num] - tree->rates->nd_t[tree->n_root->num]) ;
	      
	      r[i] /= dt[i];
	    }
	}
      
      For(i,dim) z_mean[i] = tree->rates->u_ml_l[ref_br]/(r[ref_br]*dt[ref_br]) - tree->rates->u_ml_l[i]/(r[i]*dt[i]); 
      z_mean[ref_br] = tree->rates->u_ml_l[ref_br];
      
      For(i,dim) 
	{
	  For(j,dim)
	    {
	      z_cov[i*dim+j] = 
		tree->rates->cov[ref_br*dim+ref_br] / (r[ref_br]*dt[ref_br]*r[ref_br]*dt[ref_br]) - 
		tree->rates->cov[ref_br*dim+j] / (r[ref_br]*dt[ref_br]*r[j]*dt[j]) - 
		tree->rates->cov[ref_br*dim+i] / (r[ref_br]*dt[ref_br]*r[i]*dt[i]) + 
		tree->rates->cov[i*dim+j] / (r[i]*dt[i]*r[j]*dt[j]) ; 
	    }
	}
      
      For(i,dim)
	{
	  z_cov[ref_br*dim+i] = 
	    tree->rates->cov[ref_br*dim+ref_br] / (r[ref_br]*dt[ref_br]) -
	    tree->rates->cov[ref_br*dim+i] / (r[i]*dt[i]) ;
	  
	  z_cov[i*dim+ref_br] = z_cov[ref_br*dim+i];
	}
      
      z_cov[ref_br*dim+ref_br] = tree->rates->cov[ref_br*dim+ref_br];
      
      
      /*   PhyML_Printf("\n"); */
      /*   For(i,dim)  */
      /*     { */
      /*       PhyML_Printf("[%13G]  ",z_mean[i]); */
      /*       For(j,dim) */
      /* 	{ */
      /* 	  PhyML_Printf("%13G ",z_cov[i*dim+j]); */
      /* 	} */
      /*       PhyML_Printf("\n"); */
      /*     } */
      /*   PhyML_Printf("\n"); */
      

      
      For(i,dim) is_1[i] = 0;
      is_1[ref_br] = 1;
      For(i,dim) z[i] = .0;
      z[ref_br] = tree->rates->u_cur_l[ref_br];
      Normal_Conditional(z_mean,z_cov,z,dim,is_1,1,cond_mu,cond_cov);
      
      /*   PhyML_Printf("\n. L0=%G COND_L0=%G VAR_L0=%G", */
      /* 	 tree->rates->u_cur_l[0], */
      /* 	 cond_mu[0], */
      /* 	 sqrt(cond_cov[0])); */
      
      if(cond_cov[ref_br] < 0.0 || cond_mu[ref_br] < 0.0)
	{
	  PhyML_Printf("\n. cond_cov=%f cond_mu=%f",cond_cov[ref_br],cond_mu[ref_br]);
	  PhyML_Printf("\n. l=%f ref_br=%d",tree->rates->u_cur_l[ref_br],ref_br);
	  PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run);
	  Exit("\n");
	}

      tree->rates->clock_r = Rnorm_Trunc(cond_mu[ref_br],sqrt(cond_cov[ref_br*dim+ref_br]),BL_MIN,BL_MAX,&err);
      tree->rates->clock_r /= (r[ref_br]*dt[ref_br]);
      
      Free(r);
      Free(dt);
      
      /**************************************************************/
      /**************************************************************/
      /**************************************************************/
      /**************************************************************/
      /**************************************************************/
      
      if(err)
	{
	  PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run);
	  Exit("\n");
	}
      
      
      RATES_Update_Cur_Bl(tree);
      
      /* Optional here but useful for monitoring ESS for likelihoods */
      tree->c_lnL        = Dnorm_Multi_Given_InvCov_Det(tree->rates->u_cur_l,tree->rates->u_ml_l,tree->rates->invcov,tree->rates->covdet,2*tree->n_otu-3,YES);
      tree->rates->c_lnL = RATES_Lk_Rates(tree);
      
      tree->mcmc->run++;
      MCMC_Print_Param(tree->mcmc,tree);
      
      /* !!!!!!!!!!!!!!!!!! */
      if(!(tree->mcmc->run%tree->mcmc->norm_freq)) 
	{
	  RATES_Normalise_Rates(tree);
	}
    }
}

/*********************************************************/

void RATES_Posterior_Rates(t_tree *tree)
{
/*   RATES_Posterior_Rates_Pre(tree->n_root,tree->n_root->v[0],tree); */
/*   RATES_Posterior_Rates_Pre(tree->n_root,tree->n_root->v[1],tree); */

  int node_num;
  int i;
  
  For(i,2*tree->n_otu-2)
    {
      node_num = Rand_Int(0,2*tree->n_otu-3);
      RATES_Posterior_Rates_Pre(tree->noeud[node_num]->anc,tree->noeud[node_num],tree);
    }
}

/*********************************************************/

void RATES_Posterior_Times(t_tree *tree)
{
  int node_num;

  int i;
  
  For(i,2*tree->n_otu-1)
    {
      node_num = Rand_Int(tree->n_otu,2*tree->n_otu-2);
      
      if(tree->noeud[node_num] != tree->n_root)
	{
	  if(tree->rates->met_within_gibbs)
	    MCMC_Times_Pre(tree->noeud[node_num]->anc,tree->noeud[node_num],1,tree);
	  else
	    RATES_Posterior_Times_Pre(tree->noeud[node_num]->anc,tree->noeud[node_num],tree);
	}
      else
	{
	  if(tree->rates->met_within_gibbs)
	    MCMC_Times_Pre(NULL,tree->noeud[node_num],1,tree);
	  else
	    RATES_Posterior_Time_Root(tree);
	}
    }
}

/*********************************************************/

void RATES_Posterior_Rates_Pre(t_node *a, t_node *d, t_tree *tree)
{
  phydbl like_mean, like_var;
  phydbl prior_mean, prior_var;
  phydbl post_mean, post_var, post_sd;
  phydbl dt,el,vl,ra,rd,min_r,max_r,cr,nu,cel,cvl,cel_tmp,cvl_tmp,cer,cvr;
  phydbl cur_l,new_l;
  int dim;
  short int *is_1;
  phydbl *cond_mu, *cond_cov; 
  phydbl l_opp; /* length of the branch connected to the root, opposite to the one connected to d */
  t_edge *b;
  int i;
  t_node *v2,*v3;
  phydbl T0,T1,T2,T3;
  phydbl U0,U1,U2,U3;
  phydbl V1,V2,V3;
  phydbl r_min,r_max;
  phydbl l_min,l_max;
  int err;
  
/*   if(a == tree->n_root) { tree->rates->nd_r[d->num] = 1.0; return; } */

  dim = 2*tree->n_otu-3;

/*   is_1     = (short int *)mCalloc(2*tree->n_otu-3,sizeof(short int)); */
/*   cond_mu  = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl)); */
/*   cond_cov = (phydbl *)mCalloc((2*tree->n_otu-3)*(2*tree->n_otu-3),sizeof(phydbl)); */

  is_1     = tree->rates->_2n_vect5;
  cond_mu  = tree->rates->_2n_vect1;
  cond_cov = tree->rates->_2n2n_vect1;

  b = NULL;
  if(a == tree->n_root) b = tree->e_root;
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
  
  if(T1 - T0 < tree->rates->min_dt)
    {
      PhyML_Printf("\n. T1 = %f T0 = %f",T1,T0);
      PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run);
      Exit("\n");
    }

  if(!d->tax)
    {
      if(T2 - T1 < tree->rates->min_dt)
	{
	  PhyML_Printf("\n. T2 = %f T1 = %f",T2,T1);
	  PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run);
	  Exit("\n");
	}
      
      if(T3 - T1 < tree->rates->min_dt)
	{
	  PhyML_Printf("\n. T3 = %f T1 = %f",T3,T1);
	  PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run);
	  Exit("\n");
	}
    }

  V1 = tree->rates->nu * (T1 - T0);
  V2 = tree->rates->nu * (T2 - T1);
  V3 = tree->rates->nu * (T3 - T1);

  dt     = tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num];
  el     = tree->rates->u_ml_l[b->num];
  vl     = tree->rates->cov[b->num*dim+b->num];
  ra     = tree->rates->nd_r[a->num];
  rd     = tree->rates->nd_r[d->num];
  min_r  = tree->rates->min_rate;
  max_r  = tree->rates->max_rate;
  cr     = tree->rates->clock_r;
  nu     = tree->rates->nu;
  cel    = el;   /* conditional branch length expectation */
  cvl    = vl;   /* conditional branch length variance */
  cer    = -1.0; /* conditional rate expectation */
  cvr    = -1.0; /* conditional rate variance */
  cur_l  = tree->rates->cur_l[d->num];

  l_opp  = -1.0;
  if(a == tree->n_root)
    {
      l_opp =
	(d == a->v[0])?
	(tree->rates->cur_l[a->v[1]->num]):
	(tree->rates->cur_l[a->v[0]->num]);
    }
  
  r_min = tree->rates->min_rate;
  r_max = tree->rates->max_rate;

  l_min = r_min * (T1 - T0) * tree->rates->clock_r;
  l_max = r_max * (T1 - T0) * tree->rates->clock_r;

  if(a == tree->n_root)
    {
      l_min += l_opp;
      l_max += l_opp;
    }

  /* Likelihood */

  cel=0.0;
  For(i,dim) if(i != b->num) cel += tree->rates->reg_coeff[b->num*dim+i] * (tree->rates->u_cur_l[i] - tree->rates->u_ml_l[i]);
  cel += tree->rates->u_ml_l[b->num];

  cvl = tree->rates->cond_var[b->num];

  if(a == tree->n_root)
    {
      if(!Norm_Trunc_Mean_Sd(cel,sqrt(cvl),l_opp,BL_MAX,&cel_tmp,&cvl_tmp)) return;
      cvl_tmp *= cvl_tmp;
      cvl = cvl_tmp;
      cel = cel_tmp;
    }

  if(isnan(cvl) || isnan(cel)) 
    {
      PhyML_Printf("\n. Warning: invalid expected and/or std. dev. values. Skipping this step.\n"); 
      return;
    }

  like_mean = cel;
  like_var  = cvl;
 
  /* Prior */
  if(!d->tax)
    {
      cvr = 1./((1./V1)+(1./V2)+(1./V3));
      cer = cvr*(U0/V1 + U2/V2 + U3/V3);
    }
  else
    {
      cvr = V1;
      cer = U0;
    }

  prior_mean = log(cer) + log(dt) + log(cr);
  prior_mean = exp(prior_mean);
  if(a == tree->n_root) prior_mean += l_opp;

  prior_var = log(cvr) + 2.*log(dt) + 2.*log(cr) ;
  prior_var = exp(prior_var);
  
  /* Posterior */
  post_mean = (prior_mean/prior_var + like_mean/like_var)/(1./prior_var + 1./like_var);

  post_var  = 1./(1./prior_var + 1./like_var);
  post_sd   = sqrt(post_var);

  new_l = Rnorm_Trunc(post_mean,post_sd,l_min,l_max,&err);
  
  if(err)
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run);
      PhyML_Printf("\n. cvr=%G dt=%G cr=%G",cvr,dt,cr);
      PhyML_Printf("\n. prior_mean = %G prior_var = %G",prior_mean,prior_var);
      PhyML_Printf("\n. like_mean = %G like_var = %G",like_mean,like_var);
      PhyML_Printf("\n. post_mean = %G post_var = %G",post_mean,post_var);
      PhyML_Printf("\n. clock_r = %f",tree->rates->clock_r);
      PhyML_Printf("\n. T0=%f T1=%f T2=%f T3=%f",T0,T1,T2,T3);
      PhyML_Printf("\n. U0=%f U1=%f U2=%f U3=%f",U0,U1,U2,U3);
      PhyML_Printf("\n. l_min=%f l_max=%f",l_min,l_max);
      PhyML_Printf("\n. Setting t_edge length to %f",cur_l);
      new_l = cur_l;
    }

  if(a == tree->n_root) rd = (new_l-l_opp)/(dt*cr);
  else                  rd = new_l / (dt*cr);
  
  if(isnan(rd))
    {
      PhyML_Printf("\n. rd=%G new_l=%G prior_mean=%G prior_var=%G like_mean=%G like_var=%G post_mean=%G post_var=%G dt=%f ra=%G cr=%G root=%d",
	     rd,new_l,
	     prior_mean,prior_var,
	     like_mean,like_var,
	     post_mean,post_var,
	     dt,ra,cr,(a==tree->n_root)?(1):(0));
      PhyML_Printf("\n. t(root)=%f t0=%f t1=%f",
	     tree->rates->nd_t[tree->n_root->num],
	     tree->rates->nd_t[tree->n_root->v[0]->num],
	     tree->rates->nd_t[tree->n_root->v[1]->num]);
      PhyML_Printf("\n. cvl_tmp = %f",cvl_tmp);
      PhyML_Printf("\n. Run = %d",tree->mcmc->run);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }
  
  if(rd < r_min) 
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run);
      PhyML_Printf("\n. cvr=%G dt=%G cr=%G",cvr,dt,cr);
      PhyML_Printf("\n. prior_mean = %G prior_var = %G",prior_mean,prior_var);
      PhyML_Printf("\n. like_mean = %G like_var = %G",like_mean,like_var);
      PhyML_Printf("\n. post_mean = %G post_var = %G",post_mean,post_var);
      PhyML_Printf("\n. clock_r = %f",tree->rates->clock_r);
      PhyML_Printf("\n. T0=%f T1=%f T2=%f T3=%f",T0,T1,T2,T3);
      PhyML_Printf("\n. U0=%f U1=%f U2=%f U3=%f",U0,U1,U2,U3);
      PhyML_Printf("\n. l_min=%f l_max=%f",l_min,l_max);
      PhyML_Printf("\n. Setting rd to %f",tree->rates->min_rate);
      rd = r_min;
    }

  if(rd > r_max) 
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run);
      PhyML_Printf("\n. cvr=%G dt=%G cr=%G",cvr,dt,cr);
      PhyML_Printf("\n. prior_mean = %G prior_var = %G",prior_mean,prior_var);
      PhyML_Printf("\n. like_mean = %G like_var = %G",like_mean,like_var);
      PhyML_Printf("\n. post_mean = %G post_var = %G",post_mean,post_var);
      PhyML_Printf("\n. clock_r = %f",tree->rates->clock_r);
      PhyML_Printf("\n. T0=%f T1=%f T2=%f T3=%f",T0,T1,T2,T3);
      PhyML_Printf("\n. U0=%f U1=%f U2=%f U3=%f",U0,U1,U2,U3);
      PhyML_Printf("\n. l_min=%f l_max=%f",l_min,l_max);
      PhyML_Printf("\n. Setting rd to %f",tree->rates->max_rate);
      rd = r_max;
    }

  tree->rates->nd_r[d->num] = rd;
  
  RATES_Update_Cur_Bl(tree);

  /* Optional here but useful for monitoring ESS for likelihoods */
  /*   tree->c_lnL        = Dnorm_Multi_Given_InvCov_Det(tree->rates->u_cur_l,tree->rates->u_ml_l,tree->rates->invcov,tree->rates->covdet,2*tree->n_otu-3,YES); */
  /*   tree->rates->c_lnL = RATES_Lk_Rates(tree); */

  /*   RATES_Check_Lk_Rates(tree,&err); */
  
  /*   if(err) */
  /*     { */
  /*       PhyML_Printf("\n. Small lk detected"); */
  /*       PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
  /*       Exit("\n"); */
  /*     } */


  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc,tree);  
  if(!(tree->mcmc->run%tree->mcmc->norm_freq))
    {
      RATES_Normalise_Rates(tree);
    }
}

/*********************************************************/

void RATES_Posterior_Times_Pre(t_node *a, t_node *d, t_tree *tree)
{
  /*
           T0 (a)
            |
            | l1,U1,b1
            |
           T1 (d)
           / \
 l2,u2,b2 /   \ L3,u3,b3
         /     \
       t2       t3
      (v2)     (v3)

   */

  phydbl l1,l2,l3;
  phydbl El1,El2,El3;
  phydbl u0,u1,u2,u3;
  phydbl t0,t1,t2,t3;
  phydbl t_min, t_max;
  phydbl t_max_12,t_max_13;
  phydbl bl_min, bl_max;
  phydbl t1_new;
  phydbl X,Y;
  phydbl cr;
  int    i,j;
  phydbl *mu, *cov;
  phydbl *cond_mu, *cond_cov;
  short int *is_1;
  phydbl sig11, sig1X, sig1Y, sig22, sig2X, sig2Y, sigXX, sigXY, sigYY;
  phydbl cov11,cov12,cov13,cov22,cov23,cov33;
  int dim;
  t_edge *b1, *b2, *b3;
  t_node *v2,*v3;
  phydbl *l2XY;
  phydbl l_opp;
  t_edge *buff_b;
  t_node *buff_n;
  int err;
  int num_1, num_2, num_3;

  dim = 2*tree->n_otu-3;

  if(d->tax) return;

  if(fabs(tree->rates->t_prior_min[d->num] - tree->rates->t_prior_max[d->num]) < 1.E-10) return;

  l2XY     = tree->rates->_2n_vect2;
  mu       = tree->rates->_2n_vect3;
  cov      = tree->rates->_2n2n_vect1;
  cond_mu  = tree->rates->_2n_vect1;
  cond_cov = tree->rates->_2n2n_vect2;
  is_1     = tree->rates->_2n_vect5;

  b1 = NULL;
  if(a == tree->n_root) b1 = tree->e_root;
  else For(i,3) if(d->v[i] == a) { b1 = d->b[i]; break; }

  b2 = b3 = NULL;
  v2 = v3 = NULL;
  For(i,3)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      {
	if(!v2) { v2 = d->v[i]; b2 = d->b[i]; }
	else    { v3 = d->v[i]; b3 = d->b[i]; }
      }

  t2  = tree->rates->nd_t[v2->num];
  t3  = tree->rates->nd_t[v3->num];

  buff_n = NULL;
  buff_b = NULL;
  if(t3 > t2)
    {
      buff_n = v2;
      v2     = v3;
      v3     = buff_n;

      buff_b = b2;
      b2     = b3;
      b3     = buff_b;
    }  
    
  t0 = tree->rates->nd_t[a->num];
  t1 = tree->rates->nd_t[d->num];
  t2 = tree->rates->nd_t[v2->num];
  t3 = tree->rates->nd_t[v3->num];
  u0 = tree->rates->nd_r[a->num];
  u1 = tree->rates->nd_r[d->num];
  u2 = tree->rates->nd_r[v2->num];
  u3 = tree->rates->nd_r[v3->num];
  l1 = tree->rates->cur_l[d->num];
  l2 = tree->rates->cur_l[v2->num];
  l3 = tree->rates->cur_l[v3->num];
  cr = tree->rates->clock_r;

/*   For(i,dim) is_1[i] = 0; */
/*   is_1[b1->num] = 1; */
/*   is_1[b2->num] = 1; */
/*   is_1[b3->num] = 1; */

/*   For(i,dim)     cond_mu[i]  = 0.0; */
/*   For(i,dim*dim) cond_cov[i] = 0.0; */
/*   Normal_Conditional(tree->rates->u_ml_l,tree->rates->cov,tree->rates->u_cur_l,dim,is_1,3,cond_mu,cond_cov); */
  
/*   El1 = cond_mu[b1->num]; */
/*   El2 = cond_mu[b2->num]; */
/*   El3 = cond_mu[b3->num]; */

/*   cov11 = cond_cov[b1->num*dim+b1->num]; */
/*   cov12 = cond_cov[b1->num*dim+b2->num]; */
/*   cov13 = cond_cov[b1->num*dim+b3->num]; */
/*   cov23 = cond_cov[b2->num*dim+b3->num]; */
/*   cov22 = cond_cov[b2->num*dim+b2->num]; */
/*   cov33 = cond_cov[b3->num*dim+b3->num]; */


/*   El1 = tree->rates->u_ml_l[b1->num]; */
/*   El2 = tree->rates->u_ml_l[b2->num]; */
/*   El3 = tree->rates->u_ml_l[b3->num]; */

/*   cov11 = tree->rates->cov[b1->num*dim+b1->num]; */
/*   cov12 = tree->rates->cov[b1->num*dim+b2->num]; */
/*   cov13 = tree->rates->cov[b1->num*dim+b3->num]; */
/*   cov23 = tree->rates->cov[b2->num*dim+b3->num]; */
/*   cov22 = tree->rates->cov[b2->num*dim+b2->num]; */
/*   cov33 = tree->rates->cov[b3->num*dim+b3->num]; */

/*   PhyML_Printf("\n- El1=%10f El2=%10f El3=%10f",El1,El2,El3); */
/*   PhyML_Printf("\n- cov11=%10f cov12=%10f cov13=%10f cov23=%10f cov22=%10f cov33=%10f", cov11,cov12,cov13,cov23,cov22,cov33); */

  num_1 = num_2 = num_3 = -1;
  if(b1->num < b2->num && b2->num < b3->num) { num_1 = 0; num_2 = 1; num_3 = 2; }
  if(b1->num < b3->num && b3->num < b2->num) { num_1 = 0; num_2 = 2; num_3 = 1; }
  if(b2->num < b1->num && b1->num < b3->num) { num_1 = 1; num_2 = 0; num_3 = 2; }
  if(b2->num < b3->num && b3->num < b1->num) { num_1 = 2; num_2 = 0; num_3 = 1; }
  if(b3->num < b2->num && b2->num < b1->num) { num_1 = 2; num_2 = 1; num_3 = 0; }
  if(b3->num < b1->num && b1->num < b2->num) { num_1 = 1; num_2 = 2; num_3 = 0; }

  cov11 = tree->rates->trip_cond_cov[d->num * 9 + num_1 * 3 + num_1];
  cov12 = tree->rates->trip_cond_cov[d->num * 9 + num_1 * 3 + num_2];
  cov13 = tree->rates->trip_cond_cov[d->num * 9 + num_1 * 3 + num_3];
  cov23 = tree->rates->trip_cond_cov[d->num * 9 + num_2 * 3 + num_3];
  cov22 = tree->rates->trip_cond_cov[d->num * 9 + num_2 * 3 + num_2];
  cov33 = tree->rates->trip_cond_cov[d->num * 9 + num_3 * 3 + num_3];

  El1=0.0;
  For(i,dim)
    if(i != b1->num && i != b2->num && i != b3->num)
      El1 += tree->rates->trip_reg_coeff[d->num * (6*tree->n_otu-9) + num_1 * dim +i] * (tree->rates->u_cur_l[i] - tree->rates->u_ml_l[i]);
  El1 += tree->rates->u_ml_l[b1->num];

  El2=0.0;
  For(i,dim)
    if(i != b1->num && i != b2->num && i != b3->num)
      El2 += tree->rates->trip_reg_coeff[d->num * (6*tree->n_otu-9) + num_2 * dim +i] * (tree->rates->u_cur_l[i] - tree->rates->u_ml_l[i]);
  El2 += tree->rates->u_ml_l[b2->num];

  El3=0.0;
  For(i,dim)
    if(i != b1->num && i != b2->num && i != b3->num)
      El3 += tree->rates->trip_reg_coeff[d->num * (6*tree->n_otu-9) + num_3 * dim +i] * (tree->rates->u_cur_l[i] - tree->rates->u_ml_l[i]);
  El3 += tree->rates->u_ml_l[b3->num];

/*   PhyML_Printf("\n+ El1=%10f El2=%10f El3=%10f",El1,El2,El3); */
/*   PhyML_Printf("\n+ cov11=%10f cov12=%10f cov13=%10f cov23=%10f cov22=%10f cov33=%10f", cov11,cov12,cov13,cov23,cov22,cov33); */

  if(El1 < BL_MIN) El1 = BL_MIN;
  if(El2 < BL_MIN) El2 = BL_MIN;
  if(El3 < BL_MIN) El3 = BL_MIN;

  if(El1 > BL_MAX) El1 = BL_MAX;
  if(El2 > BL_MAX) El2 = BL_MAX;
  if(El3 > BL_MAX) El3 = BL_MAX;

  t1_new = +1;

  /* Ui's are relative substitution rates */
/*   t_min =     t0 + (1./tree->rates->nu)*pow((u1-u0)/tree->rates->z_max,2); */
/*   t_max = MIN(t3 - (1./tree->rates->nu)*pow((u1-u3)/tree->rates->z_max,2), */
/* 	      t2 - (1./tree->rates->nu)*pow((u1-u2)/tree->rates->z_max,2)); */

/*   t_min = t0; */
/*   t_max = MIN(t3,t2); */
  
  t_min    = RATES_Find_Min_Dt_Bisec(u1,u0,t0,MIN(t2,t3),tree->rates->nu,tree->rates->p_max,(u1 < u0)?YES:NO);
  t_max_12 = RATES_Find_Max_Dt_Bisec(u2,u1,t0,t2,tree->rates->nu,tree->rates->p_max,(u2 < u1)?YES:NO);
  t_max_13 = RATES_Find_Max_Dt_Bisec(u3,u1,t0,t3,tree->rates->nu,tree->rates->p_max,(u3 < u1)?YES:NO);
  t_max = MIN(t_max_12,t_max_13);


  t_min = MAX(t_min,tree->rates->t_prior_min[d->num]);
  t_max = MIN(t_max,tree->rates->t_prior_max[d->num]);
    
  if(fabs(t_max - t_min) < 1.E-10) return;   
  if(t_max < t_min) 
    {
      PhyML_Printf("\n%3d t0=%10f t2=%10f t3=%10f t_max=%10f t_min=%10f u0=%10f u1=%10f u2=%10f u3=%10f",
		   d->num,t0,t2,t3,t_max,t_min,u0,u1,u2,u3);      
      return;
    }

  tree->rates->t_prior[d->num] = Uni()*(t_max - t_min) + t_min;

  bl_min = bl_max = -1.0;


  /* Ui's are now absolute substitution rates */
  u0 *= cr;
  u1 *= cr;
  u2 *= cr;
  u3 *= cr;

  l_opp = -1.;
  if(a == tree->n_root)
    {
      if(d == tree->n_root->v[0])
	l_opp = tree->rates->cur_l[tree->n_root->v[1]->num];
      else
	l_opp = tree->rates->cur_l[tree->n_root->v[0]->num];

      l1 =
	tree->rates->cur_l[tree->n_root->v[0]->num] +
	tree->rates->cur_l[tree->n_root->v[1]->num] ;
    }
   
  /* Not used but here for the sake of clarity (we want X=0 and Y=0) */
  X = l1/u1 + t0 + l2/u2 - t2;
  Y = l1/u1 + t0 + l3/u3 - t3;
  
  l2XY[0] = 0.0; /* does not matter */
  l2XY[1] = 0.0; /* constraint 1 (X=0) */
  l2XY[2] = 0.0; /* constraint 2 (Y=0) */
   
/*   if(cov11 > cov22) */
    mu[0] = El1;
/*   else */
/*     mu[0] = El2; */

  if(a == tree->n_root)
    {
      mu[1] = (El1 - l_opp)/u1 + t0 + El2/u2 - t2;
      mu[2] = (El1 - l_opp)/u1 + t0 + El3/u3 - t3;
    }
  else
    {
      mu[1] = El1/u1 + t0 + El2/u2 - t2;
      mu[2] = El1/u1 + t0 + El3/u3 - t3;
    }

  sig11 = cov11;
  sig22 = cov22;
  sig1X = (1./u1) * cov11 + (1./u2) * cov12;
  sig1Y = (1./u1) * cov11 + (1./u3) * cov13;
  sig2X = (1./u1) * cov12 + (1./u2) * cov22;
  sig2Y = (1./u1) * cov12 + (1./u3) * cov23;
  sigXX = 1./(u1*u1) * cov11 + 1./(u2*u2) * cov22 + 2./(u1*u2) * cov12;
  sigYY = 1./(u1*u1) * cov11 + 1./(u3*u3) * cov33 + 2./(u1*u3) * cov13;
  sigXY = 1./(u1*u1) * cov11 + 1./(u1*u2) * cov12 + 1./(u1*u3) * cov13 + 1./(u2*u3) * cov23;

/*   if(cov11 > cov22) */
/*     { */
      cov[0*3+0] = sig11; cov[0*3+1] = sig1X; cov[0*3+2] = sig1Y;
      cov[1*3+0] = sig1X; cov[1*3+1] = sigXX; cov[1*3+2] = sigXY;
      cov[2*3+0] = sig1Y; cov[2*3+1] = sigXY; cov[2*3+2] = sigYY;
/*     } */
/*   else */
/*     { */
/*       cov[0*3+0] = sig22; cov[0*3+1] = sig2X; cov[0*3+2] = sig2Y; */
/*       cov[1*3+0] = sig2X; cov[1*3+1] = sigXX; cov[1*3+2] = sigXY; */
/*       cov[2*3+0] = sig2Y; cov[2*3+1] = sigXY; cov[2*3+2] = sigYY; */
/*     } */
  
/*   PhyML_Printf("\n"); */
/*   For(i,3) */
/*     { */
/*       For(j,3) */
/* 	{ */
/* 	  PhyML_Printf("%f ",cov[i*3+j]); */
/* 	} */
/*       PhyML_Printf("\n"); */
/*     } */
/*   PhyML_Printf("\n"); */

  is_1[0] = is_1[1] = is_1[2] = 0;
  is_1[0] = 1;
  Normal_Conditional(mu,cov,l2XY,3,is_1,1,cond_mu,cond_cov);

  if(cond_cov[0*3+0] < 0.0)
    {
      PhyML_Printf("\n. Conditional mean=%G var=%G",cond_mu[0],cond_cov[0*3+0]);
      PhyML_Printf("\n. t0=%G t1=%f t2=%f l1=%G l2=%G El1=%G El2=%G Nu=%G",t0,t1,t2,l1,l2,El1,El2,tree->rates->nu);
      PhyML_Printf("\n. COV11=%f COV12=%f COV13=%f COV22=%f COV23=%f COV33=%f",cov11,cov12,cov13,cov22,cov23,cov33);
      PhyML_Printf("\n");
      For(i,3)
	{
	  PhyML_Printf(". mu%d=%12lf\t",i,mu[i]);
	  For(j,3)
	    {
	      PhyML_Printf("%12lf ",cov[i*3+j]);
	    }
	  PhyML_Printf("\n");
	}
      
      cond_cov[0*3+0] = 1.E-10;
    }

/*   if(cov11 > cov22) */
/*     { */
      bl_min = (t_min - t0) * u1;
      bl_max = (t_max - t0) * u1;
      
      if(a == tree->n_root)
	{
	  l1 = Rnorm_Trunc(cond_mu[0],sqrt(cond_cov[0*3+0]),l_opp+bl_min,l_opp+bl_max,&err);
	  t1_new = (l1-l_opp)/u1 + t0;
	}
      else
	{
	  l1 = Rnorm_Trunc(cond_mu[0],sqrt(cond_cov[0*3+0]),bl_min,bl_max,&err);
	  t1_new = l1/u1 + t0;
	}

      if(err)
	{
	  PhyML_Printf("\n");
	  PhyML_Printf("\n. Root ? %s",(tree->n_root==a)?("yes"):("no")); 
	  PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run);
	  PhyML_Printf("\n. t0=%f t1=%f t2=%f t3=%f",t0,t1,t2,t3);
	  PhyML_Printf("\n. t_min=%f t_max=%f",t_min,t_max);
	  PhyML_Printf("\n. bl_min=%f bl_max=%f",bl_min,bl_max);
	  PhyML_Printf("\n. cond_mu[0]=%f cond_cov[0]=%f",cond_mu[0],sqrt(cond_cov[0*3+0]));
	  PhyML_Printf("\n. El1=%f El2=%f El3=%f",El1,El2,El3);
	  PhyML_Printf("\n. l1=%f l2=%f l3=%f",l1,l2,l3);
	  PhyML_Printf("\n. u0=%f u1=%f u2=%f u3=%f",u0,u1,u2,u3);
	  PhyML_Printf("\n. COV11=%f COV22=%f",cov11,cov22);
	  PhyML_Printf("\n. Clock rate = %f",tree->rates->clock_r);
	  PhyML_Printf("\n. Setting t1_new to %f",t1);
	  t1_new = t1;
/* 	  Exit("\n"); */
	}
/*     } */
/*   else */
/*     { */
/*       bl_min = (t2 - t_max) * u2; */
/*       bl_max = (t2 - t_min) * u2; */

/*       l2 = Rnorm_Trunc(cond_mu[0],sqrt(cond_cov[0*3+0]),bl_min,bl_max,&err); */
/*       t1_new = -l2/u2 + t2; */
      
/*       if(err) */
/* 	{ */
/* 	  PhyML_Printf("\n"); */
/* 	  PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run); */
/* 	  PhyML_Printf("\n. t0=%f t1=%f t2=%f t3=%f",t0,t1,t2,t3); */
/* 	  PhyML_Printf("\n. t_min=%f t_max=%f",t_min,t_max); */
/* 	  PhyML_Printf("\n. bl_min=%f bl_max=%f",bl_min,bl_max); */
/* 	  PhyML_Printf("\n. cond_mu[0]=%f cond_cov[0]=%f",cond_mu[0],sqrt(cond_cov[0*3+0])); */
/* 	  PhyML_Printf("\n. El1=%f El2=%f El3=%f",El1,El2,El3); */
/* 	  PhyML_Printf("\n. l1=%f l2=%f l3=%f",l1,l2,l3); */
/* 	  PhyML_Printf("\n. u0=%f u1=%f u2=%f u3=%f",u0,u1,u2,u3); */
/* 	  PhyML_Printf("\n. COV11=%f COV22=%f",cov11,cov22); */
/* 	  PhyML_Printf("\n. Clock rate = %f",tree->rates->clock_r); */
/* 	  PhyML_Printf("\n. Setting t1_new to %f",t1); */
/* 	  t1_new = t1; */
/* /\* 	  Exit("\n"); *\/ */
/* 	} */
/*     } */
  
  if(t1_new < t0)
    {
      t1_new = t0+1.E-4;
      PhyML_Printf("\n");
      PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
      PhyML_Printf("\n. t0 = %f t1_new = %f t1 = %f",t0,t1_new,t1);
      PhyML_Printf("\n. t_min=%f t_max=%f",t_min,t_max);
      PhyML_Printf("\n. l1 = %f",l1);
      PhyML_Printf("\n. bl_min = %f bl_max = %f",bl_min,bl_max);
      PhyML_Printf("\n. (t1-t0)=%f (t2-t1)=%f",t1-t0,t2-t1);
      PhyML_Printf("\n. l1 = %f l2 = %f cov11=%f cov22=%f cov33=%f",l1,l2,cov11,cov22,cov33);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
/*       Exit("\n"); */
    }
  if(t1_new > MIN(t2,t3))
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
      PhyML_Printf("\n. t0 = %f t1_new = %f t1 = %f t2 = %f t3 = %f MIN(t2,t3)=%f",t0,t1_new,t1,t2,t3,MIN(t2,t3));
      PhyML_Printf("\n. t_min=%f t_max=%f",t_min,t_max);
      PhyML_Printf("\n. l2 = %f",l2);
      PhyML_Printf("\n. bl_min = %f bl_max = %f",bl_min,bl_max);
      PhyML_Printf("\n. (t1-t0)=%f (t2-t1)=%f",t1-t0,t2-t1);
      PhyML_Printf("\n. l1 = %f l2 = %f cov11=%f cov22=%f cov33=%f",l1,l2,cov11,cov22,cov33);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
/*       Exit("\n"); */
    }

  if(isnan(t1_new))
    {
      PhyML_Printf("\n. run=%d",tree->mcmc->run);
      PhyML_Printf("\n. mean=%G var=%G",cond_mu[0],cond_cov[0*3+0]);
      PhyML_Printf("\n. t1=%f l1=%G u1=%G t0=%G",t1,l1,u1,t0);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
/*       Exit("\n"); */
    }
  

  if(t1_new - t0 < tree->rates->min_dt)
    {
      PhyML_Printf("\n. t1_new = %f t0 = %f t2=%f t3=%f",t1_new,t0,t2,t3);
      PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run);
      Exit("\n");
    }

  if(t2 - t1_new  < tree->rates->min_dt)
    {
      PhyML_Printf("\n. t1_new = %f t0 = %f t2=%f t3=%f",t1_new,t0,t2,t3);
      PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run);
      Exit("\n");
    }
  if(t3 - t1_new  < tree->rates->min_dt)
    {
      PhyML_Printf("\n. t1_new = %f t0 = %f t2=%f t3=%f",t1_new,t0,t2,t3);
      PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run);
      Exit("\n");
    }

/*   tree->rates->nd_t[d->num] = tree->rates->t_prior[d->num]; */
  tree->rates->nd_t[d->num] = t1_new;
  RATES_Update_Cur_Bl(tree);
  
/*   tree->c_lnL        = Dnorm_Multi_Given_InvCov_Det(tree->rates->u_cur_l,tree->rates->u_ml_l,tree->rates->invcov,tree->rates->covdet,2*tree->n_otu-3,YES); */
/*   tree->rates->c_lnL = RATES_Lk_Rates(tree); */


  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc,tree);
  if(!(tree->mcmc->run%tree->mcmc->norm_freq)) 
    {
      RATES_Normalise_Rates(tree);
    }
}

/*********************************************************/

void RATES_Posterior_Time_Root(t_tree *tree)
{
  /* 
        Root, T0
        / \
  l1   /   \ l2
      /     \
     v1,t1   v2,t2

  */

  phydbl t1,t2;
  phydbl u0,u1,u2;
  phydbl cel,cvl;
  int i,dim;
  t_edge *b;
  t_node *root;
  phydbl t0,t0_min, t0_max;
  phydbl t_max_01, t_max_02;
  int err;
  phydbl cr;
  phydbl bl_min, bl_max;
  phydbl new_L;

  dim = 2*tree->n_otu-3;

  b = tree->e_root;
  root = tree->n_root;

  t1 = tree->rates->nd_t[root->v[0]->num];
  t2 = tree->rates->nd_t[root->v[1]->num];
  u1 = tree->rates->nd_r[root->v[0]->num];
  u2 = tree->rates->nd_r[root->v[1]->num];
  u0 = 1.0;
  cr = tree->rates->clock_r;

  t0_min = -MDBL_MAX;

/*   t0_max = MIN(t1 - (1./tree->rates->nu)*pow((u0-u1)/tree->rates->z_max,2), */
/* 	       t2 - (1./tree->rates->nu)*pow((u0-u2)/tree->rates->z_max,2)); */
/*   t0_max = MIN(t1,t2); */
  
  
/*   if(u1 > u0) t_max_01 = t1 - (1./tree->rates->nu)*pow((u1-u0)/(PointNormal(tree->rates->p_max*Pnorm(.0,.0,1.))),2); */
/*   else        t_max_01 = t1 - (1./tree->rates->nu)*pow((u1-u0)/(PointNormal(tree->rates->p_max*(1.-Pnorm(.0,.0,1.)))),2); */

/*   if(u2 > u0) t_max_02 = t2 - (1./tree->rates->nu)*pow((u2-u0)/(PointNormal(tree->rates->p_max*Pnorm(.0,.0,1.))),2); */
/*   else        t_max_02 = t2 - (1./tree->rates->nu)*pow((u2-u0)/(PointNormal(tree->rates->p_max*(1.-Pnorm(.0,.0,1.)))),2); */
  

  t_max_01 = RATES_Find_Max_Dt_Bisec(u1,u0,tree->rates->t_prior_min[root->num],t1,tree->rates->nu,tree->rates->p_max,(u1 < u0)?YES:NO);
  t_max_02 = RATES_Find_Max_Dt_Bisec(u2,u0,tree->rates->t_prior_min[root->num],t2,tree->rates->nu,tree->rates->p_max,(u2 < u0)?YES:NO);

  t0_max = MIN(t_max_01,t_max_02);

  t0_min = MAX(t0_min,tree->rates->t_prior_min[root->num]);
  t0_max = MIN(t0_max,tree->rates->t_prior_max[root->num]);


  u0 *= cr;
  u1 *= cr;
  u2 *= cr;

/*   if(fabs(t0_max - t0_min) < 1.E-10) return; */
  if(fabs(tree->rates->t_prior_min[root->num] - tree->rates->t_prior_max[root->num]) < 1.E-10) return;

  cel=0.0;
  For(i,dim) if(i != b->num) cel += tree->rates->reg_coeff[b->num*dim+i] * (tree->rates->u_cur_l[i] - tree->rates->u_ml_l[i]);
  cel += tree->rates->u_ml_l[b->num];

  cvl = tree->rates->cond_var[b->num];

/*   cel = tree->rates->u_ml_l[b->num]; */
/*   cvl = tree->rates->cov[b->num*dim+b->num]; */

/*   Et0 = (u1*t1 + u2*t2 - cel)/(u1+u2); */
/*   Vt0 = cvl/pow(u1+u2,2); */

/*   t0  = Rnorm_Trunc(Et0,sqrt(Vt0),t0_min,t0_max,&err); */

  
  bl_min = u1 * (t1 - t0_max) + u2 * (t2 - t0_max); 
  bl_max = u1 * (t1 - t0_min) + u2 * (t2 - t0_min);

  if(bl_min > bl_max) return;

  new_L = Rnorm_Trunc(cel,sqrt(cvl),bl_min,bl_max,&err);

  t0 = (u1*t1 + u2*t2 - new_L)/(u1+u2);

  if(t0 > t1 || t0 > t2)
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. Run = %d",tree->mcmc->run);
      PhyML_Printf("\n. t0=%f t1=%f t2=%f",t0,t1,t2);
      PhyML_Printf("\n. t0_min=%f t0_max=%f",t0_min,t0_max);
      PhyML_Printf("\n. new_L=%f cel=%f",new_L,cel);
      PhyML_Printf("\n. u0=%f u1=%f u2=%f",u0/cr,u1/cr,u2/cr);
      PhyML_Printf("\n. Nu = %f Clock = %f",tree->rates->nu,tree->rates->clock_r);
      PhyML_Printf("\n. Setting t0 to %f",tree->rates->nd_t[root->num]);
/*       t0 = tree->rates->nd_t[root->num];       */
      return;
    }

  if(t0 < t0_min || t0 > t0_max)
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. Run = %d",tree->mcmc->run);
      PhyML_Printf("\n. t0=%f t1=%f t2=%f",t0,t1,t2);
      PhyML_Printf("\n. t0_min=%f t0_max=%f",t0_min,t0_max);
      PhyML_Printf("\n. u0=%f u1=%f u2=%f",u0/cr,u1/cr,u2/cr);
      PhyML_Printf("\n. Nu = %f Clock = %f",tree->rates->nu,tree->rates->clock_r);
      PhyML_Printf("\n. Setting t0 to %f",tree->rates->nd_t[root->num]);
      t0 = tree->rates->nd_t[root->num];      
    }

  tree->rates->nd_t[root->num] = t0;

  RATES_Update_Cur_Bl(tree);
  
/*   tree->c_lnL        = Dnorm_Multi_Given_InvCov_Det(tree->rates->u_cur_l,tree->rates->u_ml_l,tree->rates->invcov,tree->rates->covdet,2*tree->n_otu-3,YES); */
/*   tree->rates->c_lnL = RATES_Lk_Rates(tree); */

/*   RATES_Check_Lk_Rates(tree,&err); */

/*   if(err) */
/*     { */
/*       PhyML_Printf("\n. Small lk detected"); */
/*       PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
/*       Exit("\n"); */
/*     } */

/*   PhyML_Printf("\n. T0=%f T1=%f T2=%f T0_MIN=%f T0_MAX=%f U0=%f U1=%f U2=%f L=%f ML(L)=%f NU=%f", */
/* 	       t0,t1,t2,t0_min,t0_max,u0/cr,u1/cr,u2/cr, */
/* 	       tree->e_root->l,tree->rates->u_ml_l[tree->e_root->num],tree->rates->nu); */


  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc,tree);
  if(!(tree->mcmc->run%tree->mcmc->norm_freq)) 
    {
      RATES_Normalise_Rates(tree);
    }
}

/*********************************************************/

void RATES_Update_Cur_Bl(t_tree *tree)
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

void RATES_Update_Cur_Bl_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  phydbl down_t,up_t,dt,rr,cr;
  
  down_t = tree->rates->nd_t[d->num];
  up_t   = tree->rates->nd_t[a->num];
  dt     = down_t - up_t;
  rr     = tree->rates->nd_r[d->num];
  cr     = tree->rates->clock_r;

  tree->rates->cur_l[d->num] = dt*rr*cr;
  
  if(tree->rates->cur_l[d->num] < BL_MIN) tree->rates->cur_l[d->num] = BL_MIN;

  if(b) 
    {
      b->l                         = tree->rates->cur_l[d->num];
      tree->rates->u_cur_l[b->num] = tree->rates->cur_l[d->num];
    }

  if(b && isnan(b->l))
    {
      PhyML_Printf("\n. dt=%G rr=%G cr=%G",dt,rr,cr);
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

void RATES_Bl_To_Ml(t_tree *tree)
{
  RATES_Bl_To_Ml_Pre(tree->n_root,tree->n_root->v[0],NULL,tree);
  RATES_Bl_To_Ml_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);
  tree->rates->u_ml_l[tree->e_root->num] = tree->t_edges[tree->e_root->num]->l;
  tree->rates->ml_l[tree->n_root->v[0]->num] = tree->rates->u_ml_l[tree->e_root->num] * tree->n_root_pos;
  tree->rates->ml_l[tree->n_root->v[1]->num] = tree->rates->u_ml_l[tree->e_root->num] * (1. - tree->n_root_pos);
}

/*********************************************************/

void RATES_Bl_To_Ml_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{

  if(b) 
    {
      tree->rates->u_ml_l[b->num] = b->l;
      tree->rates->ml_l[d->num]   = b->l;
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

void RATES_Get_Cov_Matrix_Rooted(phydbl *unroot_cov, t_tree *tree)
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

void RATES_Get_Cov_Matrix_Rooted_Pre(t_node *a, t_node *d, t_edge *b, phydbl *cov, t_tree *tree)
{
  int i, dim, done_left;
  t_node *n;

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

void RATES_Covariance_Mu(t_tree *tree)
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

void RATES_Variance_Mu_Pre(t_node *a, t_node *d, t_tree *tree)
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

void RATES_Fill_Lca_Table(t_tree *tree)
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
/* Get V(L_{i} | L_{-i}) for all i */
void RATES_Get_Conditional_Variances(t_tree *tree)
{
  int i,j;
  short int *is_1;
  phydbl *a;
  int dim;
  t_edge *b;
  phydbl *cond_mu, *cond_cov;

  dim      = 2*tree->n_otu-3;
  a        = tree->rates->_2n_vect1;
  is_1     = tree->rates->_2n_vect5;
  b        = NULL;
  cond_mu  = tree->rates->_2n_vect2;
  cond_cov = tree->rates->_2n2n_vect1;

  For(i,dim) a[i] = tree->rates->u_ml_l[i] * (Uni() * 0.2 + 0.9);

  For(i,dim)
    {
      b = tree->t_edges[i];

      For(j,dim) is_1[j] = 0;
      is_1[b->num]       = 1;

      For(j,dim*dim) cond_cov[j] = 0.0;
      For(j,dim)     cond_mu[j]  = 0.0;
      Normal_Conditional(tree->rates->u_ml_l,tree->rates->cov,a,dim,is_1,1,cond_mu,cond_cov);
      
      tree->rates->cond_var[b->num] = cond_cov[b->num*dim+b->num];
    }
}

/*********************************************************/

void RATES_Get_All_Reg_Coeff(t_tree *tree)
{
  int i,j;
  short int *is_1;
  phydbl *a;
  int dim;
  t_edge *b;

  dim  = 2*tree->n_otu-3;
  a    = tree->rates->_2n_vect1;
  is_1 = tree->rates->_2n_vect5;
  b    = NULL;

  For(i,dim) a[i] = tree->rates->u_ml_l[i] * (Uni() * 0.2 + 0.9);

  For(i,dim)
    {
      b = tree->t_edges[i];

      For(j,dim) is_1[j] = 0;
      is_1[b->num]       = 1;
      
      Get_Reg_Coeff(tree->rates->u_ml_l,tree->rates->cov,a,dim,is_1,1,tree->rates->reg_coeff+b->num*dim);
    }
}

/*********************************************************/

/* Get V(L_{i} | L_{-i}) for all i */
void RATES_Get_Trip_Conditional_Variances(t_tree *tree)
{
  int i,j;
  short int *is_1;
  phydbl *a;
  phydbl *cond_mu, *cond_cov;
  t_node *n;
  int n_otu;

  a        = tree->rates->_2n_vect1;
  is_1     = tree->rates->_2n_vect5;
  cond_mu  = tree->rates->_2n_vect2;
  cond_cov = tree->rates->_2n2n_vect1;
  n        = NULL;
  n_otu    = tree->n_otu;

  For(i,2*n_otu-3) a[i] = tree->rates->u_ml_l[i] * (Uni() * 0.2 + 0.9);

  For(i,2*n_otu-2)
    {
      n = tree->noeud[i];
      if(!n->tax)
	{
	  For(j,2*n_otu-3) is_1[j] = 0;
	  is_1[n->b[0]->num] = 1;
	  is_1[n->b[1]->num] = 1;
	  is_1[n->b[2]->num] = 1;
	  
	  Normal_Conditional_Unsorted(tree->rates->u_ml_l,tree->rates->cov,a,2*n_otu-3,is_1,3,cond_mu,cond_cov);
	  
	  For(j,9) tree->rates->trip_cond_cov[n->num*9+j] = cond_cov[j];
	}
    }
}

/*********************************************************/

void RATES_Get_All_Trip_Reg_Coeff(t_tree *tree)
{
  int i,j;
  short int *is_1;
  phydbl *a;
  t_node *n;
  int n_otu;

  a     = tree->rates->_2n_vect1;
  is_1  = tree->rates->_2n_vect5;
  n     = NULL;
  n_otu = tree->n_otu;

  For(i,2*n_otu-3) a[i] = tree->rates->u_ml_l[i] * (Uni() * 0.2 + 0.9);

  For(i,2*n_otu-2)
    {
      n = tree->noeud[i];
      if(!n->tax)
	{	  
	  For(j,2*n_otu-3) is_1[j] = 0;
	  is_1[n->b[0]->num] = 1;
	  is_1[n->b[1]->num] = 1;
	  is_1[n->b[2]->num] = 1;
	  
	  Get_Reg_Coeff(tree->rates->u_ml_l,tree->rates->cov,a,2*n_otu-3,is_1,3,tree->rates->trip_reg_coeff+n->num*(6*n_otu-9));
	}
    }
}

/*********************************************************/

void RATES_Check_Lk_Rates(t_tree *tree, int *err)
{
  int i;
  phydbl u,u_anc;
  phydbl t,t_anc;
  phydbl z;
  
  *err = 0;

  z = 0.0;
  For(i,2*tree->n_otu-2)
    {
      u     = tree->rates->nd_r[i];
      u_anc = tree->rates->nd_r[tree->noeud[i]->anc->num];
      t     = tree->rates->nd_t[i];
      t_anc = tree->rates->nd_t[tree->noeud[i]->anc->num];

/*       z = fabs(u - u_anc) / sqrt(tree->rates->nu * (t - t_anc)); */

/*       if(z > tree->rates->z_max)  */
/* 	{ */
/* 	  PhyML_Printf("\n. %d %d u=%f u_anc=%f t=%f t_anc=%f z=%f",i,tree->noeud[i]->anc->num,u,u_anc,t,t_anc,z); */
/* 	  PhyML_Printf("\n. %d %d %d",tree->n_root->num,tree->n_root->v[0]->num,tree->n_root->v[1]->num); */
/* 	  *err = 1; */
/* 	} */
    }
}

/*********************************************************/

phydbl RATES_Expected_Tree_Length(t_tree *tree)
{
/*   int i; */
/*   phydbl u,u_anc,t,t_anc; */
/*   phydbl exp_r; */

/*   exp_r = 0.0; */
/*   For(i,2*tree->n_otu-2) */
/*     { */
/*       u     = tree->rates->nd_r[i]; */
/*       u_anc = tree->rates->nd_r[tree->noeud[i]->anc->num]; */
/*       t     = tree->rates->nd_t[i]; */
/*       t_anc = tree->rates->nd_t[tree->noeud[i]->anc->num]; */
      
/*       exp_r += Cond_Exp(u_anc,sqrt((t-t_anc)*tree->rates->nu),tree->rates->min_rate,tree->rates->max_rate); */
/*     } */

/*   return exp_r / (phydbl)(2*tree->n_otu-2); */

  int n;
  phydbl mean;

  n    = 0;
  mean = 0.0;
  RATES_Expected_Tree_Length_Pre(tree->n_root,tree->n_root->v[0],1.0,&mean,&n,tree);
  RATES_Expected_Tree_Length_Pre(tree->n_root,tree->n_root->v[1],1.0,&mean,&n,tree);

  if(n != 2*tree->n_otu-2)
    {
      PhyML_Printf("\n. n=%d 2n-2=%d",n,2*tree->n_otu-2);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  return mean;
}

/*********************************************************/

void RATES_Expected_Tree_Length_Pre(t_node *a, t_node *d, phydbl eranc, phydbl *mean, int *n, t_tree *tree)
{
  phydbl erdes;
  int i;
  phydbl u,u_anc,t,t_anc;
  phydbl loc_mean;
  int loc_n;

  u     = tree->rates->nd_r[d->num];
  u_anc = tree->rates->nd_r[a->num];
  t     = tree->rates->nd_t[d->num];
  t_anc = tree->rates->nd_t[a->num];
  
  erdes = Norm_Trunc_Mean(eranc,sqrt((t-t_anc)*tree->rates->nu),tree->rates->min_rate,tree->rates->max_rate);

  loc_mean = *mean;
  loc_n = *n;

  loc_mean *= loc_n;
  loc_mean += erdes;
  loc_mean /= (loc_n + 1);

  *mean = loc_mean;
  *n    = loc_n + 1;


  if(d->tax) return;
  else
    For(i,3)
      if(d->v[i] != a && d->b[i] != tree->e_root)
	RATES_Expected_Tree_Length_Pre(d,d->v[i],erdes,mean,n,tree);
}

/*********************************************************/

void RATES_Normalise_Rates(t_tree *tree)
{
  phydbl expr,curr;
  int i;
  

  curr = RATES_Check_Mean_Rates(tree);

  /* Set expected mean rate to one such that clock_r is
     easy to interpret regarding the mean values of nd_r */
  expr = 1.0;
  
  For(i,2*tree->n_otu-2) tree->rates->nd_r[i] *= expr/curr;

  tree->rates->clock_r *= curr/expr;
  /* Branch lengths therefore do not change */

  RATES_Update_Cur_Bl(tree);
}

/*********************************************************/

phydbl RATES_Find_Max_Dt_Bisec(phydbl r, phydbl r_mean, phydbl ta, phydbl tc, phydbl nu, phydbl threshp, int inf)
{
  phydbl cdfr,cdf0;
  phydbl sd;
  phydbl trunc_cdf;
  phydbl ori_tc, ori_ta;
  phydbl tb;

  ori_tc = tc;
  ori_ta = ta;

/*   PhyML_Printf("\n Max %s r=%f r_mean%f",inf?"inf":"sup",r,r_mean); */
  do
    {
      tb = ta + (tc - ta)/2.;


      sd   = sqrt(nu*(ori_tc - tb));
      cdfr = Pnorm(r ,r_mean,sd);
      cdf0 = Pnorm(.0,r_mean,sd);
      
      if(inf)
	trunc_cdf = (cdfr - cdf0)/(1. - cdf0);
      else
	trunc_cdf = (1. - cdfr)/(1. - cdf0);
	
/*       PhyML_Printf("\n. ta=%15f tb=%15f tc=%15f cdf = %15f",ta,tb,tc,trunc_cdf); */

      if(trunc_cdf > threshp)
	{
	  ta = tb;
	}
      else
	{
	  tc = tb;
	}

    }while((tc - ta)/(ori_tc - ori_ta) > 0.001);

  if(tb < ori_ta)
    {
      PhyML_Printf("\n. tb < ta r=%f r_mean=%f",r,r_mean);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }
  if(tb > ori_tc)
    {
      PhyML_Printf("\n. tb > tc r=%f r_mean=%f ori_ta=%f ori_tc=%f tb=%f",r,r_mean,ori_ta,ori_tc,tb);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  return tb;
}

/*********************************************************/

phydbl RATES_Find_Min_Dt_Bisec(phydbl r, phydbl r_mean, phydbl ta, phydbl tc, phydbl nu, phydbl threshp, int inf)
{
  phydbl cdfr,cdf0;
  phydbl sd;
  phydbl trunc_cdf;
  phydbl ori_tc, ori_ta;
  phydbl tb;

  ori_tc = tc;
  ori_ta = ta;

/*   PhyML_Printf("\n Min %s r=%f r_mean=%f",inf?"inf":"sup",r,r_mean); */
  do
    {

      tb = ta + (tc - ta)/2.;


      sd   = sqrt(nu*(tb - ori_ta));
      cdfr = Pnorm(r ,r_mean,sd);
      cdf0 = Pnorm(.0,r_mean,sd);
      
      if(inf)
	trunc_cdf = (cdfr - cdf0)/(1. - cdf0);
      else
	trunc_cdf = (1. - cdfr)/(1. - cdf0);
	
/*       PhyML_Printf("\n. ta=%15f tb=%15f tc=%15f cdf = %15f",ta,tb,tc,trunc_cdf); */

      if(trunc_cdf > threshp)
	{
	  tc = tb;
	}
      else
	{
	  ta = tb;
	}

    }while((tc - ta)/(ori_tc - ori_ta) > 0.001);

  if(tb < ori_ta)
    {
      PhyML_Printf("\n. tb < ta r=%f r_mean=%f",r,r_mean);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }
  if(tb > ori_tc)
    {
      PhyML_Printf("\n. tb > tc r=%f r_mean=%f ori_ta=%f ori_tc=%f tb=%f",r,r_mean,ori_ta,ori_tc,tb);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  return tb;
}

/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/


#endif
