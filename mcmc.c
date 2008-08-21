/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/



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
#include "mcmc.h"
#include <stdlib.h>
#include <unistd.h>


/*********************************************************/

void MCMC(arbre *tree)
{
  FILE *fp;
  int pid;
  char *filename;
  phydbl u;
  int n_moves;

  filename = (char *)mCalloc(T_MAX_FILE,sizeof(char));

  pid = getpid();
  printf("\n\n. pid=%d\n\n",pid);
  
  switch(tree->rates->model)
    {
    case COMPOUND_COR :
      {
	if(tree->rates->approx == 1)
	  strcpy(filename,"cor1");
	else
	  strcpy(filename,"cor2");
	break;
      }
    case COMPOUND_NOCOR :
      {
	strcpy(filename,"uncor");
	break;
      }
    case EXPONENTIAL :
      {
	strcpy(filename,"expo");
	break;
      }
    }

  sprintf(filename+strlen(filename),".%d",pid);
  fp = fopen(filename,"w");
  
  tree->mcmc->sample_interval = 100;

  MCMC_Print_Param(fp,tree);
  MCMC_Print_Param(stdout,tree);

/*   MCMC_Randomize_Branch_Lengths(tree); */
/*   MCMC_Randomize_Lexp(tree); */
  MCMC_Randomize_Alpha(tree);
/*   MCMC_Randomize_Node_Times(tree); */


  RATES_Lk_Rates(tree);
  Lk(tree);

  n_moves = 6;
  do
    {
      tree->mcmc->run++;
      
      u = Uni();
      u = rint(u*(n_moves));


      switch((int)u)
	{
	case 0 : { MCMC_Lexp(tree);      break; }
	case 1 : { MCMC_Alpha(tree);     break; }
 	case 2 : { MCMC_Rates(tree);     break; }
	case 3 : { MCMC_Times(tree);     break; }
/* /\* 	case 4 : { MCMC_Nu(tree);        break; } *\/ */
/* 	case 5 : { MCMC_No_Change(tree); break; } */
	}
	    
      MCMC_Print_Param(fp,tree);
      MCMC_Print_Param(stdout,tree);

      if(!(tree->mcmc->run%50))
	{
	  RATES_Adjust_Clock_Rate(tree);
	  RATES_Lk_Rates(tree);
	  Lk(tree);
	}
    }
  while(tree->mcmc->run < 200000);

  fclose(fp);
  Free(filename);
}

/*********************************************************/

void MCMC_Lexp(arbre *tree)
{
  phydbl cur_lexp,new_lexp,cur_lnL,new_lnL;
  phydbl u,alpha,prior_mean_lexp,ratio;
  
  if((tree->rates->model != COMPOUND_NOCOR) &&
     (tree->rates->model != COMPOUND_COR)) return;

  cur_lnL = UNLIKELY;
  new_lnL = UNLIKELY;
  cur_lexp = -1.0;
  new_lexp = -1.0;
  prior_mean_lexp = 0.1;
  ratio = -1.0;

  cur_lnL    = tree->rates->c_lnL;
  cur_lexp   = tree->rates->lexp;
  
  u = Uni();
/*   if(tree->mcmc->run < 2000) */
/*     new_lexp   = cur_lexp  * exp(1.0*(u-0.5)); */
/*   else */
    new_lexp   = cur_lexp  * exp(H_MCMC_LEXP*(u-0.5));

  if((new_lexp  > 1.E-5) && (new_lexp  < 2.0))
    {
      tree->rates->lexp = new_lexp;
      
      new_lnL = RATES_Lk_Rates(tree);
      
      /*       printf("\n. run %4d new_lexp = %f new_lnL = %f",run+1,new_lexp,new_lnL); */
      /* 	  ratio = exp(new_lnL-cur_lnL)*(new_lexp/cur_lexp)*(new_alpha/cur_alpha); */
      
      ratio = 
	exp(new_lnL-cur_lnL)*
	(new_lexp/cur_lexp) *
	exp((1./prior_mean_lexp)*(cur_lexp-new_lexp));
      
      alpha = MIN(1.,ratio);
      
      u = Uni();
      if(u > alpha) /* Reject */
	{
	  tree->rates->lexp  = cur_lexp;
	  RATES_Lk_Rates(tree);
	}
      else
	{
	  tree->mcmc->acc_lexp++;
	}
    }
}

/*********************************************************/

void MCMC_Nu(arbre *tree)
{
  phydbl cur_nu,new_nu,cur_lnL,new_lnL;
  phydbl u,alpha,prior_mean_nu,ratio;
  
  if(tree->rates->model != EXPONENTIAL) return;

  cur_lnL = UNLIKELY;
  new_lnL = UNLIKELY;
  cur_nu = -1.0;
  new_nu = -1.0;
  prior_mean_nu = 1.0;
  ratio = -1.0;

  cur_lnL    = tree->rates->c_lnL;
  cur_nu   = tree->rates->nu;
  
  u = Uni();
  if(tree->mcmc->run < 2000)
    new_nu   = cur_nu  * exp(1.0*(u-0.5));
  else
    new_nu   = cur_nu  * exp(H_MCMC_NU*(u-0.5));

  if((new_nu  > 1.E-5) && (new_nu  < 10.0))
    {
      tree->rates->nu = new_nu;
      
      new_lnL = RATES_Lk_Rates(tree);
      
      /*       printf("\n. run %4d new_nu = %f new_lnL = %f",run+1,new_nu,new_lnL); */
      /* 	  ratio = exp(new_lnL-cur_lnL)*(new_nu/cur_nu)*(new_alpha/cur_alpha); */
      
      ratio = 
	exp(new_lnL-cur_lnL)*
	(new_nu/cur_nu) *
	exp((1./prior_mean_nu)*(cur_nu-new_nu));
      
      alpha = MIN(1.,ratio);
      
      u = Uni();
      if(u > alpha) /* Reject */
	{
	  tree->rates->nu  = cur_nu;
	  RATES_Lk_Rates(tree);
	}
      else
	{
	  tree->mcmc->acc_nu++;
	}
    }
}

/*********************************************************/

void MCMC_Alpha(arbre *tree)
{
  phydbl cur_lnL,new_lnL,cur_alpha,new_alpha;
  phydbl u,alpha,ratio;
  
  if((tree->rates->model != COMPOUND_NOCOR) &&
     (tree->rates->model != COMPOUND_COR)) return;

  cur_lnL = UNLIKELY;
  new_lnL = UNLIKELY;
  ratio = -1.0;

  cur_lnL    = tree->rates->c_lnL;
  cur_alpha  = tree->rates->alpha;
  
  u =  Uni();
  new_alpha = cur_alpha + (u * 2.*cur_alpha/10. - cur_alpha/10.);

  if(new_alpha > 1.0 && new_alpha < 7.0)
    {
      tree->rates->alpha = new_alpha;
      new_lnL = RATES_Lk_Rates(tree);      
      ratio = exp(new_lnL-cur_lnL);
      alpha = MIN(1.,ratio);
      u = Uni();
 
      if(u > alpha) /* Reject */
	{
	  tree->rates->alpha = cur_alpha;
	  RATES_Lk_Rates(tree);
	}
    }
}

/*********************************************************/

void MCMC_Times(arbre *tree)
{
  MCMC_Times_Pre(tree->n_root,tree->n_root->v[0],NULL,tree);
  MCMC_Times_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);
}

/*********************************************************/

void MCMC_Rates(arbre *tree)
{
/*   phydbl new_lnL_rate, new_lnL_data; */
/*   phydbl cur_lnL_rate, cur_lnL_data; */
/*   phydbl ratio,alpha,u; */
 
/*   cur_lnL_data = tree->c_lnL; */
/*   cur_lnL_rate = tree->rates->c_lnL; */
  
/*   Record_Br_Len(NULL,tree); */

  MCMC_Rates_Pre(tree->n_root,tree->n_root->v[0],NULL,tree);
  MCMC_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);

/*   new_lnL_rate = RATES_Lk_Rates(tree); */
/*   new_lnL_data = Return_Lk(tree); */
  
/*   ratio = */
/*     (new_lnL_data + new_lnL_rate ) - */
/*     (cur_lnL_data + cur_lnL_rate ); */
  
/*   ratio = exp(ratio); */
  
/*   alpha = MIN(1.,ratio); */
  
/*   u = Uni(); */
  
/*   if(u > alpha) /\* Reject *\/ */
/*     { */
/*       Restore_Br_Len(NULL,tree); */
      
/*       RATES_Lk_Rates(tree); */
/*       Lk(tree); */
      
      
/*       if((fabs(cur_lnL_data - tree->c_lnL) > 1.E-3) || (fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-0)) */
/* 	{ */
/* 	  printf("\n. lexp=%f alpha=%f cur_lnL_data = %f vs %f ; cur_lnL_rates = %f vs %f", */
/* 		 tree->rates->lexp, */
/* 		 tree->rates->alpha, */
/* 		 cur_lnL_data,tree->c_lnL, */
/* 		 cur_lnL_rate,tree->rates->c_lnL); */
	  
/* 	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
/* 	  Warn_And_Exit(""); */
/* 	} */
/*       if(fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-3) */
/* 	{ */
/* 	  printf("\n. WARNING: numerical precision issue detected (diff=%G). Reseting the likelihood.\n",cur_lnL_rate - tree->rates->c_lnL); */
/* 	  RATES_Lk_Rates(tree); */
/* 	  Lk(tree); */
/* 	} */
/*     } */
/*   else */
/*     { */
/*       tree->mcmc->acc_rates++; */
/*     } */

}

/*********************************************************/

void MCMC_One_Rate(edge *b, arbre *tree)
{
  phydbl u;
  phydbl new_mu, cur_mu, new_lnL_data, cur_lnL_data, new_lnL_rate, cur_lnL_rate, cur_l, ratio, dt, alpha;
  

  if(b != tree->e_root)
    {
      cur_l = b->l;
      
      dt = fabs(tree->rates->t[b->left->num] - tree->rates->t[b->rght->num]);
      if(dt < MIN_DT) dt = MIN_DT;
      cur_mu = b->l / (dt*tree->rates->clock_r);
      
      cur_lnL_data = tree->c_lnL;
      cur_lnL_rate = tree->rates->c_lnL;
      /*      cur_lnL_rate = RATES_Lk_Rates(tree); */
      
      u = Uni();
/*       if(tree->mcmc->run < 10000) new_mu = cur_mu * exp(1.*(u-0.5)); */
/*       else new_mu = cur_mu * exp(0.5*(u-0.5)); */
      new_mu = cur_mu * exp(0.1*(u-0.5));

      b->l = dt * new_mu * tree->rates->clock_r; 
      new_lnL_rate = RATES_Lk_Change_One_Rate(b,new_mu,tree); /* Must be called before Lk_At_Given_Edge because
								 the rate change (branch length change) is performed
								 by this function */
      new_lnL_data = Lk_At_Given_Edge(b,tree);
      
      ratio = 
	(new_lnL_data + new_lnL_rate + log(new_mu)) -  
	(cur_lnL_data + cur_lnL_rate + log(cur_mu));
      
      ratio = exp(ratio);
      
      alpha = MIN(1.,ratio);
      
      u = Uni();
      
      if(u > alpha) /* Reject */
	{
	  b->l = cur_l; 
	  
	  RATES_Lk_Rates(tree);
	  Lk(tree);
	  
	  if((fabs(cur_lnL_data - tree->c_lnL) > 1.E-3) || (fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-0))
	    {
	      printf("\n. cur_lnL_data = %f vs %f ; cur_lnL_rates = %f vs %f",
		     cur_lnL_data,tree->c_lnL,
		     cur_lnL_rate,tree->rates->c_lnL);
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	  if(fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-3)
	    {
	      printf("\n. WARNING: numerical precision issue detected (diff=%G). Reseting the likelihood.\n",cur_lnL_rate - tree->rates->c_lnL);
	      RATES_Lk_Rates(tree);
	    }
	}
      else
	{
	  RATES_Lk_Rates(tree);
	  Lk(tree);
	}
    }
}

/*********************************************************/

void MCMC_Rates_Pre(node *a, node *d, edge *b, arbre *tree)
{
  phydbl u;
  phydbl new_mu, cur_mu, new_lnL_data, cur_lnL_data, new_lnL_rate, cur_lnL_rate, cur_l, ratio, dt, alpha;
  int i;
  
  if(b)
    {
      cur_l = b->l;

      dt = fabs(tree->rates->t[a->num] - tree->rates->t[d->num]);
      if(dt < MIN_DT) dt = MIN_DT;
      cur_mu = b->l / (dt*tree->rates->clock_r);
      
      cur_lnL_data = tree->c_lnL;
      cur_lnL_rate = tree->rates->c_lnL;

      u = Uni();

      new_mu = cur_mu * exp(H_MCMC_RATES*(u-0.5));
/*       new_mu =  u * 2*cur_mu/100. + cur_mu - cur_mu/100.; */
/*       new_mu =  u * 2*cur_mu/100. + cur_mu - cur_mu/100.; */

      b->l = dt * new_mu * tree->rates->clock_r; 


      new_lnL_rate = RATES_Lk_Change_One_Rate(b,new_mu,tree); /* Must be called before Lk_At_Given_Edge because
								 the rate change (branch length change) is performed
								 by this function */
      new_lnL_data = Lk_At_Given_Edge(b,tree);
      
      ratio =
	(new_lnL_data + new_lnL_rate + log(new_mu)) -
	(cur_lnL_data + cur_lnL_rate + log(cur_mu));
      
      ratio = exp(ratio);
	
      alpha = MIN(1.,ratio);
      
      u = Uni();
      
      if(u > alpha) /* Reject */
	{
	  b->l = cur_l;

	  RATES_Lk_Change_One_Rate(b,cur_mu,tree); /* Must be called before Lk_At_Given_Edge because
						      the rate change (branch length change) is performed
						      by this function */
/* 	  RATES_Lk_Rates(tree); */
	  Lk_At_Given_Edge(b,tree);
	  

	  if((fabs(cur_lnL_data - tree->c_lnL) > 1.E-3) || (fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-0))
	    {
	      printf("\n. lexp=%f alpha=%f b->l = (%f %f) dt=(%f %f); cur_lnL_data = %f vs %f ; cur_lnL_rates = %f vs %f",
		     tree->rates->lexp,
		     tree->rates->alpha,
		     b->l,cur_l,
		     dt,fabs(tree->rates->t[a->num] - tree->rates->t[d->num]),
		     cur_lnL_data,tree->c_lnL,
		     cur_lnL_rate,tree->rates->c_lnL);
	      
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	  if(fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-3)
	    {
	      printf("\n. WARNING: numerical precision issue detected (diff=%G). Reseting the likelihood.\n",cur_lnL_rate - tree->rates->c_lnL);
	      RATES_Lk_Rates(tree);
	      Lk(tree);
	    }
	}
      else
	{
	  tree->mcmc->acc_rates++;
	  /* 	  printf("\n. Accept %f->%f (%f %f)", */
	  /* 		 cur_lnL_rate,new_lnL_rate, */
	  /* 		 cur_mu,new_mu); */
	  /* 	  Lk(tree); */
	  /* 	  printf("\n. Accept change no branch %d",b->num); */
	}
    }

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      Update_P_Lk(tree,d->b[i],d);
	      MCMC_Rates_Pre(d,d->v[i],d->b[i],tree);
	    }
	}
      if(b) { Update_P_Lk(tree,b,d); }
      else  { Update_P_Lk(tree,tree->e_root,d); }
    }
}

/*********************************************************/

void MCMC_No_Change(arbre *tree)
{
  phydbl u;
  int i;
  node *n;
  edge *b;
  phydbl mu_anc,dt_anc;
  phydbl new_mu, cur_mu, new_lnL_data, cur_lnL_data, new_lnL_rate, cur_lnL_rate, ratio, dt, alpha, cur_l;

  
  n = NULL;
  b = NULL;

  cur_lnL_data = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL;

  u = Uni();
  u *= tree->n_otu-3;
  u += tree->n_otu;

  n = tree->noeud[(int)u];

  if(n->tax)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  if(n == tree->n_root)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  u = Uni();
  For(i,3)
    {
      if((n->v[i] != n->anc) && (n->b[i] != tree->e_root))
	{
	  if(u < 0.5) b=n->b[i];
	  else        b=n->b[i];
	}
    }
  
  if(b == tree->e_root)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  cur_l = b->l;
  
      
  dt = fabs(tree->rates->t[b->left->num] - tree->rates->t[b->rght->num]);
  if(dt < MIN_DT) dt = MIN_DT;
  cur_mu = b->l / (dt*tree->rates->clock_r);
     
  mu_anc = -1.0;

  dt_anc = fabs(tree->rates->t[n->num] - tree->rates->t[n->anc->num]);
  if((n == tree->n_root->v[0]) || (n == tree->n_root->v[1]))
    dt_anc = fabs(tree->rates->t[n->num] - tree->rates->t[tree->n_root->num]);

  For(i,3)
    {
      if(n->v[i] == n->anc)
	{
	  mu_anc = n->b[i]->l / (dt_anc*tree->rates->clock_r);
	  break;
	}
      if(n->b[i] == tree->e_root)
	{
	  if(n == tree->n_root->v[0])
	    mu_anc = tree->n_root->l[0] / (dt_anc*tree->rates->clock_r);
	  if(n == tree->n_root->v[1])
	    mu_anc = tree->n_root->l[1] / (dt_anc*tree->rates->clock_r);
	  else
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	  break;
	}
    }

  new_mu = mu_anc;
  
  if(new_mu < 0.0)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }


  b->l = dt * new_mu * tree->rates->clock_r;

/*   Record_Br_Len(NULL,tree); */
/*   For(i,2*tree->n_otu-3) */
/*     { */
/*       if(tree->t_edges[i] != tree->e_root) */
/* 	{ */
/* 	  u = Uni(); */
/* 	  if(u < 0.3) */
/* 	    { */
/* 	      tree->t_edges[i]->l =  */
/* 		fabs(tree->rates->t[tree->t_edges[i]->left->num] - tree->rates->t[tree->t_edges[i]->rght->num]) * */
/* 		tree->rates->clock_r; */
/* 	    } */
/* 	} */
/*     } */

  new_lnL_rate = RATES_Lk_Rates(tree);
  new_lnL_data = Return_Lk(tree);
      
  ratio = 
    (new_lnL_data + new_lnL_rate) -  
    (cur_lnL_data + cur_lnL_rate);
  
  ratio = exp(ratio);
	
  alpha = MIN(1.,ratio);
      
  u = Uni();
      
  if(u > alpha) /* Reject */
    {
/*       printf("\n. Reject"); */
/*       Restore_Br_Len(NULL,tree); */

      b->l = cur_l;
      
      RATES_Lk_Rates(tree);
      Lk(tree);
	  
      if((fabs(cur_lnL_data - tree->c_lnL) > 1.E-3) || (fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-0))
	{
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      if(fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-3)
	{
	  printf("\n. WARNING: numerical precision issue detected (diff=%G). Reseting the likelihood.\n",cur_lnL_rate - tree->rates->c_lnL);
	  RATES_Lk_Rates(tree);
	}
    }
  else
    {
/*       printf("\n. Accept"); */
    }
}

/*********************************************************/

void MCMC_Rates_Root_Branch(arbre *tree)
{
  phydbl dt0, dt1, dens;
  phydbl cur_l0, cur_l1;
  phydbl cur_mu0, cur_mu1;
  phydbl new_mu0,new_mu1;
  phydbl u, ratio, alpha;
  phydbl cur_lnL_data, cur_lnL_rate, new_lnL_rate, new_lnL_data;
  
  cur_l0 = tree->n_root->l[0];
  cur_l1 = tree->n_root->l[1];

  cur_lnL_data = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL;
/*   cur_lnL_rate = RATES_Lk_Rates(tree); */
  
  dt0 = fabs(tree->rates->t[tree->n_root->v[0]->num] - tree->rates->t[tree->n_root->num]);
  if(dt0 < MIN_DT) dt0 = MIN_DT;

  dt1 = fabs(tree->rates->t[tree->n_root->v[1]->num] - tree->rates->t[tree->n_root->num]);
  if(dt1 < MIN_DT) dt1 = MIN_DT;

  cur_mu0 = tree->n_root->l[0] / (dt0 * tree->rates->clock_r);
  cur_mu1 = tree->n_root->l[1] / (dt1 * tree->rates->clock_r);

  u = Uni();
  new_mu0 = cur_mu0 * exp(H_MCMC_RATES * (u-0.5));
  
  u = Uni();
  new_mu1 = cur_mu1 * exp(H_MCMC_RATES * (u-0.5));

  tree->n_root->l[0] = new_mu0 * dt0 * tree->rates->clock_r;
  tree->n_root->l[1] = new_mu1 * dt1 * tree->rates->clock_r;
  tree->e_root->l = tree->n_root->l[0] + tree->n_root->l[1];

  RATES_Update_Triplet(tree->e_root->left,tree);
  RATES_Update_Triplet(tree->e_root->rght,tree);

  dens = RATES_Lk_Rates_Core(new_mu0,new_mu1,dt0,dt1,tree);
  dens *= RATES_Dmu(new_mu0,dt0,tree->rates->alpha,1./tree->rates->alpha,tree->rates->lexp,0);

  tree->rates->c_lnL -= tree->rates->triplet[tree->n_root->num];
  tree->rates->c_lnL += log(dens);
  tree->rates->triplet[tree->n_root->num] = log(dens);
  new_lnL_rate = tree->rates->c_lnL;

/*   new_lnL_rate = RATES_Lk_Rates(tree); */

  new_lnL_data = Lk_At_Given_Edge(tree->e_root,tree);

  ratio = 
    (new_lnL_data + new_lnL_rate + log(new_mu0) + log(new_mu1)) -
     (cur_lnL_data + cur_lnL_rate + log(cur_mu0) + log(cur_mu1));

  ratio = exp(ratio);
  
  alpha = MIN(1.,ratio);
  
  u = Uni();
  
  if(u > alpha) /* Reject */
    {
      tree->n_root->l[0] = cur_l0;
      tree->n_root->l[1] = cur_l1;
      tree->e_root->l = tree->n_root->l[0] + tree->n_root->l[1];

      RATES_Update_Triplet(tree->e_root->left,tree);
      RATES_Update_Triplet(tree->e_root->rght,tree);

      dens = RATES_Lk_Rates_Core(cur_mu0,cur_mu1,dt0,dt1,tree);
      dens *= RATES_Dmu(cur_mu0,dt0,tree->rates->alpha,1./tree->rates->alpha,tree->rates->lexp,0);
      
      tree->rates->c_lnL -= tree->rates->triplet[tree->n_root->num];
      tree->rates->c_lnL += log(dens);
      tree->rates->triplet[tree->n_root->num] = log(dens);

/*       RATES_Lk_Rates(tree); */

      Lk_At_Given_Edge(tree->e_root,tree);
      
      if(fabs(tree->c_lnL - cur_lnL_data) > 1.E-3)
	{
	  printf("\n. c_lnL = %f cur_lnL_data = %f",tree->c_lnL,cur_lnL_data);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      if(fabs(tree->rates->c_lnL - cur_lnL_rate) > 1.E-3)
	{
	  printf("\n. run=%d",tree->mcmc->run);
	  printf("\n. lnL=[%f %f]",tree->rates->c_lnL,cur_lnL_rate);
	  printf("\n. l0=[%f %f]",tree->n_root->l[0],cur_l0);
	  printf("\n. l1=[%f %f]",tree->n_root->l[1],cur_l1);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
    }
  else /* Accept */
    {
      tree->mcmc->acc_rates++;
/*       tree->rates->triplet[tree->n_root->num] = log(dens); */
/*       tree->rates->c_lnL = new_lnL_rate; */
    }
}

/*********************************************************/
/* TO DO  Node times at each extremity of the root edge */
void MCMC_Times_Pre(node *a, node *d, edge *b, arbre *tree)
{
  phydbl u;
  phydbl cur_dt0,cur_dt1,cur_dt2;
  phydbl new_dt0,new_dt1,new_dt2;
  phydbl t_inf,t_sup;
  phydbl cur_t, new_t;
  phydbl cur_lnL_data,cur_lnL_rate;
  phydbl new_lnL_rate;
  phydbl ratio,alpha;
  int    i,dir0,dir1,dir2;
  phydbl cur_l0, cur_l1, cur_l2;
  phydbl cur_lnL_times, new_lnL_times;

  if(d->tax) return; /* Won't change time at tip */

  if(b)
    {
      cur_lnL_data  = tree->c_lnL;
      cur_lnL_rate  = tree->rates->c_lnL;
      cur_t         = tree->rates->t[d->num];

      cur_lnL_times = RATES_Yule(tree);
      
      dir0=dir1=dir2=-1;
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      if(dir1 < 0) dir1 = i;
	      else if(dir2 < 0) dir2 = i;
	    }
	  else dir0 = i;
	}
      
      cur_dt0 = fabs(tree->rates->t[a->num] - tree->rates->t[d->num]);
      cur_dt1 = fabs(tree->rates->t[d->num] - tree->rates->t[d->v[dir1]->num]);
      cur_dt2 = fabs(tree->rates->t[d->num] - tree->rates->t[d->v[dir2]->num]);

      t_inf = MIN(tree->rates->t[d->v[dir1]->num],tree->rates->t[d->v[dir2]->num]);
      t_sup = tree->rates->t[a->num];


      if(t_inf-t_sup < 2*MIN_DT)
	{
	  printf("\n. t_inf = %f t_sup = %f",t_inf,t_sup);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      u = Uni();
      new_t = u * (t_inf-t_sup-2*MIN_DT) + t_sup + MIN_DT;
      tree->rates->t[d->num] = new_t;

      new_dt0 = fabs(tree->rates->t[a->num] - tree->rates->t[d->num]);
      if(new_dt0 < MIN_DT) new_dt0 = MIN_DT;
      new_dt1 = fabs(tree->rates->t[d->num] - tree->rates->t[d->v[dir1]->num]);
      if(new_dt1 < MIN_DT) new_dt1 = MIN_DT;
      new_dt2 = fabs(tree->rates->t[d->num] - tree->rates->t[d->v[dir2]->num]);
      if(new_dt2 < MIN_DT) new_dt2 = MIN_DT;

      cur_l0 = d->b[dir0]->l;
      cur_l1 = d->b[dir1]->l;
      cur_l2 = d->b[dir2]->l;

/*       d->b[dir0]->l = log(d->b[dir0]->l) + log(new_dt0) - log(cur_dt0); */
/*       d->b[dir1]->l = log(d->b[dir1]->l) + log(new_dt1) - log(cur_dt1); */
/*       d->b[dir2]->l = log(d->b[dir2]->l) + log(new_dt2) - log(cur_dt2); */

/*       d->b[dir0]->l = exp(d->b[dir0]->l); */
/*       d->b[dir1]->l = exp(d->b[dir1]->l); */
/*       d->b[dir2]->l = exp(d->b[dir2]->l); */
  
/*       Update_PMat_At_Given_Edge(d->b[dir1],tree); */
/*       Update_PMat_At_Given_Edge(d->b[dir2],tree); */
/*       Update_P_Lk(tree,b,d); */

/*       new_lnL_data  = Lk_At_Given_Edge(b,tree); */
      new_lnL_rate  = RATES_Lk_Change_One_Time(d,new_t,tree); 
      new_lnL_times = RATES_Yule(tree);

/*       ratio = */
/* 	(new_lnL_data + new_lnL_rate + new_lnL_times) - */
/* 	(cur_lnL_data + cur_lnL_rate + cur_lnL_times); */
/*       ratio = */
/* 	(new_lnL_data + new_lnL_rate) - */
/* 	(cur_lnL_data + cur_lnL_rate); */
      ratio =
	(new_lnL_rate + new_lnL_times + log(cur_dt0) + log(cur_dt1) + log(cur_dt2)) -
	(cur_lnL_rate + cur_lnL_times + log(new_dt0) + log(new_dt1) + log(new_dt2));
           
      ratio = exp(ratio);
      
      alpha = MIN(1.,ratio);
      
      u = Uni();
      
      if(u > alpha) /* Reject */
	{
/* 	  d->b[dir0]->l = cur_l0; */
/* 	  d->b[dir1]->l = cur_l1; */
/* 	  d->b[dir2]->l = cur_l2; */

/* 	  Update_PMat_At_Given_Edge(d->b[dir1],tree); */
/* 	  Update_PMat_At_Given_Edge(d->b[dir2],tree); */
/* 	  Update_P_Lk(tree,b,d); */
	  
/* 	  Lk_At_Given_Edge(b,tree); */

	  RATES_Lk_Change_One_Time(d,cur_t,tree);

/* 	  printf("\n. Reject"); */
	  
	  if((fabs(cur_lnL_data - tree->c_lnL) > 1.E-3) || (fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-0))
	    {
	      printf("\n. lexp = %f alpha = %f",
		     tree->rates->lexp,
		     tree->rates->alpha);

	      printf("\n. l0=%f l1=%f l2=%f \n. cur_l0=%f cur_l1=%f cur_l2=%f",
		     d->b[dir0]->l,
		     d->b[dir1]->l,
		     d->b[dir2]->l,
		     cur_l0,
		     cur_l1,
		     cur_l2);
		   
	      printf("\n.  dt0 = %f %f dt1 = %f %f dt2 =%f %f",
		     cur_dt0,new_dt0,
		     cur_dt1,new_dt1,
		     cur_dt2,new_dt2);

	      printf("\n. b->l = %f cur_t = %f ; cur_lnL_data = %f vs %f ; cur_lnL_rates = %f vs %f",
		     b->l,cur_t,
		     cur_lnL_data,tree->c_lnL,
		     cur_lnL_rate,tree->rates->c_lnL);
	      
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }

	  if(fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-3)
	    {
	      printf("\n. WARNING: numerical precision issue detected (diff=%G). Reseting the likelihood.\n",cur_lnL_rate - tree->rates->c_lnL);
	      RATES_Lk_Rates(tree);
	    }
	}
      else
	{
	  tree->mcmc->acc_times++;
/* 	  printf("\n. Accept"); */
	}
    }

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      Update_P_Lk(tree,d->b[i],d);
	      MCMC_Times_Pre(d,d->v[i],d->b[i],tree);
	    }
	}
      if(b) Update_P_Lk(tree,b,d);
      else Update_P_Lk(tree,tree->e_root,d);
    }
}

/*********************************************************/

void MCMC_One_Node_Time(node *d, arbre *tree)
{
  phydbl u;
  phydbl cur_dt0,cur_dt1,cur_dt2;
  phydbl new_dt0,new_dt1,new_dt2;
  phydbl t_inf,t_sup;
  phydbl cur_t, new_t;
  phydbl cur_lnL_data,cur_lnL_rate;
  phydbl new_lnL_data,new_lnL_rate;
  phydbl ratio,alpha;
  int    i,dir0,dir1,dir2;
  phydbl cur_l0, cur_l1, cur_l2;
  phydbl cur_lnL_times, new_lnL_times;

  if(d->tax) return; /* Won't change time at tip */

  cur_lnL_data  = tree->c_lnL;
  cur_lnL_rate  = tree->rates->c_lnL;
  cur_t         = tree->rates->t[d->num];
  
  cur_lnL_times = RATES_Yule(tree);
  
  dir0=dir1=dir2=-1;
  For(i,3)
    {
      if((d->v[i] != d->anc) && (d->b[i] != tree->e_root))
	{
	  if(dir1 < 0) dir1 = i;
	  else if(dir2 < 0) dir2 = i;
	}
      else dir0 = i;
    }
  
  cur_dt0 = fabs(tree->rates->t[d->anc->num] - tree->rates->t[d->num]);
  cur_dt1 = fabs(tree->rates->t[d->num] - tree->rates->t[d->v[dir1]->num]);
  cur_dt2 = fabs(tree->rates->t[d->num] - tree->rates->t[d->v[dir2]->num]);
  
  t_inf = MIN(tree->rates->t[d->v[dir1]->num],tree->rates->t[d->v[dir2]->num]);
  t_sup = tree->rates->t[d->anc->num];
  
  
  if(t_inf-t_sup < 2*MIN_DT)
    {
      printf("\n. d->anc = %d",d->anc->num);
      printf("\n. t_inf = %f t_sup = %f",t_inf,t_sup);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  
  u = Uni();
  new_t = u * (t_inf-t_sup-2*MIN_DT) + t_sup + MIN_DT;
  tree->rates->t[d->num] = new_t;
  
  new_dt0 = fabs(tree->rates->t[d->anc->num] - tree->rates->t[d->num]);
  if(new_dt0 < MIN_DT) new_dt0 = MIN_DT;
  new_dt1 = fabs(tree->rates->t[d->num] - tree->rates->t[d->v[dir1]->num]);
  if(new_dt1 < MIN_DT) new_dt1 = MIN_DT;
  new_dt2 = fabs(tree->rates->t[d->num] - tree->rates->t[d->v[dir2]->num]);
  if(new_dt2 < MIN_DT) new_dt2 = MIN_DT;
  
  cur_l0 = d->b[dir0]->l;
  cur_l1 = d->b[dir1]->l;
  cur_l2 = d->b[dir2]->l;
  
  d->b[dir0]->l = log(d->b[dir0]->l) + log(new_dt0) - log(cur_dt0);
  d->b[dir1]->l = log(d->b[dir1]->l) + log(new_dt1) - log(cur_dt1);
  d->b[dir2]->l = log(d->b[dir2]->l) + log(new_dt2) - log(cur_dt2);
  
  d->b[dir0]->l = exp(d->b[dir0]->l);
  d->b[dir1]->l = exp(d->b[dir1]->l);
  d->b[dir2]->l = exp(d->b[dir2]->l);
  
  new_lnL_data  = Return_Lk(tree);
  new_lnL_rate  = RATES_Lk_Change_One_Time(d,new_t,tree); 
  new_lnL_times = RATES_Yule(tree);
  
  ratio =
    (new_lnL_data + new_lnL_rate + new_lnL_times) -
    (cur_lnL_data + cur_lnL_rate + cur_lnL_times);
  
  ratio = exp(ratio);
  
  alpha = MIN(1.,ratio);
  
  u = Uni();
  
  if(u > alpha) /* Reject */
    {
      d->b[dir0]->l = cur_l0;
      d->b[dir1]->l = cur_l1;
      d->b[dir2]->l = cur_l2;
      
      Lk(tree);
      RATES_Lk_Change_One_Time(d,cur_t,tree);
      
      if((fabs(cur_lnL_data - tree->c_lnL) > 1.E-3) || (fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-0))
	{
	  printf("\n. lexp = %f alpha = %f",
		 tree->rates->lexp,
		 tree->rates->alpha);
	  
	  printf("\n. l0=%f l1=%f l2=%f \n. cur_l0=%f cur_l1=%f cur_l2=%f",
		 d->b[dir0]->l,
		 d->b[dir1]->l,
		 d->b[dir2]->l,
		 cur_l0,
		 cur_l1,
		 cur_l2);
	  
	  printf("\n.  dt0 = %f %f dt1 = %f %f dt2 =%f %f",
		 cur_dt0,new_dt0,
		 cur_dt1,new_dt1,
		 cur_dt2,new_dt2);
	  	  
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      
      if(fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-3)
	{
	  printf("\n. WARNING: numerical precision issue detected (diff=%G). Reseting the likelihood.\n",cur_lnL_rate - tree->rates->c_lnL);
	  RATES_Lk_Rates(tree);
	}
    }	
}


/*********************************************************/

void MCMC_Print_Param(FILE *fp, arbre *tree)
{
  int i;

  if(!(tree->mcmc->run%tree->mcmc->sample_interval)) 
    {
      if(tree->mcmc->run == 0)
	{
	  fprintf(fp,"\n");
	  fprintf(fp,"Run\t");
	  fprintf(fp,"ClockRate\t");
	  fprintf(fp,"LnLSeq\t");
	  fprintf(fp,"LnLRate\t");
	  fprintf(fp,"Lexp\t");
	  fprintf(fp,"Alpha\t");
	  fprintf(fp,"Nu\t");
	  if(fp != stdout) for(i=tree->n_otu;i<2*tree->n_otu-1;i++) fprintf(fp,"T%d\t",i);
/*  	  if(fp != stdout) For(i,2*tree->n_otu-3) fprintf(fp,"B%d\t",i); */
/* 	  For(i,5) fprintf(fp,"B%d ",i); */
	}

      fprintf(fp,"\n");

/*       tree->rates->approx = 1; */
/*       fprintf(fp,"RATES_lnL = [1] %f",RATES_Lk_Rates(tree)); */
/*       tree->rates->approx = 2; */
/*       fprintf(fp," [2] %f\t",RATES_Lk_Rates(tree)); */
	      

      fprintf(fp,"%6d\t",tree->mcmc->run);
      fprintf(fp,"%4.2f\t",RATES_Check_Mean_Rates(tree));
      fprintf(fp,"%15.2f\t",tree->c_lnL);
      fprintf(fp,"%15.2f\t",tree->rates->c_lnL);
      fprintf(fp,"%15f\t",tree->rates->lexp);
      fprintf(fp,"%4.2f\t",tree->rates->alpha);
      fprintf(fp,"%15f\t",tree->rates->nu);
      if(fp != stdout) for(i=tree->n_otu;i<2*tree->n_otu-1;i++) fprintf(fp,"%8f\t",tree->rates->t[i]-tree->rates->true_t[i]);
/*       if(fp != stdout) For(i,2*tree->n_otu-3) fprintf(fp,"%8f\t",tree->t_edges[i]->l); */

      if(fp == stdout) printf("%2.5f\t%2.5f\t%2.5f",
			      (phydbl)tree->mcmc->acc_lexp/tree->mcmc->run,
			      (phydbl)tree->mcmc->acc_rates/((2*tree->n_otu-2)*tree->mcmc->run),
			      (phydbl)tree->mcmc->acc_times/tree->mcmc->run);



      fflush(NULL);
    }
}

/*********************************************************/

tmcmc *MCMC_Make_MCMC_Struct()
{
  tmcmc *mcmc;
  mcmc = (tmcmc *)mCalloc(1,sizeof(tmcmc));
  return(mcmc);
}

/*********************************************************/

void MCMC_Init_MCMC_Struct(tmcmc *mcmc)
{
  mcmc->acc_lexp  = 0;
  mcmc->acc_rates = 0;
  mcmc->acc_times = 0;
  mcmc->acc_nu    = 0;

  mcmc->run = 0;
  mcmc->sample_interval = 500;
}

/*********************************************************/

void MCMC_Randomize_Branch_Lengths(arbre *tree)
{
  int i;
  phydbl u;

  For(i,2*tree->n_otu-3)
    {
      if(tree->t_edges[i] != tree->e_root)
	{
	  u = Uni();
	  tree->t_edges[i]->l *= -log(u);
	}
      else
	{
	  printf("\n. Didn't randomize root edge.");
	}
    }
}

/*********************************************************/

void MCMC_Randomize_Lexp(arbre *tree)
{
  phydbl u;

  tree->rates->lexp = -1.0;
  do
    {
      u = Uni();
      tree->rates->lexp = -log(u);
    }
  while((tree->rates->lexp < 1.E-5) && (tree->rates->lexp > 2.0));
}

/*********************************************************/

void MCMC_Randomize_Nu(arbre *tree)
{
  phydbl u;
  do
    {
      u = Uni();
      tree->rates->nu = -log(u);
    }
  while((tree->rates->nu < 1.E-5) && (tree->rates->nu > 10.0));
}

/*********************************************************/

void MCMC_Randomize_Alpha(arbre *tree)
{
  phydbl u;

  u = Uni();
  tree->rates->alpha = u*6.0+1.0;
}

/*********************************************************/

void MCMC_Randomize_Node_Times(arbre *tree)
{
  MCMC_Randomize_Node_Times_Pre(tree->n_root,tree->n_root->v[0],tree);
  MCMC_Randomize_Node_Times_Pre(tree->n_root,tree->n_root->v[1],tree);
}

/*********************************************************/

void MCMC_Randomize_Node_Times_Pre(node *a, node *d, arbre *tree)
{
  int i;

  if(d->tax) return;
  else
    {
      node *v1, *v2; /* the two sons of d */
      phydbl t_sup, t_inf;
      phydbl u;

      v1 = v2 = NULL;
      For(i,3) if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	{
	  if(!v1) v1 = d->v[i]; 
	  else    v2 = d->v[i];
	}
	  
      t_inf = MIN(tree->rates->t[v1->num],tree->rates->t[v2->num]);
      t_sup = tree->rates->t[a->num];

      u = Uni();
      u *= (t_inf - t_sup);
      u += t_sup;

      if((d != tree->n_root->v[0]) && (d != tree->n_root->v[1])) tree->rates->t[d->num] = u;
      else printf("\n. Didn't randomize the root");
    }

  
  For(i,3)
    {
      if((d->v[i] != a) && (d->b[i] != tree->e_root))
	{
	  MCMC_Randomize_Node_Times_Pre(d,d->v[i],tree);	      
	}
    }
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
