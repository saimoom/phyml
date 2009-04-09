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
#include "numeric.h"
#include <stdlib.h>
#include <unistd.h>


/*********************************************************/

void MCMC(char *filename, arbre *tree)
{
  int n_moves,i;

  tree->mcmc = (tmcmc *)MCMC_Make_MCMC_Struct(tree);
  MCMC_Init_MCMC_Struct(filename,tree->mcmc);
  
  if(filename) tree->mcmc->out_fp = fopen(tree->mcmc->out_filename,"w");
  else tree->mcmc->out_fp = stderr;

  For(i,2*tree->n_otu-2) tree->rates->nd_r[i] = tree->rates->true_r[i];
  For(i,2*tree->n_otu-2) tree->rates->nd_t[i] = tree->rates->true_t[i];
  RATES_Update_Cur_Bl(tree);

  MCMC_Print_Param(tree->mcmc->out_fp,tree);
  MCMC_Print_Param(stdout,tree);

  MCMC_Randomize_Rates(tree);
  /* MCMC_Randomize_Nu(tree); */
/*   MCMC_Randomize_Lexp(tree); */
/*   MCMC_Randomize_Jumps(tree); */
/*   MCMC_Randomize_Alpha(tree); */
  MCMC_Randomize_Node_Times(tree);

  RATES_Lk_Rates(tree);

  if(tree->rates->lk_approx == NORMAL)
    tree->c_lnL = Dnorm_Multi_Given_InvCov_Det(tree->rates->u_cur_l,tree->rates->u_ml_l,tree->rates->invcov,tree->rates->covdet,2*tree->n_otu-3,YES);
  else Lk(tree);

  tree->mcmc->sample_interval = 2*tree->n_otu-2;
  n_moves = 11;
  do
    {            

/*       MCMC_Lexp(tree);         */
/*       MCMC_Alpha(tree);        */
      /* MCMC_Nu(tree); */
      MCMC_Rates_Local(tree);
/*       MCMC_Rates_Global(tree); */
      MCMC_Times_Local(tree);
/*       MCMC_Times_Global(tree); */
/*       MCMC_Stick_Rates(tree);  */
/*       MCMC_Mixing_Step(tree);  */
/*       MCMC_Jumps_Local(tree);  */
    }
  while(tree->mcmc->run < 1E+6);

  fclose(tree->mcmc->out_fp);
  MCMC_Free_MCMC(tree->mcmc);
}

/*********************************************************/

void MCMC_Lexp(arbre *tree)
{
  phydbl cur_lexp,new_lexp;
  phydbl cur_lnL_rates,new_lnL_rates;
/*   phydbl cur_lnL_jps,new_lnL_jps; */
  phydbl u,alpha,prior_mean_lexp,ratio;
  
  if((tree->rates->model != COMPOUND_NOCOR) &&
     (tree->rates->model != COMPOUND_COR)) return;

  new_lnL_rates   = UNLIKELY;
  cur_lexp        = -1.0;
  new_lexp        = -1.0;
  prior_mean_lexp = 0.03;
  ratio           = -1.0;

  cur_lnL_rates = tree->rates->c_lnL;
/*   cur_lnL_jps   = RATES_Lk_Jumps(tree); */

  cur_lexp = tree->rates->lexp;
  
  if(cur_lexp > 2.0) printf("\n. cur_lexp = %f iter=%d",cur_lexp,tree->mcmc->run);

  u = Uni();
  new_lexp = cur_lexp * exp(H_MCMC_LEXP*(u-0.5));

  if((new_lexp  > 1.E-5) && (new_lexp  < 2.0))
    {
      tree->rates->lexp = new_lexp;
      
      new_lnL_rates = RATES_Lk_Rates(tree);
/*       new_lnL_jps   = RATES_Lk_Jumps(tree); */

/*       ratio = (new_lnL_rates + new_lnL_jps) - (cur_lnL_rates + cur_lnL_jps) + log(new_lexp/cur_lexp); */
      ratio = (new_lnL_rates) - (cur_lnL_rates) + log(new_lexp/cur_lexp);

      ratio = exp(ratio);

/*       ratio =  */
/* 	exp(new_lnL-cur_lnL)* */
/* 	(new_lexp/cur_lexp) * */
/* 	exp((1./prior_mean_lexp)*(cur_lexp-new_lexp)); */
      
      alpha = MIN(1.,ratio);
      
      u = Uni();
      if(u > alpha) /* Reject */
	{
	  tree->rates->lexp = cur_lexp;
	  RATES_Lk_Rates(tree);
	}
      else
	{
	  tree->mcmc->acc_lexp++;
	}
    }

  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc->out_fp,tree);
}

/*********************************************************/

void MCMC_Nu(arbre *tree)
{
  phydbl cur_nu,new_nu,cur_lnL,new_lnL;
  phydbl u,alpha,prior_mean_nu,ratio;
  
  if(
     (tree->rates->model == COMPOUND_COR)   || 
     (tree->rates->model == COMPOUND_NOCOR) ||
     (tree->rates->model == GAMMA) ||
     (tree->rates->model == EXPONENTIAL)
     ) return;

  cur_lnL       = UNLIKELY;
  new_lnL       = UNLIKELY;
  cur_nu        = -1.0;
  new_nu        = -1.0;
  prior_mean_nu =  1.0;
  ratio         = -1.0;

  cur_lnL = tree->rates->c_lnL;
  cur_nu  = tree->rates->nu;
  
  u = Uni();
  new_nu = cur_nu  * exp(H_MCMC_NU*(u-0.5));

  if((new_nu > 1.E-7) && (new_nu  < 1.E-0))
    {
      tree->rates->nu = new_nu;
      
      new_lnL = RATES_Lk_Rates(tree);
            
      ratio = (new_lnL-cur_lnL) + log(new_nu/cur_nu);

      ratio = exp(ratio);

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

  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc->out_fp,tree);

}

/*********************************************************/

void MCMC_Alpha(arbre *tree)
{
  phydbl cur_lnL,new_lnL,cur_alpha,new_alpha;
  phydbl u,alpha,ratio;
  
  if((tree->rates->model != COMPOUND_NOCOR) &&
     (tree->rates->model != COMPOUND_COR)   &&
     (tree->rates->model != GAMMA)) return;

  cur_lnL = UNLIKELY;
  new_lnL = UNLIKELY;
  ratio = -1.0;

  cur_lnL    = tree->rates->c_lnL;
  cur_alpha  = tree->rates->alpha;
  
  u =  Uni();
  new_alpha = cur_alpha + (u * 2.*cur_alpha/10. - cur_alpha/10.);

  if(new_alpha > 0.1 && new_alpha < 10.0)
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
  
  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc->out_fp,tree);

}

/*********************************************************/

void MCMC_Times_Local(arbre *tree)
{
  int local;
  int node_num;

  local = 1;
  
  node_num = Rand_Int(0,2*tree->n_otu-3);
  MCMC_Times_Pre(tree->noeud[node_num]->anc,tree->noeud[node_num],local,tree);

/*   MCMC_Times_Pre(tree->n_root,tree->n_root->v[0],local,tree); */
/*   MCMC_Times_Pre(tree->n_root,tree->n_root->v[1],local,tree); */
}

/*********************************************************/

void MCMC_Rates_Local(arbre *tree)
{
  int local;
  int node_num;

  local = 1;

  node_num = Rand_Int(0,2*tree->n_otu-3);
  MCMC_Rates_Pre(tree->noeud[node_num]->anc,tree->noeud[node_num],local,tree);

/*   MCMC_Rates_Pre(tree->n_root,tree->n_root->v[0],local,tree); */
/*   MCMC_Rates_Pre(tree->n_root,tree->n_root->v[1],local,tree); */

}

/*********************************************************/

void MCMC_Stick_Rates(arbre *tree)
{
  tree->both_sides = 1;
  Lk(tree);

  MCMC_Stick_Rates_Pre(tree->n_root,tree->n_root->v[0],tree);
  MCMC_Stick_Rates_Pre(tree->n_root,tree->n_root->v[1],tree);
}

/*********************************************************/

void MCMC_Rates_Pre(node *a, node *d, int local, arbre *tree)
{
  phydbl u;
  phydbl new_lnL_data, cur_lnL_data, new_lnL_rate, cur_lnL_rate;
  phydbl ratio, alpha;
  phydbl new_mu, cur_mu;
  int i;
  edge *b;

  b = NULL;
  
  cur_mu       = tree->rates->nd_r[d->num];
  cur_lnL_data = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL;

  if(a == tree->n_root) b = tree->e_root;
  else For(i,3) if(d->v[i] == a) {b = d->b[i]; break;}
  
  u = Uni();
  new_mu = cur_mu * exp(H_MCMC_RATES*(u-0.5));
  
  if(new_mu < tree->rates->min_rate)
    {
      new_mu = tree->rates->max_rate - fmod(new_mu,tree->rates->max_rate-tree->rates->min_rate);
    }
  if(new_mu > tree->rates->max_rate)
    {
      new_mu = tree->rates->min_rate + fmod(new_mu,tree->rates->max_rate-tree->rates->min_rate);
    }
  if(new_mu < tree->rates->min_rate || new_mu > tree->rates->max_rate)
    {
      PhyML_Printf("\n Problem with setting new rate.\n");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }       


  if(local)
    {
      tree->rates->nd_r[d->num] = new_mu;
      RATES_Update_Cur_Bl(tree);

      if(tree->rates->lk_approx == NORMAL)
	new_lnL_data = Dnorm_Multi_Given_InvCov_Det(tree->rates->u_cur_l,tree->rates->u_ml_l,tree->rates->invcov,tree->rates->covdet,2*tree->n_otu-3,YES);
      else
	new_lnL_data = Return_Lk(tree);

      tree->c_lnL  = new_lnL_data;
      
      new_lnL_rate = RATES_Lk_Rates(tree);
	    
      ratio =
	(new_lnL_data + new_lnL_rate + log(new_mu)) -
	(cur_lnL_data + cur_lnL_rate + log(cur_mu));

/*       ratio = */
/*       	(new_lnL_rate + log(new_mu)) - */
/*       	(cur_lnL_rate + log(cur_mu)); */

      ratio = exp(ratio);	
      alpha = MIN(1.,ratio);
      
      u = Uni();
      
      if(u > alpha) /* Reject */
	{
	  tree->rates->nd_r[d->num] = cur_mu;
	  RATES_Update_Cur_Bl(tree);

	  if(tree->rates->lk_approx == NORMAL)
	    new_lnL_data = Dnorm_Multi_Given_InvCov_Det(tree->rates->u_cur_l,tree->rates->u_ml_l,tree->rates->invcov,tree->rates->covdet,2*tree->n_otu-3,YES);
	  else
	    new_lnL_data = Return_Lk(tree);

	  tree->c_lnL  = new_lnL_data;
	  
	  RATES_Lk_Rates(tree);

	  if((tree->mcmc->run > 10) && ((fabs(cur_lnL_data - tree->c_lnL) > 1.E-3) || (fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-0)))
	    {
	      printf("\n. Run=%d a=%d d=%d lexp=%f alpha=%f mu(d)=%f n(d)=%d mu(a)=%f n(a)=%d dt=(%f); cur_lnL_data = %f vs %f ; cur_lnL_rates = %f vs %f",
		     tree->mcmc->run,
		     d->anc->num,d->num,
		     tree->rates->lexp,
		     tree->rates->alpha,
		     tree->rates->nd_r[d->num],tree->rates->n_jps[d->num],
		     tree->rates->nd_r[d->anc->num],tree->rates->n_jps[d->anc->num],
		     fabs(tree->rates->nd_t[a->num] - tree->rates->nd_t[d->num]),
		     cur_lnL_data,tree->c_lnL,
		     cur_lnL_rate,tree->rates->c_lnL);
	      
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	  
	  if(fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-3)
	    {
	      printf("\n. WARNING: numerical precision issue detected (diff=%G). Reseting the likelihood.\n",cur_lnL_rate - tree->rates->c_lnL);

	      if(tree->rates->lk_approx == NORMAL)
		tree->c_lnL = Dnorm_Multi_Given_InvCov_Det(tree->rates->u_cur_l,tree->rates->u_ml_l,tree->rates->invcov,tree->rates->covdet,2*tree->n_otu-3,YES);
	      else
		tree->c_lnL = Return_Lk(tree);

	      RATES_Lk_Rates(tree);
	    }
	}
      else
	{
	  tree->mcmc->acc_rates++;
	}
    }
      
  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc->out_fp,tree);

/*   if(d->tax) return; */
/*   else */
/*     { */
/*       For(i,3) */
/* 	{ */
/* 	  if((d->v[i] != a) && (d->b[i] != tree->e_root)) */
/* 	    { */
/* /\* 	      Update_P_Lk(tree,d->b[i],d); *\/ */
/* 	      MCMC_Rates_Pre(d,d->v[i],local,tree); */
/* 	    } */
/* 	} */
/* /\*       if(a != tree->n_root) { Update_P_Lk(tree,b,d); } *\/ */
/* /\*       else                  { Update_P_Lk(tree,tree->e_root,d); } *\/ */
/*     } */
}

/*********************************************************/

void MCMC_Times_Pre(node *a, node *d, int local, arbre *tree)
{
  phydbl u;
  phydbl t_min,t_max;
  phydbl cur_t, new_t;
  phydbl cur_lnL_times, new_lnL_times;
  phydbl cur_lnL_rate, new_lnL_rate;
  phydbl cur_lnL_data, new_lnL_data;
  phydbl ratio,alpha;
  int    i;
  phydbl t0,t1,t2,t3;
  node *v2,*v3, *buff_n;
  phydbl u0,u1,u2,u3;

  if(d->tax) return; /* Won't change time at tip */

  RATES_Record_Times(tree);
  RATES_Record_Rates(tree);

  cur_lnL_data  = tree->c_lnL;
  cur_lnL_rate  = tree->rates->c_lnL;
  cur_t         = tree->rates->nd_t[d->num];
  cur_lnL_times = RATES_Yule(tree);
  new_lnL_data  = cur_lnL_data;
  new_lnL_times = cur_lnL_times;

  buff_n = v2 = v3 = NULL;
  For(i,3)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      {
	if(!v2) { v2 = d->v[i]; }
	else    { v3 = d->v[i]; }
      }
  
  t2 = tree->rates->nd_t[v2->num];
  t3 = tree->rates->nd_t[v3->num];

  if(t3 > t2)
    {
      buff_n = v2;
      v2     = v3;
      v3     = buff_n;
    }
  
  t0 = tree->rates->nd_t[a->num];
  t1 = tree->rates->nd_t[d->num];
  t2 = tree->rates->nd_t[v2->num];
  t3 = tree->rates->nd_t[v3->num];
  u0 = tree->rates->nd_r[a->num]  * tree->rates->clock_r;
  u1 = tree->rates->nd_r[d->num]  * tree->rates->clock_r;
  u2 = tree->rates->nd_r[v2->num] * tree->rates->clock_r;
  u3 = tree->rates->nd_r[v3->num] * tree->rates->clock_r;

  t_min = t0 + (1./tree->rates->nu)*pow((u1-u0)/1.96,2);
  t_max = t3 - (1./tree->rates->nu)*pow((u1-u3)/1.96,2);
  /* t_min = t0; */
  /* t_max = t3; */
  
  if(t_max < t_min)
    {
      PhyML_Printf("\n. t_max = %f t_min=%f",t_max,t_min);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }


  u = Uni();
    
  new_t = u*(t_max-t_min)+t_min;  

  if(new_t < t_min)
    {
      new_t = t_max - fmod(fabs(new_t),fabs(t_max-t_min));
    }
  if(new_t > t_max)
    {
      new_t = t_min + fmod(fabs(new_t),fabs(t_max-t_min));
    }
  if(new_t < t_min || new_t > t_max)
    {
      PhyML_Printf("\n. Problem with setting new time.\n");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }


  tree->rates->nd_t[d->num] = new_t;
  
  if(local)
    {
      new_lnL_rate  = RATES_Lk_Rates(tree);

      if(tree->rates->lk_approx == NORMAL)
	{
	  RATES_Update_Cur_Bl(tree);
	  new_lnL_data = Dnorm_Multi_Given_InvCov_Det(tree->rates->u_cur_l,tree->rates->u_ml_l,tree->rates->invcov,tree->rates->covdet,2*tree->n_otu-3,YES);
	}
      else
	{
	  new_lnL_data = Return_Lk(tree);
	}

      tree->c_lnL = new_lnL_data;
      
      /* ratio = */
      /* 	(new_lnL_data + new_lnL_rate) - */
      /* 	(cur_lnL_data + cur_lnL_rate); */
      
      ratio = new_lnL_data - cur_lnL_data;
      
      ratio = exp(ratio);
      
      alpha = MIN(1.,ratio);
      
      u = Uni();
      
      
      if(u > alpha) /* Reject */
	{
	  RATES_Reset_Times(tree);
	  RATES_Reset_Rates(tree);
	  
	  RATES_Lk_Rates(tree);

	  if(tree->rates->lk_approx == NORMAL)
	    {
	      RATES_Update_Cur_Bl(tree);
	      new_lnL_data = Dnorm_Multi_Given_InvCov_Det(tree->rates->u_cur_l,tree->rates->u_ml_l,tree->rates->invcov,tree->rates->covdet,2*tree->n_otu-3,YES);
	    }
	  else
	    {
	      new_lnL_data = Return_Lk(tree);
	    }

	  tree->c_lnL = new_lnL_data;
	  
	  if((tree->mcmc->run > 10) && ((fabs(cur_lnL_data - tree->c_lnL) > 1.E-3) || (fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-0)))
	    {
	      printf("\n. Run=%d",tree->mcmc->run);
	      
	      printf("\n. lexp = %f alpha = %f",
		     tree->rates->lexp,
		     tree->rates->alpha);
	      	      
	      printf("\n. cur_t = %f ; cur_lnL_data = %f vs %f ; cur_lnL_rates = %f vs %f",
		     cur_t,
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
	}
    }

  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc->out_fp,tree);
      
/*   if(d->tax) return; */
/*   else */
/*     { */
/*       For(i,3) */
/* 	{ */
/* 	  if((d->v[i] != a) && (d->b[i] != tree->e_root)) */
/* 	    { */
/* 	      MCMC_Times_Pre(d,d->v[i],local,tree); */
/* 	    } */
/* 	} */
/*     } */
}

/*********************************************************/

void MCMC_Jumps_Local(arbre *tree)
{
  int local;
  local = 1;
  MCMC_Jumps_Pre(tree->n_root,tree->n_root->v[0],local,tree);
  MCMC_Jumps_Pre(tree->n_root,tree->n_root->v[1],local,tree);
}

/*********************************************************/

void MCMC_Jumps_Pre(node *a, node *d, int local, arbre *tree)
{
  phydbl u;
  phydbl new_lnL_rate, cur_lnL_rate;
/*   phydbl new_lnL_jps, cur_lnL_jps; */
  phydbl ratio, alpha;
  int new_jps, cur_jps;
  int i;

  cur_jps      = tree->rates->n_jps[d->num];
  cur_lnL_rate = tree->rates->c_lnL;
/*   cur_lnL_jps  = RATES_Lk_Jumps(tree); */
  
  u = Uni();

  new_jps = cur_jps + (int)((u-0.5)*4.);

  if(local && new_jps > -1 && new_jps < 30)
    {
      tree->rates->n_jps[d->num] = new_jps;

      new_lnL_rate = RATES_Lk_Rates(tree);
/*       new_lnL_jps  = RATES_Lk_Jumps(tree); */

      ratio = (new_lnL_rate) - (cur_lnL_rate);
/*       ratio = (new_lnL_rate + new_lnL_jps) - (cur_lnL_rate + cur_lnL_jps); */
      ratio = exp(ratio);
      alpha = MIN(1.,ratio);
      
      u = Uni();
      
      if(u > alpha) /* Reject */
	{
	  tree->rates->n_jps[d->num] = cur_jps;

	  RATES_Lk_Rates(tree);
	  
	  if(fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-0)
	    {
	      printf("\n. lexp=%f alpha=%f dt=(%f); cur_lnL_rates = %f vs %f",
		     tree->rates->lexp,
		     tree->rates->alpha,
		     fabs(tree->rates->nd_t[a->num] - tree->rates->nd_t[d->num]),
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
    }

  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc->out_fp,tree);

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      MCMC_Jumps_Pre(d,d->v[i],local,tree);
	    }
	}
    }
}

/*********************************************************/

void MCMC_Stick_Rates_Pre(node *a, node *d, arbre *tree)
{
  phydbl u;
  phydbl new_lnL_data, cur_lnL_data, new_lnL_rate, cur_lnL_rate;
  phydbl dta,dtd;
  phydbl ratio, alpha, hr;
  phydbl new_mu, cur_mu;
  int i;
  edge *b;

  b = NULL;

  if(a != tree->n_root)
    {      
      cur_mu       = tree->rates->nd_r[d->num];
      cur_lnL_data = tree->c_lnL;
      cur_lnL_rate = tree->rates->c_lnL;
      
      dta = fabs(tree->rates->nd_t[a->num] - tree->rates->nd_t[a->anc->num]);
      dtd = fabs(tree->rates->nd_t[d->num] - tree->rates->nd_t[d->anc->num]);

      For(i,3) if(d->v[i] == a) {b = d->b[i]; break;}

      u = Uni();
      
      if(u < exp(-tree->rates->lexp * (dta+dtd)))
	{
	  new_mu = tree->rates->nd_r[a->num];
      
/* 	  new_lnL_rate = RATES_Lk_Rates(tree); */
	  new_lnL_rate = RATES_Lk_Change_One_Rate(d,new_mu,tree);
	  new_lnL_data = Lk_At_Given_Edge(b,tree);
	  
	  hr = (1. - exp(-tree->rates->lexp * (dta+dtd))) / exp(-tree->rates->lexp * (dta+dtd));
	  	  
	  ratio =
	    (new_lnL_data + new_lnL_rate) -
	    (cur_lnL_data + cur_lnL_rate) +
	    log(hr);
	  
	  ratio = exp(ratio);	
	  alpha = MIN(1.,ratio);
	  
	  u = Uni();
	  
	  if(u > alpha) /* Reject */
	    {
/* 	      RATES_Lk_Rates(tree); */
	      RATES_Lk_Change_One_Rate(d,cur_mu,tree);
	      Lk_At_Given_Edge(b,tree);
	      
	      if((fabs(cur_lnL_data - tree->c_lnL) > 1.E-3) || (fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-0))
		{
		  printf("\n. lexp=%f alpha=%f b->l = (%f) dt=(%f); cur_lnL_data = %f vs %f ; cur_lnL_rates = %f vs %f",
			 tree->rates->lexp,
			 tree->rates->alpha,
			 b->l,
			 fabs(tree->rates->nd_t[a->num] - tree->rates->nd_t[d->num]),
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
	}
    }

  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc->out_fp,tree);

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      Update_P_Lk(tree,d->b[i],d);
	      MCMC_Stick_Rates_Pre(d,d->v[i],tree);
	    }
	}
      if(a != tree->n_root) { Update_P_Lk(tree,b,d); }
      else                  { Update_P_Lk(tree,tree->e_root,d); }
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
	  fprintf(fp,"TreeSize\t");
	  fprintf(fp,"LnLSeq\t");
	  fprintf(fp,"LnLRate\t");
	  fprintf(fp,"RootPos[%f]\t",tree->n_root_pos);
	  fprintf(fp,"Nu\t");

	  if(fp != stdout) 
	    {
	      for(i=tree->n_otu;i<2*tree->n_otu-1;i++)
		{
		  if(i == tree->n_root->num)
		    {
		      fprintf(fp,"XXT%d[%f]\t",i,tree->rates->true_t[i]);
		    }
		  else
		    {
		      phydbl min_t = +1E+5;
		      int j;
		      
		      For(j,3) 
			if(tree->noeud[i]->v[j] != tree->noeud[i]->anc && tree->noeud[i]->b[j] != tree->e_root)
			  {
			    if(tree->rates->true_t[tree->noeud[i]->v[j]->num] < min_t)
			      {
				min_t = tree->rates->true_t[tree->noeud[i]->v[j]->num];
			      }
			  }

		      if(tree->noeud[i] == tree->n_root->v[0] || tree->noeud[i] == tree->n_root->v[1])
			{
			  fprintf(fp,"**T%d[%.1f<%.1f<%.1f]\t",i,
				  -100.,
				  tree->rates->true_t[i],
				  min_t);
			}
		      else
			{
			  fprintf(fp,"  T%d[%.1f<%.1f<%.1f]\t",
				  i,
				  tree->rates->true_t[tree->noeud[i]->anc->num],
				  tree->rates->true_t[i],
				  min_t);
			}
		    }
		}
	    }

	  if(fp != stdout)
	    for(i=0;i<2*tree->n_otu-2;i++) 
	      if(
		 (tree->noeud[i] == tree->n_root->v[0] && tree->noeud[i]->tax) ||
		 (tree->noeud[i] == tree->n_root->v[1] && tree->noeud[i]->tax)
		 )
		fprintf(fp,"00R%d[%f]\t",i,tree->rates->true_r[i]);
	      else if((tree->noeud[i] == tree->n_root->v[0]) || (tree->noeud[i] == tree->n_root->v[1]))
		fprintf(fp,"11R%d[%f]\t",i,tree->rates->true_r[i]);
	      else fprintf(fp,"  R%d[%f]\t",i,tree->rates->true_r[i]);

	  if(fp != stdout) 
	    for(i=0;i<2*tree->n_otu-3;i++) 
	      {
		if(tree->t_edges[i] == tree->e_root) fprintf(fp,"**L%d[%f,%G]\t",i,tree->rates->u_ml_l[i],tree->rates->cov[i*(2*tree->n_otu-3)+i]);
		else fprintf(fp,"  L%d[%f,%G]\t",i,tree->rates->u_ml_l[i],tree->rates->cov[i*(2*tree->n_otu-3)+i]);
	      }	  

	  if(fp != stdout) fprintf(fp,"AccRate\t");
	}

      fprintf(fp,"\n");
      fprintf(fp,"%6d\t",tree->mcmc->run);
/*       fprintf(fp,"%4.2f\t",RATES_Check_Mean_Rates(tree)); */
      fprintf(fp,"%4.2f\t",Get_Tree_Size(tree));

      fprintf(fp,"%15lf\t",tree->c_lnL);

      fprintf(fp,"%15lf\t",tree->rates->c_lnL);

      fprintf(fp,"%15lf\t",tree->rates->cur_l[tree->n_root->v[0]->num] / tree->rates->u_cur_l[tree->e_root->num]-tree->n_root_pos);
      fprintf(fp,"%15lf\t",tree->rates->nu);
      if(fp != stdout) for(i=tree->n_otu;i<2*tree->n_otu-1;i++) fprintf(fp,"%8f\t",tree->rates->nd_t[i]-tree->rates->true_t[i]);
      if(fp != stdout) for(i=0;i<2*tree->n_otu-2;i++) fprintf(fp,"%8f\t",tree->rates->nd_r[i]);
      if(fp != stdout) for(i=0;i<2*tree->n_otu-3;i++) fprintf(fp,"%8f\t",tree->rates->u_cur_l[i]);
      if(fp != stdout) 
	{
	  if(tree->mcmc->run)
	    fprintf(fp,"%8f\t",(phydbl)tree->mcmc->acc_rates/tree->mcmc->run);
	  else
	    fprintf(fp,"%8f\t",0.0);
	}
      fflush(NULL);
    }
}

/*********************************************************/

tmcmc *MCMC_Make_MCMC_Struct(arbre *tree)
{
  tmcmc *mcmc;
  
  mcmc               = (tmcmc *)mCalloc(1,sizeof(tmcmc));
  mcmc->dt_prop      = (phydbl *)mCalloc(tree->n_otu-2,sizeof(phydbl));
  mcmc->p_no_jump    = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));
  mcmc->t_rate_jumps = (phydbl *)mCalloc(10*tree->n_otu,sizeof(phydbl));
  mcmc->t_rank       = (int *)mCalloc(tree->n_otu-1,sizeof(int));
  mcmc->r_path       = (phydbl *)mCalloc(tree->n_otu-2,sizeof(phydbl));
  mcmc->out_filename = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  return(mcmc);
}

/*********************************************************/

void MCMC_Free_MCMC(tmcmc *mcmc)
{
  Free(mcmc->dt_prop);
  Free(mcmc->p_no_jump);
  Free(mcmc->t_rate_jumps);
  Free(mcmc->t_rank);
  Free(mcmc->r_path);
  Free(mcmc->out_filename);
  Free(mcmc);
}

/*********************************************************/

void MCMC_Init_MCMC_Struct(char *filename, tmcmc *mcmc)
{
  int pid;

  mcmc->acc_lexp        = 0;
  mcmc->acc_rates       = 0;
  mcmc->acc_times       = 0;
  mcmc->acc_nu          = 0;
  mcmc->run             = 0;
  mcmc->sample_interval = 100;  
  mcmc->n_rate_jumps    = 0;

  if(filename)
    {
      strcpy(mcmc->out_filename,filename);
      pid = getpid();
      sprintf(mcmc->out_filename+strlen(mcmc->out_filename),".%d",pid);
    }
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

void MCMC_Randomize_Rates(arbre *tree)
{
  int i;
  phydbl u;

  For(i,2*tree->n_otu-2)
    {
      u = Uni();
      tree->rates->nd_r[i] = -log(u);
      if(tree->rates->nd_r[i] < tree->rates->min_rate) tree->rates->nd_r[i] = tree->rates->min_rate; 
      if(tree->rates->nd_r[i] > tree->rates->max_rate) tree->rates->nd_r[i] = tree->rates->max_rate; 
    }
}

/*********************************************************/

void MCMC_Randomize_Jumps(arbre *tree)
{
  int i;
  phydbl u;

  For(i,2*tree->n_otu-2)
    {
      u = Uni();
      u *= 4.;
      tree->rates->n_jps[i] = (int)(u)+1;
      if(tree->rates->n_jps[i] > 10) tree->rates->n_jps[i] = 10;
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
  while((tree->rates->lexp < 1.E-5) || (tree->rates->lexp > 2.0));
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
  while((tree->rates->nu < 1.E-5) || (tree->rates->nu > 1.E-1));
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
      
      t_inf = MIN(tree->rates->nd_t[v1->num],tree->rates->nd_t[v2->num]);
      t_sup = tree->rates->nd_t[a->num];
      
      u = Uni();
      u *= (t_inf - t_sup);
      u += t_sup;
      
      tree->rates->nd_t[d->num] = u;
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

void MCMC_Mixing_Step(arbre *tree)
{
  phydbl new_lnL_time, new_lnL_rate, new_lnL_data;
  phydbl cur_lnL_time, cur_lnL_rate, cur_lnL_data;
  phydbl ratio,alpha,u,hr;
  int i;
  phydbl multiplier;

  cur_lnL_rate = tree->rates->c_lnL;
  cur_lnL_time = RATES_Yule(tree);
  cur_lnL_data = tree->c_lnL;
  
  new_lnL_data = cur_lnL_data;

  RATES_Record_Times(tree);

  u = Uni();
  multiplier = 1.0*exp(u-0.5);

  For(i,2*tree->n_otu-1) tree->rates->nd_t[i] *= multiplier;
  tree->rates->clock_r /= multiplier;
  
  RATES_Update_Cur_Bl(tree);
  
  hr = pow(multiplier,-(tree->n_otu-1));
  
  new_lnL_rate = RATES_Lk_Rates(tree);
  new_lnL_time = RATES_Yule(tree);
  /* new_lnL_data = Return_Lk(tree); */
  
  ratio = (new_lnL_rate + new_lnL_time) - (cur_lnL_rate + cur_lnL_time) + log(hr);
  ratio = exp(ratio);
  alpha = MIN(1.,ratio);
  
  u = Uni();
  
  if(u > alpha) /* Reject */
    {
      RATES_Reset_Times(tree);
      tree->rates->clock_r *= multiplier;

      RATES_Lk_Rates(tree);
      
      if(fabs(cur_lnL_rate - tree->rates->c_lnL) > 1.E-0)
	{
	  printf("\n. lexp=%f alpha=%f ; cur_lnL_rates = %f vs %f",
		 tree->rates->lexp,
		 tree->rates->alpha,
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
/*       printf("\n. Accept global times"); */
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
