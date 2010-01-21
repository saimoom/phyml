/*

PhyML:  a program that  computes maximum likelihood phyLOGenies from
DNA or AA homoLOGous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/



#include "mcmc.h"

/*********************************************************/

void MCMC(t_tree *tree)
{
  MCMC_Print_Param(tree->mcmc,tree);

  if(tree->mcmc->randomize)
    {
/*       MCMC_Randomize_Nu(tree); */
      MCMC_Randomize_Node_Times(tree);
      MCMC_Randomize_Rates(tree);
/*       MCMC_Randomize_Clock_Rate(tree); /\* Clock Rate must be the last parameter to be randomized *\/ */
    }

  int i;
  	

  RATES_Update_Cur_Bl(tree);
  RATES_Lk_Rates(tree);
  Lk(tree);

  do
    {
      tree->rates->c_lnL = RATES_Lk_Rates(tree);      
      tree->c_lnL        = Lk(tree);

/*       MCMC_Nu(tree); */
/*       MCMC_Clock_Rate(tree); */
      MCMC_Times_Local(tree);
      MCMC_Rates_Local(tree);
    }
  while(tree->mcmc->run < tree->mcmc->n_tot_run);
}

/*********************************************************/

void MCMC_Nu(t_tree *tree)
{
  phydbl cur_nu,new_nu,cur_lnL_rate,new_lnL_rate,cur_lnL;
  phydbl u,alpha,prior_mean_nu,ratio;
  phydbl min_nu,max_nu;

  if(
     (tree->rates->model == COMPOUND_COR)   || 
     (tree->rates->model == COMPOUND_NOCOR) ||
     (tree->rates->model == GAMMA) ||
     (tree->rates->model == EXPONENTIAL)
     ) return;

  cur_nu        = -1.0;
  new_nu        = -1.0;
  prior_mean_nu =  1.0;
  ratio         = -1.0;

  cur_lnL      = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL;
  cur_nu       = tree->rates->nu;
  
  min_nu = tree->rates->min_nu;
  max_nu = tree->rates->max_nu;

  u = Uni();
/*   new_nu = u*(max_nu - min_nu) + min_nu; */
  new_nu = cur_nu * EXP(tree->mcmc->h_nu*(u-0.5));

/*  "Mirroring" messes up the sampling. Can be seen by checking the prior densities obtained */
/*   if(new_nu < min_nu) */
/*     { */
/*       new_nu = max_nu - fmod(new_nu,max_nu-min_nu); */
/*     } */
/*   if(new_nu > max_nu) */
/*     { */
/*       new_nu = min_nu + fmod(new_nu,max_nu-min_nu); */
/*     } */
  if(new_nu > max_nu || new_nu < min_nu)
    {      
      return;
      PhyML_Printf("\n. max_nu = %G min_nu = %G",max_nu,min_nu);
      PhyML_Printf("\n. Problem with setting autocorrelation parameter.\n");
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  tree->rates->nu = new_nu;

  new_lnL_rate = RATES_Lk_Rates(tree);

  ratio =
/*     0.0; */
    (new_lnL_rate - cur_lnL_rate) +
    (LOG(new_nu)  - LOG(cur_nu));


  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  u = Uni();
  if(u > alpha) /* Reject */
    {
      tree->rates->nu = cur_nu;
      RATES_Lk_Rates(tree);

      tree->rates->nu    = cur_nu;
      tree->rates->c_lnL = cur_lnL_rate;
    }

  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc,tree);

  if(!(tree->mcmc->run%tree->mcmc->norm_freq))
    {
      RATES_Normalise_Rates(tree);
      RATES_Lk_Rates(tree);
    }
}

/*********************************************************/

void MCMC_Clock_Rate(t_tree *tree)
{
  phydbl cur_cr,new_cr,cur_lnL,new_lnL,cur_lnL_rate;
  phydbl u,alpha,prior_mean_cr,ratio;
  
  cur_lnL       = UNLIKELY;
  new_lnL       = UNLIKELY;
  cur_cr        = -1.0;
  new_cr        = -1.0;
  prior_mean_cr =  1.0;
  ratio         = -1.0;

  cur_lnL      = tree->c_lnL;
  cur_cr       = tree->rates->clock_r;  
  cur_lnL_rate = tree->rates->c_lnL;


  u = Uni();
  new_cr = cur_cr * EXP(tree->mcmc->h_clock*(u-0.5));
/*   new_cr = tree->rates->min_clock + u*(tree->rates->max_clock - tree->rates->min_clock); */

/*  "Mirroring" messes up the sampling. Can be seen by checking the prior densities obtained */
/*   if(new_cr < tree->rates->min_clock) */
/*     { */
/*       new_cr = tree->rates->max_clock - fmod(new_cr,tree->rates->max_clock-tree->rates->min_clock); */
/*     } */
/*   if(new_cr > tree->rates->max_clock) */
/*     { */
/*       new_cr = tree->rates->min_clock + fmod(new_cr,tree->rates->max_clock-tree->rates->min_clock); */
/*     } */
  if(new_cr > tree->rates->max_clock || new_cr < tree->rates->min_clock)
    {
      return;
      PhyML_Printf("\n Problem with setting new clock rate.\n");
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  tree->rates->clock_r = new_cr;
  RATES_Update_Cur_Bl(tree);

  new_lnL = Lk(tree);

  ratio =
/*     0.0; */
    (new_lnL - cur_lnL) +
    (LOG(new_cr) - LOG(cur_cr));


  ratio = EXP(ratio);
  
  alpha = MIN(1.,ratio);
  u = Uni();
  if(u > alpha) /* Reject */
    {
      tree->rates->clock_r = cur_cr;
      RATES_Update_Cur_Bl(tree);
      tree->c_lnL = cur_lnL;      
    }

  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc,tree);
  if(!(tree->mcmc->run%tree->mcmc->norm_freq)) 
    {
      RATES_Normalise_Rates(tree);
      RATES_Lk_Rates(tree);
    }
}

/*********************************************************/

void MCMC_Times_Local(t_tree *tree)
{
  int local;
  int node_num;
  int i;

  local = 1;
  
  For(i,2*tree->n_otu-1)
    {
      node_num = Rand_Int(tree->n_otu,2*tree->n_otu-2);

      if(tree->noeud[node_num] != tree->n_root)
	MCMC_Times_Pre(tree->noeud[node_num]->anc,tree->noeud[node_num],local,tree);
      else
      	MCMC_Time_Root(tree);
    }

  RATES_Lk_Rates(tree); /* Rates likelihood needs to be updated here */

/*   MCMC_Times_Pre(tree->n_root,tree->n_root->v[0],local,tree); */
/*   MCMC_Times_Pre(tree->n_root,tree->n_root->v[1],local,tree); */
}

/*********************************************************/

void MCMC_Rates_Local(t_tree *tree)
{
  int node_num;
  int i;

  
  For(i,2*tree->n_otu-2)
    {
      node_num = Rand_Int(0,2*tree->n_otu-3);
      MCMC_One_Rate(tree->noeud[node_num]->anc,tree->noeud[node_num],tree);
    }
}

/*********************************************************/

void MCMC_Stick_Rates(t_tree *tree)
{
  tree->both_sides = 1;
  Lk(tree);

  MCMC_Stick_Rates_Pre(tree->n_root,tree->n_root->v[0],tree);
  MCMC_Stick_Rates_Pre(tree->n_root,tree->n_root->v[1],tree);
}

/*********************************************************/

void MCMC_One_Rate(t_node *a, t_node *d, t_tree *tree)
{
  phydbl u;
  phydbl new_lnL_data, cur_lnL_data, new_lnL_rate, cur_lnL_rate;
  phydbl ratio, alpha;
  phydbl new_mu, cur_mu;
  phydbl r_min, r_max;
  phydbl nu;
  t_node *v2,*v3;
  phydbl T0,T1,T2,T3;
  phydbl U0,U1,U2,U3;
  phydbl V1,V2,V3;
  int i,err;
 
/*   if(a == tree->n_root) { tree->rates->nd_r[d->num] = 1.0; return; } */

  nu           = tree->rates->nu;
  cur_mu       = tree->rates->nd_r[d->num];
  cur_lnL_data = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL;
  r_min        = tree->rates->min_rate;
  r_max        = tree->rates->max_rate;


  
  u = Uni();

/*   new_mu = cur_mu * EXP(tree->mcmc->h_rates*(u-0.5)); */
/*   if(new_mu < r_min || new_mu > r_max) return; */

  new_mu = Uni() * (r_max - r_min) + r_min;  


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


  V1 = tree->rates->nu * (T1 - T0);
  V2 = tree->rates->nu * (T2 - T1);
  V3 = tree->rates->nu * (T3 - T1);

  if(!d->tax)
    {
      cur_lnL_rate = 
	Log_Dnorm_Trunc(U1,U0,SQRT(V1),r_min,r_max,&err) +
	Log_Dnorm_Trunc(U2,U1,SQRT(V2),r_min,r_max,&err) +
	Log_Dnorm_Trunc(U3,U1,SQRT(V3),r_min,r_max,&err);
      
      new_lnL_rate = 
	Log_Dnorm_Trunc(new_mu,U0,SQRT(V1),r_min,r_max,&err) +
	Log_Dnorm_Trunc(U2,new_mu,SQRT(V2),r_min,r_max,&err) +
	Log_Dnorm_Trunc(U3,new_mu,SQRT(V3),r_min,r_max,&err);
    }
  else
    {
      cur_lnL_rate = Log_Dnorm_Trunc(U1,U0,SQRT(V1),r_min,r_max,&err);
      new_lnL_rate = Log_Dnorm_Trunc(new_mu,U0,SQRT(V1),r_min,r_max,&err);
    }


  tree->rates->nd_r[d->num] = new_mu;  
  RATES_Update_Cur_Bl(tree);


/*   new_lnL_data = Lk(tree); */
/*   new_lnL_rate = RATES_Lk_Rates(tree); */



  ratio =
/*     (new_lnL_data - cur_lnL_data) + */
    (new_lnL_rate - cur_lnL_rate) ;
/*     (LOG(new_mu) - LOG(cur_mu)); */

/*   if(!(tree->mcmc->run%100000)) */
/*     printf("\n. %s red=%20f",(d->tax)?("EXT"):("INT"),(new_lnL_rate - cur_lnL_rate)); */

/* /\* !!!!!!!!!!!!!!!!!!!!!! *\/ */
/*   new_lnL_rate = RATES_Lk_Rates(tree); */

/*   if(!(tree->mcmc->run%100000)) */
/*     printf("\tfull=%20f new_mu=%20f",(new_lnL_rate - tmp),new_mu); */


  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  u = Uni();
  
  if(u > alpha) /* Reject */
    {
      tree->rates->nd_r[d->num] = cur_mu;
      RATES_Update_Cur_Bl(tree);
      tree->c_lnL        = cur_lnL_data;
      tree->rates->c_lnL = cur_lnL_rate;
    }
  else
    {
      tree->mcmc->acc_rates++;
    }

  /* !!!!!!!!!!!!!!!!!!!!!! */
  tree->rates->c_lnL = RATES_Lk_Rates(tree);

  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc,tree);
  if(!(tree->mcmc->run%tree->mcmc->norm_freq))
    {
      RATES_Normalise_Rates(tree);
      RATES_Lk_Rates(tree);
    }
}

/*********************************************************/

void MCMC_Times_Pre(t_node *a, t_node *d, int local, t_tree *tree)
{
  phydbl u;
  phydbl t_min,t_max;
/*   phydbl t_max_12, t_max_13; */
  phydbl t1_cur, t1_new;
  phydbl cur_lnL_times, new_lnL_times;
  phydbl cur_lnL_rate, new_lnL_rate;
  phydbl cur_lnL_data, new_lnL_data;
  phydbl ratio,alpha;
  int    i;
  phydbl t0,t1,t2,t3;
  t_node *v2,*v3, *buff_n;
  phydbl u0,u1,u2,u3;


  if(d->tax) return; /* Won't change time at tip */

  if(FABS(tree->rates->t_prior_min[d->num] - tree->rates->t_prior_max[d->num]) < 1.E-10) return;

  RATES_Record_Times(tree);
  RATES_Record_Rates(tree);

  cur_lnL_data  = tree->c_lnL;
  cur_lnL_rate  = tree->rates->c_lnL;
  t1_cur        = tree->rates->nd_t[d->num];
  cur_lnL_times = RATES_Yule(tree);
  
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
  u0 = tree->rates->nd_r[a->num];
  u1 = tree->rates->nd_r[d->num];
  u2 = tree->rates->nd_r[v2->num];
  u3 = tree->rates->nd_r[v3->num];

  t_min = t0;
  t_max = MIN(t2,t3);

  t_min = MAX(t_min,tree->rates->t_prior_min[d->num]);
  t_max = MIN(t_max,tree->rates->t_prior_max[d->num]);

  t_min += tree->rates->min_dt;
  t_max -= tree->rates->min_dt;

  if(t_min > t_max) return;

  tree->rates->t_prior[d->num] = Uni() * (t_max - t_min) + t_min;
  

  u = Uni();

  t1_new = u*(t_max-t_min)+t_min;
/*   t1_new = t1_cur * EXP(tree->mcmc->h_times*(u-0.5)); */

/*  "Mirroring" messes up the sampling. Can be seen by checking the prior densities obtained */
/*   if(t1_new < t_min) */
/*     { */
/*       t1_new = t_max - fmod(FABS(t1_new),FABS(t_max-t_min)); */
/*     } */

/*   if(t1_new > t_max) */
/*     { */
/*       t1_new = t_min + fmod(FABS(t1_new),FABS(t_max-t_min)); */
/*     } */

  if(t1_new < t_min || t1_new > t_max)
    {
      return;
      PhyML_Printf("\n. Problem with setting new time.\n");
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  
  if(local)
    {
      tree->rates->nd_t[d->num] = t1_new;

      RATES_Update_Cur_Bl(tree);
      new_lnL_data = Lk(tree);

      ratio =
	0.0;
/* 	(new_lnL_data - cur_lnL_data) + */
/* 	(LOG(FABS(t1_new)) - LOG(FABS(t1_cur))); */


      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      u = Uni();
      
      if(u > alpha) /* Reject */
	{
	  RATES_Reset_Times(tree);
 	  RATES_Update_Cur_Bl(tree);
	  tree->c_lnL        = cur_lnL_data;
	  tree->rates->c_lnL = cur_lnL_rate;
	}
      else
	{
	  tree->mcmc->acc_times++;
	}
    }

  if(t1_new < t0)
    {
      t1_new = t0+1.E-4;
      PhyML_Printf("\n");
      PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
      PhyML_Printf("\n. t0 = %f t1_new = %f t1 = %f",t0,t1_new,t1);
      PhyML_Printf("\n. t_min=%f t_max=%f",t_min,t_max);
      PhyML_Printf("\n. (t1-t0)=%f (t2-t1)=%f",t1-t0,t2-t1);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
/*       Exit("\n"); */
    }
  if(t1_new > MIN(t2,t3))
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
      PhyML_Printf("\n. t0 = %f t1_new = %f t1 = %f t2 = %f t3 = %f MIN(t2,t3)=%f",t0,t1_new,t1,t2,t3,MIN(t2,t3));
      PhyML_Printf("\n. t_min=%f t_max=%f",t_min,t_max);
      PhyML_Printf("\n. (t1-t0)=%f (t2-t1)=%f",t1-t0,t2-t1);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
/*       Exit("\n"); */
    }

  if(isnan(t1_new))
    {
      PhyML_Printf("\n. run=%d",tree->mcmc->run);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
/*       Exit("\n"); */
    }
  
  /* !!!!!!!!!!!!!!!!!!!!!! */
  tree->rates->c_lnL = RATES_Lk_Rates(tree);

  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc,tree);
  if(!(tree->mcmc->run%tree->mcmc->norm_freq))
    {
      RATES_Normalise_Rates(tree);
      RATES_Lk_Rates(tree);
    }
}

/*********************************************************/

void MCMC_Time_Root(t_tree *tree)
{
  phydbl u;
  phydbl t_min,t_max;
/*   phydbl t_max_02, t_max_03; */
  phydbl cur_t, new_t;
  phydbl cur_lnL_times, new_lnL_times;
  phydbl cur_lnL_rate, new_lnL_rate;
  phydbl cur_lnL_data, new_lnL_data;
  phydbl ratio,alpha;
  phydbl t0,t2,t3;
  t_node *v2,*v3, *buff_n;
  phydbl u0,u2,u3;
  t_node *root;

  RATES_Record_Times(tree);
  
  root = tree->n_root;

  if(FABS(tree->rates->t_prior_min[root->num] - tree->rates->t_prior_max[root->num]) < 1.E-10) return;

  cur_lnL_data  = tree->c_lnL;
  cur_lnL_rate  = tree->rates->c_lnL;
  cur_t         = tree->rates->nd_t[root->num];
  cur_lnL_times = RATES_Yule(tree);
  new_lnL_data  = cur_lnL_data;
  new_lnL_times = cur_lnL_times;
  new_lnL_rate  = cur_lnL_rate;

  v2 = root->v[0];
  v3 = root->v[1];

  t2 = tree->rates->nd_t[v2->num];
  t3 = tree->rates->nd_t[v3->num];

  if(t3 > t2)
    {
      buff_n = v2;
      v2     = v3;
      v3     = buff_n;
    }
  
  t0 = tree->rates->nd_t[root->num];
  t2 = tree->rates->nd_t[v2->num];
  t3 = tree->rates->nd_t[v3->num];
  u2 = tree->rates->nd_r[v2->num];
  u3 = tree->rates->nd_r[v3->num];
  u0 = 1.0;

  t_min = tree->rates->t_prior_min[root->num];

  t_max = MIN(t2,t3);

  t_min = MAX(t_min,tree->rates->t_prior_min[root->num]);
  t_max = MIN(t_max,tree->rates->t_prior_max[root->num]);
  
  if(FABS(t_min - t_max) < 1.E-10) return;

  tree->rates->t_prior[root->num] = Uni() * (t_max - t_min) + t_min;
  
  if(t_max < t_min) return;
  
  u0 *= tree->rates->clock_r;
  u2 *= tree->rates->clock_r;
  u3 *= tree->rates->clock_r;

  u = Uni();
/*   new_t = u*(t_max-t_min)+t_min; */
  new_t = cur_t * EXP(tree->mcmc->h_times*(u-0.5));

/*  "Mirroring" messes up the sampling. Can be seen by checking the prior densities obtained */
/*   if(new_t < t_min) */
/*     { */
/*       new_t = t_max - fmod(FABS(new_t),FABS(t_max-t_min)); */
/*     } */
/*   if(new_t > t_max) */
/*     { */
/*       new_t = t_min + fmod(FABS(new_t),FABS(t_max-t_min)); */
/*     } */
  if(new_t < t_min || new_t > t_max)
    {
      return;
      PhyML_Printf("\n. Problem with setting new time.\n");
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  tree->rates->nd_t[root->num] = new_t;
  
  RATES_Update_Cur_Bl(tree);
  
  new_lnL_data = Lk(tree);

  ratio = 
/*     0.0; */
    (new_lnL_data - cur_lnL_data) +
    (LOG(FABS(new_t)) - LOG(FABS(cur_t)));

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  if(u > alpha) /* Reject */
    {
      RATES_Reset_Times(tree);
      RATES_Update_Cur_Bl(tree);
      tree->c_lnL        = cur_lnL_data;
      tree->rates->c_lnL = cur_lnL_rate;
    }
  else
    {
      tree->mcmc->acc_times++;
    }
    
  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc,tree);

  if(!(tree->mcmc->run%tree->mcmc->norm_freq))
    {
      RATES_Normalise_Rates(tree);
      RATES_Lk_Rates(tree);
    }
}

/*********************************************************/

void MCMC_Jumps_Local(t_tree *tree)
{
  int local;
  local = 1;
  MCMC_Jumps_Pre(tree->n_root,tree->n_root->v[0],local,tree);
  MCMC_Jumps_Pre(tree->n_root,tree->n_root->v[1],local,tree);
}

/*********************************************************/

void MCMC_Jumps_Pre(t_node *a, t_node *d, int local, t_tree *tree)
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
      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      
      u = Uni();
      
      if(u > alpha) /* Reject */
	{
	  tree->rates->n_jps[d->num] = cur_jps;

	  RATES_Lk_Rates(tree);
	  
	  if(FABS(cur_lnL_rate - tree->rates->c_lnL) > 1.E-0)
	    {
	      PhyML_Printf("\n. lexp=%f alpha=%f dt=(%f); cur_lnL_rates = %f vs %f",
		     tree->rates->lexp,
		     tree->rates->alpha,
		     FABS(tree->rates->nd_t[a->num] - tree->rates->nd_t[d->num]),
		     cur_lnL_rate,tree->rates->c_lnL);
	      
	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	  
	  if(FABS(cur_lnL_rate - tree->rates->c_lnL) > 1.E-3)
	    {
	      PhyML_Printf("\n. WARNING: numerical precision issue detected (diff=%G). Reseting the likelihood.\n",cur_lnL_rate - tree->rates->c_lnL);
	      RATES_Lk_Rates(tree);
	    }
	}
    }

  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc,tree);

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

void MCMC_Stick_Rates_Pre(t_node *a, t_node *d, t_tree *tree)
{
  phydbl u;
  phydbl new_lnL_data, cur_lnL_data, new_lnL_rate, cur_lnL_rate;
  phydbl dta,dtd;
  phydbl ratio, alpha, hr;
  phydbl new_mu, cur_mu;
  int i;
  t_edge *b;

  b = NULL;

  if(a != tree->n_root)
    {      
      cur_mu       = tree->rates->nd_r[d->num];
      cur_lnL_data = tree->c_lnL;
      cur_lnL_rate = tree->rates->c_lnL;
      
      dta = FABS(tree->rates->nd_t[a->num] - tree->rates->nd_t[a->anc->num]);
      dtd = FABS(tree->rates->nd_t[d->num] - tree->rates->nd_t[d->anc->num]);

      For(i,3) if(d->v[i] == a) {b = d->b[i]; break;}

      u = Uni();
      
      if(u < EXP(-tree->rates->lexp * (dta+dtd)))
	{
	  new_mu = tree->rates->nd_r[a->num];
      
/* 	  new_lnL_rate = RATES_Lk_Rates(tree); */
	  new_lnL_rate = RATES_Lk_Change_One_Rate(d,new_mu,tree);
	  new_lnL_data = Lk_At_Given_Edge(b,tree);
	  
	  hr = (1. - EXP(-tree->rates->lexp * (dta+dtd))) / EXP(-tree->rates->lexp * (dta+dtd));
	  	  
	  ratio =
	    (new_lnL_data + new_lnL_rate) -
	    (cur_lnL_data + cur_lnL_rate) +
	    LOG(hr);
	  
	  ratio = EXP(ratio);	
	  alpha = MIN(1.,ratio);
	  
	  u = Uni();
	  
	  if(u > alpha) /* Reject */
	    {
/* 	      RATES_Lk_Rates(tree); */
	      RATES_Lk_Change_One_Rate(d,cur_mu,tree);
	      Lk_At_Given_Edge(b,tree);
	      
	      if((FABS(cur_lnL_data - tree->c_lnL) > 1.E-3) || (FABS(cur_lnL_rate - tree->rates->c_lnL) > 1.E-0))
		{
		  PhyML_Printf("\n. lexp=%f alpha=%f b->l = (%f) dt=(%f); cur_lnL_data = %f vs %f ; cur_lnL_rates = %f vs %f",
			 tree->rates->lexp,
			 tree->rates->alpha,
			 b->l,
			 FABS(tree->rates->nd_t[a->num] - tree->rates->nd_t[d->num]),
			 cur_lnL_data,tree->c_lnL,
			 cur_lnL_rate,tree->rates->c_lnL);
		  
		  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
		  Warn_And_Exit("");
		}
	      
	      if(FABS(cur_lnL_rate - tree->rates->c_lnL) > 1.E-3)
		{
		  PhyML_Printf("\n. WARNING: numerical precision issue detected (diff=%G). Reseting the likelihood.\n",cur_lnL_rate - tree->rates->c_lnL);
		  RATES_Lk_Rates(tree);
		  Lk(tree);
		}
	    }
	}
    }

  tree->mcmc->run++;
  MCMC_Print_Param(tree->mcmc,tree);

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


void MCMC_Print_Param(tmcmc *mcmc, t_tree *tree)
{
  int i;
  FILE *fp;

  fp = mcmc->out_fp_stats;
  
  if(!(mcmc->run%(mcmc->n_tot_run/10))) 
    { 
      PhyML_Fprintf(stdout,"."); 
      fflush(NULL); 
    }

  MCMC_Print_Means(mcmc,tree);
  MCMC_Print_Last(mcmc,tree);

/*   if(!(mcmc->run%10)) PhyML_Printf("\r. [%10d/%10d]",tree->mcmc->run,tree->mcmc->n_tot_run); */

  if(!(mcmc->run%mcmc->sample_interval)) 
    {      
      if(tree->mcmc->run == 0)
	{
	  time(&(mcmc->t_beg));

	  PhyML_Fprintf(fp,"\n");
	  PhyML_Fprintf(fp,"Run\t");
	  PhyML_Fprintf(fp,"Time\t");
	  PhyML_Fprintf(fp,"TreeSize\t");
	  PhyML_Fprintf(fp,"ExpTreeSize[%f]\t",RATES_Expected_Tree_Length(tree));
	  PhyML_Fprintf(fp,"LnLSeq\t");
	  PhyML_Fprintf(fp,"LnLRate\t");
/* 	  PhyML_Fprintf(fp,"RootPos[%f]\t",tree->n_root_pos); */
	  PhyML_Fprintf(fp,"Nu\t");
	  PhyML_Fprintf(fp,"Clock\t");

	  if(fp != stdout)
	    {
	      for(i=tree->n_otu;i<2*tree->n_otu-1;i++)
		{
		  PhyML_Fprintf(fp,"T%d[%f]\t",i,tree->rates->true_t[i]);
		}
	    }

/* 	  if(fp != stdout) */
/* 	    { */
/* 	      for(i=tree->n_otu;i<2*tree->n_otu-1;i++) */
/* 		{ */
/* 		  PhyML_Fprintf(fp,"pT%d\t",i); */
/* 		} */
/* 	    } */

	  if(fp != stdout)
	    for(i=0;i<2*tree->n_otu-2;i++)
	      if((tree->noeud[i] == tree->n_root->v[0] && tree->noeud[i]->tax) ||
		 (tree->noeud[i] == tree->n_root->v[1] && tree->noeud[i]->tax))
		PhyML_Fprintf(fp,"00R%d[%f]\t",i,tree->rates->true_r[i]);
	      else if((tree->noeud[i] == tree->n_root->v[0]) || (tree->noeud[i] == tree->n_root->v[1]))
		PhyML_Fprintf(fp,"11R%d[%f]\t",i,tree->rates->true_r[i]);
	      else PhyML_Fprintf(fp,"  R%d[%.2f,%.2f]\t",i,tree->rates->true_r[i],tree->rates->true_r[tree->noeud[i]->anc->num]);

	  if(fp != stdout)
	    for(i=0;i<2*tree->n_otu-3;i++)
	      {
		if(tree->t_edges[i] == tree->e_root)
		  PhyML_Fprintf(fp,"**L%d[%f,%d,%d]\t",i,
				tree->rates->u_ml_l[i],
				tree->t_edges[i]->left->num,
				tree->t_edges[i]->rght->num);
		else
		  PhyML_Fprintf(fp,"  L%d[%f,%d,%d]\t",i,
				tree->rates->u_ml_l[i],
				tree->t_edges[i]->left->num,
				tree->t_edges[i]->rght->num);
	      }

/* 	  if(fp != stdout) PhyML_Fprintf(fp,"AccRate\t"); */

	  PhyML_Fprintf(mcmc->out_fp_trees,"#NEXUS\n");
	  PhyML_Fprintf(mcmc->out_fp_trees,"BEGIN TREES;\n");

	}

      PhyML_Fprintf(fp,"\n");
      PhyML_Fprintf(fp,"%6d\t",tree->mcmc->run);
      time(&mcmc->t_cur);
      PhyML_Fprintf(fp,"%6d\t",(int)(mcmc->t_cur-mcmc->t_beg));
      
      RATES_Update_Cur_Bl(tree);
/*       PhyML_Fprintf(fp,"%4.2f\t",Get_Tree_Size(tree)-tree->rates->true_tree_size); */
      PhyML_Fprintf(fp,"%G\t",RATES_Check_Mean_Rates(tree));
      PhyML_Fprintf(fp,"%G\t",RATES_Expected_Tree_Length(tree));

      PhyML_Fprintf(fp,"%G\t",tree->c_lnL);

      PhyML_Fprintf(fp,"%G\t",tree->rates->c_lnL);

/*       PhyML_Fprintf(fp,"%15lf\t",tree->rates->cur_l[tree->n_root->v[0]->num] / tree->rates->u_cur_l[tree->e_root->num]-tree->n_root_pos); */
      PhyML_Fprintf(fp,"%G\t",tree->rates->nu);
      PhyML_Fprintf(fp,"%G\t",tree->rates->clock_r);
/*       if(fp != stdout) for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(fp,"%8f\t",tree->rates->t_prior[i] - tree->rates->true_t[i]); */
      if(fp != stdout) for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(fp,"%G\t",tree->rates->nd_t[i] - tree->rates->true_t[i]);
/*       if(fp != stdout) for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(fp,"%.1f\t",tree->rates->nd_t[i]); */
/*       if(fp != stdout) for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(fp,"%G\t",tree->rates->t_prior[i] - tree->rates->true_t[i]); */
 /*      if(fp != stdout) for(i=0;i<2*tree->n_otu-2;i++) PhyML_Fprintf(fp,"%.20lf\t",tree->rates->nd_r[i] - tree->rates->true_r[i]); */
      if(fp != stdout) for(i=0;i<2*tree->n_otu-2;i++) PhyML_Fprintf(fp,"%.20lf\t",tree->rates->nd_r[i]);
      if(fp != stdout) for(i=0;i<2*tree->n_otu-3;i++) PhyML_Fprintf(fp,"%G\t",tree->rates->u_cur_l[i]-tree->rates->u_ml_l[i]);
/*       if(fp != stdout)  */
/* 	{ */
/* 	  if(tree->mcmc->run) */
/* 	    PhyML_Fprintf(fp,"%8f\t",(phydbl)tree->mcmc->acc_rates/tree->mcmc->run); */
/* 	  else */
/* 	    PhyML_Fprintf(fp,"%8f\t",0.0); */
/* 	} */
      fflush(NULL);


      // TREES
      char *s;     
      Branch_Lengths_To_Time_Lengths(tree);
      s = Write_Tree(tree);
      PhyML_Fprintf(mcmc->out_fp_trees,"TREE %8d [%f] = [&R] %s\n",mcmc->run,tree->c_lnL,s);
      Free(s);
      RATES_Update_Cur_Bl(tree);
    }
}

/*********************************************************/

void MCMC_Print_Means(tmcmc *mcmc, t_tree *tree)
{

  if(!(mcmc->run%mcmc->sample_interval)) 
    {
      int i;
      char *s;

      s = (char *)mCalloc(T_MAX_FILE,sizeof(char));

      strcpy(s,tree->mcmc->out_filename);
      strcat(s,".means");
      
      fclose(mcmc->out_fp_means);

      mcmc->out_fp_means = fopen(s,"w");
      
      for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(mcmc->out_fp_means,"T%d\t",i);	  

      PhyML_Fprintf(mcmc->out_fp_means,"\n");      

      for(i=tree->n_otu;i<2*tree->n_otu-1;i++) tree->rates->t_mean[i] *= (phydbl)(mcmc->run / mcmc->sample_interval);

      for(i=tree->n_otu;i<2*tree->n_otu-1;i++)
	{
	  tree->rates->t_mean[i] += tree->rates->nd_t[i];
	  tree->rates->t_mean[i] /= (phydbl)(mcmc->run / mcmc->sample_interval + 1);

/* 	  PhyML_Fprintf(tree->mcmc->out_fp_means,"%d\t",tree->mcmc->run / tree->mcmc->sample_interval);	   */
	  PhyML_Fprintf(tree->mcmc->out_fp_means,"%.1f\t",tree->rates->t_mean[i]);
	}

      PhyML_Fprintf(tree->mcmc->out_fp_means,"\n");
      fflush(NULL);

      Free(s);
    }
}

/*********************************************************/

void MCMC_Print_Last(tmcmc *mcmc, t_tree *tree)
{

  if(!(mcmc->run%mcmc->sample_interval)) 
    {
      int i;
      char *s;

      s = (char *)mCalloc(T_MAX_FILE,sizeof(char));

      strcpy(s,tree->mcmc->out_filename);
      strcat(s,".lasts");
      
      fclose(mcmc->out_fp_last);

      mcmc->out_fp_last = fopen(s,"w");

/*       rewind(mcmc->out_fp_last); */

      PhyML_Fprintf(tree->mcmc->out_fp_last,"Time\t");

      for(i=tree->n_otu;i<2*tree->n_otu-1;i++)
	PhyML_Fprintf(tree->mcmc->out_fp_last,"T%d\t",i);

      PhyML_Fprintf(tree->mcmc->out_fp_last,"\n");

      if(mcmc->run)
	{
	  time(&(mcmc->t_cur));
	  PhyML_Fprintf(tree->mcmc->out_fp_last,"%d\t",(int)(mcmc->t_cur-mcmc->t_beg));
/* 	  PhyML_Fprintf(tree->mcmc->out_fp_last,"%d\t",(int)(mcmc->t_beg)); */
	}

      for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(tree->mcmc->out_fp_last,"%.1f\t",tree->rates->nd_t[i]);

      PhyML_Fprintf(tree->mcmc->out_fp_last,"\n");
      fflush(NULL);

      Free(s);
  }
}


/*********************************************************/
tmcmc *MCMC_Make_MCMC_Struct(t_tree *tree)
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

void MCMC_Init_MCMC_Struct(char *filename, tmcmc *mcmc, t_tree *tree)
{
  int pid;

  mcmc->acc_lexp        = 0;
  mcmc->acc_rates       = 0;
  mcmc->acc_times       = 0;
  mcmc->acc_nu          = 0;
  mcmc->run             = 0;
  mcmc->sample_interval = 10000;
  mcmc->n_rate_jumps    = 0;
  mcmc->n_tot_run       = 1.E+6;
  mcmc->randomize       = 1;
/*   mcmc->norm_freq       = 20*tree->n_otu; */
  mcmc->norm_freq       = 1E+9;

  mcmc->h_times         = 0.3;
  mcmc->h_rates         = 0.1;
  mcmc->h_nu            = 1.0;
  mcmc->h_clock         = 0.1;

  if(filename)
    {
      strcpy(mcmc->out_filename,filename);
      pid = getpid();
      sprintf(mcmc->out_filename+strlen(mcmc->out_filename),".%d",pid);
    }

  if(filename) 
    {
      char *s;

      s = (char *)mCalloc(T_MAX_NAME,sizeof(char));

      tree->mcmc->out_fp_stats = fopen(tree->mcmc->out_filename,"w");

      strcpy(s,tree->mcmc->out_filename);
      strcat(s,".trees");
      tree->mcmc->out_fp_trees = fopen(s,"w");

      strcpy(s,tree->mcmc->out_filename);
      strcat(s,".means");
      tree->mcmc->out_fp_means = fopen(s,"w");

      strcpy(s,tree->mcmc->out_filename);
      strcat(s,".lasts");
      tree->mcmc->out_fp_last  = fopen(s,"w");

      Free(s);
    }
  else 
    {
      tree->mcmc->out_fp_stats = stderr;
      tree->mcmc->out_fp_trees = stderr;
      tree->mcmc->out_fp_means = stderr;
      tree->mcmc->out_fp_last  = stderr;
    }
}

/*********************************************************/

void MCMC_Close_MCMC(tmcmc *mcmc)
{
  fclose(mcmc->out_fp_trees);
  fclose(mcmc->out_fp_stats);
  fclose(mcmc->out_fp_means);
  fclose(mcmc->out_fp_last);
}

/*********************************************************/

void MCMC_Randomize_Branch_Lengths(t_tree *tree)
{
  int i;
  phydbl u;

  For(i,2*tree->n_otu-3)
    {
      if(tree->t_edges[i] != tree->e_root)
	{
	  u = Uni();
	  tree->t_edges[i]->l *= -LOG(u);
	}
      else
	{
	  PhyML_Printf("\n. Didn't randomize root edge.");
	}
    }
}

/*********************************************************/

void MCMC_Randomize_Rates(t_tree *tree)
{

  /* Should be called once t_node times have been determined */

  int i;
  phydbl u;
  phydbl r_min,r_max;

  For(i,2*tree->n_otu-2) tree->rates->nd_r[i] = 1.0;

  r_min = 0.9;
  r_max = 1.1;
  u     = 0.0;

/*   For(i,2*tree->n_otu-2) */
/*     { */
/*       u = Uni(); */
/*       u = u * (r_max-r_min) + r_min; */
/*       tree->rates->nd_r[i] = u; */

/*       if(tree->rates->nd_r[i] < tree->rates->min_rate) tree->rates->nd_r[i] = tree->rates->min_rate;  */
/*       if(tree->rates->nd_r[i] > tree->rates->max_rate) tree->rates->nd_r[i] = tree->rates->max_rate;  */
/*     } */

  MCMC_Randomize_Rates_Pre(tree->n_root,tree->n_root->v[0],tree);
  MCMC_Randomize_Rates_Pre(tree->n_root,tree->n_root->v[1],tree);
}

/*********************************************************/

void MCMC_Randomize_Rates_Pre(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  phydbl mean_r, var_r;
  phydbl min_r, max_r;
  int err;

/*   mean_r = tree->rates->nd_r[a->num]; */
/*   var_r  = tree->rates->nu * (tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num]); */
  mean_r = 1.0;
  var_r  = 0.5;
  min_r  = tree->rates->min_rate;
  max_r  = tree->rates->max_rate;
  
  tree->rates->nd_r[d->num] = Rnorm_Trunc(mean_r,SQRT(var_r),min_r,max_r,&err);

  if(d->tax) return;
  else
    {
      For(i,3)
	if(d->v[i] != a && d->b[i] != tree->e_root)
	  MCMC_Randomize_Rates_Pre(d,d->v[i],tree);
    }
}

/*********************************************************/

void MCMC_Randomize_Jumps(t_tree *tree)
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

void MCMC_Randomize_LEXP(t_tree *tree)
{
  phydbl u;

  tree->rates->lexp = -1.0;
  do
    {
      u = Uni();
      tree->rates->lexp = -LOG(u);
    }
  while((tree->rates->lexp < 1.E-5) || (tree->rates->lexp > 2.0));
}

/*********************************************************/

void MCMC_Randomize_Nu(t_tree *tree)
{
  phydbl min_nu,max_nu;
  phydbl u;

  min_nu = tree->rates->min_nu;
  max_nu = tree->rates->max_nu;

  u = Uni();
  tree->rates->nu = (max_nu - min_nu) * u + min_nu;
}

/*********************************************************/

void MCMC_Randomize_Clock_Rate(t_tree *tree)
{
  phydbl u,t,t_anc;
  phydbl ml_size,cur_size;
  int i;


  ml_size = 0.0;
  For(i,2*tree->n_otu-3) ml_size += tree->rates->u_ml_l[i];


  cur_size = 0.0;
  For(i,2*tree->n_otu-2) 
    {
      u = tree->rates->nd_r[i];
      t = tree->rates->nd_t[i];
      t_anc = tree->rates->nd_t[tree->noeud[i]->anc->num];

      cur_size += u * (t - t_anc);
    }
  tree->rates->clock_r = ml_size / cur_size;
}

/*********************************************************/

void MCMC_Randomize_Alpha(t_tree *tree)
{
  phydbl u;

  u = Uni();
  tree->rates->alpha = u*6.0+1.0;
}

/*********************************************************/

void MCMC_Randomize_Node_Times(t_tree *tree)
{
  phydbl t_sup, t_inf;
  phydbl u;
  phydbl min_prior;
  int i;

  min_prior = 0.0;
  
  MCMC_Randomize_Node_Times_Bottom_Up(tree->n_root,tree->n_root->v[0],tree);
  MCMC_Randomize_Node_Times_Bottom_Up(tree->n_root,tree->n_root->v[1],tree);

  if(tree->rates->t_has_prior[tree->n_root->num])
    {
      t_inf = tree->rates->t_prior_min[tree->n_root->num];
      t_sup = tree->rates->t_prior_max[tree->n_root->num];
      
      u = Uni();
      u *= (t_sup - t_inf);
      u += t_inf;
      
      tree->rates->nd_t[tree->n_root->num] = u;
    }
  else
    {
      For(i,2*tree->n_otu-2)
	{
	  if(tree->rates->t_has_prior[i] && tree->rates->t_prior_min[i] < min_prior)
	    {
	      min_prior = tree->rates->t_prior_min[i];
	    }
	}
      tree->rates->nd_t[tree->n_root->num] = 10. * min_prior;
    }
  
  MCMC_Randomize_Node_Times_Top_Down(tree->n_root,tree->n_root->v[0],tree);
  MCMC_Randomize_Node_Times_Top_Down(tree->n_root,tree->n_root->v[1],tree);

  For(i,2*tree->n_otu-1)
    {
      if(!tree->rates->t_has_prior[i])
	{
	  tree->rates->t_prior_min[i] = 1.E+2 * tree->rates->nd_t[tree->n_root->num];
	  tree->rates->t_prior_max[i] = 0.0;
	}
    }
}

/*********************************************************/

void MCMC_Randomize_Node_Times_Bottom_Up(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  phydbl t_inf,t_sup;
  phydbl u;

  if(d->tax) 
    {
      if(tree->rates->t_has_prior[d->num])
	{
	  t_inf = tree->rates->t_prior_min[d->num];
	  t_sup = tree->rates->t_prior_max[d->num];
	  
	  u = Uni();
	  u *= (t_sup - t_inf);
	  u += t_inf;
	  
	  tree->rates->nd_t[d->num] = u;
	}
      return;
    }
  else 
    {
      t_node *v1, *v2; /* the two sons of d */

      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      MCMC_Randomize_Node_Times_Bottom_Up(d,d->v[i],tree);	      
	    }
	}
      
      v1 = v2 = NULL;
      For(i,3) if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	{
	  if(!v1) v1 = d->v[i]; 
	  else    v2 = d->v[i];
	}
      
      if(tree->rates->t_has_prior[d->num])
	{
	  t_sup = MIN(tree->rates->t_prior_max[d->num],
		      MIN(tree->rates->t_prior_max[v1->num],tree->rates->t_prior_max[v2->num]));

	  tree->rates->t_prior_max[d->num] = t_sup;

	  if(tree->rates->t_prior_max[d->num] < tree->rates->t_prior_min[d->num])
	    {
	      PhyML_Printf("\n. prior_min=%f prior_max=%f",tree->rates->t_prior_min[d->num],tree->rates->t_prior_max[d->num]);
	      PhyML_Printf("\n. Inconsistency in the prior settings detected at t_node %d",d->num);
	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	      Warn_And_Exit("\n");
	    }
	}
      else
	{
	  tree->rates->t_prior_max[d->num] = 
	    MIN(tree->rates->t_prior_max[v1->num],
		tree->rates->t_prior_max[v2->num]);
	}
    }
}

/*********************************************************/

void MCMC_Randomize_Node_Times_Top_Down(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;      
      t_node *v1, *v2; /* the two sons of d */
      phydbl u;
      phydbl t_inf, t_sup;

      
      v1 = v2 = NULL;
      For(i,3) if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	{
	  if(!v1) v1 = d->v[i]; 
	  else    v2 = d->v[i];
	}

      t_inf = MAX(tree->rates->nd_t[a->num],tree->rates->t_prior_min[d->num]);
      t_sup = tree->rates->t_prior_max[d->num];
	  	  
      if(t_inf > t_sup)
	{
	  PhyML_Printf("\n. t_inf=%f t_sup=%f.",t_inf,t_sup);
	  PhyML_Printf("\n. prior_min=%f prior_max=%f.",tree->rates->t_prior_min[d->num],tree->rates->t_prior_max[d->num]);
	  PhyML_Printf("\n. Inconsistency in the prior settings detected at t_node %d",d->num);
	  PhyML_Printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
	  Warn_And_Exit("\n");
	}
	  
      u = Uni();
      u *= (t_sup - t_inf);
      u += t_inf;
      
      tree->rates->nd_t[d->num] = u;
      
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      MCMC_Randomize_Node_Times_Top_Down(d,d->v[i],tree);
	    }
	}
    }
}

/*********************************************************/

void MCMC_Mixing_Step(t_tree *tree)
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
  multiplier = 1.0*EXP(u-0.5);

  For(i,2*tree->n_otu-1) tree->rates->nd_t[i] *= multiplier;
  tree->rates->clock_r /= multiplier;
  
  RATES_Update_Cur_Bl(tree);
  
  hr = POW(multiplier,-(tree->n_otu-1));
  
  new_lnL_rate = RATES_Lk_Rates(tree);
  new_lnL_time = RATES_Yule(tree);
  /* new_lnL_data = Lk(tree); */
  
  ratio = (new_lnL_rate + new_lnL_time) - (cur_lnL_rate + cur_lnL_time) + LOG(hr);
  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  u = Uni();
  
  if(u > alpha) /* Reject */
    {
      RATES_Reset_Times(tree);
      tree->rates->clock_r *= multiplier;

/*       RATES_Lk_Rates(tree); */
      tree->rates->c_lnL = cur_lnL_rate;
      
      if(FABS(cur_lnL_rate - tree->rates->c_lnL) > 1.E-0)
	{
	  PhyML_Printf("\n. lexp=%f alpha=%f ; cur_lnL_rates = %f vs %f",
		 tree->rates->lexp,
		 tree->rates->alpha,
		 cur_lnL_rate,tree->rates->c_lnL);

	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      if(FABS(cur_lnL_rate - tree->rates->c_lnL) > 1.E-3)
	{
	  PhyML_Printf("\n. WARNING: numerical precision issue detected (diff=%G run=%d). Reseting the likelihood.\n",cur_lnL_rate - tree->rates->c_lnL,tree->mcmc->run);
	  RATES_Lk_Rates(tree);
	  Lk(tree);
	}
    }
  else
    {
/*       PhyML_Printf("\n. Accept global times"); */
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
