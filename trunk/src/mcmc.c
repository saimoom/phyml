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
  int move,n_moves;
  phydbl u;
  int i;

  if(tree->mcmc->randomize)
    {
      MCMC_Randomize_Nu(tree);
      MCMC_Randomize_Node_Times(tree);
      MCMC_Randomize_Rates(tree);
      MCMC_Randomize_Clock_Rate(tree); /* Clock Rate must be the last parameter to be randomized */
    }

  MCMC_Print_Param(tree->mcmc,tree);
  
  RATES_Update_Cur_Bl(tree);
  
  n_moves = 5;
  do
    {
      if(tree->mcmc->run > tree->mcmc->n_tot_run / 3) MCMC_Adjust_Tuning_Parameter(tree->mcmc);

      tree->rates->c_lnL = RATES_Lk_Rates(tree);
      tree->c_lnL        = Lk(tree);
      
      u = Uni();

      For(move,tree->mcmc->n_moves) if(u < tree->mcmc->move_weight[move]) break;

      switch(move)
	{
	case 0:
	  {
	    MCMC_Clock_Rate(tree);
	    break;
	  }
	case 1:
	  {
	    /* MCMC_Tree_Height(tree); */
	    MCMC_Swing(tree);
	    break;
	  }
	case 2:
	  {
	    MCMC_SubTree_Height(tree);
	    break;
	  }
	case 3:
	  {
	    MCMC_Nu(tree);
	    break;
	  }
	case 4:
	  {
	    MCMC_Times_Local(tree);
	    break;
	  }
	case 5:
	  {
	    MCMC_Rates_Local(tree);
	    break;
	  }
	default:
	  {
	    PhyML_Printf("\n. u=%f move=%d",u,move);
	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Exit("\n");
	  }
	}
      RATES_Check_Node_Times(tree);
    }
  while(tree->mcmc->run < tree->mcmc->n_tot_run);
}

/*********************************************************/

void MCMC_Nu(t_tree *tree)
{
  phydbl cur_nu,new_nu,cur_lnL_rate,new_lnL_rate;
  phydbl u,alpha,ratio;
  phydbl min_nu,max_nu;
  phydbl K,mult;

  if(
     (tree->rates->model == COMPOUND_COR)   || 
     (tree->rates->model == COMPOUND_NOCOR) ||
     (tree->rates->model == GAMMA) ||
     (tree->rates->model == EXPONENTIAL)
     ) return;

  cur_nu        = -1.0;
  new_nu        = -1.0;
  ratio         = -1.0;

  K            = tree->mcmc->h_nu;
  cur_lnL_rate = tree->rates->c_lnL;
  cur_nu       = tree->rates->nu;
  
  min_nu = tree->rates->min_nu;
  max_nu = tree->rates->max_nu;

  u = Uni();
  /* mult = u*(K-1./K)+1./K; */
  mult = EXP(K*(u-0.5));
  new_nu = cur_nu * mult;

  if(new_nu > max_nu || new_nu < min_nu) 
    {
      /* printf("\n. Extreme value of nu: %G [%G;%G] K=%f\n",new_nu,min_nu,max_nu,K); */
      tree->mcmc->run_nu++;
      return;
    }

  tree->rates->nu = new_nu;

  new_lnL_rate = RATES_Lk_Rates(tree);

  ratio = LOG(mult);
  if(tree->mcmc->use_data == YES) ratio += (new_lnL_rate - cur_lnL_rate);


  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  phydbl eps = 0.01;
  u = Uni();
  if(u > alpha) /* Reject */
    {
      tree->rates->nu    = cur_nu;
      tree->rates->c_lnL = cur_lnL_rate;
    }
  else
    {
      tree->mcmc->acc_nu++;
    }

  tree->mcmc->run++;
  tree->mcmc->run_nu++;
  MCMC_Print_Param(tree->mcmc,tree);

  if(!(tree->mcmc->run%tree->mcmc->norm_freq))
    {
      RATES_Normalise_Rates(tree);
      RATES_Lk_Rates(tree);
      Lk(tree);
    }
}

/*********************************************************/

void MCMC_Clock_Rate(t_tree *tree)
{
  phydbl cur_cr,new_cr,cur_lnL,new_lnL;
  phydbl u,alpha,ratio;
  phydbl K,mult;

  
  cur_lnL       = UNLIKELY;
  new_lnL       = UNLIKELY;
  cur_cr        = -1.0;
  new_cr        = -1.0;
  ratio         = -1.0;

  K            = tree->mcmc->h_clock;
  cur_lnL      = tree->c_lnL;
  cur_cr       = tree->rates->clock_r;  

  u = Uni();
  /* mult = u*(K-1./K)+1./K; */
  mult = EXP(K*(u-0.5));

  new_cr = cur_cr * mult;

  if(new_cr > tree->rates->max_clock || new_cr < tree->rates->min_clock) 
    {
      /* printf("\n. Extreme value of clock_r: %G [%G;%G]\n",new_cr,tree->rates->min_clock,tree->rates->max_clock); */
      tree->mcmc->run_clock++;
      return;
    }

  tree->rates->clock_r = new_cr;
  RATES_Update_Cur_Bl(tree);
  new_lnL = Lk(tree);

  ratio = LOG(mult);
  if(tree->mcmc->use_data == YES) ratio += (new_lnL - cur_lnL);

  ratio = EXP(ratio);
  
  alpha = MIN(1.,ratio);
  u = Uni();
  if(u > alpha) /* Reject */
    {
      tree->rates->clock_r = cur_cr;
      RATES_Update_Cur_Bl(tree);
      tree->c_lnL = cur_lnL;
    }
  else
    {
      tree->mcmc->acc_clock++;
    }

  tree->mcmc->run++;
  tree->mcmc->run_clock++;
  MCMC_Print_Param(tree->mcmc,tree);
  
  if(!(tree->mcmc->run%tree->mcmc->norm_freq))
    {
      RATES_Normalise_Rates(tree);
      RATES_Lk_Rates(tree);
      Lk(tree);
    }
}

/*********************************************************/

void MCMC_Times_Local(t_tree *tree)
{
  int node_num;
  int i;
  
  For(i,2*tree->n_otu-1)
    {
      node_num = Rand_Int(tree->n_otu,2*tree->n_otu-2);
      
      if(tree->noeud[node_num] != tree->n_root)
	MCMC_Times_Pre(tree->noeud[node_num]->anc,tree->noeud[node_num],tree);

      /* !!!!!!! */
      break;

    }

  RATES_Lk_Rates(tree); /* Rates likelihood needs to be updated here */
}

/*********************************************************/

void MCMC_Rates_Local(t_tree *tree)
{
  int node_num;
  int i;

  RATES_Lk_Rates(tree);

  For(i,2*tree->n_otu-2)
    {
      node_num = Rand_Int(0,2*tree->n_otu-3);
      MCMC_One_Rate(tree->noeud[node_num]->anc,tree->noeud[node_num],tree);
      
      /* !!!!!!! */
      break;
    }
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
  /* t_node *v2,*v3; */
  /* phydbl T0,T1,T2,T3; */
  /* phydbl U0,U1,U2,U3; */
  /* phydbl V1,V2,V3; */
  phydbl tmp;
  phydbl K,mult;

  if(tree->rates->ml_l[d->num] < BL_MIN*1.1) 
    {
      tree->rates->nd_r[d->num] = tree->rates->min_rate;
      return;
    }


  nu           = tree->rates->nu;
  cur_mu       = tree->rates->nd_r[d->num];
  cur_lnL_data = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL;
  r_min        = tree->rates->min_rate;
  r_max        = tree->rates->max_rate;
  K            = tree->mcmc->h_rates;
  tmp = cur_lnL_rate;
  
  u = Uni();
  /* mult = u*(K-1./K)+1./K; */
  mult = EXP(K*(u-0.5));
  new_mu = cur_mu * mult;


  if(new_mu < r_min || new_mu > r_max) 
    {
      tree->mcmc->run_rates++;
      return;
    }

  tree->rates->nd_r[d->num] = new_mu;  
  RATES_Update_Cur_Bl(tree);

  new_lnL_data = Lk(tree);
  new_lnL_rate = RATES_Lk_Rates(tree);

  ratio = LOG(mult);
  ratio += (new_lnL_rate - cur_lnL_rate);
  if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);

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


  tree->mcmc->run++;
  tree->mcmc->run_rates++;
  MCMC_Print_Param(tree->mcmc,tree);

  if(!(tree->mcmc->run%tree->mcmc->norm_freq))
    {
      RATES_Normalise_Rates(tree);
      RATES_Lk_Rates(tree);
      Lk(tree);
    }
}

/*********************************************************/

void MCMC_Times_Pre(t_node *a, t_node *d, t_tree *tree)
{
  phydbl u;
  phydbl t_min,t_max;
  phydbl t1_cur, t1_new;
  phydbl cur_lnL_data, new_lnL_data;
  phydbl ratio,alpha;
  int    i;
  phydbl t0,t2,t3;
  t_node *v2,*v3;
  phydbl K,mult;

  if(d->tax) return; /* Won't change time at tip */

  if(FABS(tree->rates->t_prior_min[d->num] - tree->rates->t_prior_max[d->num]) < 1.E-10) return;



  K             = tree->mcmc->h_times;
  cur_lnL_data  = tree->c_lnL;
  t1_cur        = tree->rates->nd_t[d->num];
    
  v2 = v3 = NULL;

  if(a)
    {
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  {
	    if(!v2) { v2 = d->v[i]; }
	    else    { v3 = d->v[i]; }
	  }
    }
  else
    {
      v2 = tree->n_root->v[0];
      v3 = tree->n_root->v[1];
    }
  
  t0 = (a)?(tree->rates->nd_t[a->num]):(tree->rates->t_prior_min[d->num]);
  t2 = tree->rates->nd_t[v2->num];
  t3 = tree->rates->nd_t[v3->num];

  t_min = MAX(t0,tree->rates->t_prior_min[d->num]);
  t_max = MIN(MIN(t2,t3),tree->rates->t_prior_max[d->num]);

  t_min += tree->rates->min_dt;
  t_max -= tree->rates->min_dt;

  if(t_min > t_max) 
    {
      PhyML_Printf("\n. t_min = %f t_max = %f",t_min,t_max);
      PhyML_Printf("\n. prior_min = %f prior_max = %f",tree->rates->t_prior_min[d->num],tree->rates->t_prior_max[d->num]);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      /* Exit("\n"); */
    }

  u = Uni();
  /* mult = u*(K-1./K)+1./K; */
  mult = EXP(K*(u-0.5));
  t1_new = t1_cur * mult;

  if(t1_new < t_min || t1_new > t_max) 
    {
      tree->mcmc->run_times++;
      return;
    }

  tree->rates->nd_t[d->num] = t1_new;
  RATES_Update_Cur_Bl(tree);
  new_lnL_data = Lk(tree);

  ratio = LOG(mult);
  if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  if(u > alpha) /* Reject */
    {
      tree->rates->nd_t[d->num] = t1_cur;
      RATES_Update_Cur_Bl(tree);
      tree->c_lnL = cur_lnL_data;
    }
  else
    {
      tree->mcmc->acc_times++;
    }    

  if(t1_new < t0)
    {
      t1_new = t0+1.E-4;
      PhyML_Printf("\n");
      PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
      PhyML_Printf("\n. t0 = %f t1_new = %f",t0,t1_new);
      PhyML_Printf("\n. t_min=%f t_max=%f",t_min,t_max);
      PhyML_Printf("\n. (t1-t0)=%f (t2-t1)=%f",t1_cur-t0,t2-t1_cur);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
/*       Exit("\n"); */
    }
  if(t1_new > MIN(t2,t3))
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
      PhyML_Printf("\n. t0 = %f t1_new = %f t1 = %f t2 = %f t3 = %f MIN(t2,t3)=%f",t0,t1_new,t1_cur,t2,t3,MIN(t2,t3));
      PhyML_Printf("\n. t_min=%f t_max=%f",t_min,t_max);
      PhyML_Printf("\n. (t1-t0)=%f (t2-t1)=%f",t1_cur-t0,t2-t1_cur);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
/*       Exit("\n"); */
    }

  if(isnan(t1_new))
    {
      PhyML_Printf("\n. run=%d",tree->mcmc->run);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
/*       Exit("\n"); */
    }
  
  tree->mcmc->run++;
  tree->mcmc->run_times++;
  MCMC_Print_Param(tree->mcmc,tree);

  if(!(tree->mcmc->run%tree->mcmc->norm_freq))
    {
      RATES_Normalise_Rates(tree);
      RATES_Lk_Rates(tree);
      Lk(tree);
    }
}

/*********************************************************/

void MCMC_Tree_Height(t_tree *tree)
{
  int i;
  phydbl K,mult,u,alpha,ratio,beta;
  phydbl cur_lnL_data,new_lnL_data;
  phydbl add;
  phydbl floor;

  RATES_Record_Times(tree);

  K = tree->mcmc->h_tree_height;
  cur_lnL_data  = tree->c_lnL;

  u = Uni();
  mult = EXP(K*(u-0.5));
  
  floor = 0.0;
  For(i,tree->n_otu) 
    if(tree->rates->t_prior_max[i] < floor) 
      floor = tree->rates->t_prior_max[i];

  Scale_Subtree_Height(tree->n_root,mult,floor,tree);

  For(i,2*tree->n_otu-1)
    {
      if(tree->rates->nd_t[i] > tree->rates->t_prior_max[i] ||
  	 tree->rates->nd_t[i] < tree->rates->t_prior_min[i])
  	{
  	  RATES_Reset_Times(tree);
	  tree->mcmc->run_tree_height++;
  	  return;
  	}
    }
  
  RATES_Update_Cur_Bl(tree);
  new_lnL_data = Lk(tree);

  /* Yes, the Hastings ratio is actually mult */
  ratio = LOG(mult);
  if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  if(u > alpha)
    {
      RATES_Reset_Times(tree);
      RATES_Update_Cur_Bl(tree);
      tree->c_lnL = cur_lnL_data;
    }
  else
    {
      tree->mcmc->acc_tree_height++;
    }

  tree->mcmc->run++;
  tree->mcmc->run_tree_height++;
  MCMC_Print_Param(tree->mcmc,tree);
  
  if(!(tree->mcmc->run%tree->mcmc->norm_freq))
    {
      RATES_Normalise_Rates(tree);
      RATES_Lk_Rates(tree);
      Lk(tree);
    }
}

/*********************************************************/

void MCMC_SubTree_Height(t_tree *tree)
{
  int i;
  phydbl K,mult,u,alpha,ratio,beta;
  phydbl cur_lnL_data,new_lnL_data;
  phydbl add;
  phydbl floor;
  int target;

  RATES_Record_Times(tree);

  K = tree->mcmc->h_subtree_height;
  cur_lnL_data  = tree->c_lnL;

  u = Uni();
  /* mult = u*(K-1./K)+1./K; */
  mult = EXP(K*(u-0.5));
  
  target = Rand_Int(tree->n_otu,2*tree->n_otu-2);

  floor = 0.0;
  For(i,tree->n_otu) 
    if((tree->rates->nd_t[i] > tree->rates->nd_t[target]) &&
       (tree->rates->t_prior_max[i] < floor)) 
      floor = tree->rates->t_prior_max[i];

  Scale_Subtree_Height(tree->noeud[target],mult,floor,tree);

  For(i,2*tree->n_otu-1)
    {
      if(tree->rates->nd_t[i] > tree->rates->t_prior_max[i] ||
  	 tree->rates->nd_t[i] < tree->rates->t_prior_min[i])
  	{
  	  RATES_Reset_Times(tree);
	  tree->mcmc->run_tree_height++;
  	  return;
  	}
    }
  
  /* tree->rates->clock_r /= mult; */

  RATES_Update_Cur_Bl(tree);
  new_lnL_data = Lk(tree);

  /* Yes, the Hastings ratio is actually mult */
  ratio = LOG(mult);
  if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  if(u > alpha)
    {
      /* tree->rates->clock_r *= mult; */
      RATES_Reset_Times(tree);
      RATES_Update_Cur_Bl(tree);
      tree->c_lnL = cur_lnL_data;
    }
  else
    {
      tree->mcmc->acc_subtree_height++;
    }

  tree->mcmc->run++;
  tree->mcmc->run_subtree_height++;
  MCMC_Print_Param(tree->mcmc,tree);
  
  if(!(tree->mcmc->run%tree->mcmc->norm_freq))
    {
      RATES_Normalise_Rates(tree);
      RATES_Lk_Rates(tree);
      Lk(tree);
    }
}

/*********************************************************/

void MCMC_Swing(t_tree *tree)
{
  int i;
  phydbl K,mult,u,alpha,ratio,beta;
  phydbl cur_lnL_data,new_lnL_data;
  phydbl add;
  phydbl floor;

  RATES_Record_Times(tree);

  K = tree->mcmc->h_tree_height;
  cur_lnL_data  = tree->c_lnL;

  u = Uni();
  mult = EXP(K*(u-0.5));
  
  floor = 0.0;
  For(i,tree->n_otu) 
    if(tree->rates->t_prior_max[i] < floor) 
      floor = tree->rates->t_prior_max[i];

  Scale_Subtree_Height(tree->n_root,mult,floor,tree);

  For(i,2*tree->n_otu-1)
    {
      if(tree->rates->nd_t[i] > tree->rates->t_prior_max[i] ||
  	 tree->rates->nd_t[i] < tree->rates->t_prior_min[i])
  	{
  	  RATES_Reset_Times(tree);
	  tree->mcmc->run_tree_height++;
  	  return;
  	}
    }
  
  tree->rates->clock_r /= mult;

  RATES_Update_Cur_Bl(tree);
  new_lnL_data = Lk(tree);

  /* Yes, the Hastings ratio is actually 0.0 */
  ratio = 0.0;
  if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  if(u > alpha)
    {
      tree->rates->clock_r *= mult;
      RATES_Reset_Times(tree);
      RATES_Update_Cur_Bl(tree);
      tree->c_lnL = cur_lnL_data;
    }
  else
    {
      tree->mcmc->acc_tree_height++;
    }

  tree->mcmc->run++;
  tree->mcmc->run_tree_height++;
  MCMC_Print_Param(tree->mcmc,tree);
  
  if(!(tree->mcmc->run%tree->mcmc->norm_freq))
    {
      RATES_Normalise_Rates(tree);
      RATES_Lk_Rates(tree);
      Lk(tree);
    }
}

/*********************************************************/

void MCMC_Print_Param(tmcmc *mcmc, t_tree *tree)
{
  int i;
  FILE *fp;
  phydbl eps = 0.01;

  if(tree->mcmc->run > mcmc->n_tot_run) return;

  fp = mcmc->out_fp_stats;
  
  if(tree->mcmc->run == 0)
    {
      PhyML_Fprintf(stdout," ["); 
      fflush(NULL); 
    }

  if(!(mcmc->run%(mcmc->n_tot_run/10))) 
    { 
      PhyML_Fprintf(stdout,"."); 
      fflush(NULL); 
    }

  if(tree->mcmc->run == mcmc->n_tot_run)
    {
      PhyML_Fprintf(stdout,"]"); 
      fflush(NULL); 
    }


/*   MCMC_Print_Means(mcmc,tree); */
/*   MCMC_Print_Last(mcmc,tree); */

  if(!(mcmc->run%mcmc->sample_interval)) 
    {      
      if(tree->mcmc->run == 0)
	{
	  time(&(mcmc->t_beg));

	  PhyML_Fprintf(fp,"\n");
	  PhyML_Fprintf(fp,"Run\t");
	  PhyML_Fprintf(fp,"Time\t");
	  PhyML_Fprintf(fp,"MeanRate\t");

	  PhyML_Fprintf(fp,"AccClock\t");
	  PhyML_Fprintf(fp,"AccNu\t");
	  PhyML_Fprintf(fp,"AccTreeHeight\t");
	  PhyML_Fprintf(fp,"AccSubTreeHeight\t");
	  PhyML_Fprintf(fp,"AccTimes\t");
	  PhyML_Fprintf(fp,"AccRates\t");

	  PhyML_Fprintf(fp,"HClock\t");
	  PhyML_Fprintf(fp,"HNu\t");
	  PhyML_Fprintf(fp,"HTreeHeight\t");
	  PhyML_Fprintf(fp,"HSubTreeHeight\t");
	  PhyML_Fprintf(fp,"HTimes\t");
	  PhyML_Fprintf(fp,"HRates\t");

	  PhyML_Fprintf(fp,"LnLSeq\t");
	  PhyML_Fprintf(fp,"LnLSeq.Approx\t");
	  PhyML_Fprintf(fp,"LnLRate\t");
	  PhyML_Fprintf(fp,"Nu\t");
	  PhyML_Fprintf(fp,"Clock\t");

	  if(fp != stdout)
	    {
	      for(i=tree->n_otu;i<2*tree->n_otu-1;i++)
		{
		  PhyML_Fprintf(fp,"T%d%s\t",i,tree->noeud[i] == tree->n_root?"*":"");
		}
	    }


	  if(fp != stdout)
	    {
	      For(i,2*tree->n_otu-3)
		{
		  if(tree->t_edges[i] == tree->e_root)
		    PhyML_Fprintf(fp,"*L%d[%f]\t",i,tree->rates->u_ml_l[i]);
		  else
		    PhyML_Fprintf(fp," L%d[%f]\t",i,tree->rates->u_ml_l[i]);		    
		}
	    }


	  PhyML_Fprintf(mcmc->out_fp_trees,"#NEXUS\n");
	  PhyML_Fprintf(mcmc->out_fp_trees,"BEGIN TREES;\n");
	  PhyML_Fprintf(mcmc->out_fp_trees,"\tTRANSLATE\n");
	  For(i,tree->n_otu-1) PhyML_Fprintf(mcmc->out_fp_trees,"\t%3d\t%s,\n",tree->noeud[i]->num,tree->noeud[i]->name);
	  PhyML_Fprintf(mcmc->out_fp_trees,"\t%3d\t%s;\n",tree->noeud[i]->num,tree->noeud[i]->name);
	  tree->write_tax_names = NO;
	}

      PhyML_Fprintf(fp,"\n");
      PhyML_Fprintf(fp,"%6d\t",tree->mcmc->run);
      time(&mcmc->t_cur);
      PhyML_Fprintf(fp,"%6d\t",(int)(mcmc->t_cur-mcmc->t_beg));
      
      RATES_Update_Cur_Bl(tree);
      PhyML_Fprintf(fp,"%f\t",RATES_Check_Mean_Rates(tree));
            
      PhyML_Fprintf(fp,"%f\t",(phydbl)(tree->mcmc->acc_clock+eps)/(tree->mcmc->run_clock+eps));
      PhyML_Fprintf(fp,"%f\t",(phydbl)(tree->mcmc->acc_nu+eps)/(tree->mcmc->run_nu+eps));
      PhyML_Fprintf(fp,"%f\t",(phydbl)(tree->mcmc->acc_tree_height+eps)/(tree->mcmc->run_tree_height+eps));
      PhyML_Fprintf(fp,"%f\t",(phydbl)(tree->mcmc->acc_subtree_height+eps)/(tree->mcmc->run_subtree_height+eps));
      PhyML_Fprintf(fp,"%f\t",(phydbl)(tree->mcmc->acc_times+eps)/(tree->mcmc->run_times+eps));
      PhyML_Fprintf(fp,"%f\t",(phydbl)(tree->mcmc->acc_rates+eps)/(tree->mcmc->run_rates+eps));


      PhyML_Fprintf(fp,"%f\t",(phydbl)(tree->mcmc->h_clock));
      PhyML_Fprintf(fp,"%f\t",(phydbl)(tree->mcmc->h_nu));
      PhyML_Fprintf(fp,"%f\t",(phydbl)(tree->mcmc->h_tree_height));
      PhyML_Fprintf(fp,"%f\t",(phydbl)(tree->mcmc->h_subtree_height));
      PhyML_Fprintf(fp,"%f\t",(phydbl)(tree->mcmc->h_times));
      PhyML_Fprintf(fp,"%f\t",(phydbl)(tree->mcmc->h_rates));



      PhyML_Fprintf(fp,"%.1f\t",tree->c_lnL);

      /* !!!!!!!!!!!!!!!!!!!!!!! */
      phydbl tmp=tree->c_lnL;
      tree->rates->lk_approx = NORMAL; Lk(tree);
      PhyML_Fprintf(fp,"%.1f\t",tree->c_lnL);
      /* tree->rates->lk_approx = EXACT; tree->c_lnL = tmp; */

      PhyML_Fprintf(fp,"%.1f\t",tree->rates->c_lnL);
      PhyML_Fprintf(fp,"%G\t",tree->rates->nu);
      PhyML_Fprintf(fp,"%G\t",tree->rates->clock_r);
      for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(fp,"%.1f\t",tree->rates->nd_t[i] - tree->rates->true_t[i]);
      /* for(i=0;i<2*tree->n_otu-2;i++) PhyML_Fprintf(fp,"%f\t",tree->rates->nd_r[i]); */
      /* if(fp != stdout) for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(fp,"%G\t",tree->rates->t_prior[i]); */
      For(i,2*tree->n_otu-3) PhyML_Fprintf(fp,"%f\t",tree->t_edges[i]->l);
      fflush(NULL);


      // TREES
      char *s;
      Branch_Lengths_To_Time_Lengths(tree);
      /* tree->bl_ndigits = 7; */
      tree->bl_ndigits = 0;
      s = Write_Tree(tree);
      tree->bl_ndigits = 7;
      PhyML_Fprintf(mcmc->out_fp_trees,"TREE %8d [%f] = [&R] %s\n",mcmc->run,tree->c_lnL,s);
      Free(s);
      RATES_Update_Cur_Bl(tree);
    }

  if(tree->mcmc->run == mcmc->n_tot_run)
    {
      PhyML_Fprintf(mcmc->out_fp_trees,"END;\n");
      fflush(NULL); 
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
      
      PhyML_Fprintf(mcmc->out_fp_means,"#");
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

      PhyML_Fprintf(mcmc->out_fp_last,"#");
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
  mcmc->move_weight  = (phydbl *)mCalloc(N_MAX_MOVES,sizeof(phydbl));

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
  Free(mcmc->move_weight);
  Free(mcmc);
}

/*********************************************************/

void MCMC_Init_MCMC_Struct(char *filename, tmcmc *mcmc, t_tree *tree)
{
  int pid;
  int i;
  phydbl sum;

  mcmc->use_data           = tree->io->use_data;
  mcmc->adjust_tuning      = YES;
  mcmc->acc_tree_height    = 0;
  mcmc->acc_subtree_height = 0;
  mcmc->acc_clock          = 0;
  mcmc->acc_rates          = 0;
  mcmc->acc_times          = 0;
  mcmc->run_nu             = 0;
  mcmc->run_tree_height    = 0;
  mcmc->run_subtree_height = 0;
  mcmc->run_clock          = 0;
  mcmc->run_rates          = 0;
  mcmc->run_times          = 0;
  mcmc->run                = 0;
  mcmc->sample_interval    = 1E+3;
  mcmc->n_rate_jumps       = 0;
  mcmc->n_tot_run          = 1E+6;
  mcmc->randomize          = 1;
  mcmc->norm_freq          = 20*tree->n_otu;
/*   mcmc->norm_freq       = 1E+9; */

  /* mcmc->h_times         = 1.2; */
  /* mcmc->h_rates         = 1.2; */
  /* mcmc->h_nu            = 1.1; */
  /* mcmc->h_clock         = 1.1; */

  mcmc->h_times          = 2.00;
  mcmc->h_rates          = 4.00;
  mcmc->h_nu             = 1.00;
  mcmc->h_clock          = 1.00;
  mcmc->h_tree_height    = 2.00;
  mcmc->h_subtree_height = 1.00;

  mcmc->n_moves          = 6;
  
  mcmc->move_weight[0] = (phydbl)tree->n_otu; /* Clock */
  mcmc->move_weight[1] = (phydbl)tree->n_otu; /* Tree height */
  mcmc->move_weight[2] = (phydbl)tree->n_otu; /* SubTree height */
  mcmc->move_weight[3] = (phydbl)tree->n_otu; /* Nu */
  mcmc->move_weight[4] = 2.*(phydbl)tree->n_otu; /* Node heights */
  mcmc->move_weight[5] = 2.*(phydbl)tree->n_otu; /* Edge rates */

  sum = 0.0;
  For(i,N_MAX_MOVES) sum += mcmc->move_weight[i];
  For(i,N_MAX_MOVES) mcmc->move_weight[i] /= sum;
  for(i=1;i<N_MAX_MOVES;i++) mcmc->move_weight[i] += mcmc->move_weight[i-1];

  For(i,mcmc->n_moves)
    {
      printf("\n. ## %f",mcmc->move_weight[i]);
    }

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

/*       strcpy(s,tree->mcmc->out_filename); */
/*       strcat(s,".means"); */
/*       tree->mcmc->out_fp_means = fopen(s,"w"); */

/*       strcpy(s,tree->mcmc->out_filename); */
/*       strcat(s,".lasts"); */
/*       tree->mcmc->out_fp_last  = fopen(s,"w"); */

      Free(s);
    }
  else 
    {
      tree->mcmc->out_fp_stats = stderr;
      tree->mcmc->out_fp_trees = stderr;
      /* tree->mcmc->out_fp_means = stderr; */
      /* tree->mcmc->out_fp_last  = stderr; */
    }
}

/*********************************************************/

void MCMC_Copy_MCMC_Struct(tmcmc *ori, tmcmc *cpy, char *filename, t_tree *tree)
{
  int pid;
  int i;
  
  cpy->use_data           = ori->use_data        ;
  cpy->adjust_tuning      = ori->adjust_tuning   ;
  cpy->acc_tree_height    = ori->acc_tree_height ;
  cpy->acc_subtree_height = ori->acc_subtree_height ;
  cpy->acc_clock          = ori->acc_clock       ;
  cpy->acc_rates          = ori->acc_rates       ;
  cpy->acc_times          = ori->acc_times       ;
  cpy->run_nu             = ori->run_nu          ;
  cpy->run_tree_height    = ori->run_tree_height ;
  cpy->run_subtree_height = ori->run_subtree_height ;
  cpy->run_clock          = ori->run_clock       ;
  cpy->run_rates          = ori->run_rates       ;
  cpy->run_times          = ori->run_times       ;
  cpy->run                = ori->run             ;
  cpy->sample_interval    = ori->sample_interval ;
  cpy->n_rate_jumps       = ori->n_rate_jumps    ;
  cpy->n_tot_run          = ori->n_tot_run       ;
  cpy->randomize          = ori->randomize       ;
  cpy->norm_freq          = ori->norm_freq       ;

  cpy->h_times            = ori->h_times         ;
  cpy->h_rates            = ori->h_rates         ;
  cpy->h_nu               = ori->h_nu            ;
  cpy->h_clock            = ori->h_clock         ;
  cpy->h_tree_height      = ori->h_tree_height   ;
  cpy->h_subtree_height   = ori->h_subtree_height   ;

  cpy->n_moves            = ori->n_moves         ;

  For(i,cpy->n_moves) cpy->move_weight[i] = ori->move_weight[i];

  if(filename) 
    {
      char *s;

      s = (char *)mCalloc(T_MAX_NAME,sizeof(char));

      strcpy(cpy->out_filename,filename);
      pid = getpid();
      sprintf(cpy->out_filename+strlen(cpy->out_filename),".%d",pid);

      cpy->out_fp_stats = fopen(cpy->out_filename,"w");

      strcpy(s,cpy->out_filename);
      strcat(s,".trees");
      cpy->out_fp_trees = fopen(s,"w");

      Free(s);
    }
  else 
    {
      cpy->out_fp_stats = stderr;
      cpy->out_fp_trees = stderr;
    }
}


/*********************************************************/

void MCMC_Close_MCMC(tmcmc *mcmc)
{
  fclose(mcmc->out_fp_trees);
  fclose(mcmc->out_fp_stats);
  /* fclose(mcmc->out_fp_means); */
  /* fclose(mcmc->out_fp_last); */
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
  int iter;

  t_inf = tree->rates->t_prior_min[tree->n_root->num];
  /* t_sup = MIN(tree->rates->t_prior_max[tree->n_root->num], */
  /* 	      MIN(tree->rates->nd_t[tree->n_root->v[0]->num], */
  /* 		  tree->rates->nd_t[tree->n_root->v[1]->num])); */
  t_sup = tree->rates->t_prior_max[tree->n_root->num];

  u = Uni();
  u *= (t_sup - t_inf);
  u += t_inf;
  
  tree->rates->nd_t[tree->n_root->num] = u;


  MCMC_Randomize_Node_Times_Top_Down(tree->n_root,tree->n_root->v[0],tree);
  MCMC_Randomize_Node_Times_Top_Down(tree->n_root,tree->n_root->v[1],tree);
  
  For(iter,10)
    {
      MCMC_Randomize_Node_Times_Bottom_Up(tree->n_root,tree->n_root->v[0],tree);
      MCMC_Randomize_Node_Times_Bottom_Up(tree->n_root,tree->n_root->v[1],tree);
    }

  /* TIMES_Print_Node_Times(tree->n_root,tree->n_root->v[0],tree); */
  /* TIMES_Print_Node_Times(tree->n_root,tree->n_root->v[1],tree); */

  RATES_Check_Node_Times(tree);    
}

/*********************************************************/

void MCMC_Randomize_Node_Times_Bottom_Up(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;      
      phydbl u;
      phydbl t_inf, t_sup;
      t_node *v1, *v2;


      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      MCMC_Randomize_Node_Times_Bottom_Up(d,d->v[i],tree);
	    }
	}

      v1 = v2 = NULL;
      For(i,3)
      	{
      	  if(d->v[i] != a && d->b[i] != tree->e_root)
      	    {
      	      if(!v1) v1 = d->v[i];
      	      else    v2 = d->v[i];
      	    }
      	}

      t_sup = MIN(tree->rates->nd_t[v1->num],tree->rates->nd_t[v2->num]);
      t_inf = tree->rates->nd_t[a->num],tree->rates->nd_t[v2->num];

      u = Uni();
      u *= (t_sup - t_inf);
      u += t_inf;
      
      if(u > tree->rates->t_prior_min[d->num] && u < tree->rates->t_prior_max[d->num])
	tree->rates->nd_t[d->num] = u;
    }  
}

/*********************************************************/

void MCMC_Randomize_Node_Times_Top_Down(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;      
      phydbl u;
      phydbl t_inf, t_sup;
      /* t_node *v1, *v2; */

      t_inf = MAX(tree->rates->nd_t[a->num],tree->rates->t_prior_min[d->num]);
      t_sup = tree->rates->t_prior_max[d->num];

      /* t_inf = tree->rates->t_prior_min[d->num]/3.; */
      /* t_sup = MIN(tree->rates->nd_t[v1->num],tree->rates->nd_t[v2->num]); */

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
      
      /* v1 = v2 = NULL; */
      /* For(i,3) */
      /* 	{ */
      /* 	  if(d->v[i] != a && d->b[i] != tree->e_root) */
      /* 	    { */
      /* 	      if(!v1) v1 = d->v[i]; */
      /* 	      else    v2 = d->v[i]; */
      /* 	    } */
      /* 	} */

    }
}

/*********************************************************/


void MCMC_Adjust_Tuning_Parameter(tmcmc *mcmc)
{
  if(mcmc->adjust_tuning)
    {
      phydbl eps=0.01;
      phydbl scale=1.2;
      phydbl rate;
      phydbl rate_inf,rate_sup;

      /* printf("\n. run=%d sample_interval=%d",mcmc->run,mcmc->sample_interval); */
      /* printf("\n. %f %f %f %f %f", */
      /* 	     (phydbl)(mcmc->acc_nu+eps)/(mcmc->run_nu+eps), */
      /* 	     (phydbl)(mcmc->acc_clock+eps)/(mcmc->run_clock+eps), */
      /* 	     (phydbl)(mcmc->acc_tree_height+eps)/(mcmc->run_tree_height+eps), */
      /* 	     (phydbl)(mcmc->acc_rates+eps)/(mcmc->run_rates+eps), */
      /* 	     (phydbl)(mcmc->acc_times+eps)/(mcmc->run_times+eps)); */
      
      rate_inf = 0.1;
      rate_sup = 0.5;
      
      if(mcmc->run_nu > 100)
	{
	  /* printf("\n. nu %d %d [%f]",mcmc->acc_nu,mcmc->run_nu,mcmc->h_nu); */
	  rate = (phydbl)(mcmc->acc_nu+eps)/(mcmc->run_nu+eps); 
	  if(rate < rate_inf)
	    {
	      mcmc->h_nu /= scale;
	    }
	  else if(rate > rate_sup)
	    {
	      mcmc->h_nu *= scale;
	    }

	  mcmc->run_nu = 0;
	  mcmc->acc_nu = 0;
	}
      
      if(mcmc->run_clock > 100)
	{
	  /* printf("\n. clock %d %d [%f]",mcmc->acc_clock,mcmc->run_clock,mcmc->h_clock); */
	  rate = (phydbl)(mcmc->acc_clock+eps)/(mcmc->run_clock+eps); 
	  if(rate < rate_inf)
	    {
	      mcmc->h_clock /= scale;
	    }
	  else if(rate > rate_sup)
	    {
	      mcmc->h_clock *= scale;
	    }
	  mcmc->run_clock = 0;
	  mcmc->acc_clock = 0;
	}
      
      if(mcmc->run_tree_height > 100)
	{
	  /* printf("\n. height %d %d [%f]",mcmc->acc_tree_height,mcmc->run_tree_height,mcmc->h_tree_height); */
	  rate = (phydbl)(mcmc->acc_tree_height+eps)/(mcmc->run_tree_height+eps); 
	  if(rate < rate_inf)
	    {
	      mcmc->h_tree_height /= scale;
	    }
	  else if(rate > rate_sup)
	    {
	      mcmc->h_tree_height *= scale;
	    }
	  mcmc->run_tree_height = 0;
	  mcmc->acc_tree_height = 0;
	}

      if(mcmc->run_subtree_height > 100)
	{
	  /* printf("\n. subheight %d %d [%f]",mcmc->acc_subtree_height,mcmc->run_subtree_height,mcmc->h_subtree_height); */
	  rate = (phydbl)(mcmc->acc_subtree_height+eps)/(mcmc->run_subtree_height+eps); 
	  if(rate < rate_inf)
	    {
	      mcmc->h_subtree_height /= scale;
	    }
	  else if(rate > rate_sup)
	    {
	      mcmc->h_subtree_height *= scale;
	    }
	  mcmc->run_subtree_height = 0;
	  mcmc->acc_subtree_height = 0;
	}

      
      if(mcmc->run_rates > 100)
	{
	  /* printf("\n. rates %d %d [%f]",mcmc->acc_rates,mcmc->run_rates,mcmc->h_tree_height); */
	  rate = (phydbl)(mcmc->acc_rates+eps)/(mcmc->run_rates+eps); 
	  if(rate < rate_inf)
	    {
	      mcmc->h_rates /= scale;
	    }
	  else if(rate > rate_sup)
	    {
	      mcmc->h_rates *= scale;
	    }
	  mcmc->run_rates = 0;
	  mcmc->acc_rates = 0;
	}
      
      
      if(mcmc->run_times > 100)
	{
	  /* printf("\n. times %d %d [%f]",mcmc->acc_times,mcmc->run_times,mcmc->h_times); */
	  rate = (phydbl)(mcmc->acc_times+eps)/(mcmc->run_times+eps); 
	  if(rate < rate_inf)
	    {
	      mcmc->h_times /= scale;
	    }
	  else if(rate > rate_sup)
	    {
	      mcmc->h_times *= scale;
	    }
	  mcmc->run_times = 0;
	  mcmc->acc_times = 0;	  
	}
    }
}

/*********************************************************/

void MCMC_One_Length(t_edge *b, t_tree *tree)
{
  phydbl u;
  phydbl new_lnL_data, cur_lnL_data;
  phydbl ratio, alpha;
  phydbl new_l, cur_l;
  phydbl K,mult;

  if(tree->rates->u_ml_l[b->num] < 1.1*BL_MIN) return;

  cur_l        = b->l;
  cur_lnL_data = tree->c_lnL;
  K            = 1.;
  
  u = Uni();
  mult = EXP(K*(u-0.5));
  new_l = cur_l * mult;

  if(new_l < BL_MIN || new_l > BL_MAX) return;

  b->l = new_l;  
  new_lnL_data = Lk_At_Given_Edge(b,tree);
  /* new_lnL_data = Lk(tree); */

  ratio =
    (new_lnL_data - cur_lnL_data) +
    (LOG(mult));

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  u = Uni();
  
  if(u > alpha) /* Reject */
    {
      b->l = cur_l;
      tree->c_lnL = cur_lnL_data;
    }
}

/*********************************************************/

void MCMC_Br_Lens(t_tree *tree)
{
  MCMC_Br_Lens_Pre(tree->noeud[0],
		   tree->noeud[0]->v[0],
		   tree->noeud[0]->b[0],tree);
}

/*********************************************************/

void MCMC_Br_Lens_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  int i;
  
  if(!a->tax) Update_P_Lk(tree,b,a);

  MCMC_One_Length(b,tree);

  if(d->tax) return;
  else 
    {
      For(i,3) if(d->v[i] != a)
	{
	  MCMC_Br_Lens_Pre(d,d->v[i],d->b[i],tree);
	}
    }
  For(i,3) if((d->v[i] == a) && !(d->v[i]->tax)) Update_P_Lk(tree,d->b[i],d);
  
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
