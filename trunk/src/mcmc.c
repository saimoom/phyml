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
  int move;
  phydbl u;
  int i;
  int change;
  int node_num;

  if(tree->mcmc->randomize)
    {
      MCMC_Randomize_Nu(tree);
      MCMC_Randomize_Node_Times(tree);
      MCMC_Randomize_Rates(tree);
      MCMC_Randomize_Clock_Rate(tree); /* Clock Rate must be the last parameter to be randomized */
    }

  tree->rates->lk_approx = NORMAL;
  tree->both_sides       = NO;

  MCMC_Print_Param(tree->mcmc,tree);

  Update_Ancestors(tree->n_root,tree->n_root->v[0],tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);

  RATES_Update_Cur_Bl(tree);
  RATES_Lk_Rates(tree);
  Lk(tree);
  
  do
    {
      u = Uni();

      For(move,tree->mcmc->n_moves) if(u < tree->mcmc->move_weight[move]) break;
      
      /* Clock rate */
      if(move == tree->mcmc->num_move_clock_r)
	{
	  MCMC_Single_Param_Generic(&(tree->rates->clock_r),tree->rates->min_clock,tree->rates->max_clock,move,
				    NULL,&(tree->c_lnL),
				    NULL,Wrap_Lk,MCMC_MOVE_SCALE,NULL,tree,NULL);
	}


      /* Nu */
      else if(move == tree->mcmc->num_move_nu)
	{
	  MCMC_Single_Param_Generic(&(tree->rates->nu),tree->rates->min_nu,tree->rates->max_nu,move,
				    NULL,&(tree->rates->c_lnL),
				    NULL,Wrap_Lk_Rates,MCMC_MOVE_SCALE,NULL,tree,NULL);
	}


      /* Tree height */
      else if(move == tree->mcmc->num_move_tree_height)
	{
	  MCMC_Tree_Height(tree);
	  /* 	    MCMC_Swing(tree); */
	}


      /* Subtree height */
      else if(move == tree->mcmc->num_move_subtree_height)
	{
	  MCMC_SubTree_Height(tree);
	}


      /* Ts/tv ratio */
      else if(move == tree->mcmc->num_move_kappa)
	{
	  change = NO;
	  
	  if(tree->rates->lk_approx == NORMAL)
	    {
	      tree->rates->lk_approx = EXACT;
	      if(tree->mcmc->use_data == YES) Lk(tree);
	      change = YES;
	    }
	  
	  MCMC_Single_Param_Generic(&(tree->mod->kappa),0.0,100.,move,
				    NULL,&(tree->c_lnL),
				    NULL,Wrap_Lk,MCMC_MOVE_SCALE,NULL,tree,NULL);
	  
	  if(change == YES)
	    {
	      tree->rates->lk_approx = NORMAL;
	      if(tree->mcmc->use_data == YES) Lk(tree);
	    }
	}

      /* Times */
      else if(move >= tree->mcmc->num_move_nd_t) 
	{
	  node_num = move - tree->mcmc->num_move_nd_t + tree->n_otu;
	  MCMC_One_Time(tree->noeud[node_num]->anc,tree->noeud[node_num],NO,tree);
	}

      
      /* Rates */
      else if(move >= tree->mcmc->num_move_nd_r) 
	{

	  tree->both_sides = YES;
	  Lk(tree);
	  MCMC_One_Rate(tree->n_root,tree->n_root->v[0],YES,tree);
	  MCMC_One_Rate(tree->n_root,tree->n_root->v[1],YES,tree);

/* 	  RATES_Posterior_One_Rate(tree->n_root,tree->n_root->v[0],YES,tree); */
/* 	  RATES_Posterior_One_Rate(tree->n_root,tree->n_root->v[1],YES,tree); */
/* 	  RATES_Lk_Rates(tree); */
/* 	  if(tree->mcmc->use_data == YES) Lk(tree); */

/* /\* 	  MCMC_Sim_Rate(tree->n_root,tree->n_root->v[0],tree); *\/ */
/* /\* 	  MCMC_Sim_Rate(tree->n_root,tree->n_root->v[1],tree); *\/ */
/* /\* 	  RATES_Lk_Rates(tree); *\/ */
/* /\* 	  if(tree->mcmc->use_data == YES) Lk(tree); *\/ */

	}

      tree->mcmc->run++;
      MCMC_Adjust_Tuning_Parameter(tree->mcmc);
      RATES_Update_Norm_Fact(tree);
      RATES_Lk_Rates(tree);
      if(tree->mcmc->use_data == YES) Lk(tree);
      MCMC_Print_Param(tree->mcmc,tree);

/*       RATES_Update_Mean_Br_Len(tree->mcmc->run,tree); */
/*       RATES_Update_Cov_Br_Len(tree->mcmc->run,tree); */

/*       if(RATES_Check_Node_Times(tree)) */
/* 	{ */
/* 	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
/* 	  Warn_And_Exit(""); */
/* 	} */

    }
  while(tree->mcmc->run < tree->mcmc->chain_len);
}

/*********************************************************/


void MCMC_Single_Param_Generic(phydbl *val, 
			       phydbl lim_inf, 
			       phydbl lim_sup, 
			       int move_num,
			       phydbl *lnPrior,
			       phydbl *lnLike,
			       phydbl (*prior_func)(t_tree *), 
			       phydbl (*like_func)(t_edge *,t_tree *,supert_tree *),
			       int move_type,
			       t_edge *branch, t_tree *tree, supert_tree *stree)
{
  phydbl cur_val,new_val,new_lnLike,new_lnPrior,cur_lnLike,cur_lnPrior;
  phydbl u,alpha,ratio;
  phydbl min_nu,max_nu;
  phydbl K,mult;
  phydbl new_lnval, cur_lnval;

  Record_Br_Len(NULL,tree);

  cur_val       = *val;
  new_val       = -1.0;
  ratio         =  0.0;
  mult          = -1.0;
  K             = tree->mcmc->tune_move[move_num];
  cur_lnLike    = (lnLike)?(*lnLike):(UNLIKELY);
  cur_lnPrior   = (lnPrior)?(*lnPrior):(UNLIKELY);
  cur_lnval     = LOG(*val);

  u = Uni();

  if(move_type == MCMC_MOVE_SCALE)
    {
      /* Thorne's move. log(HR) = + log(mult) */
        mult    = EXP(K*(u-0.5));
        new_val = cur_val * mult;
    }
  else if(move_type == MCMC_MOVE_LOG_RANDWALK)
    {
      /* Random walk on the log scale. Watch out for the change of variable. */
      new_lnval = u * (2.*K) + cur_lnval - K;
      new_val   = EXP(new_lnval);
    }
  else if(move_type == MCMC_MOVE_RANDWALK)
    {
      /* Random walk on the original scale. log(HR) = 0 */
      new_val = u*(2.*K) + cur_val - K;
    }


  if(new_val > lim_sup || new_val < lim_inf)
    {
      tree->mcmc->run_move[move_num]++;
      return;
    }

  *val = new_val;
  
  RATES_Update_Cur_Bl(tree);

  ratio = 0.0;

  if(move_type == MCMC_MOVE_LOG_RANDWALK) ratio += (new_lnval - cur_lnval); /* Change of variable */
  if(move_type == MCMC_MOVE_SCALE)        ratio += LOG(mult); /* Hastings ratio */
  if(move_type == MCMC_MOVE_RANDWALK)     ratio += .0; /* Hastings ratio */ 

  if(prior_func) /* Prior ratio */
    { 
      new_lnPrior = (*prior_func)(tree); 
      ratio += (new_lnPrior - cur_lnPrior); 
    }
  if(like_func && tree->mcmc->use_data == YES)  /* Likelihood ratio */
    { 
      new_lnLike  = (*like_func)(branch,tree,stree);  
      ratio += (new_lnLike - cur_lnLike);  
    }

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  u = Uni();
  if(u > alpha) /* Reject */
    {
      *val     = cur_val;
      if(lnPrior) *lnPrior = cur_lnPrior;
      if(lnLike)  *lnLike  = cur_lnLike;
      Restore_Br_Len(NULL,tree);
    }
  else /* Accept */
    {
      tree->mcmc->acc_move[move_num]++;
      if(lnPrior) *lnPrior = new_lnPrior;
      if(lnLike)  *lnLike  = new_lnLike;
    }

  RATES_Update_Cur_Bl(tree);

  tree->mcmc->run_move[move_num]++;
}

/*********************************************************/

void MCMC_One_Rate(t_node *a, t_node *d, int traversal, t_tree *tree)
{
  t_edge *b;
  int i;

  b = NULL;
  if(a == tree->n_root) b = tree->e_root;
  else
    For(i,3) if(d->v[i] == a) { b = d->b[i]; break; }

  MCMC_Single_Param_Generic(&(tree->rates->nd_r[d->num]),
			    tree->rates->min_rate,
			    tree->rates->max_rate,
			    tree->mcmc->num_move_nd_r,
			    &(tree->rates->c_lnL),&(tree->c_lnL),
			    RATES_Lk_Rates,Wrap_Lk_At_Given_Edge,MCMC_MOVE_RANDWALK,b,tree,NULL);
/* 			    RATES_Lk_Rates,Wrap_Lk,MCMC_MOVE_RANDWALK,NULL,tree,NULL); */
  

/*   phydbl u; */
/*   phydbl new_lnL_data, cur_lnL_data, new_lnL_rate, cur_lnL_rate; */
/*   phydbl ratio, alpha; */
/*   phydbl new_mu, cur_mu; */
/*   phydbl r_min, r_max; */
/*   phydbl nu; */
/*   phydbl tmp; */
/*   phydbl K,mult; */
/*   int err; */
/*   phydbl new_logr, cur_logr; */

/*   nu           = tree->rates->nu; */
/*   cur_mu       = tree->rates->nd_r[d->num]; */
/*   cur_lnL_data = tree->c_lnL; */
/*   cur_lnL_rate = tree->rates->c_lnL; */
/*   r_min        = tree->rates->min_rate; */
/*   r_max        = tree->rates->max_rate; */
/*   K            = tree->mcmc->tune_move[d->num]; */
/*   tmp          = cur_lnL_rate; */
/*   new_lnL_data = UNLIKELY; */
/*   cur_logr     = LOG(cur_mu); */


/*   u = Uni(); */

/*   /\* !!!!!!!!!!!!!!! *\/ */
/* /\*   mult = EXP(K*(u-0.5)); *\/ */
/* /\*   mult = u*(K-1./K)+1./K; *\/ */
/* /\*   new_logr = cur_logr * mult; *\/ */

  
/*   /\* new_logr ~ U[cur_log-K; cur_logr+K]. Hastings ratio = 1 *\/ */
/*   new_logr = u * (2.*K) + cur_logr - K; */

/*   new_mu = EXP(new_logr); */

/*   if(new_mu < r_min || new_mu > r_max) */
/*     { */
/*       tree->mcmc->run_move[d->num]++; */
/*       return; */
/*     } */

/*   tree->rates->nd_r[d->num] = new_mu; */

/*   RATES_Update_Cur_Bl(tree); */

/*   if(tree->mcmc->use_data) new_lnL_data = Lk(tree); */
/*   new_lnL_rate = RATES_Lk_Rates(tree); */

/*   ratio = 0.0; */
/*   ratio += new_logr - cur_logr; */
/* /\*   ratio += LOG(mult); *\/ */
/* /\*   ratio -= LOG(mult); *\/ */
/*   ratio += (new_lnL_rate - cur_lnL_rate); */
/*   if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data); */

/*   ratio = EXP(ratio); */
/*   alpha = MIN(1.,ratio); */
  
/*   u = Uni(); */
  
/*   if(u > alpha) /\* Reject *\/ */
/*     { */
/*       tree->rates->nd_r[d->num] = cur_mu; */
/*       RATES_Update_Cur_Bl(tree); */
/*       tree->c_lnL        = cur_lnL_data; */
/*       tree->rates->c_lnL = cur_lnL_rate; */
/*     } */
/*   else */
/*     { */
/*       tree->mcmc->acc_move[d->num]++; */
/*     } */

/*   tree->mcmc->run_move[d->num]++; */






  if(traversal == YES)
    {
      if(d->tax == YES) return;
      else
	{
	  For(i,3)
	    if(d->v[i] != a)
	      {
		if(tree->rates->lk_approx == EXACT) Update_P_Lk(tree,d->b[i],d);
		MCMC_One_Rate(d,d->v[i],YES,tree);
	      }
	}
      if(tree->rates->lk_approx == EXACT) Update_P_Lk(tree,b,d);
    }
}

/*********************************************************/

void MCMC_One_Time(t_node *a, t_node *d, int traversal, t_tree *tree)
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
  int move_num;


  if(d->tax) return; /* Won't change time at tip */
  if(d == tree->n_root) return; /* Won't change time at root. Use tree height move instead */


  if(FABS(tree->rates->t_prior_min[d->num] - tree->rates->t_prior_max[d->num]) < 1.E-10) return;

/*   move_num       = d->num-tree->n_otu+tree->mcmc->num_move_nd_t; */
  move_num       = tree->mcmc->num_move_nd_t;
  K              = tree->mcmc->tune_move[move_num];
  cur_lnL_data   = tree->c_lnL;
  t1_cur         = tree->rates->nd_t[d->num];
  new_lnL_data   = UNLIKELY;

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
/*   mult = EXP(K*(u-0.5)); */
/*   t1_new = t1_cur * mult; */
  /* !!!!!!!!!!!!!!!!!!! */
  t1_new = u*(t_max-t_min) + t_min;

  if(t1_new < t_min || t1_new > t_max) 
    {
      tree->mcmc->run_move[move_num]++;
      return;
    }

  tree->rates->nd_t[d->num] = t1_new;
  RATES_Update_Cur_Bl(tree);
  if(tree->mcmc->use_data) new_lnL_data = Lk(tree);

  ratio = 0.0;
/*   ratio += LOG(mult); */
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
      tree->mcmc->acc_move[move_num]++;
      RATES_Lk_Rates(tree);
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
  
  tree->mcmc->run_move[move_num]++;

  if(traversal == YES)
    {
      if(d->tax == YES) return;
      else
	For(i,3)
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    MCMC_One_Time(d,d->v[i],YES,tree);
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
  int n_nodes;

  RATES_Record_Times(tree);

  K = tree->mcmc->tune_move[tree->mcmc->num_move_tree_height];
  cur_lnL_data = tree->c_lnL;
  new_lnL_data = UNLIKELY;

  u = Uni();
  mult = EXP(K*(u-0.5));
  
  floor = tree->rates->t_floor[tree->n_root->num];

  Scale_Subtree_Height(tree->n_root,mult,floor,&n_nodes,tree);

  For(i,2*tree->n_otu-1)
    {
      if(tree->rates->nd_t[i] > tree->rates->t_prior_max[i] ||
  	 tree->rates->nd_t[i] < tree->rates->t_prior_min[i])
  	{
  	  RATES_Reset_Times(tree);
	  tree->mcmc->run_move[tree->mcmc->num_move_tree_height]++;
  	  return;
  	}
    }
  
  RATES_Update_Cur_Bl(tree);
  if(tree->mcmc->use_data) new_lnL_data = Lk(tree);

  /* The Hastings ratio is actually mult^(n) when changing the absolute
     node heights. When considering the relative heights, this ratio combined
     to the Jacobian for the change of variable ends up to mult. 
  */
  ratio = 0.0;
  ratio += LOG(mult);
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
      tree->mcmc->acc_move[tree->mcmc->num_move_tree_height]++;
      RATES_Lk_Rates(tree);
    }

  tree->mcmc->run_move[tree->mcmc->num_move_tree_height]++;
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
  int n_nodes;

  RATES_Record_Times(tree);

  K = tree->mcmc->tune_move[tree->mcmc->num_move_subtree_height];
  cur_lnL_data = tree->c_lnL;
  new_lnL_data = UNLIKELY;
  
  u = Uni();
  mult = EXP(K*(u-0.5));

  target = Rand_Int(tree->n_otu,2*tree->n_otu-3);

  floor = tree->rates->t_floor[target];

  if(tree->noeud[target] == tree->n_root)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }


  if(!Scale_Subtree_Height(tree->noeud[target],mult,floor,&n_nodes,tree))
    {
      RATES_Reset_Times(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_subtree_height]++;
      return;
    }

  For(i,2*tree->n_otu-1)
    {
      if(tree->rates->nd_t[i] > tree->rates->t_prior_max[i] ||
  	 tree->rates->nd_t[i] < tree->rates->t_prior_min[i])
  	{
  	  RATES_Reset_Times(tree);
	  tree->mcmc->run_move[tree->mcmc->num_move_subtree_height]++;
  	  return;
  	}
    }
  
  RATES_Update_Cur_Bl(tree);
  if(tree->mcmc->use_data) new_lnL_data = Lk(tree);

  /* The Hastings ratio here is mult^(n_nodes) and the ratio of the prior joint densities
     of the modified node heigths given the unchanged one is 1. This is didifferent from the 
     case where all the nodes, including the root node, are scaled. 
  */
  ratio = 0.0;
  ratio += (n_nodes) * LOG(mult);
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
      RATES_Lk_Rates(tree);
      tree->mcmc->acc_move[tree->mcmc->num_move_subtree_height]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_subtree_height]++;
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

  K = tree->mcmc->tune_move[tree->mcmc->num_move_tree_height];
  cur_lnL_data  = tree->c_lnL;
  new_lnL_data  = UNLIKELY;

  u = Uni();
  mult = EXP(K*(u-0.5));
  
  floor = tree->rates->t_floor[tree->n_root->num];

  Scale_Subtree_Height(tree->n_root,mult,floor,tree);

  For(i,2*tree->n_otu-1)
    {
      if(tree->rates->nd_t[i] > tree->rates->t_prior_max[i] ||
  	 tree->rates->nd_t[i] < tree->rates->t_prior_min[i])
  	{
  	  RATES_Reset_Times(tree);
	  tree->mcmc->run_move[tree->mcmc->num_move_tree_height]++;
  	  return;
  	}
    }
  
  tree->rates->clock_r /= mult;
  
  if(tree->rates->clock_r > tree->rates->max_clock || 
     tree->rates->clock_r < tree->rates->min_clock)
    {
      tree->rates->clock_r *= mult;
      RATES_Reset_Times(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_tree_height]++;
      return;      
    }

  RATES_Update_Cur_Bl(tree);
  if(tree->mcmc->use_data) new_lnL_data = Lk(tree);

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
      tree->mcmc->acc_move[tree->mcmc->num_move_tree_height]++;
      RATES_Lk_Rates(tree);
    }

  tree->mcmc->run_move[tree->mcmc->num_move_tree_height]++;
  
}

/*********************************************************/

void MCMC_Print_Param(t_mcmc *mcmc, t_tree *tree)
{
  int i;
  FILE *fp;
  phydbl eps = 0.01;
  char *s;

  if(tree->mcmc->run > mcmc->chain_len) return;

  

  s = (char *)mCalloc(100,sizeof(char));

  fp = mcmc->out_fp_stats;
  
  if(tree->mcmc->run == 0)
    {
      PhyML_Fprintf(stdout," [");
      fflush(NULL);
    }

  if(!(mcmc->run%(mcmc->chain_len/10)))
    {
      PhyML_Fprintf(stdout,".");
      fflush(NULL);
    }

  if(tree->mcmc->run == mcmc->chain_len)
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
/* 	  PhyML_Fprintf(fp,"Time\t"); */
/* 	  PhyML_Fprintf(fp,"MeanRate\t"); */
	  PhyML_Fprintf(fp,"NormFact\t");

	  For(i,mcmc->n_moves)
	    {
	      strcpy(s,"Acc.");
	      PhyML_Fprintf(fp,"%s\t",strcat(s,mcmc->move_name[i]));
	    }

	  For(i,mcmc->n_moves)
	    {
	      strcpy(s,"Tune.");
	      PhyML_Fprintf(fp,"%s\t",strcat(s,mcmc->move_name[i]));
	    }

	  PhyML_Fprintf(fp,"LnL.Seq.%s\t",(tree->rates->lk_approx == NORMAL)?"Exact":"Approx");
	  PhyML_Fprintf(fp,"LnL.Seq.%s\t",(tree->rates->lk_approx == NORMAL)?"Approx":"Exact");
	  PhyML_Fprintf(fp,"LnPrior.Rate\t");
	  PhyML_Fprintf(fp,"Clock.Rate\t");
	  PhyML_Fprintf(fp,"Average.Rate\t");
	  PhyML_Fprintf(fp,"Nu\t");
	  PhyML_Fprintf(fp,"TsTv\t");

	  if(fp != stdout)
	    {
	      for(i=tree->n_otu;i<2*tree->n_otu-1;i++)
	  	{
/* 	  	  PhyML_Fprintf(fp,"%sT%d[%d,%d]\t", */
/* 				tree->noeud[i] == tree->n_root?"*":"", */
/* 	  			i,tree->noeud[i]->rank,tree->noeud[i]->rank_max); */
/* 	  	  PhyML_Fprintf(fp,"T%d[%7.1f;%7.1f]\t",i,tree->rates->t_prior_min[i],tree->rates->t_prior_max[i]); */
	  	  PhyML_Fprintf(fp,"T%d\t",i);
	  	}
	    }


	  if(fp != stdout)
	    {
	      For(i,2*tree->n_otu-2)
	  	{
	  	  if(tree->noeud[i] == tree->n_root->v[0])
	  	    PhyML_Fprintf(fp,"0R%d\t",i);
	  	  else if(tree->noeud[i] == tree->n_root->v[1])
	  	    PhyML_Fprintf(fp,"1R%d\t",i);
	  	  else
	  	    PhyML_Fprintf(fp," R%d\t",i);

/* 		  PhyML_Fprintf(fp," R%d[%f]\t",i,tree->rates->mean_l[i]); */
	  	}
	    }


/* 	  if(fp != stdout) */
/* 	    { */
/* 	      For(i,2*tree->n_otu-3) */
/* 	  	{ */
/* 		  PhyML_Fprintf(fp," L[%f]%d\t",i,tree->rates->u_ml_l[i]); */
/* 	  	} */
/* 	    } */


	  PhyML_Fprintf(mcmc->out_fp_trees,"#NEXUS\n");
	  PhyML_Fprintf(mcmc->out_fp_trees,"BEGIN TREES;\n");
	  PhyML_Fprintf(mcmc->out_fp_trees,"\tTRANSLATE\n");
	  For(i,tree->n_otu-1) PhyML_Fprintf(mcmc->out_fp_trees,"\t%3d\t%s,\n",tree->noeud[i]->num,tree->noeud[i]->name);
	  PhyML_Fprintf(mcmc->out_fp_trees,"\t%3d\t%s;\n",tree->noeud[i]->num,tree->noeud[i]->name);
	  tree->write_tax_names = NO;
	}

      PhyML_Fprintf(fp,"\n");

      PhyML_Fprintf(fp,"%6d\t",tree->mcmc->run);

/*       time(&mcmc->t_cur); */
/*       PhyML_Fprintf(fp,"%6d\t",(int)(mcmc->t_cur-mcmc->t_beg)); */
      
/*       RATES_Update_Cur_Bl(tree); */
/*       PhyML_Fprintf(fp,"%f\t",RATES_Check_Mean_Rates(tree)); */

      PhyML_Fprintf(fp,"%f\t",tree->rates->norm_fact);

      For(i,tree->mcmc->n_moves)
      	{
	  PhyML_Fprintf(fp,"%f\t",
			(phydbl)((tree->mcmc->acc_move[i]+eps)/(tree->mcmc->run_move[i]+eps)));
	}

      For(i,tree->mcmc->n_moves) PhyML_Fprintf(fp,"%f\t",(phydbl)(tree->mcmc->tune_move[i]));


      tree->rates->lk_approx = (tree->rates->lk_approx == NORMAL)?EXACT:NORMAL;
      Lk(tree);
      PhyML_Fprintf(fp,"%.1f\t",(tree->mcmc->use_data)?tree->c_lnL:0.0);
      tree->rates->lk_approx = (tree->rates->lk_approx == NORMAL)?EXACT:NORMAL;
      Lk(tree);
      PhyML_Fprintf(fp,"%.1f\t",(tree->mcmc->use_data)?tree->c_lnL:0.0);

/*       PhyML_Fprintf(fp,"0\t0\t"); */


      PhyML_Fprintf(fp,"%G\t",tree->rates->c_lnL);
      PhyML_Fprintf(fp,"%G\t",tree->rates->clock_r);
      PhyML_Fprintf(fp,"%G\t",RATES_Average_Substitution_Rate(tree));
      PhyML_Fprintf(fp,"%G\t",tree->rates->nu);
      PhyML_Fprintf(fp,"%G\t",tree->mod->kappa);
      for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(fp,"%.1f\t",tree->rates->nd_t[i]);
      for(i=0;i<2*tree->n_otu-2;i++) PhyML_Fprintf(fp,"%.4f\t",tree->rates->nd_r[i]);
      /* if(fp != stdout) for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(fp,"%G\t",tree->rates->t_prior[i]); */
/*       For(i,2*tree->n_otu-3) PhyML_Fprintf(fp,"%f\t",EXP(tree->t_edges[i]->l)); */
/*       For(i,2*tree->n_otu-3) PhyML_Fprintf(fp,"%f\t",tree->t_edges[i]->l); */
      fflush(NULL);


      // TREES
/*       char *s_tree; */
/* /\*       Branch_Lengths_To_Time_Lengths(tree); *\/ */
/*       Branch_Lengths_To_Rate_Lengths(tree); */
/*       tree->bl_ndigits = 5; */
/* /\*       tree->bl_ndigits = 0; *\/ */
/*       s_tree = Write_Tree(tree); */
/*       tree->bl_ndigits = 7; */
/*       PhyML_Fprintf(mcmc->out_fp_trees,"TREE %8d [%f] = [&R] %s\n",mcmc->run,tree->c_lnL,s_tree); */
/*       Free(s_tree); */
/*       RATES_Update_Cur_Bl(tree); */
    }

  if(tree->mcmc->run == mcmc->chain_len)
    {
      PhyML_Fprintf(mcmc->out_fp_trees,"END;\n");
      fflush(NULL); 
    }
  
  Free(s);
  
}

/*********************************************************/

void MCMC_Print_Means(t_mcmc *mcmc, t_tree *tree)
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

void MCMC_Print_Last(t_mcmc *mcmc, t_tree *tree)
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
t_mcmc *MCMC_Make_MCMC_Struct()
{
  t_mcmc *mcmc;
  int i;

  mcmc               = (t_mcmc *)mCalloc(1,sizeof(t_mcmc));
  mcmc->out_filename = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  mcmc->move_weight  = (phydbl *)mCalloc(N_MAX_MOVES,sizeof(phydbl));

  return(mcmc);
}

/*********************************************************/

void MCMC_Free_MCMC(t_mcmc *mcmc)
{
  int i;

  Free(mcmc->out_filename);
  Free(mcmc->move_weight);
  Free(mcmc->acc_move);
  Free(mcmc->run_move);
  Free(mcmc->tune_move);
  For(i,mcmc->n_moves) Free(mcmc->move_name[i]);
  Free(mcmc->move_name);
  Free(mcmc);
}

/*********************************************************/

void MCMC_Init_MCMC_Struct(char *filename, t_mcmc *mcmc)
{
  int pid;
  int i;
  phydbl sum;

  mcmc->use_data         = YES;
  mcmc->adjust_tuning    = YES;
  mcmc->run              = 0;
  mcmc->sample_interval  = 1E+3;
  mcmc->chain_len        = 1E+6;
  mcmc->chain_len_burnin = 1E+5;
  mcmc->randomize        = 1;
  mcmc->norm_freq        = 1000;

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

      mcmc->out_fp_stats = fopen(mcmc->out_filename,"w");

      strcpy(s,mcmc->out_filename);
      strcat(s,".trees");
      mcmc->out_fp_trees = fopen(s,"w");

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
      mcmc->out_fp_stats = stderr;
      mcmc->out_fp_trees = stderr;
      /* tree->mcmc->out_fp_means = stderr; */
      /* tree->mcmc->out_fp_last  = stderr; */
    }
}

/*********************************************************/

void MCMC_Copy_MCMC_Struct(t_mcmc *ori, t_mcmc *cpy, char *filename)
{
  int pid;
  int i;
  
  cpy->use_data           = ori->use_data        ;
  cpy->adjust_tuning      = ori->adjust_tuning   ;
  cpy->sample_interval    = ori->sample_interval ;
  cpy->chain_len          = ori->chain_len       ;
  cpy->randomize          = ori->randomize       ;
  cpy->norm_freq          = ori->norm_freq       ;
  cpy->n_moves            = ori->n_moves         ;

  For(i,cpy->n_moves) 
    {
      cpy->move_weight[i] = ori->move_weight[i];
      cpy->run_move[i]   = ori->run_move[i];
      cpy->acc_move[i]   = ori->acc_move[i];
      cpy->tune_move[i]  = ori->tune_move[i];
      strcpy(cpy->move_name[i],ori->move_name[i]);
    }

  
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

void MCMC_Close_MCMC(t_mcmc *mcmc)
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

  if(tree->mod->log_l == NO)
    For(i,2*tree->n_otu-3) tree->t_edges[i]->l = Rexp(10.);
  else
    For(i,2*tree->n_otu-3) tree->t_edges[i]->l = -4* Uni();
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
  max_nu = 1.0;

  u = Uni();
  tree->rates->nu = (1.0 - min_nu) * u + min_nu;
}

/*********************************************************/

void MCMC_Randomize_Clock_Rate(t_tree *tree)
{
  phydbl u;
  u = Uni();
  tree->rates->clock_r = u * (tree->rates->max_clock - tree->rates->min_clock) + tree->rates->min_clock;
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
  int i;
  phydbl dt,min_dt;
  int min_node;

  t_inf = tree->rates->t_prior_min[tree->n_root->num];
  t_sup = tree->rates->t_prior_max[tree->n_root->num];

  u = Uni();
  u *= (t_sup - t_inf);
  u += t_inf;
  
  tree->rates->nd_t[tree->n_root->num] = u;


  MCMC_Randomize_Node_Times_Top_Down(tree->n_root,tree->n_root->v[0],tree);
  MCMC_Randomize_Node_Times_Top_Down(tree->n_root,tree->n_root->v[1],tree);
  
  
  iter = 0;
  do
    {
      min_dt = 1E+5;
      For(i,2*tree->n_otu-2) 
	{
	  dt = tree->rates->nd_t[i] - tree->rates->nd_t[tree->noeud[i]->anc->num];
	  if(dt < min_dt)
	    {
	      min_dt = dt;
	      min_node = i;
	    }
	}

      if(min_dt > -.1 * tree->rates->nd_t[tree->n_root->num]/(phydbl)(tree->n_otu-1)) break;

      MCMC_Randomize_Node_Times_Bottom_Up(tree->n_root,tree->n_root->v[0],tree);
      MCMC_Randomize_Node_Times_Bottom_Up(tree->n_root,tree->n_root->v[1],tree);
      
      iter++;
    }
  while(iter < 200);

  if(iter == 200)
    {      
      PhyML_Printf("\n. min_dt = %f",min_dt);
      PhyML_Printf("\n. min->t=%f min->anc->t=%f",tree->rates->nd_t[min_node],tree->rates->nd_t[tree->noeud[min_node]->anc->num]);
      PhyML_Printf("\n. up=%f down=%f",tree->rates->t_prior_min[min_node],tree->rates->t_floor[tree->noeud[min_node]->anc->num]);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }


/*   PhyML_Printf("\n. Needed %d iterations to randomize node heights.",iter); */

/*   TIMES_Print_Node_Times(tree->n_root,tree->n_root->v[0],tree); */
/*   TIMES_Print_Node_Times(tree->n_root,tree->n_root->v[1],tree); */

  if(RATES_Check_Node_Times(tree))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
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
      t_inf = tree->rates->nd_t[a->num];

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

      t_inf = MAX(tree->rates->nd_t[a->num],tree->rates->t_prior_min[d->num]);
      t_sup = tree->rates->t_prior_max[d->num];

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


void MCMC_Adjust_Tuning_Parameter(t_mcmc *mcmc)
{
  if(mcmc->adjust_tuning)
    {


      phydbl eps=0.01;
      phydbl scale=1.1;
      phydbl rate;
      phydbl rate_inf,rate_sup;
      int i;

      rate_inf = 0.1;
      rate_sup = 0.5;
      
      For(i,mcmc->n_moves)
	{
	  if(mcmc->run_move[i] > 100)
	    {
	      /* PhyML_Printf("\n. %s acc=%d run=%d tune=%f", */
	      /* 		   mcmc->move_name[i], */
	      /* 		   mcmc->acc_move[i], */
	      /* 		   mcmc->run_move[i], */
	      /* 		   mcmc->tune_move[i]); */

	      rate = (phydbl)(mcmc->acc_move[i]+eps)/(mcmc->run_move[i]+eps); 
	      if(rate < rate_inf)
		{
		  mcmc->tune_move[i] /= scale;
		}
	      else if(rate > rate_sup)
		{
		  mcmc->tune_move[i] *= scale;
		}
	      
	      mcmc->run_move[i] = 0;
	      mcmc->acc_move[i] = 0;
	    }
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


  cur_l        = b->l;
  cur_lnL_data = tree->c_lnL;
  K            = 0.1;
  
  u = Uni();
  mult = EXP(K*(u-0.5));
  /* mult = u*(K-1./K)+1./K; */
  new_l = cur_l * mult;

  if(new_l < tree->mod->l_min || new_l > tree->mod->l_max) return;

  b->l = new_l;
  new_lnL_data = Lk_At_Given_Edge(b,tree);
/*   tree->both_sides = NO; */
/*   new_lnL_data = Lk(tree); */


  ratio =
    (new_lnL_data - cur_lnL_data) +
    (LOG(mult));


  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  u = Uni();
  
  if(u > alpha) /* Reject */
    {
      b->l = cur_l;
      Update_PMat_At_Given_Edge(b,tree);
      tree->c_lnL = cur_lnL_data;
    }

}

/*********************************************************/

void MCMC_Scale_Br_Lens(t_tree *tree)
{
  phydbl u;
  phydbl new_lnL_data, cur_lnL_data;
  phydbl ratio, alpha;
  phydbl K,mult;
  int i;

  Record_Br_Len(NULL,tree);

  cur_lnL_data = tree->c_lnL;
  K            = 1.2;
  
  u = Uni();
  mult = u*(K-1./K)+1./K;

  For(i,2*tree->n_otu-3) 
    {
      tree->t_edges[i]->l *= mult;
      if(tree->t_edges[i]->l < tree->mod->l_min || 
	 tree->t_edges[i]->l > tree->mod->l_max) return;
    }

  tree->both_sides = NO;
  new_lnL_data = Lk(tree);

  ratio =
    (new_lnL_data - cur_lnL_data) +
    (2*tree->n_otu-5) * (LOG(mult));

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  u = Uni();
  
  if(u > alpha) /* Reject */
    {
      Restore_Br_Len(NULL,tree);
      tree->c_lnL = cur_lnL_data;
    }
}

/*********************************************************/

void MCMC_Br_Lens(t_tree *tree)
{
  MCMC_Br_Lens_Pre(tree->noeud[0],
  		   tree->noeud[0]->v[0],
  		   tree->noeud[0]->b[0],tree);

  /* int i; */
  /* For(i,2*tree->n_otu-3) */
  /*   { */
  /*     MCMC_One_Length(tree->t_edges[Rand_Int(0,2*tree->n_otu-4)],acc,run,tree); */
  /*   } */
}

/*********************************************************/

void MCMC_Br_Lens_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  int i;


  if(a == tree->n_root || d == tree->n_root)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  MCMC_One_Length(b,tree);
  if(d->tax) return;
  else 
    {
      For(i,3) 
	if(d->v[i] != a)
	  {
	    Update_P_Lk(tree,d->b[i],d);
	    MCMC_Br_Lens_Pre(d,d->v[i],d->b[i],tree);
	  }
      Update_P_Lk(tree,b,d);
    }  
}

/*********************************************************/

void MCMC_Nu(t_tree *tree)
{
  phydbl cur_nu,new_nu,cur_lnL_rate,new_lnL_rate;
  phydbl u,alpha,ratio;
  phydbl min_nu,max_nu;
  phydbl K,mult;

  cur_nu        = -1.0;
  new_nu        = -1.0;
  ratio         = -1.0;

  K = tree->mcmc->tune_move[tree->mcmc->num_move_nu];

  cur_lnL_rate = tree->rates->c_lnL;
  cur_nu       = tree->rates->nu;
  
  min_nu = tree->rates->min_nu;
  max_nu = tree->rates->max_nu;

  u = Uni();
  mult = EXP(K*(u-0.5));
  new_nu = cur_nu * mult;

  if(new_nu > max_nu || new_nu < min_nu) return;

  tree->rates->nu = new_nu;

  new_lnL_rate = RATES_Lk_Rates(tree);

  ratio = LOG(mult);
  ratio +=
    (new_lnL_rate - cur_lnL_rate);

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  u = Uni();
  if(u > alpha) /* Reject */
    {
      tree->rates->nu    = cur_nu;
      tree->rates->c_lnL = cur_lnL_rate;
    }
}

/*********************************************************/


/* Only works when simulating from prior */
void MCMC_Sim_Rate(t_node *a, t_node *d, t_tree *tree)
{
  int err;

  tree->rates->nd_r[d->num] = Rnorm_Trunc(tree->rates->nd_r[a->num],
					  SQRT(tree->rates->nu * (tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num])),
					  tree->rates->min_rate,
					  tree->rates->max_rate,
					  &err);
  
  if(d->tax) return;
  else
    {
      int i;

      For(i,3)
	if(d->v[i] != a && d->b[i] != tree->e_root)
	  MCMC_Sim_Rate(d,d->v[i],tree);
    }
}

/*********************************************************/

void MCMC_Complete_MCMC(t_mcmc *mcmc, t_tree *tree)
{
  int i;
  phydbl sum;

  mcmc->n_moves = 0;
  mcmc->num_move_nd_r = mcmc->n_moves;
  mcmc->n_moves += 2*tree->n_otu-2;
  mcmc->num_move_nd_t = mcmc->n_moves;
  mcmc->n_moves += tree->n_otu-1;

  mcmc->num_move_nu             = mcmc->n_moves++;
  mcmc->num_move_clock_r        = mcmc->n_moves++;
  mcmc->num_move_tree_height    = mcmc->n_moves++;
  mcmc->num_move_subtree_height = mcmc->n_moves++;
  mcmc->num_move_kappa          = mcmc->n_moves++;


  mcmc->run_move = (int *)mCalloc(mcmc->n_moves,sizeof(int));
  mcmc->acc_move = (int *)mCalloc(mcmc->n_moves,sizeof(int));

  mcmc->move_name = (char **)mCalloc(mcmc->n_moves,sizeof(char *));
  For(i,mcmc->n_moves) mcmc->move_name[i] = (char *)mCalloc(50,sizeof(char));

  For(i,2*tree->n_otu-2) strcpy(mcmc->move_name[i],"rate");  
  for(i=2*tree->n_otu-2;i<tree->n_otu+1+2*tree->n_otu-2;i++) strcpy(mcmc->move_name[i],"time");
  strcpy(mcmc->move_name[mcmc->num_move_nu],"nu");
  strcpy(mcmc->move_name[mcmc->num_move_clock_r],"clock");
  strcpy(mcmc->move_name[mcmc->num_move_tree_height],"tree_height");
  strcpy(mcmc->move_name[mcmc->num_move_subtree_height],"subtree_height");
  strcpy(mcmc->move_name[mcmc->num_move_kappa],"kappa");
  
  mcmc->tune_move = (phydbl *)mCalloc(mcmc->n_moves,sizeof(phydbl));
  For(i,2*tree->n_otu-2) mcmc->tune_move[i] = 0.05; /* Rates */
  for(i= 2*tree->n_otu-2; i < tree->n_otu+1+2*tree->n_otu-2; i++) mcmc->tune_move[i] = 1.0;  /* Times */
  mcmc->tune_move[mcmc->num_move_clock_r]         = 0.1;
  mcmc->tune_move[mcmc->num_move_tree_height]     = 0.1;
  mcmc->tune_move[mcmc->num_move_subtree_height]  = 0.5;
  mcmc->tune_move[mcmc->num_move_nu]              = 0.1;
  mcmc->tune_move[mcmc->num_move_kappa]           = 0.5;   

  
  mcmc->move_weight = (phydbl *)mCalloc(mcmc->n_moves,sizeof(phydbl));

  For(i,2*tree->n_otu-2) mcmc->move_weight[i] = (phydbl)(1./(2.*tree->n_otu-2)); /* Rates */
  for(i= 2*tree->n_otu-2; i < tree->n_otu+1+2*tree->n_otu-2; i++) mcmc->move_weight[i] = (phydbl)(1./(tree->n_otu-2));  /* Times */
  mcmc->move_weight[mcmc->num_move_clock_r]         = 1.0;
  mcmc->move_weight[mcmc->num_move_tree_height]     = 1.0;
  mcmc->move_weight[mcmc->num_move_subtree_height]  = 1.0;
  mcmc->move_weight[mcmc->num_move_nu]              = 1.0;
  mcmc->move_weight[mcmc->num_move_kappa]           = 0.5;   
 
  sum = 0.0;
  For(i,mcmc->n_moves) sum += mcmc->move_weight[i];
  For(i,mcmc->n_moves) mcmc->move_weight[i] /= sum;
  for(i=1;i<mcmc->n_moves;i++) mcmc->move_weight[i] += mcmc->move_weight[i-1];

}

/*********************************************************/
/*********************************************************/
