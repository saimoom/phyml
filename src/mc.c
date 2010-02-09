/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/


/* Routines for molecular clock trees and molecular dating */


#include "mc.h"

/*********************************************************/
#ifdef MC

int MC_main(int argc, char **argv)
{
  align **data;
  calign *cdata;
  option *io;
  t_tree *tree;
  int n_otu, num_data_set;
  int num_tree,tree_line_number,num_rand_tree;
  matrix *mat;
  model *mod;
  time_t t_beg,t_end;
  phydbl best_lnL,most_likely_size,tree_size;
  int r_seed;
  char *most_likely_tree;

  
#ifdef MPI
  int rc;
  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    PhyML_Printf("\n. Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  MPI_Comm_size(MPI_COMM_WORLD,&Global_numTask);
  MPI_Comm_rank(MPI_COMM_WORLD,&Global_myRank);
#endif

#ifdef QUIET
  setvbuf(stdout,NULL,_IOFBF,2048);
#endif

  tree             = NULL;
  mod              = NULL;
  data             = NULL;
  most_likely_tree = NULL;
  best_lnL         = UNLIKELY;
  most_likely_size = -1.0;
  tree_size        = -1.0;

  io = (option *)Get_Input(argc,argv);
  r_seed = (io->r_seed < 0)?(time(NULL)):(io->r_seed);
/*   r_seed = 1248325214; */
/*   r_seed = 1254029112; */
/*   r_seed = 1254094805; */
/*   r_seed = 1254097617; */
/*   r_seed = 1254103645; */
/*   r_seed = 1262843891; */
/*   r_seed = 1262925254; */
/*   r_seed = 1263111043; */
/*   r_seed = 1263158433; */
/*   r_seed = 1263168633; */
/*   r_seed = 1263372871; */
/*   r_seed = 1263375598; */
/*   r_seed = 1263781194; */
/*   r_seed = 1263894098; */
/*   r_seed = 1264014238; */
/*   r_seed = 1264121265; */
/*   r_seed = 1264246773; */
/*   r_seed = 1264249283; */
/*   r_seed = 1264319131; */
/*   r_seed = 1264735233; */
  /* !!!!!!!!!!!!!!!!!!!!!!!! */

  srand(r_seed); rand();
  PhyML_Printf("\n. Seed: %d\n",r_seed);
  PhyML_Printf("\n. Pid: %d\n",getpid());
  Make_Model_Complete(io->mod);
  mod = io->mod;
  if(io->in_tree == 2) Test_Multiple_Data_Set_Format(io);
  else io->n_trees = 1;

  io->colalias = 1;  /* Do not compress sites if you're using Evolve function */

  mat = NULL;
  tree_line_number = 0;


  if((io->n_data_sets > 1) && (io->n_trees > 1))
    {
      io->n_data_sets = MIN(io->n_trees,io->n_data_sets);
      io->n_trees     = MIN(io->n_trees,io->n_data_sets);
    }

  For(num_data_set,io->n_data_sets)
    {
      n_otu = 0;
      best_lnL = UNLIKELY;
      data = Get_Seq(io);

      if(data)
	{
	  if(io->n_data_sets > 1) PhyML_Printf("\n. Data set [#%d]\n",num_data_set+1);
	  PhyML_Printf("\n. Compressing sequences...\n");
	  cdata = Compact_Data(data,io);
	  Free_Seq(data,cdata->n_otu);
	  Check_Ambiguities(cdata,io->mod->io->datatype,io->mod->state_len);

	  for(num_tree=(io->n_trees == 1)?(0):(num_data_set);num_tree < io->n_trees;num_tree++)
	    {
	      if(!io->mod->s_opt->random_input_tree) io->mod->s_opt->n_rand_starts = 1;

	      For(num_rand_tree,io->mod->s_opt->n_rand_starts)
		{
		  if((io->mod->s_opt->random_input_tree) && (io->mod->s_opt->topo_search != NNI_MOVE))
		    PhyML_Printf("\n. [Random start %3d/%3d]\n",num_rand_tree+1,io->mod->s_opt->n_rand_starts);

		  Init_Model(cdata,mod,io);

		  /* A BioNJ tree is built here */
		  if(!io->in_tree) tree = Dist_And_BioNJ(cdata,mod,io);
		  /* A user-given tree is used here instead of BioNJ */
		  else             tree = Read_User_Tree(cdata,mod,io);

		  if(!tree) continue;

		  time(&t_beg);
		  time(&(tree->t_beg));

		  int n_otu;
		  int i;


/* 		  n_otu = 15; */
/* 		  tree = Generate_Random_Tree_From_Scratch(n_otu,1); */

		  tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
		  RATES_Init_Rate_Struct(tree->rates,tree->n_otu);

		  Update_Ancestors(tree->n_root,tree->n_root->v[0],tree);
		  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
		  tree->n_root->anc = NULL;
		  
		  RATES_Fill_Lca_Table(tree);

		  tree->mod         = mod;
		  tree->io          = io;
		  tree->data        = cdata;
		  tree->both_sides  = 1;
		  tree->n_pattern   = tree->data->crunch_len/tree->mod->state_len;

/*  		  For(i,tree->n_otu) strcpy(tree->noeud[i]->name,cdata->c_seq[i]->name); */

		  Fill_Dir_Table(tree);
		  Update_Dirs(tree);
		  Make_Tree_4_Pars(tree,cdata,cdata->init_len);
		  Make_Tree_4_Lk(tree,cdata,cdata->init_len);

/* 		  Evolve(tree->data,tree->mod,tree); /\* do not forget to make sure sequences are not compressed! *\/ */
		  Init_Ui_Tips(tree);
		  Init_P_Pars_Tips(tree);
		  if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(tree);
		  else Init_P_Lk_Tips_Int(tree);

		  For(i,2*tree->n_otu-1) tree->rates->true_t[i] = tree->rates->nd_t[i];
		  For(i,2*tree->n_otu-2) tree->rates->true_r[i] = tree->rates->nd_r[i];

/* 		  tree->data->format = 1; */
/* 		  Print_CSeq(stdout,tree->data); */
		  char *s;
/* 		  Branch_Lengths_To_Time_Lengths(tree); */
		  s = Write_Tree(tree);
		  PhyML_Fprintf(stdout,"TREE %8d [%f] = %s\n",0,0.0,s);
		  Free(s);

		  Read_Clade_Priors(io->clade_list_file,tree);


		  /************************************/
/* 		  IMPORTANCE SAMPLING STUFF */
		  t_node *buff;

		  buff = tree->n_root;
		  tree->n_root = NULL;
		  PhyML_Printf("\n");
		  Round_Optimize(tree,tree->data,ROUND_MAX);
		  tree->n_root = buff;

		  Record_Br_Len(NULL,tree);
		  PhyML_Printf("\n");
		  PhyML_Printf("\n. Computing Hessian...\n");
		  tree->rates->bl_from_rt = 0;
		  phydbl *cov;
		  cov = Hessian(tree);
		  For(i,(2*tree->n_otu-3)*(2*tree->n_otu-3)) cov[i] *= -1.0;
		  if(!Matinv(cov,2*tree->n_otu-3,2*tree->n_otu-3))
		    {
		      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		      Exit("\n");      
		    }
		  Free(tree->rates->cov); 
		  tree->rates->cov = cov;
		  For(i,(2*tree->n_otu-3)*(2*tree->n_otu-3)) tree->rates->invcov[i] = tree->rates->cov[i];
		  if(!Matinv(tree->rates->invcov,2*tree->n_otu-3,2*tree->n_otu-3))
		    {
		      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		      Exit("\n");      
		    }
		  tree->rates->covdet = Matrix_Det(tree->rates->cov,2*tree->n_otu-3,YES);
		  Restore_Br_Len(NULL,tree);

		  tree->rates->true_tree_size = Get_Tree_Size(tree);
		  RATES_Bl_To_Ml(tree);
		  RATES_Get_Conditional_Variances(tree);
		  RATES_Get_All_Reg_Coeff(tree);
		  RATES_Get_Trip_Conditional_Variances(tree);
		  RATES_Get_All_Trip_Reg_Coeff(tree);


		  Lk(tree);
		  PhyML_Printf("\n. p(data|model) [exact] ~ %f",tree->c_lnL);
		  For(i,2*tree->n_otu-3) tree->rates->u_cur_l[i] = tree->t_edges[i]->l;
		  tree->c_lnL = Dnorm_Multi_Given_InvCov_Det(tree->rates->u_cur_l,tree->rates->u_ml_l,tree->rates->invcov,tree->rates->covdet,2*tree->n_otu-3,YES);
		  PhyML_Printf("\n. p(data|model) [normal approx] ~ %f",tree->c_lnL);
/* 		  PhyML_Printf("\n. p(model) [normal approx] ~ %f",RATES_Lk_Rates(tree)); */

		  tree->rates->bl_from_rt = 1;

/* 		  PhyML_Printf("\n. Burnin...\n"); */
/* 		  tree->mcmc = (tmcmc *)MCMC_Make_MCMC_Struct(tree); */
/* 		  MCMC_Init_MCMC_Struct("burnin",tree->mcmc,tree); */
/* 		  tree->rates->lk_approx = NORMAL; */
/* 		  tree->mcmc->n_tot_run  = 1E+6; */
/* 		  MCMC(tree); */
/* 		  MCMC_Close_MCMC(tree->mcmc); */
/* 		  MCMC_Free_MCMC(tree->mcmc); */

/* 		  PhyML_Printf("\n. Gibbs sampling...\n"); */
/* 		  tree->mcmc = (tmcmc *)MCMC_Make_MCMC_Struct(tree); */
/* 		  MCMC_Init_MCMC_Struct("gibbs",tree->mcmc,tree); */
/* 		  tree->rates->lk_approx = NORMAL; */
/* 		  RATES_Lk_Rates(tree); */
/* 		  MCMC_Print_Param(tree->mcmc,tree); */
		  		  
/* 		  time(&t_beg); */
/* 		  tree->mcmc->n_tot_run = 1E+8; */
/* 		  phydbl u; */
/* 		  int acc,n_trials; */
/* 		  acc = n_trials = 0; */
/* 		  do */
/* 		    { */
/* 		      tree->rates->c_lnL = RATES_Lk_Rates(tree); */
/* 		      tree->c_lnL        = Lk(tree); */

/* 		      MCMC_Nu(tree); */
/* 		      RATES_Posterior_Clock_Rate(tree); */
		      
/* 		      RATES_Posterior_Times(tree); */
/* 		      RATES_Posterior_Rates(&acc,&n_trials,tree); */
/* 		    } */
/* 		  while(tree->mcmc->run < tree->mcmc->n_tot_run); */
/* 		  time(&t_end); */
/* 		  Print_Time_Info(t_beg,t_end); */
/* 		  MCMC_Close_MCMC(tree->mcmc); */
/* 		  MCMC_Free_MCMC(tree->mcmc); */

/* 		  PhyML_Printf("\n. End of Gibbs sampling...\n"); */
/* 		  system("sleep 1s"); */
/* 		  END OF IMPORTANCE SAMPLING STUFF */


		  /************************************/
		  /************************************/
		  /************************************/
		  /************************************/

		  /* THORNE NORMAL */
		  tree->both_sides        = 1;
		  tree->rates->model      = THORNE;
		  tree->rates->bl_from_rt = 1;

		  PhyML_Printf("\n. Burnin...\n");
		  tree->mcmc = (tmcmc *)MCMC_Make_MCMC_Struct(tree);
		  MCMC_Init_MCMC_Struct("burnin",tree->mcmc,tree);
		  tree->rates->lk_approx = NORMAL;
		  tree->mcmc->n_tot_run  = 1E+6;
		  MCMC(tree);
		  MCMC_Close_MCMC(tree->mcmc);
		  MCMC_Free_MCMC(tree->mcmc);

		  tree->mcmc = (tmcmc *)MCMC_Make_MCMC_Struct(tree);
		  MCMC_Init_MCMC_Struct("thorne.normal",tree->mcmc,tree);
		  tree->rates->lk_approx = NORMAL;
		  
		  tree->mcmc->n_tot_run = 1E+8;
		  tree->mcmc->randomize = NO;
		  time(&t_beg);
		  PhyML_Printf("\n. Thorne (approx)...\n");
		  MCMC(tree);
		  PhyML_Printf("\n. End of Thorne (approx)...\n");
		  system("sleep 3s");
		  time(&t_end);
		  Print_Time_Info(t_beg,t_end);
		  MCMC_Close_MCMC(tree->mcmc);
		  MCMC_Free_MCMC(tree->mcmc);



		  /* THORNE EXACT */
		  tree->both_sides        = 0;
		  tree->rates->model      = THORNE;
		  tree->rates->bl_from_rt = 1;

		  PhyML_Printf("\n. Burnin...\n");
		  tree->mcmc = (tmcmc *)MCMC_Make_MCMC_Struct(tree);
		  MCMC_Init_MCMC_Struct("burnin",tree->mcmc,tree);
		  tree->rates->lk_approx = NORMAL;
		  tree->mcmc->n_tot_run  = 1E+6;
		  MCMC(tree);
		  MCMC_Close_MCMC(tree->mcmc);
		  MCMC_Free_MCMC(tree->mcmc);

		  tree->mcmc = (tmcmc *)MCMC_Make_MCMC_Struct(tree);
		  MCMC_Init_MCMC_Struct("thorne.exact",tree->mcmc,tree);
		  tree->rates->lk_approx = EXACT;
		  
		  tree->mcmc->n_tot_run       = 5E+5;
		  tree->mcmc->sample_interval = 1E+3;

		  tree->mcmc->randomize = 0;
		  time(&t_beg);
		  PhyML_Printf("\n. Thorne (exact)...\n");
		  MCMC(tree);
		  PhyML_Printf("\n. End of Thorne (exact)...\n");
		  system("sleep 3s");
		  time(&t_end);
		  Print_Time_Info(t_beg,t_end);
		  MCMC_Close_MCMC(tree->mcmc);
		  MCMC_Free_MCMC(tree->mcmc);
		  Exit("\n");

		  /* END OF COMPOUND POISSON STUFF */
		  /************************************/
		  /************************************/
		  /************************************/
		  /************************************/

		  /* Print the tree estimated using the current random (or BioNJ) starting tree */
		  if(io->mod->s_opt->n_rand_starts > 1)
		    {
		      Br_Len_Involving_Invar(tree);
		      Print_Tree(io->fp_out_trees,tree);
		      fflush(NULL);
		    }

		  /* Record most likely tree in a string of characters */
		  if(tree->c_lnL > best_lnL)
		    {
		      best_lnL = tree->c_lnL;
		      Br_Len_Involving_Invar(tree);
		      most_likely_tree = Write_Tree(tree);
		      most_likely_size = Get_Tree_Size(tree);
		    }

/* 		  JF(tree); */

		  time(&t_end);
		  Print_Fp_Out(io->fp_out_stats,t_beg,t_end,tree,
			       io,num_data_set+1,
			       (tree->mod->s_opt->n_rand_starts > 1)?
			       (num_rand_tree):(num_tree));
		  
		  if(tree->io->print_site_lnl) Print_Site_Lk(tree,io->fp_out_lk);

		  /* Start from BioNJ tree */
		  if((num_rand_tree == io->mod->s_opt->n_rand_starts-1) && (tree->mod->s_opt->random_input_tree))
		    {
		      /* Do one more iteration in the loop, but don't randomize the tree */
		      num_rand_tree--;
		      tree->mod->s_opt->random_input_tree = 0;
		    }
		  
		  Free_Spr_List(tree);
		  Free_One_Spr(tree->best_spr);
		  if(tree->mat) Free_Mat(tree->mat);
		  Free_Triplet(tree->triplet_struct);
		  Free_Tree_Pars(tree);
		  Free_Tree_Lk(tree);
		  Free_Tree(tree);
		}

	      /* Launch bootstrap analysis */
	      if(tree->mod->bootstrap) 
		{
		  PhyML_Printf("\n. Launch bootstrap analysis on the most likely tree...\n");

                  #ifdef MPI
		  MPI_Bcast (most_likely_tree, strlen(most_likely_tree)+1, MPI_CHAR, 0, MPI_COMM_WORLD);
		  PhyML_Printf("\n. The bootstrap analysis will use %d CPUs.",Global_numTask);
		  #endif

		  Bootstrap_From_String(most_likely_tree,cdata,mod,io);
		}
	      else if(tree->io->ratio_test) 
		{
		  /* Launch aLRT */
		  PhyML_Printf("\n. Compute aLRT branch supports on the most likely tree...\n");
		  aLRT_From_String(most_likely_tree,cdata,mod,io);
		}

	      /* Print the most likely tree in the output file */
	      PhyML_Printf("\n. Printing the most likely tree in file '%s'...\n",io->out_tree_file);
	      if(io->n_data_sets == 1) rewind(io->fp_out_tree);

	      PhyML_Fprintf(io->fp_out_tree,"%s\n",most_likely_tree);

	      if(io->n_trees > 1 && io->n_data_sets > 1) break;
	    }
	  Free_Cseq(cdata);
	}
    }


  if(io->mod->s_opt->n_rand_starts > 1) PhyML_Printf("\n\n. Best log likelihood : %f\n",best_lnL);

  Free_Model(mod);

  if(io->fp_in_align)     fclose(io->fp_in_align);
  if(io->fp_in_tree)    fclose(io->fp_in_tree);
  if(io->fp_out_lk)     fclose(io->fp_out_lk);
  if(io->fp_out_tree)   fclose(io->fp_out_tree);
  if(io->fp_out_trees)  fclose(io->fp_out_trees);
  if(io->fp_out_stats)  fclose(io->fp_out_stats);

  Free(most_likely_tree);
  Free_Input(io);

  time(&t_end);
  Print_Time_Info(t_beg,t_end);

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}

#endif

/*********************************************************/

void MC_Least_Square_Node_Times(t_edge *e_root, t_tree *tree)
{

  /* Solve A.x = b, where x are the t_node time estimated
     under the least square criterion.

     A is a n x n matrix, with n being the number of
     nodes in a rooted tree (i.e. 2*n_otu-1).

   */

  phydbl *A, *b, *x;
  int n;
  int i,j;
  t_node *root;
  
  
  PhyML_Printf("\n. Making the tree molecular clock like.");
  
  n = 2*tree->n_otu-1;
  
  A = (phydbl *)mCalloc(n*n,sizeof(phydbl));
  b = (phydbl *)mCalloc(n,  sizeof(phydbl));
  x = (phydbl *)mCalloc(n,  sizeof(phydbl));
  
  
  if(!tree->n_root && e_root) Add_Root(e_root,tree);
  else if(!e_root)            Add_Root(tree->t_edges[0],tree);
  
  root = tree->n_root;

  MC_Least_Square_Node_Times_Pre(root,root->v[0],A,b,n,tree);
  MC_Least_Square_Node_Times_Pre(root,root->v[1],A,b,n,tree);
  
  b[root->num] = tree->e_root->l/2.;
  
  A[root->num * n + root->num]       = 1.0;
  A[root->num * n + root->v[0]->num] = -.5;
  A[root->num * n + root->v[1]->num] = -.5;
    
  if(!Matinv(A, n, n))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }

  For(i,n) x[i] = .0;
  For(i,n) For(j,n) x[i] += A[i*n+j] * b[j];

  For(i,n-1) { tree->rates->nd_t[tree->noeud[i]->num] = x[i]; }
  tree->rates->nd_t[root->num] = x[n-1];
  tree->n_root->l[0] = tree->rates->nd_t[root->num] - tree->rates->nd_t[root->v[0]->num];
  tree->n_root->l[1] = tree->rates->nd_t[root->num] - tree->rates->nd_t[root->v[1]->num];


  /* Rescale the t_node times such that the time at the root
     is -100. This constraint implies that the clock rate
     is fixed to the actual tree length divided by the tree
     length measured in term of differences of t_node times */

  phydbl scale_f,time_tree_length,tree_length;

  scale_f = -100./tree->rates->nd_t[root->num];
  For(i,2*tree->n_otu-1) tree->rates->nd_t[i] *= scale_f;
  For(i,2*tree->n_otu-1) if(tree->rates->nd_t[i] > .0) tree->rates->nd_t[i] = .0;


  time_tree_length = 0.0;
  For(i,2*tree->n_otu-3)
    if(tree->t_edges[i] != tree->e_root)
      time_tree_length +=
	FABS(tree->rates->nd_t[tree->t_edges[i]->left->num] -
	     tree->rates->nd_t[tree->t_edges[i]->rght->num]);
  time_tree_length += FABS(tree->rates->nd_t[root->num] - tree->rates->nd_t[root->v[0]->num]);
  time_tree_length += FABS(tree->rates->nd_t[root->num] - tree->rates->nd_t[root->v[1]->num]);
  
  tree_length = 0.0;
  For(i,2*tree->n_otu-3) tree_length += tree->t_edges[i]->l;

  tree->rates->clock_r = tree_length / time_tree_length;

  Free(A);
  Free(b);
  Free(x);

}

/*********************************************************/

void MC_Least_Square_Node_Times_Pre(t_node *a, t_node *d, phydbl *A, phydbl *b, int n, t_tree *tree)
{
  if(d->tax)
    {
      A[d->num * n + d->num] = 1.;
      
      /* Set the time stamp at tip nodes to 0.0 */
/*       PhyML_Printf("\n. Tip t_node date set to 0"); */
      b[d->num] = 0.0;
      return;
    }
  else
    {
      int i;
      
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  MC_Least_Square_Node_Times_Pre(d,d->v[i],A,b,n,tree);
      
      A[d->num * n + d->num] = 1.;
      b[d->num] = .0;
      For(i,3)
	{
	  A[d->num * n + d->v[i]->num] = -1./3.;
	  if(d->v[i] != a) b[d->num] += d->b[i]->l;
	  else             b[d->num] -= d->b[i]->l;
	}
      b[d->num] /= 3.;
    }
}

/*********************************************************/

/* Adjust t_node times in order to have correct time stamp ranking with
 respect to the tree topology */

void MC_Adjust_Node_Times(t_tree *tree)
{
  MC_Adjust_Node_Times_Pre(tree->n_root->v[0],tree->n_root->v[1],tree);
  MC_Adjust_Node_Times_Pre(tree->n_root->v[1],tree->n_root->v[0],tree);

  if(tree->rates->nd_t[tree->n_root->num] > MIN(tree->rates->nd_t[tree->n_root->v[0]->num],
						tree->rates->nd_t[tree->n_root->v[1]->num]))
    {
      tree->rates->nd_t[tree->n_root->num] = MIN(tree->rates->nd_t[tree->n_root->v[0]->num],
						 tree->rates->nd_t[tree->n_root->v[1]->num]);
    }
}

/*********************************************************/

void MC_Adjust_Node_Times_Pre(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;
      phydbl min_height;

      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  {
	    MC_Adjust_Node_Times_Pre(d,d->v[i],tree);
	  }

      min_height = 0.0;
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      if(tree->rates->nd_t[d->v[i]->num] < min_height)
		{
		  min_height = tree->rates->nd_t[d->v[i]->num];
		}
	    }
	}

      if(tree->rates->nd_t[d->num] > min_height) tree->rates->nd_t[d->num] = min_height;

      if(tree->rates->nd_t[d->num] < -100.) tree->rates->nd_t[d->num] = -100.;

    }
}

/*********************************************************/

  /* Multiply each time stamp at each internal 
     t_node by  'tree->time_stamp_mult'.
   */

void MC_Mult_Time_Stamps(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-2) tree->rates->nd_t[tree->noeud[i]->num] *= FABS(tree->mod->s_opt->tree_size_mult);
  tree->rates->nd_t[tree->n_root->num] *= FABS(tree->mod->s_opt->tree_size_mult);
}

/*********************************************************/

/* Divide each time stamp at each internal 
   t_node by  'tree->time_stamp_mult'.
*/
void MC_Div_Time_Stamps(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-2) tree->rates->nd_t[tree->noeud[i]->num] /= FABS(tree->mod->s_opt->tree_size_mult);
  tree->rates->nd_t[tree->n_root->num] /= FABS(tree->mod->s_opt->tree_size_mult);
}

/*********************************************************/

void MC_Bl_From_T(t_tree *tree)
{
  phydbl mean_rate, branch_rate;

  /* Branch lengths are deduced from time stamps */
  MC_Bl_From_T_Post(tree->n_root,tree->n_root->v[0],NULL,tree);
  MC_Bl_From_T_Post(tree->n_root,tree->n_root->v[1],NULL,tree);

  mean_rate   = tree->rates->clock_r;
  branch_rate = tree->rates->br_r[tree->e_root->num];

  if(tree->rates->use_rates)
    tree->e_root->l = 
      mean_rate * branch_rate * (tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->e_root->left->num]) + 
      mean_rate * branch_rate * (tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->e_root->rght->num]);
  else
    tree->e_root->l = 
      (tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->e_root->left->num]) + 
      (tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->e_root->rght->num]);

  /* Actual formula =>  tree->e_root->l = 
     (tree->n_root->t - tree->e_root->left->t) + 
     (tree->n_root->t - tree->e_root->rght->t); */
  
  tree->n_root_pos = (tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->e_root->left->num])/tree->e_root->l;

}

/*********************************************************/

void MC_Bl_From_T_Post(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{

  if(b)
    {
      phydbl mean_rate, branch_rate;

      mean_rate   = tree->rates->clock_r;
      branch_rate = tree->rates->br_r[b->num];

      if(tree->rates->use_rates)
	{
	  b->l = (tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num]) * mean_rate * branch_rate;
	}
      else
	{
	  b->l = (tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num]);
	}

      if(b->l < 0.0)
	{
	  PhyML_Printf("\n. Correction failed.");
	  PhyML_Printf("\n. d->t = %f a->t = %f",tree->rates->nd_t[d->num],tree->rates->nd_t[a->num]);
	  PhyML_Printf("\n. a->num=%d d->num=%d",a->num,d->num);
	  Warn_And_Exit("\n");
	}
    }

  if(d->tax) return;
  else
    {
      int i;
      For(i,3) if((d->v[i] != a) && (d->b[i] != tree->e_root)) MC_Bl_From_T_Post(d,d->v[i],d->b[i],tree);
    }
}

/*********************************************************/



void MC_Print_Node_Times(t_node *a, t_node *d, t_tree *tree)
{
  t_edge *b;
  int i;
  
  b = NULL;
  For(i,3) if((d->v[i]) && (d->v[i] == a)) {b = d->b[i]; break;}

  PhyML_Printf("\n. (%3d %3d) a->t = %f d->t = %f (#=%f) b->l = %f",
	 a->num,d->num,
	 tree->rates->nd_t[a->num],
	 tree->rates->nd_t[d->num],
	 tree->rates->nd_t[a->num]-tree->rates->nd_t[d->num],
	 (b)?(b->l):(-1.0));
  if(d->tax) return;
  else
    {
      int i;
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  MC_Print_Node_Times(d,d->v[i],tree);
    }
}

/*********************************************************/

void MC_Compute_Rates_And_Times_Least_Square_Adjustments(t_tree *tree)
{
  MC_Compute_Rates_And_Times_Least_Square_Adjustments_Post(tree->n_root,tree->n_root->v[0],NULL,tree);
  MC_Compute_Rates_And_Times_Least_Square_Adjustments_Post(tree->n_root,tree->n_root->v[1],NULL,tree);
}

/*********************************************************/

void MC_Compute_Rates_And_Times_Least_Square_Adjustments_Post(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  int i;

  if(d->tax) return;

  if(b)
    {
      phydbl t0, t1, t2, t3;
      phydbl mu1, mu2, mu3;
      phydbl K;

      t0 = tree->rates->nd_t[a->num];
      t1 = tree->rates->nd_t[d->num];

      mu1 = tree->rates->br_r[b->num];
      mu2 = -1.;
      mu3 = -1.;

      t2 = t3 = -1.;
      For(i,3)
	if(d->v[i] != a)
	  {
	    if(t2 < 0) 
	      {
		t2  = tree->rates->nd_t[d->v[i]->num];
		mu2 = tree->rates->br_r[d->b[i]->num];
	      }
	    if(t3 < 0) 
	      {
		t3  = tree->rates->nd_t[d->v[i]->num];
		mu3 = tree->rates->br_r[d->b[i]->num];
	      }
	  }

      PhyML_Printf("\n. t0=%f, t1=%f, t2=%f, t3=%f, mu1=%f, mu2=%f, mu3=%f",
	     t0,t1,t2,t3,
	     mu1,mu2,mu3);

      if(t2 < 1.E-10 && t3 < 1.E-10)
	{
	  K = t0 / t1;
	}
      else
	{
	  K = 
	    (POW(mu2,2)*t2)/(POW(t1-t2,2)) +
	    (POW(mu1,2)*t0)/(POW(t0-t1,2)) +
	    (POW(mu3,2)*t3)/(POW(t1-t3,2)) +
	    (POW(mu1,2)*t0)/(POW(t0-t1,2)) +
	    mu2 * mu1 * (t0 + t2) / ((t1-t2)*(t0-t1)) +
	    mu3 * mu1 * (t0 + t3) / ((t1-t3)*(t0-t1)) ;

	  K /= 
	    t1 *
	    (
	    (POW(mu2,2))/(POW(t1-t2,2))  +
	    (POW(mu1,2))/(POW(t0-t1,2))  +
	    (POW(mu3,2))/(POW(t1-t3,2))  +
	    (POW(mu1,2))/(POW(t0-t1,2))  +
	    2.*mu2*mu1/((t1-t2)*(t0-t1)) +
	    2.*mu3*mu1/((t1-t3)*(t0-t1)) 
	    );
	}
      PhyML_Printf("\n. K = %f",K); 
    }

  For(i,3) 
    if(d->v[i] != a)
      MC_Compute_Rates_And_Times_Least_Square_Adjustments_Post(d,d->v[i],d->b[i],tree);

}

/*********************************************************/

int MC_Check_MC(t_tree *tree)
{
  int i;
  phydbl dist,eps;

  eps = 1.E-06;

  Dist_To_Root(tree->n_root,tree);
  dist = tree->noeud[0]->dist_to_root;
  
  For(i,tree->n_otu)
    {
      if(FABS(tree->noeud[i]->dist_to_root - dist) > eps)
	{
	  return 0;
	}
    }
  return 1;
}

/*********************************************************/

/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
