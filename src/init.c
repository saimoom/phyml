/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "init.h"

void Init_Eigen_Struct(eigen *this)
{
}

void Init_Scalar_Dbl(scalar_dbl *p)
{
  p->v = -1.;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Scalar_Int(scalar_int *p)
{
  p->v = -1.;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Vect_Dbl(int len, vect_dbl *p)
{
  p->len = len;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Vect_Int(int len, vect_int *p)
{
  p->len = len;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Tree(t_tree *tree, int n_otu)
{
  tree->n_otu                     = n_otu;
  tree->mat                       = NULL;
  tree->n_root                    = NULL;
  tree->e_root                    = NULL;
  tree->ps_tree                   = NULL;
  tree->short_l                   = NULL;
  tree->mutmap                    = NULL;
  tree->next                      = NULL;
  tree->prev                      = NULL;
  tree->child                     = NULL;
  tree->parent                    = NULL;
  tree->is_mixt_tree              = NO;
  
  tree->tree_num                  = 0;
  tree->depth_curr_path           = 0;
  tree->has_bip                   = NO;
  tree->n_moves                   = 0;
  tree->n_improvements            = 0;
  tree->bl_from_node_stamps       = 0;
  tree->lock_topo                 = 0;
  tree->ps_page_number            = 0;
  tree->init_lnL                  = UNLIKELY;
  tree->best_lnL                  = UNLIKELY;
  tree->old_lnL                   = UNLIKELY;
  tree->c_lnL                     = UNLIKELY;
  tree->sum_min_sum_scale         = .0;
  tree->n_swap                    = 0;
  tree->best_pars                 = 1E+5;
  tree->n_pattern                 = -1;
  tree->n_root_pos                = -1.;
  tree->print_labels              = 1;
  tree->print_boot_val            = 0;
  tree->print_alrt_val            = 0;
  tree->num_curr_branch_available = 0;
  tree->tip_order_score           = .0;
  tree->write_tax_names           = YES;
  tree->update_alias_subpatt      = NO;
  tree->bl_ndigits                = 8;
  tree->n_short_l                 = 100;
  tree->norm_scale                = 0.0;
  tree->br_len_recorded           = NO;
  tree->max_spr_depth             = 0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Init_Edge_Light(t_edge *b, int num)
{
  b->num                  = num;
  b->bip_score            = 0;
  b->dist_btw_edges       = .0;
  b->topo_dist_btw_edges  = 0;
  b->has_zero_br_len      = NO;
  b->n_jumps              = 0;
  b->gamma_prior_mean     = 1.E-0;
  b->gamma_prior_var      = 1.E-1;
  b->does_exist           = YES;
  b->l->v                    = -1.;
  b->bin_cod_num          = -1.;

  b->next                 = NULL;
  b->prev                 = NULL;
  b->child                = NULL;
  b->parent               = NULL;

  b->p_lk_left            = NULL;
  b->p_lk_rght            = NULL;
  b->p_lk_loc_left        = NULL;
  b->p_lk_loc_rght        = NULL;
  b->Pij_rr               = NULL;
  b->labels               = NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Init_Node_Light(t_node *n, int num)
{
  n->num                    = num;
  n->tax                    = -1;
  n->dist_to_root           = .0;
  n->common                 = 1;
  n->ext_node               = NULL;
  n->name                   = NULL;
  n->ori_name               = NULL;
  n->y_rank                 = 0.;
  n->y_rank_ori             = 0.;
  n->y_rank_max             = 0.;
  n->y_rank_min             = 0.;
  n->anc                    = NULL;
  n->rank                   = 0;
  n->match_node             = NULL;
  n->id_rank                = 0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Init_NNI(nni *a_nni)
{
  a_nni->left         = NULL;
  a_nni->rght         = NULL;
  a_nni->b            = NULL;
  a_nni->init_l       = -1.;
  a_nni->init_lk      = .0;
  a_nni->score        = +1.0;
  a_nni->best_l       = -1.;
  a_nni->swap_node_v1 = NULL;
  a_nni->swap_node_v2 = NULL;
  a_nni->swap_node_v3 = NULL;
  a_nni->swap_node_v4 = NULL;
  a_nni->lk0          = UNLIKELY;
  a_nni->lk1          = UNLIKELY;
  a_nni->lk2          = UNLIKELY;
  a_nni->l0           = -1.0;
  a_nni->l1           = -1.0;
  a_nni->l2           = -1.0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Nexus_Format(nexcom **com)
{

  /*****************************/

  strcpy(com[0]->name,"dimensions");
  com[0]->nparm = 2;
  com[0]->nxt_token_t = NEXUS_PARM;
  com[0]->cur_token_t = NEXUS_COM;

  com[0]->parm[0] = Make_Nexus_Parm();
  strcpy(com[0]->parm[0]->name,"ntax");
  com[0]->parm[0]->fp = Read_Nexus_Dimensions;
  com[0]->parm[0]->com = com[0];
  com[0]->parm[0]->nxt_token_t = NEXUS_EQUAL;
  com[0]->parm[0]->cur_token_t = NEXUS_PARM;

  com[0]->parm[1] = Make_Nexus_Parm();
  strcpy(com[0]->parm[1]->name,"nchar");
  com[0]->parm[1]->fp = Read_Nexus_Dimensions;
  com[0]->parm[1]->com = com[0];
  com[0]->parm[1]->nxt_token_t = NEXUS_EQUAL;
  com[0]->parm[1]->cur_token_t = NEXUS_PARM;

  /*****************************/

  strcpy(com[1]->name,"format");
  com[1]->nparm = 11;
  com[1]->nxt_token_t = NEXUS_PARM;
  com[1]->cur_token_t = NEXUS_COM;

  com[1]->parm[0] = Make_Nexus_Parm();
  strcpy(com[1]->parm[0]->name,"datatype");
  com[1]->parm[0]->fp = Read_Nexus_Format;
  com[1]->parm[0]->com = com[1];
  com[1]->parm[0]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[0]->cur_token_t = NEXUS_PARM;

  com[1]->parm[1] = Make_Nexus_Parm();
  strcpy(com[1]->parm[1]->name,"respectcase");
  com[1]->parm[1]->fp = Read_Nexus_Format;
  com[1]->parm[1]->com = com[1];
  com[1]->parm[1]->nxt_token_t = NEXUS_PARM;
  com[1]->parm[1]->cur_token_t = NEXUS_VALUE;

  com[1]->parm[2] = Make_Nexus_Parm();
  strcpy(com[1]->parm[2]->name,"missing");
  com[1]->parm[2]->fp = Read_Nexus_Format;
  com[1]->parm[2]->com = com[1];
  com[1]->parm[2]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[2]->cur_token_t = NEXUS_PARM;

  com[1]->parm[3] = Make_Nexus_Parm();
  strcpy(com[1]->parm[3]->name,"gap");
  com[1]->parm[3]->fp = Read_Nexus_Format;
  com[1]->parm[3]->com = com[1];
  com[1]->parm[3]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[3]->cur_token_t = NEXUS_PARM;

  com[1]->parm[4] = Make_Nexus_Parm();
  strcpy(com[1]->parm[4]->name,"symbols");
  com[1]->parm[4]->fp = Read_Nexus_Format;
  com[1]->parm[4]->com = com[1];
  com[1]->parm[4]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[4]->cur_token_t = NEXUS_PARM;

  com[1]->parm[5] = Make_Nexus_Parm();
  strcpy(com[1]->parm[5]->name,"equate");
  com[1]->parm[5]->fp = Read_Nexus_Format;
  com[1]->parm[5]->com = com[1];
  com[1]->parm[5]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[5]->cur_token_t = NEXUS_PARM;

  com[1]->parm[6] = Make_Nexus_Parm();
  strcpy(com[1]->parm[6]->name,"matchchar");
  com[1]->parm[6]->fp = Read_Nexus_Format;
  com[1]->parm[6]->com = com[1];
  com[1]->parm[6]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[6]->cur_token_t = NEXUS_PARM;

  com[1]->parm[7] = Make_Nexus_Parm();
  strcpy(com[1]->parm[7]->name,"transpose");
  com[1]->parm[7]->fp = Read_Nexus_Format;
  com[1]->parm[7]->com = com[1];
  com[1]->parm[7]->nxt_token_t = NEXUS_PARM;
  com[1]->parm[7]->cur_token_t = NEXUS_VALUE;

  com[1]->parm[8] = Make_Nexus_Parm();
  strcpy(com[1]->parm[8]->name,"interleave");
  com[1]->parm[8]->fp = Read_Nexus_Format;
  com[1]->parm[8]->com = com[1];
  com[1]->parm[8]->nxt_token_t = NEXUS_PARM;
  com[1]->parm[8]->cur_token_t = NEXUS_VALUE;

  com[1]->parm[9] = Make_Nexus_Parm();
  strcpy(com[1]->parm[9]->name,"items");
  com[1]->parm[9]->fp = Read_Nexus_Format;
  com[1]->parm[9]->com = com[1];
  com[1]->parm[9]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[9]->cur_token_t = NEXUS_PARM;

  com[1]->parm[10] = Make_Nexus_Parm();
  strcpy(com[1]->parm[10]->name,"statesformat");
  com[1]->parm[10]->fp = Read_Nexus_Format;
  com[1]->parm[10]->com = com[1];
  com[1]->parm[10]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[10]->cur_token_t = NEXUS_PARM;

  /*****************************/

  strcpy(com[2]->name,"eliminate");
  com[2]->nparm = 0;
  com[2]->nxt_token_t = NEXUS_VALUE;
  com[2]->cur_token_t = NEXUS_COM;

  /*****************************/

  strcpy(com[3]->name,"taxlabels");
  com[3]->nparm = 0;
  com[3]->nxt_token_t = -1;
  com[3]->cur_token_t = -1;
 
 /*****************************/

  strcpy(com[4]->name,"charstatelabels");
  com[4]->nparm = 0;
  com[4]->nxt_token_t = -1;
  com[4]->cur_token_t = -1;

  /*****************************/

  strcpy(com[5]->name,"charlabels");
  com[5]->nparm = 0;
  com[5]->nxt_token_t = -1;
  com[5]->cur_token_t = -1;

  /*****************************/

  strcpy(com[6]->name,"statelabels");
  com[6]->nparm = 0;
  com[6]->nxt_token_t = -1;
  com[6]->cur_token_t = -1;

  /*****************************/

  strcpy(com[7]->name,"matrix");
  com[7]->nparm = 1;
  com[7]->nxt_token_t = NEXUS_COM;
  com[7]->cur_token_t = NEXUS_VALUE; /* This will allow us to skip directly 
					to the matrix reading function */

  com[7]->parm[0] = Make_Nexus_Parm();
  strcpy(com[7]->parm[0]->name,"matrix");
  com[7]->parm[0]->fp = Read_Nexus_Matrix;
  com[7]->parm[0]->com = com[7];
  com[7]->parm[0]->nxt_token_t = NEXUS_COM;
  com[7]->parm[0]->cur_token_t = -1; 

  /*****************************/

  strcpy(com[8]->name,"begin");
  com[8]->nparm = 3;

  com[8]->nxt_token_t = NEXUS_PARM;
  com[8]->cur_token_t = NEXUS_COM;

  com[8]->parm[0] = Make_Nexus_Parm();
  strcpy(com[8]->parm[0]->name,"data");
  com[8]->parm[0]->fp = Read_Nexus_Begin;
  com[8]->parm[0]->com = com[8];
  com[8]->parm[0]->nxt_token_t = NEXUS_COM;
  com[8]->parm[0]->cur_token_t = NEXUS_PARM;


  com[8]->parm[1] = Make_Nexus_Parm();
  strcpy(com[8]->parm[1]->name,"trees");
  com[8]->parm[1]->fp = Read_Nexus_Begin;
  com[8]->parm[1]->com = com[8];
  com[8]->parm[1]->nxt_token_t = NEXUS_COM;
  com[8]->parm[1]->cur_token_t = NEXUS_PARM;


  com[8]->parm[2] = Make_Nexus_Parm();
  strcpy(com[8]->parm[2]->name,"taxa");
  com[8]->parm[2]->fp = Read_Nexus_Taxa;
  com[8]->parm[2]->com = com[8];
  com[8]->parm[2]->nxt_token_t = NEXUS_COM;
  com[8]->parm[2]->cur_token_t = NEXUS_VALUE; 


  /*****************************/

  strcpy(com[9]->name,"end");
  com[9]->nparm = 0;
  com[9]->nxt_token_t = -1;
  com[9]->cur_token_t = -1;
 
  /*****************************/

  strcpy(com[10]->name,"translate");
  com[10]->nparm = 1;
  com[10]->nxt_token_t = NEXUS_COM;
  com[10]->cur_token_t = NEXUS_VALUE;

  com[10]->parm[0] = Make_Nexus_Parm();
  strcpy(com[10]->parm[0]->name,"translate");
  com[10]->parm[0]->fp = Read_Nexus_Translate;
  com[10]->parm[0]->com = com[10];
  com[10]->parm[0]->nxt_token_t = NEXUS_COM;
  com[10]->parm[0]->cur_token_t = -1; 

  /*****************************/

  strcpy(com[11]->name,"tree");
  com[11]->nparm = 1;
  com[11]->nxt_token_t = NEXUS_COM;
  com[11]->cur_token_t = NEXUS_VALUE;

  com[11]->parm[0] = Make_Nexus_Parm();
  strcpy(com[11]->parm[0]->name,"tree");
  com[11]->parm[0]->fp = Read_Nexus_Tree;
  com[11]->parm[0]->com = com[11];
  com[11]->parm[0]->nxt_token_t = -1;
  com[11]->parm[0]->cur_token_t = -1; 


  /*****************************/

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Init_Mat(matrix *mat, calign *data)
{
  int i;

  mat->n_otu = data->n_otu;
  mat->r = mat->n_otu;
  mat->curr_int = mat->n_otu;
  mat->method = 1;

  For(i,data->n_otu)
    {
      strcpy(mat->name[i],data->c_seq[i]->name);
      mat->on_off[i] = 1;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Set_Defaults_Input(option* io)
{
  io->fp_in_align                = NULL;
  io->fp_in_tree                 = NULL;
  io->fp_in_constraint_tree      = NULL;
  io->fp_out_tree                = NULL;
  io->fp_out_trees               = NULL;
  io->fp_out_boot_tree           = NULL;
  io->fp_out_boot_stats          = NULL;
  io->fp_out_stats               = NULL;
  io->long_tax_names             = NULL;
  io->short_tax_names            = NULL;
  io->lon                        = NULL;
  io->lat                        = NULL;
  io->z_scores                   = NULL;
  io->cstr_tree                  = NULL;
  io->next                       = NULL;
  io->prev                       = NULL;

  io->tree                       = NULL;
  io->mod                        = NULL;
  strcpy(io->nt_or_cd,"nucleotides");
  io->n_data_sets                = 1;
  io->interleaved                = 1;
  io->in_tree                    = 0;
  io->out_tree_file_open_mode    = 1;
  io->out_stats_file_open_mode   = 1;
  io->init_len                   = -1;
  io->n_otu                      = -1;
  io->n_data_set_asked           = -1;
  io->print_boot_trees           = 1;
  io->n_part                     = 1;
  io->ratio_test		 = 4;
  io->multigene                  = 0;
  io->config_multigene           = 0;
  io->curr_interface             = 0;
  io->r_seed                     = -1;
  io->collapse_boot              = 0;
  io->random_boot_seq_order      = 1;
  io->print_trace                = 0;
  io->print_site_lnl             = 0;
  io->m4_model                   = NO;
  io->rm_ambigu                  = 0;
  io->append_run_ID              = 0;
  io->quiet                      = 0;
  io->datatype                   = NT;
  io->colalias                   = YES;
  io->data_file_format           = PHYLIP;
  io->tree_file_format           = PHYLIP;
  io->boot_prog_every            = 20;
  io->mem_question               = YES;
  io->do_alias_subpatt           = NO;
  io->lk_approx                  = EXACT;
  io->codpos                     = -1;
  io->mutmap                     = NO;
  io->state_len                  = 1;

  MCMC_Init_MCMC_Struct(NULL,io,io->mcmc);
  RATES_Init_Rate_Struct(io->rates,NULL,-1);
  io->rates->model               = GUINDON;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Set_Defaults_Model(t_mod *mod)
{
  strcpy(mod->modelname,"HKY85");
  strcpy(mod->custom_mod_string,"000000");
  mod->next                    = NULL;
  mod->prev                    = NULL;
  mod->child                   = NULL;
  mod->parent                  = NULL;
  mod->r_mat                   = NULL;
  mod->e_frq                   = NULL;
  mod->whichmodel              = HKY85;
  mod->ras->n_catg                  = 4;
  mod->mod_num                 = 0;

  mod->kappa->v                = 4.0;
  mod->ras->alpha->v                = 1.0;
  mod->lambda->v               = 1.0;
  mod->pinvar->v               = 0.0;

  mod->kappa_old->v            = 4.0;
  mod->ras->alpha_old->v            = 1.0;
  mod->lambda_old->v           = 1.0;
  mod->pinvar_old->v           = 0.0;

  mod->bootstrap               = 0;
  mod->ras->invar                   = NO;
  mod->ns                      = 4;
  mod->use_m4mod               = NO;
  mod->ras->gamma_median            = 0;
  mod->m4mod                   = NULL;
  
  /* mod->r_mat->n_diff_rr        = 0; */
  /* mod->r_mat->rr->v            = NULL; */
  /* mod->r_mat->rr_val->v        = NULL; */
  /* mod->r_mat->n_rr_per_cat->v  = NULL; */
  mod->io                      = NULL;
  mod->log_l                   = NO;
  mod->ras->free_mixt_rates         = NO;
  mod->gamma_mgf_bl            = NO;
  mod->br_len_multiplier->v    = 1.0;
  
#ifndef PHYTIME
  mod->l_min = 1.E-8;
  mod->l_max = 100.0;
#else
  mod->l_min = 1.E-8;
  mod->l_max = 2.0;
#endif
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Set_Defaults_Optimiz(optimiz *s_opt)
{
  s_opt->print                = YES;
  s_opt->last_opt             = YES;
  s_opt->opt_subst_param      = YES;
  s_opt->opt_alpha            = YES;
  s_opt->opt_kappa            = YES;
  s_opt->opt_bl               = YES;
  s_opt->opt_lambda           = NO;
  s_opt->opt_pinvar           = NO;
  s_opt->opt_cov_delta        = NO;
  s_opt->opt_cov_alpha        = NO;
  s_opt->opt_cov_free_rates   = NO;
  s_opt->opt_rr               = NO;
  s_opt->init_lk              = UNLIKELY;
  s_opt->n_it_max             = 1000;
  s_opt->opt_topo             = YES;
  s_opt->topo_search          = NNI_MOVE;
  s_opt->random_input_tree    = 0;
  s_opt->n_rand_starts        = 5;
  s_opt->brent_it_max         = 500;
  s_opt->steph_spr            = YES;
  s_opt->user_state_freq      = NO;
  s_opt->min_diff_lk_local    = 1.E-04;
  s_opt->min_diff_lk_global   = 1.E-03;
  s_opt->min_diff_lk_move     = 1.E-02;
  s_opt->p_moves_to_examine   = 0.15;
  s_opt->fast_nni             = NO;
  s_opt->greedy               = NO;
  s_opt->general_pars         = NO;
  s_opt->tree_size_mult       = 1;
  s_opt->opt_five_branch      = YES;

  s_opt->pars_thresh          = 5;

  s_opt->hybrid_thresh        = NO;
  s_opt->quickdirty           = NO;
  s_opt->spr_pars             = YES;
  s_opt->spr_lnL              = NO;
  s_opt->min_depth_path       = 0;
  s_opt->max_depth_path       = 20;
  s_opt->deepest_path         = 20;
  s_opt->max_delta_lnL_spr    = 50.;
  s_opt->br_len_in_spr        = 10;
  s_opt->opt_free_mixt_rates  = YES;
  s_opt->constrained_br_len   = NO;
  s_opt->opt_gamma_br_len     = NO;

  s_opt->wim_n_rgrft          = -1;
  s_opt->wim_n_globl          = -1;
  s_opt->wim_max_dist         = -1;
  s_opt->wim_n_optim          = -1;
  s_opt->wim_n_best           = -1;
  s_opt->wim_inside_opt       =  0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Init_Node(xml_node *parent, xml_node *new_node, char *name)
{
  if(name) strcpy(new_node->name,name);

  new_node->parent = parent ? parent : NULL;
  new_node->next   = NULL;
  new_node->prev   = NULL;
  new_node->child  = NULL;
  new_node->ds     = NULL;

  if(parent)
    { 
      if(!parent->child)
	{
	  parent->child = new_node;
	}
      else
	{
	  xml_node *node = parent->child;
	  while(node->next) node = node->next;
	  node->next = new_node;
	  new_node->prev = node;
	}
    }

  new_node->attr = NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
