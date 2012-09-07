/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "make.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Tree_4_Pars(t_tree *tree, calign *cdata, int n_site)
{
  int i;
  
  tree->site_pars = (int *)mCalloc(tree->n_pattern, sizeof(int));
  tree->step_mat = (int *)mCalloc(tree->mod->ns * tree->mod->ns, sizeof(int));
  For(i,2*tree->n_otu-3) Make_Edge_Pars(tree->a_edges[i],tree);
  Init_Ui_Tips(tree);
  Init_P_Pars_Tips(tree); /* Must be called after Init_Ui_Tips is called */
  Get_Step_Mat(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_All_Edges_Lk(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  For(i,3) if((a->v[i]) && (a->v[i] == d)) Make_Edge_Lk(a->b[i],tree);
  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    Make_All_Edges_Lk(d,d->v[i],tree);
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_New_Edge_Label(t_edge *b)
{
  int i;

  b->labels = (char **)realloc(b->labels,(b->n_labels+BLOCK_LABELS)*sizeof(char *));

  if(!b->labels)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  else
    {
      for(i=b->n_labels;i<b->n_labels+BLOCK_LABELS;i++) b->labels[i] = (char *)mCalloc(T_MAX_LABEL,sizeof(char));
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_edge *Make_Edge_Light(t_node *a, t_node *d, int num)
{
  t_edge *b;

  b = (t_edge *)mCalloc(1,sizeof(t_edge));

  b->l = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(b->l);

  b->l_old = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(b->l_old);

  Init_Edge_Light(b,num);

  if(a && b)
    {
      b->left = a;  b->rght = d;
      if(a->tax) {b->rght = a; b->left = d;} /* root */
      /* a tip is necessary on the right side of the t_edge */

      (b->left == a)?
	(Make_Edge_Dirs(b,a,d,NULL)):
	(Make_Edge_Dirs(b,d,a,NULL));

      b->l->v             = a->l[b->l_r];
      if(a->tax) b->l->v  = a->l[b->r_l];
      b->l_old->v            = b->l->v;
    }
  else
    {
      b->left = NULL;
      b->rght = NULL;
    }

  return b;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Edge_Pars(t_edge *b, t_tree *tree)
{
/*   int site; */

  b->pars_l = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
  b->pars_r = (int *)mCalloc(tree->data->crunch_len,sizeof(int));


  b->ui_l = (unsigned int *)mCalloc(tree->data->crunch_len,sizeof(unsigned int));
  b->ui_r = (unsigned int *)mCalloc(tree->data->crunch_len,sizeof(unsigned int));


  b->p_pars_l = (int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(int ));
  b->p_pars_r = (int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(int ));

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Make_Edge_Lk(t_edge *b, t_tree *tree)
{
  int ns;
  
  ns = -1;

  ns = tree->mod->ns;

  b->l_old->v = b->l->v;

  b->div_post_pred_left = (short int *)mCalloc(ns,sizeof(short int));
  b->div_post_pred_rght = (short int *)mCalloc(ns,sizeof(short int));

  b->Pij_rr = (phydbl *)mCalloc(tree->mod->ras->n_catg*tree->mod->ns*tree->mod->ns,sizeof(phydbl));
  
  b->sum_scale_left_cat = (int *)mCalloc(tree->mod->ras->n_catg,sizeof(int));  
  b->sum_scale_rght_cat = (int *)mCalloc(tree->mod->ras->n_catg,sizeof(int));

  if(!b->left->tax)
    b->sum_scale_left = (int *)mCalloc(tree->data->crunch_len*tree->mod->ras->n_catg,sizeof(int));
  else
    b->sum_scale_left = NULL;
  
  if(!b->rght->tax)
    b->sum_scale_rght = (int *)mCalloc(tree->data->crunch_len*tree->mod->ras->n_catg,sizeof(int));
  else
    b->sum_scale_rght = NULL;
  
  
  if((!b->left->tax) || (tree->mod->s_opt->greedy))
    {
      b->p_lk_left = (phydbl *)mCalloc(tree->data->crunch_len*tree->mod->ras->n_catg*tree->mod->ns,sizeof(phydbl));
      b->p_lk_tip_l = NULL;
    }
  else if(b->left->tax)
    {
      b->p_lk_left   = NULL;      
      b->p_lk_tip_l  = (short int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(short int ));
    }  
  
  if((!b->rght->tax) || (tree->mod->s_opt->greedy))
    {
      b->p_lk_rght = (phydbl *)mCalloc(tree->data->crunch_len*tree->mod->ras->n_catg*tree->mod->ns,sizeof(phydbl));
      b->p_lk_tip_r = NULL;
    }
  else if(b->rght->tax)
    {
      b->p_lk_rght = NULL;      
      b->p_lk_tip_r  = (short int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(short int));
    }

  b->patt_id_left  = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
  b->patt_id_rght  = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
  b->p_lk_loc_left = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
  b->p_lk_loc_rght = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Make_Edge_NNI(t_edge *b)
{
  b->nni    = Make_NNI();
  b->nni->b = b;
  b->nni->left = b->left;
  b->nni->rght = b->rght;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


nni *Make_NNI()
{
  nni *a_nni;
  a_nni = (nni *)mCalloc(1,sizeof(nni ));
  Init_NNI(a_nni);
  return a_nni;
}


t_node *Make_Node_Light(int num)
{
  t_node *n;
  
  n           = (t_node *)mCalloc(1,sizeof(t_node));
  n->v        = (t_node **)mCalloc(3,sizeof(t_node *));
  n->l        = (phydbl *)mCalloc(3,sizeof(phydbl));
  n->b        = (t_edge **)mCalloc(3,sizeof(t_edge *));
  n->score    = (phydbl *)mCalloc(3,sizeof(phydbl));
  n->s_ingrp  = (int *)mCalloc(3,sizeof(int));
  n->s_outgrp = (int *)mCalloc(3,sizeof(int));

  Init_Node_Light(n,num);

  return n;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



void Make_Node_Lk(t_node *n)
{
/*   n->n_ex_nodes = (int *)mCalloc(2,sizeof(int)); */
  return;
}

nexcom **Make_Nexus_Com()
{
  nexcom **com;
  int i;

  com = (nexcom **)mCalloc(N_MAX_NEX_COM,sizeof(nexcom *));
  
  For(i,N_MAX_NEX_COM)
    {
      com[i]       = (nexcom *)mCalloc(1,sizeof(nexcom));
      com[i]->name = (char *)mCalloc(T_MAX_NEX_COM,sizeof(char));
      com[i]->parm = (nexparm **)mCalloc(N_MAX_NEX_PARM,sizeof(nexparm *));
    }

  return com;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


nexparm *Make_Nexus_Parm()
{
  nexparm *parm;

  parm        = (nexparm *)mCalloc(1,sizeof(nexparm));  
  parm->name  = (char *)mCalloc(T_MAX_TOKEN,sizeof(char ));
  parm->value = (char *)mCalloc(T_MAX_TOKEN,sizeof(char ));

  return parm;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

matrix *Make_Mat(int n_otu)
{
  matrix *mat;
  int i;

  mat = (matrix *)mCalloc(1,sizeof(matrix));

  mat->n_otu = n_otu;

  mat->P        = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->Q        = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->dist     = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->on_off   = (int *)mCalloc(n_otu,sizeof(int));
  mat->name     = (char **)mCalloc(n_otu,sizeof(char *));
  mat->tip_node = (t_node **)mCalloc(n_otu,sizeof(t_node *));


  For(i,n_otu)
    {
      mat->P[i]    = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
      mat->Q[i]    = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
      mat->dist[i] = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
      mat->name[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    }

  return mat;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_tree *Make_Tree_From_Scratch(int n_otu, calign *data)
{
  t_tree *tree;

  tree = Make_Tree(n_otu);
  Init_Tree(tree,n_otu);
  Make_All_Tree_Nodes(tree);
  Make_All_Tree_Edges(tree);
  Make_Tree_Path(tree);
  if(data)
    {
      Copy_Tax_Names_To_Tip_Labels(tree,data);
      tree->data = data;
    }
  return tree;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


t_tree *Make_Tree(int n_otu)
{
  t_tree *tree;
  tree = (t_tree *)mCalloc(1,sizeof(t_tree ));
  tree->t_dir = (short int *)mCalloc((2*n_otu-2)*(2*n_otu-2),sizeof(int));
  return tree;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Make_Tree_Path(t_tree *tree)
{
  tree->curr_path = (t_node **)mCalloc(tree->n_otu,sizeof(t_node *));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Make_All_Tree_Nodes(t_tree *tree)
{
  int i;

  tree->a_nodes = (t_node **)mCalloc(2*tree->n_otu-1,sizeof(t_node *));

  For(i,2*tree->n_otu-1)
    {
      tree->a_nodes[i] = (t_node *)Make_Node_Light(i);
      if(i < tree->n_otu) tree->a_nodes[i]->tax = 1;
      else                tree->a_nodes[i]->tax = 0;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Make_All_Tree_Edges(t_tree *tree)
{
  int i;
  tree->a_edges      = (t_edge **)mCalloc(2*tree->n_otu-2,sizeof(t_edge *));
  For(i,2*tree->n_otu-2) tree->a_edges[i] = (t_edge *)Make_Edge_Light(NULL,NULL,i);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

calign *Make_Cseq(int n_otu, int crunch_len, int state_len, int init_len, char **sp_names)
{
  calign *cdata;
  int j;

  cdata           = (calign *)mCalloc(1,sizeof(calign));
  cdata->n_otu    = n_otu;
  cdata->c_seq    = (align **)mCalloc(n_otu,sizeof(align *));
  cdata->b_frq    = (phydbl *)mCalloc(T_MAX_ALPHABET,sizeof(phydbl));
  cdata->wght     = (int *)mCalloc(crunch_len,sizeof(int));
  cdata->ambigu   = (short int *)mCalloc(crunch_len,sizeof(short int));
  cdata->invar    = (short int *)mCalloc(crunch_len,sizeof(short int));
  cdata->sitepatt = (int *)mCalloc(init_len,sizeof(int ));
  cdata->format   = 0;

  cdata->crunch_len = crunch_len;
  cdata->init_len   = init_len;
  cdata->obs_pinvar = .0;

  For(j,n_otu)
    {
      cdata->c_seq[j]            = (align *)mCalloc(1,sizeof(align));
      cdata->c_seq[j]->name      = (char *)mCalloc((int)(strlen(sp_names[j])+1),sizeof(char));
      strcpy(cdata->c_seq[j]->name,sp_names[j]);
      cdata->c_seq[j]->state     = (char *)mCalloc(crunch_len*state_len,sizeof(char));
      cdata->c_seq[j]->d_state   = (int *)mCalloc(crunch_len*state_len,sizeof(int));
      cdata->c_seq[j]->is_ambigu = (short int *)mCalloc(crunch_len,sizeof(short int));
    }

  return cdata;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


t_treelist *Make_Treelist(int list_size)
{
  t_treelist *tlist;

  tlist = (t_treelist *)mCalloc(1,sizeof(t_treelist));
  tlist->list_size = list_size;
  tlist->tree = (t_tree **)mCalloc(list_size,sizeof(t_tree *));

  return tlist;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_opt *Make_Optimiz()
{
  t_opt *s_opt;
  s_opt = (t_opt *)mCalloc(1,sizeof(t_opt));
  return s_opt;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Custom_Model(t_mod *mod)
{
  if(!mod->r_mat)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  mod->r_mat->rr->v            = (phydbl *)mCalloc(mod->ns*(mod->ns-1)/2,sizeof(phydbl));
  mod->r_mat->rr_val->v        = (phydbl *)mCalloc(mod->ns*(mod->ns-1)/2,sizeof(phydbl));
  mod->r_mat->rr_num->v        = (int *)mCalloc(mod->ns*(mod->ns-1)/2,sizeof(int *));
  mod->r_mat->n_rr_per_cat->v  = (int *)mCalloc(mod->ns*(mod->ns-1)/2,sizeof(int));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


t_mod *Make_Model_Basic()
{
  t_mod *mod;

  mod                         = (t_mod *)mCalloc(1,sizeof(t_mod));
  mod->modelname              = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  mod->custom_mod_string      = (char *)mCalloc(T_MAX_OPTION,sizeof(char));
  
  mod->ras = Make_RAS_Basic();
 
  mod->Pij_rr                 = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,mod->Pij_rr);

  mod->kappa                  = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(mod->kappa);

  mod->lambda                 = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(mod->lambda);

  mod->ras->alpha             = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(mod->ras->alpha);

  mod->ras->pinvar            = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(mod->ras->pinvar);

  mod->ras->alpha_old         = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(mod->ras->alpha_old);

  mod->kappa_old              = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(mod->kappa_old);

  mod->lambda_old             = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(mod->lambda_old);

  mod->ras->pinvar_old        = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(mod->ras->pinvar_old);

  mod->br_len_multiplier      = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(mod->br_len_multiplier);

  mod->mr                     = (scalar_dbl *)mCalloc(1,sizeof(scalar_dbl));
  Init_Scalar_Dbl(mod->mr);

  mod->user_b_freq            = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,mod->user_b_freq);
  mod->user_b_freq->v         = (phydbl *)mCalloc(T_MAX_OPTION,sizeof(phydbl));

  return mod;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/*! Call only when the values of mod->ns & ras->n_catg is set to its final value */

void Make_Model_Complete(t_mod *mod)
{
  if(mod->use_m4mod == YES)
    {
      M4_Make_Complete(mod->m4mod->n_h,mod->m4mod->n_o,mod->m4mod);
      mod->ns = mod->m4mod->n_o * mod->m4mod->n_h;
    }
  
  mod->Pij_rr->v = (phydbl *)mCalloc(mod->ras->n_catg*mod->ns*mod->ns,sizeof(phydbl));
  mod->eigen     = (eigen *)Make_Eigen_Struct(mod->ns);
  // If r_mat (e_frq) are not NULL, then they have been created elsewhere and affected.

  if(!mod->r_mat) mod->r_mat = (t_rmat *)Make_Rmat(mod->ns);
  if(!mod->e_frq) mod->e_frq = (t_efrq *)Make_Efrq(mod->ns);

  Make_RAS_Complete(mod->ras);
  
  mod->user_b_freq->len = mod->ns;
  
  if(mod->whichmodel < 0)
    {
      PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  if(mod->whichmodel == CUSTOM)
    {
      Make_Custom_Model(mod);
      Translate_Custom_Mod_String(mod);
    }

  if(mod->whichmodel == GTR)
    {
      Make_Custom_Model(mod);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_ras *Make_RAS_Basic()
{
  t_ras *ras;
  ras = (t_ras *)mCalloc(1,sizeof(t_ras));

  ras->gamma_r_proba          = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,ras->gamma_r_proba);
  ras->gamma_r_proba->v = NULL;

  ras->gamma_r_proba_unscaled = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,ras->gamma_r_proba_unscaled);
  ras->gamma_r_proba_unscaled->v = NULL;

  ras->gamma_rr               = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,ras->gamma_rr);
  ras->gamma_rr->v = NULL;

  ras->gamma_rr_unscaled      = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,ras->gamma_rr_unscaled);
  ras->gamma_rr_unscaled->v = NULL;
  
  ras->next = NULL;

  return(ras);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/*! Call only when the value of ras->n_catg is set to its final value */
void Make_RAS_Complete(t_ras *ras)
{
  if(!ras->gamma_r_proba->v)
    {
      ras->gamma_r_proba->v          = (phydbl *)mCalloc(ras->n_catg,sizeof(phydbl));
      ras->gamma_r_proba_unscaled->v = (phydbl *)mCalloc(ras->n_catg,sizeof(phydbl));
      ras->gamma_rr->v               = (phydbl *)mCalloc(ras->n_catg,sizeof(phydbl));
      ras->gamma_rr_unscaled->v      = (phydbl *)mCalloc(ras->n_catg,sizeof(phydbl));
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_efrq *Make_Efrq(int ns)
{
  t_efrq *e_frq;

  e_frq = (t_efrq *)mCalloc(1,sizeof(t_efrq));

  e_frq->pi               = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  e_frq->pi->v            = (phydbl *)mCalloc(ns,sizeof(phydbl));
  e_frq->pi->len          = ns;

  e_frq->pi_unscaled      = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  e_frq->pi_unscaled->v   = (phydbl *)mCalloc(ns,sizeof(phydbl));
  e_frq->pi_unscaled->len = ns;
  
  e_frq->mixt_weight = 1.; // Should be in a separate init function...
  e_frq->next = NULL; // same
  return e_frq;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_rmat *Make_Rmat(int ns)
{
  t_rmat *r_mat;

  r_mat                         = (t_rmat *)mCalloc(1,sizeof(t_rmat));
  
  r_mat->qmat                   = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,r_mat->qmat);

  r_mat->qmat_buff              = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,r_mat->qmat_buff);

  r_mat->rr                     = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,r_mat->rr);

  r_mat->rr_val                 = (vect_dbl *)mCalloc(1,sizeof(vect_dbl));
  Init_Vect_Dbl(0,r_mat->rr_val);

  r_mat->rr_num                 = (vect_int *)mCalloc(1,sizeof(vect_int));
  Init_Vect_Int(0,r_mat->rr_num);

  r_mat->n_rr_per_cat           = (vect_int *)mCalloc(1,sizeof(vect_int));
  Init_Vect_Int(0,r_mat->n_rr_per_cat);

  r_mat->qmat->v                = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));
  r_mat->qmat_buff->v           = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));

  r_mat->mixt_weight = 1.;
  
  r_mat->next = NULL;

  return(r_mat);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


option *Make_Input()
{
  int i;
  option* io                            = (option *)mCalloc(1,sizeof(option));

  io->in_align_file                     = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->in_tree_file                      = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->in_constraint_tree_file           = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_tree_file                     = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_trees_file                    = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_boot_tree_file                = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_boot_stats_file               = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_stats_file                    = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_lk_file                       = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_ps_file                       = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_trace_file                    = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->nt_or_cd                          = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->run_id_string                     = (char *)mCalloc(T_MAX_OPTION,sizeof(char));
  io->clade_list_file                   = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->aa_rate_mat_file                  = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->alphabet                          = (char **)mCalloc(T_MAX_ALPHABET,sizeof(char *));
  For(i,T_MAX_ALPHABET) io->alphabet[i] = (char *)mCalloc(T_MAX_STATE,sizeof(char ));
  io->treelist                          = (t_treelist *)mCalloc(1,sizeof(t_treelist));
  io->mcmc                              = (t_mcmc *)MCMC_Make_MCMC_Struct();
  io->rates                             = (t_rate *)RATES_Make_Rate_Struct(-1);
  return io;
}

eigen *Make_Eigen_Struct(int ns)
{
  eigen *eig;

  eig              = (eigen *)mCalloc(1,sizeof(eigen));
  eig->size        = ns;
  eig->space       = (phydbl *)mCalloc(2*ns,sizeof(phydbl));
  eig->space_int   = (int *)mCalloc(2*ns,sizeof(int));
  eig->e_val       = (phydbl *)mCalloc(ns,sizeof(phydbl));
  eig->e_val_im    = (phydbl *)mCalloc(ns,sizeof(phydbl));
  eig->r_e_vect    = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));
  eig->r_e_vect_im = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));
  eig->l_e_vect    = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));
  eig->q           = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));

  Init_Eigen_Struct(eig);

  return eig;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

triplet *Make_Triplet_Struct(t_mod *mod)
{
  int i,j,k;
  triplet *triplet_struct;

  triplet_struct                  = (triplet *)mCalloc(1,sizeof(triplet));
  triplet_struct->size            = mod->ns;
  triplet_struct->pi_bc           = (phydbl *)mCalloc(mod->ns,sizeof(phydbl ));
  triplet_struct->pi_cd           = (phydbl *)mCalloc(mod->ns,sizeof(phydbl ));
  triplet_struct->pi_bd           = (phydbl *)mCalloc(mod->ns,sizeof(phydbl ));
  triplet_struct->F_bc            = (phydbl *)mCalloc(mod->ns*mod->ns*mod->ras->n_catg,sizeof(phydbl));
  triplet_struct->F_cd            = (phydbl *)mCalloc(mod->ns*mod->ns*mod->ras->n_catg,sizeof(phydbl));
  triplet_struct->F_bd            = (phydbl *)mCalloc(mod->ns*mod->ns,sizeof(phydbl));
  triplet_struct->core            = (phydbl ****)mCalloc(mod->ras->n_catg,sizeof(phydbl ***));
  triplet_struct->p_one_site      = (phydbl ***)mCalloc(mod->ns,sizeof(phydbl **));
  triplet_struct->sum_p_one_site  = (phydbl ***)mCalloc(mod->ns,sizeof(phydbl **));
  triplet_struct->eigen_struct    = (eigen *)Make_Eigen_Struct(mod->ns);
  triplet_struct->mod             = mod;

  For(k,mod->ras->n_catg)
    {
      triplet_struct->core[k]                = (phydbl ***)mCalloc(mod->ns,sizeof(phydbl **));
      For(i,mod->ns)
	{
	  triplet_struct->core[k][i]         = (phydbl **)mCalloc(mod->ns,sizeof(phydbl *));
	  For(j,mod->ns)
	    triplet_struct->core[k][i][j]    = (phydbl  *)mCalloc(mod->ns,sizeof(phydbl ));
	}
    }

  For(i,mod->ns)
    {
      triplet_struct->p_one_site[i]          = (phydbl **)mCalloc(mod->ns,sizeof(phydbl *));
      For(j,mod->ns)
	triplet_struct->p_one_site[i][j]     = (phydbl  *)mCalloc(mod->ns,sizeof(phydbl ));
    }

  For(i,mod->ns)
    {
      triplet_struct->sum_p_one_site[i]      = (phydbl **)mCalloc(mod->ns,sizeof(phydbl *));
      For(j,mod->ns)
	triplet_struct->sum_p_one_site[i][j] = (phydbl  *)mCalloc(mod->ns,sizeof(phydbl ));
    }
  return triplet_struct;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Short_L(t_tree *tree)
{
  if(!tree->short_l)
    tree->short_l = (phydbl *)mCalloc(tree->n_short_l,sizeof(phydbl));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

xml_attr *XML_Make_Attribute(xml_attr *prev, char *attr_name, char *attr_value)
{
  xml_attr *new_attr;

  new_attr = (xml_attr *)mCalloc(1,sizeof(xml_attr));
  
  new_attr->prev = prev;
  new_attr->next = NULL;
  if(prev != NULL) prev->next = new_attr;
  
  new_attr->name = (char *)mCalloc(strlen(attr_name)+1,sizeof(char));
  strcpy(new_attr->name,attr_name);

  new_attr->value = (char *)mCalloc(strlen(attr_value)+1,sizeof(char));
  strcpy(new_attr->value,attr_value);

  return new_attr;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Make_Node_Id(xml_node *n, char *id)
{
  n->id = (char *)mCalloc(strlen(id)+1,sizeof(char));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Make_Node_Value(xml_node *n, char *val)
{
  n->value = (char *)mCalloc(strlen(val)+1,sizeof(char));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


xml_node *XML_Make_Node(char *name)
{
  
  xml_node *new_node = (xml_node *)mCalloc(1,sizeof(xml_node));

  if(name)
    {
      new_node->name = (char *)mCalloc(strlen(name)+1,sizeof(char));
    }
  
  new_node->ds = (t_ds *)mCalloc(1,sizeof(t_ds));

  return new_node;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Best_Spr(t_tree *tree)
{
  tree->best_spr = Make_One_Spr(tree);
  Init_One_Spr(tree->best_spr);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Make_Spr_List(t_tree *tree)
{
  int i;

  tree->size_spr_list = 2*tree->n_otu-3;
  tree->spr_list = (t_spr **)mCalloc(2*tree->n_otu-2,sizeof(t_spr *));

  For(i,2*tree->n_otu-2)
    {
      tree->spr_list[i] = Make_One_Spr(tree);
      Init_One_Spr(tree->spr_list[i]);
    }

  tree->perform_spr_right_away = NO;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_spr *Make_One_Spr(t_tree *tree)
{
  t_spr *a_spr;
  a_spr         = (t_spr *)mCalloc(1,sizeof(t_spr));
  a_spr->path   = (t_node **)mCalloc(tree->n_otu,sizeof(t_node *));
  return a_spr;
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
