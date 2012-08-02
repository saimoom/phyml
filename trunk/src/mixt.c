/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "mixt.h"


void MIXT_Connect_Edges_To_Next_Prev_Child_Parent(t_tree *tree)
{
  int i;
  t_edge *b;

  For(i,2*tree->n_otu-3)
    {
      b = tree->t_edges[i];

      if(tree->next)   b->next   = tree->next->t_edges[i];
      if(tree->prev)   b->prev   = tree->prev->t_edges[i];
      if(tree->child)  b->child  = tree->child->t_edges[i];
      if(tree->parent) b->parent = tree->parent->t_edges[i];
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Connect_Nodes_To_Next_Prev_Child_Parent(t_tree *tree)
{
  int i;
  t_node *n;

  For(i,2*tree->n_otu-2)
    {
      n = tree->t_nodes[i];

      if(tree->next)   n->next   = tree->next->t_nodes[i];
      if(tree->prev)   n->prev   = tree->prev->t_nodes[i];
      if(tree->child)  n->child  = tree->child->t_nodes[i];
      if(tree->parent) n->parent = tree->parent->t_nodes[i];
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Connect_Sprs_To_Next_Prev_Child_Parent(t_tree *tree)
{
  int i;

  For(i,2*tree->n_otu-2)
    {
      if(tree->next)   tree->spr_list[i]->next   = tree->next->spr_list[i];
      if(tree->prev)   tree->spr_list[i]->prev   = tree->prev->spr_list[i];
      if(tree->child)  tree->spr_list[i]->child  = tree->child->spr_list[i];
      if(tree->parent) tree->spr_list[i]->parent = tree->parent->spr_list[i];
    }    
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Turn_Branches_OnOff(int onoff, t_tree *tree)
{
  int i;

  For(i,2*tree->n_otu-3) tree->t_edges[i]->l->onoff = onoff;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl *MIXT_Get_Lengths_Of_This_Edge(t_edge *mixt_b)
{
  phydbl *lens;
  t_edge *b;
  int n_lens;

  lens = NULL;
  n_lens = 0;

  b = mixt_b;
  do
    {
      if(b->child) b = b->child;
      
      if(!lens) lens = (phydbl *)mCalloc(1,sizeof(phydbl));
      else      lens = (phydbl *)realloc(lens,(n_lens+1)*sizeof(phydbl));

      lens[n_lens] = b->l->v;

      n_lens++;
      b = b->next;
    }
  while(b);
  
  return(lens);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Set_Lengths_Of_This_Edge(phydbl *lens, t_edge *mixt_b)
{
  t_edge *b;
  int n_lens;

  n_lens = 0;

  b = mixt_b;
  do
    {
      if(b->child) b = b->child; 
      b->l->v = lens[n_lens];
      n_lens++;
      b = b->next;
    }
  while(b);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Post_Order_Lk(t_node *mixt_a, t_node *mixt_d, t_tree *mixt_tree)
{
  t_tree *tree;
  t_node *a,*d;

  tree = mixt_tree;
  a    = mixt_a;
  d    = mixt_d;

  do
    {
      if(tree->child)
        {
          tree = tree->child;
          a    = a->child;
          d    = d->child;
        }

      Post_Order_Lk(a,d,tree);

      tree = tree->next;
      a    = a->next;
      d    = d->next;
    }
  while(tree);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Pre_Order_Lk(t_node *mixt_a, t_node *mixt_d, t_tree *mixt_tree)
{
  t_tree *tree;
  t_node *a,*d;

  tree = mixt_tree;
  a    = mixt_a;
  d    = mixt_d;

  do
    {
      if(tree->child)
        {
          tree = tree->child;
          a    = a->child;
          d    = d->child;
        }

      Pre_Order_Lk(a,d,tree);

      tree = tree->next;
      a    = a->next;
      d    = d->next;
    }
  while(tree);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Lk(t_edge *mixt_b, t_tree *mixt_tree)
{

  t_tree *tree;
  t_edge *b;
  phydbl sum_lnL;

  tree = mixt_tree;
  b    = mixt_b;

  /*! Get all the likelihoods, without scaling any of them */
  do
    {
      if(tree->child)
        {
          tree = tree->child;
          b    = b->child;
        }
      
      if(tree->mod->ras->invar == NO) 
        {
          tree->apply_lk_scaling = NO;
          Lk(b,tree);
        }

      tree = tree->next;
      b    = b->next;
    }
  while(tree);
  
  
  /*! Apply scaling factors */
  do
    {
      int site;
      int class;
      phydbl *sum_scale_left_cat,*sum_scale_rght_cat;
      phydbl sum,tmp;
      phydbl logbig,logsmall;
      int exponent;
      phydbl site_lk_cat,site_lk,log_site_lk;
      int num_prec_issue,fact_sum_scale;

      logbig   = LOG((phydbl)BIG);
      logsmall = LOG((phydbl)SMALL);

      sum_scale_left_cat = (phydbl *)mCalloc(mixt_tree->mod->ras->n_catg,sizeof(phydbl));
      sum_scale_rght_cat = (phydbl *)mCalloc(mixt_tree->mod->ras->n_catg,sizeof(phydbl));

      mixt_tree->c_lnL = .0;

      For(site,mixt_tree->n_patterns)
        {

          tree  = mixt_tree->child;
          b     = mixt_b->child;
          class = 0;

          do
            {
              
              sum_scale_left_cat[class] =
                (b->sum_scale_left)?
                (b->sum_scale_left[site]):
                (0.0);                  

              sum_scale_rght_cat[class] =
                (b->sum_scale_rght)?
                (b->sum_scale_rght[site]):
                (0.0);
               
              sum = sum_scale_left_cat[class] + sum_scale_rght_cat[class];
              
              if(sum < .0)
                {
                  PhyML_Printf("\n== sum = %G",sum);
                  PhyML_Printf("\n== Err in file %s at line %d\n\n",__FILE__,__LINE__);
                  Warn_And_Exit("\n");
                }
              
              tmp = sum + (logbig - LOG(tree->site_lk_cat[0]))/(phydbl)LOG2;
              if(tmp < max_sum_scale) max_sum_scale = tmp; /* min of the maxs */
              
              tmp = sum + (logsmall - LOG(tree->site_lk_cat[0]))/(phydbl)LOG2;
              if(tmp > min_sum_scale) min_sum_scale = tmp; /* max of the mins */
              
              class++;

              tree = tree->next;
              b    = b->next;
            }
          while(tree->is_mixt_tree == NO);

          tree = NULL; /*! For debugging purpose */
          
          if(min_sum_scale > max_sum_scale) min_sum_scale = max_sum_scale;
          
          fact_sum_scale = (int)((max_sum_scale + min_sum_scale) / 2);
          
          For(class,mixt_tree->mod->ras->n_catg)
            {
              exponent = -(sum_scale_left_cat[class]+sum_scale_rght_cat[class])+fact_sum_scale;
              site_lk_cat = mixt_tree->site_lk_cat[class];     
              Rate_Correction(exponent,&site_lk_cat,mixt_tree);
              mixt_tree->site_lk_cat[class] = site_lk_cat;
            }
          
          site_lk = .0;
          For(class,mixt_tree->mod->ras->n_catg) 
            site_lk += 
            mixt_tree->site_lk_cat[class] * 
            mixt_tree->mod->ras->gamma_r_proba->v[class];

          /* Scaling for invariants */
          if(mixt_tree->mod->ras->invar == YES)
            {
              num_prec_issue = NO;
              inv_site_lk = Invariant_Lk(&fact_sum_scale,site,&num_prec_issue,mixt_tree);  
              
              if(num_prec_issue == YES) // inv_site_lk >> site_lk
                {
                  site_lk = inv_site_lk * mixt_tree->mod->pinvar->v;
                }
              else
                {
                  site_lk = site_lk * (1. - mixt_tree->mod->pinvar->v) + inv_site_lk * mixt_tree->mod->pinvar->v;
                }
            }
          
          log_site_lk = LOG(site_lk) - (phydbl)LOG2 * fact_sum_scale;
          
          For(class,mixt_tree->mod->ras->n_catg) 
            mixt_tree->log_site_lk_cat[class][site] = 
            LOG(mixt_tree->site_lk_cat[class]) - 
            (phydbl)LOG2 * fact_sum_scale;

            if(isinf(log_site_lk) || isnan(log_site_lk))
              {
                PhyML_Printf("\n== Site = %d",site);
                PhyML_Printf("\n== Invar = %d",mixt_tree->data->invar[site]);
                PhyML_Printf("\n== Mixt = %d",mixt_tree->is_mixt_tree);
                PhyML_Printf("\n== Lk = %G LOG(Lk) = %f < %G",site_lk,log_site_lk,-BIG);
                For(class,tree->mod->ras->n_catg) PhyML_Printf("\n== rr=%f p=%f",mixt_tree->mod->ras->gamma_rr->v[class],mixt_tree->mod->ras->gamma_r_proba->v[class]);
                PhyML_Printf("\n== Pinv = %G",mixt_tree->mod->pinvar->v);
                PhyML_Printf("\n== Bl mult = %G",mixt_tree->mod->br_len_multiplier->v);
                PhyML_Printf("\n== Err in file %s at line %d",__FILE__,__LINE__);
                Exit("\n");
              }

            mixt_tree->cur_site_lk[site] = EXP(log_site_lk);

            /* Multiply log likelihood by the number of times this site pattern is found in the data */
            mixt_tree->c_lnL_sorted[site] = mixt_tree->data->wght[site]*log_site_lk;
            
            mixt_tree->c_lnL += mixt_tree->data->wght[site]*log_site_lk;
            /*   tree->sum_min_sum_scale += (int)tree->data->wght[site]*min_sum_scale; */
            
        }

      
      Free(sum_scale_left_cat);
      Free(sum_scale_rght_cat);
      
      mixt_tree = mixt_tree->next;
    }
  while(mixt_tree);
  
  while(mixt_tree->prev) { mixt_tree = mixt_tree->prev; }

  mixt_tree->c_lnL = sum_lnL;

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

