/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "mixt.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

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

phydbl MIXT_Lk(t_edge *mixt_b, t_tree *mixt_tree)
{
  t_tree *tree,*cpy_mixt_tree;
  t_edge *b;
  phydbl sum_lnL;
  int site, class;
  phydbl *sum_scale_left_cat,*sum_scale_rght_cat;
  phydbl sum,tmp;
  phydbl logbig,logsmall;
  int exponent;
  phydbl site_lk_cat,site_lk,log_site_lk,inv_site_lk;
  int num_prec_issue,fact_sum_scale;
  phydbl max_sum_scale,min_sum_scale;

  tree          = mixt_tree;
  b             = mixt_b;
  cpy_mixt_tree = mixt_tree;


  /*! Get all the likelihoods, without scaling any of them */
  do
    {
      if(tree->child)
        {
          tree = tree->child;
          b    = b?b->child:NULL;
        }
              
      if(!(tree->mod->ras->invar == YES && mixt_tree->is_mixt_tree == YES)) 
        {
          tree->apply_lk_scaling = NO;
          Lk(b,tree);
        }

      /* printf("\n. lk = %f",tree->c_lnL); */

      tree = tree->next;
      b    = b?b->next:NULL;
    }
  while(tree);
  
  
  /*! Apply scaling factors */

  do /*! Consider each element of the data partition */
    {

      logbig   = LOG((phydbl)BIG);
      logsmall = LOG((phydbl)SMALL);

      if(!mixt_b) mixt_b = mixt_tree->t_nodes[0]->b[0];

      sum_scale_left_cat = (phydbl *)mCalloc(mixt_tree->mod->ras->n_catg,sizeof(phydbl));
      sum_scale_rght_cat = (phydbl *)mCalloc(mixt_tree->mod->ras->n_catg,sizeof(phydbl));


      mixt_tree->c_lnL = .0;

      For(site,mixt_tree->n_pattern)
        {
          tree  = mixt_tree->child;
          b     = mixt_b->child;
          class = 0;

          max_sum_scale =  (phydbl)BIG;
          min_sum_scale = -(phydbl)BIG;

          do
            {
              if(tree->mod->ras->invar == YES) 
                {
                  tree = tree->next;
                  b    = b->next;
                  continue; /*! Skip class with null rate */
                }

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
          while(tree && tree->is_mixt_tree == NO);

          tree = NULL; /*! For debugging purpose */
          
          if(min_sum_scale > max_sum_scale) min_sum_scale = max_sum_scale;
          
          fact_sum_scale = (int)((max_sum_scale + min_sum_scale) / 2);
          



          /*! Populate the mixt_tree->site_lk_cat[class] table after
            scaling */

          tree  = mixt_tree->child;
          b     = mixt_b->child;
          class = 0;

          do
            {
              if(tree->mod->ras->invar == YES) 
                {
                  tree = tree->next;
                  b    = b->next;
                  continue; /*! Skip class with null rate */
                }

              exponent = -(sum_scale_left_cat[class]+sum_scale_rght_cat[class])+fact_sum_scale;
              site_lk_cat = tree->cur_site_lk[site];
              Rate_Correction(exponent,&site_lk_cat,mixt_tree);
              mixt_tree->site_lk_cat[class] = site_lk_cat;

              class++;

              tree = tree->next;
              b    = b->next;
            }
          while(tree && tree->is_mixt_tree == NO);

          
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
                For(class,mixt_tree->mod->ras->n_catg) PhyML_Printf("\n== rr=%f p=%f",mixt_tree->mod->ras->gamma_rr->v[class],mixt_tree->mod->ras->gamma_r_proba->v[class]);
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
  
  mixt_tree = cpy_mixt_tree;

  sum_lnL = .0;
  do
    {
      sum_lnL += mixt_tree->c_lnL;
      mixt_tree = mixt_tree->next;
    }
  while(mixt_tree);

  mixt_tree = cpy_mixt_tree;

  do
    {
      mixt_tree->c_lnL = sum_lnL;
      mixt_tree = mixt_tree->next;
    }
  while(mixt_tree);

  mixt_tree = cpy_mixt_tree;
  
  return mixt_tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Update_P_Lk(t_tree *mixt_tree, t_edge *mixt_b, t_node *mixt_d)
{
  t_tree *tree;
  t_edge *b;
  t_node *d;

  tree = mixt_tree;
  b    = mixt_b;
  d    = mixt_d;

  do
    {
      if(tree->child) 
        {
          tree = tree->child;
          b    = b->child;
          d    = d->child;
        }

      Update_P_Lk(tree,b,d);

      tree = tree->next;
      b    = b->next;
      d    = d->next;

    }
  while(tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Update_PMat_At_Given_Edge(t_edge *mixt_b, t_tree *mixt_tree)
{
  t_tree *tree;
  t_edge *b;

  tree = mixt_tree;
  b    = mixt_b;

  do
    {
      if(tree->child) 
        {
          tree = tree->child;
          b    = b->child;
        }

      Update_PMat_At_Given_Edge(b,tree);

      tree = tree->next;
      b    = b->next;
    }
  while(tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int *MIXT_Get_Number_Of_Classes_In_All_Mixtures(t_tree *mixt_tree)
{
  int *n_catg;
  t_tree *tree;
  int class;
  
  n_catg = NULL;
  tree = mixt_tree;
  class = 0;
  do
    {
      if(!class) n_catg = (int *)mCalloc(1,sizeof(int));
      else       n_catg = (int *)realloc(n_catg,(class+1)*sizeof(int));

      n_catg[class] = tree->mod->ras->n_catg;
      class++;
      tree = tree->next;
    }
  while(tree);

  return(n_catg);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_tree **MIXT_Record_All_Mixtures(t_tree *mixt_tree)
{
  t_tree **tree_list;
  int n_trees;
  t_tree *tree;

  tree_list = NULL;
  n_trees   = 0;
  tree      = mixt_tree;
  do
    {
      if(!tree_list) tree_list = (t_tree **)mCalloc(1,sizeof(t_tree *));
      else           tree_list = (t_tree **)realloc(tree_list,(n_trees+1)*sizeof(t_tree *));
      
      tree_list[n_trees] = tree;
      n_trees++;

      if(tree->child) tree = tree->child;
      else            tree = tree->next;
    }
  while(tree);

  tree_list = (t_tree **)realloc(tree_list,(n_trees+1)*sizeof(t_tree *));
  tree_list[n_trees] = NULL;

  return(tree_list);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Break_All_Mixtures(int c_max, t_tree *mixt_tree)
{
  t_tree *tree;
  int c;

  if(mixt_tree->is_mixt_tree == NO) return;
  
  c = 0;
  tree = mixt_tree;
  do
    {
      if(tree->is_mixt_tree) 
        {
          c = 0;
          tree = tree->child;
        }

      if(c == (c_max-1)           && 
         tree->is_mixt_tree == NO &&
         tree->next != NULL       && 
         tree->next->is_mixt_tree == NO) tree->next = NULL;
      
      tree = tree->next;
      c++;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Reconnect_All_Mixtures(t_tree **tree_list, t_tree *mixt_tree)
{
  t_tree *tree;
  int n_trees;

  if(mixt_tree->is_mixt_tree == NO) return;

  tree = mixt_tree;
  n_trees = 0;
  do
    {
      tree->next = tree_list[n_trees+1];      
      tree = tree->next;
    }
  while(tree);
  

}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int *MIXT_Record_Has_Invariants(t_tree *mixt_tree)
{
  int *has_invariants;
  t_tree *tree;
  int n_trees;

  has_invariants = NULL;
  tree = mixt_tree;
  n_trees = 0;
  do
    {
      if(!n_trees) has_invariants = (int *)mCalloc(1,sizeof(int));
      else         has_invariants = (int *)realloc(has_invariants,(n_trees+1)*sizeof(int));
      has_invariants[n_trees] = (tree->mod->ras->invar == YES)?1:0;
      n_trees++;
      if(tree->child) tree = tree->child;
      else            tree = tree->next;
    }
  while(tree);

  return(has_invariants);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Reset_Has_Invariants(int *has_invariants, t_tree *mixt_tree)
{
  t_tree *tree;
  int n_trees;

  tree = mixt_tree;
  n_trees = 0;
  do
    {
      tree->mod->ras->invar = has_invariants[n_trees];
      n_trees++;
      if(tree->child) tree = tree->child;
      else            tree = tree->next;
    }
  while(tree);

  Free(has_invariants);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Check_Invar_Struct_In_Each_Partition_Elem(t_tree *mixt_tree)
{
  if(mixt_tree->is_mixt_tree == NO) return;
  else
    {
      t_tree *tree;
      int n_inv;
      
      n_inv = 0;
      tree = mixt_tree;
      do
        {
          if(tree->child) 
            {
              tree  = tree->child;
              n_inv = 0;
            }

          if(tree->mod->ras->invar == YES) n_inv++; 

          if(n_inv > 1)
            {
              PhyML_Printf("\n== Found %d classes of the mixture for file '%s' set to",n_inv,tree->parent->io->in_align_file);
              PhyML_Printf("\n== invariable. Only one such class per mixture is allowed.");
              PhyML_Printf("\n== Err in file %s at line %d\n\n",__FILE__,__LINE__);
              Warn_And_Exit("\n");
            }

          if(tree->parent->mod->ras->invar == NO && tree->mod->ras->invar == YES)
            {
              PhyML_Printf("\n== Unexpected settings for 'siterates' in a partition element (file '%s')",tree->parent->io->in_align_file);
              PhyML_Printf("\n== Err in file %s at line %d\n\n",__FILE__,__LINE__);
              Warn_And_Exit("\n");
            }
          
          tree = tree->next;
        }
      while(tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Check_RAS_Struct_In_Each_Partition_Elem(t_tree *mixt_tree)
{
  if(mixt_tree->is_mixt_tree == NO) return;
  else
    {
      t_tree *tree;
      int n_classes;
      
      n_classes = 0;
      tree = mixt_tree;
      do
        {
          if(tree->child) 
            {
              tree  = tree->child;
              n_classes = 0;
            }

          if(tree && tree->mod->ras->invar == NO) n_classes++;

          if((tree->next && tree->next->is_mixt_tree == YES) || (!tree->next)) /*! current tree is the last element of this mixture */
            {
              if(n_classes != tree->parent->mod->ras->n_catg)
                {                  
                  PhyML_Printf("\n== %d tree%s found for file '%s' while the corresponding",
                               n_classes,
                               (n_classes>1)?"s\0":"\0",
                               tree->parent->io->in_align_file);
                  PhyML_Printf("\n== 'siterates' element defined %d class%s.",
                               tree->parent->mod->ras->n_catg,
                               (tree->parent->mod->ras->n_catg>1)?"es\0":"\0");
                  PhyML_Printf("\n== Err in file %s at line %d\n\n",__FILE__,__LINE__);
                  Warn_And_Exit("\n");
                }
            }

          tree = tree->next;
        }
      while(tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



void MIXT_Prune_Subtree(t_node *mixt_a, t_node *mixt_d, t_edge **mixt_target, t_edge **mixt_residual, t_tree *mixt_tree)
{
  t_node *a,*d;
  t_edge **target, **residual;
  t_tree *tree;
  
  tree     = mixt_tree;
  a        = mixt_a;
  d        = mixt_d;
  target   = mixt_target;
  residual = mixt_residual;

  do
    {      
      Prune_Subtree(a,d,target,residual,tree);
      
      if(tree->child) 
        {
          tree        = tree->child;
          a           = a->child;
          d           = d->child;
          (*target)   = (*target)->child;
          (*residual) = (*residual)->child;
        }
      else            
        {
          tree        = tree->next;
          a           = a->next;
          d           = d->next;
          (*target)   = (*target)->next;
          (*residual) = (*residual)->next;
        }
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Graft_Subtree(t_edge *mixt_target, t_node *mixt_link, t_edge *mixt_residual, t_tree *mixt_tree)
{
  t_edge *target,*residual;
  t_node *link;
  t_tree *tree;

  tree     = mixt_tree;
  target   = mixt_target;
  residual = mixt_residual;
  link     = mixt_link;
  
  do
    {
      Graft_Subtree(target,link,residual,tree);
      
      if(tree->child)
        {
          tree     = tree->child;
          target   = target->child;
          residual = residual->child;
          link     = link->child;
        }
      else
        {
          tree     = tree->next;
          target   = target->next;
          residual = residual->next;
          link     = link->next;          
        }
    }

  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Fast_Br_Len(t_edge *mixt_b, t_tree *mixt_tree, int approx)
{
  Exit("\n. MIXT_Fast_Br_Len not implemented yet");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Br_Len_Brent(t_edge *mixt_b, t_tree *mixt_tree, int n_iter_max, int quickdirty)
{
  Exit("\n. MIXT_Br_Len_Brent not implemented yet");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Prepare_Tree_For_Lk(t_tree *mixt_tree)
{
  t_tree *tree;
  
  tree = mixt_tree;
  do
    {
      if(tree->child) tree = tree->child;
      Prepare_Tree_For_Lk(tree);      
      tree = tree->next;
    }
  while(tree && tree->is_mixt_tree == NO);

  if(tree) Prepare_Tree_For_Lk(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Br_Len_Involving_Invar(t_tree *mixt_tree)
{
  int i;
  t_edge *b;

  For(i,2*mixt_tree->n_otu-3) 
    {
      b = mixt_tree->t_edges[i];
      b->l->onoff = ON;
      do
	{
	  if(b->l->onoff == ON)
	    {
	      b->l->v *= (1.-mixt_tree->mod->pinvar->v);
	      b->l->onoff = OFF;
	    }
	  b = b->next;
	}
      while(b);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Br_Len_Not_Involving_Invar(t_tree *mixt_tree)
{
  int i;
  t_edge *b;


  For(i,2*mixt_tree->n_otu-3) 
    {
      b = mixt_tree->t_edges[i];
      b->l->onoff = ON;
      do
	{
	  if(b->l->onoff == ON)
	    {
	      b->l->v /= (1.-mixt_tree->mod->pinvar->v);
	      b->l->onoff = OFF;
	    }
	  b = b->next;
	}
      while(b);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl MIXT_Unscale_Br_Len_Multiplier_Tree(t_tree *mixt_tree)
{
  int i;
  t_edge *b;


  For(i,2*mixt_tree->n_otu-3) 
    {
      b = mixt_tree->t_edges[i];
      b->l->onoff = ON;
      do
	{
	  if(b->l->onoff == ON)
	    {
	      b->l->v /= mixt_tree->mod->br_len_multiplier->v;
	      b->l->onoff = OFF;
	    }
	  b = b->next;
	}
      while(b);
    }
  return(-1);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl MIXT_Rescale_Br_Len_Multiplier_Tree(t_tree *mixt_tree)
{
  int i;
  t_edge *b;


  For(i,2*mixt_tree->n_otu-3) 
    {
      b = mixt_tree->t_edges[i];
      b->l->onoff = ON;
      do
	{
	  if(b->l->onoff == ON)
	    {
	      b->l->v *= mixt_tree->mod->br_len_multiplier->v;
	      b->l->onoff = OFF;
	    }
	  b = b->next;
	}
      while(b);
    }
  return(-1);
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

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

