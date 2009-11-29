/*

PhyML:  a program that  computes maximum likelihood phyLOGenies from
DNA or AA homoLOGous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/


/* Routines for molecular clock trees and molecular dating */

#ifdef TIPORDER

#include "tiporder.h"

/*********************************************************/

int TIPORDER_main(int argc, char **argv)
{
  t_tree **list_tree,*ref_tree;
  FILE *fp_ref_tree,*fp_list_tree,*ps_tree;
  int i,j,k;
  int n_trees;
  option *ref_io,*list_io;

  ref_io  = (option *)Make_Input();
  list_io = (option *)Make_Input();

  fp_ref_tree  = (FILE *)fopen(argv[1],"r");
  fp_list_tree = (FILE *)fopen(argv[2],"r");
  ps_tree      = (FILE *)fopen(argv[3],"w");

  ref_io->fp_in_tree = fp_ref_tree;
  Read_Tree_File(ref_io);
  fclose(ref_io->fp_in_tree);
  ref_tree = ref_io->tree;
  Translate_Tax_Names(ref_io->tax_table,ref_tree);

/*   Add_Root(ref_tree->t_edges[0], ref_tree); */
/*   ref_tree->rates = RATES_Make_Rate_Struct(ref_tree->n_otu); */
/*   RATES_Init_Rate_Struct(ref_tree->rates, ref_tree->n_otu); */
/*   MC_Least_Square_Node_Times(ref_tree->t_edges[0],ref_tree); */
/*   MC_Adjust_Node_Times(ref_tree); */
/*   RATES_Update_Cur_Bl(ref_tree); */

  ref_tree->ps_tree = DR_Make_Tdraw_Struct(ref_tree);
  DR_Get_Tree_Coord(ref_tree);

  list_io->fp_in_tree = fp_list_tree;
  Read_Tree_File(list_io);
  fclose(list_io->fp_in_tree);
  list_tree = list_io->treelist->tree;
  n_trees = list_io->treelist->list_size;
  For(i,n_trees) Translate_Tax_Names(list_io->tax_table,list_tree[i]);
  

/*   list_tree = (t_tree **)mCalloc(1,sizeof(t_tree *)); */
/*   n_trees = 0; */
/*   do */
/*     { */
/*       list_tree = (t_tree **)realloc(list_tree,(n_trees+1)*sizeof(t_tree *)); */
/*       list_tree[n_trees] = Read_Tree_File_Phylip(fp_list_tree);       */
/*       PhyML_Printf("\n. Read %d trees",n_trees); */
/*       n_trees++; */
/*     } */
/*   while(list_tree[n_trees-1]); */
/*   n_trees--; */


  For(i,n_trees)
    {
/*       Add_Root(list_tree[i]->t_edges[0],list_tree[i]); */
/*       list_tree[i]->rates = RATES_Make_Rate_Struct(list_tree[i]->n_otu); */
/*       RATES_Init_Rate_Struct(list_tree[i]->rates, list_tree[i]->n_otu); */
/*       MC_Least_Square_Node_Times(list_tree[i]->t_edges[0],list_tree[i]); */
/*       MC_Adjust_Node_Times(list_tree[i]); */
/*       RATES_Update_Cur_Bl(list_tree[i]); */
/*       Dist_To_Root(list_tree[i]->n_root,list_tree[i]); */
    }

  For(i,n_trees)
    {
      For(j,ref_tree->n_otu) 
	{
	  For(k,ref_tree->n_otu) 
	    {
	      if(!strcmp(ref_tree->noeud[j]->name,list_tree[i]->noeud[k]->name))
		{
		  list_tree[i]->noeud[k]->ext_node = ref_tree->noeud[j];
		  break;
		}
	    }
	  if(k == ref_tree->n_otu)
	    {
	      PhyML_Printf("\n. Could not find matching tips for \"%s\" (tree %d)",ref_tree->noeud[j]->name,i);
	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	}
    }

  PhyML_Printf("\n. Rendering trees",n_trees); fflush(NULL);


  DR_Print_Postscript_Header(1,ps_tree);
  For(i,n_trees)
    {
      list_tree[i]->ps_tree = DR_Make_Tdraw_Struct(list_tree[i]);
      DR_Init_Tdraw_Struct(list_tree[i]->ps_tree);
      DR_Get_Tree_Box_Width(list_tree[i]->ps_tree,list_tree[i]);
      if(!list_tree[i]->n_root) Add_Root(list_tree[i]->t_edges[0],list_tree[i]);
      Dist_To_Root(list_tree[i]->n_root,list_tree[i]);
      list_tree[i]->ps_tree->max_dist_to_root = DR_Get_Max_Dist_To_Root(list_tree[i]);
      For(j,ref_tree->n_otu) list_tree[i]->ps_tree->ycoord[j] = ref_tree->ps_tree->ycoord[list_tree[i]->noeud[j]->ext_node->num];
      For(j,ref_tree->n_otu) list_tree[i]->ps_tree->xcoord[j] = ref_tree->ps_tree->xcoord[list_tree[i]->noeud[j]->ext_node->num];
      DR_Get_X_Coord(YES,list_tree[i]->ps_tree,list_tree[i]);
      DR_Get_Y_Coord(YES,list_tree[i]->ps_tree,list_tree[i]);
      if(!i) DR_Print_Tree_Postscript(1,YES,ps_tree,list_tree[i]);
      else   DR_Print_Tree_Postscript(1, NO,ps_tree,list_tree[i]);
    }
  DR_Print_Postscript_EOF(ps_tree);
    

  PhyML_Printf("\n. Getting ancestors"); fflush(NULL);
  For(i,n_trees)   
    {
      Update_Ancestors(list_tree[i]->n_root,list_tree[i]->n_root->v[0],list_tree[i]);
      Update_Ancestors(list_tree[i]->n_root,list_tree[i]->n_root->v[1],list_tree[i]);
      list_tree[i]->n_root->anc = NULL;
    }

  PhyML_Printf("\n. Getting bipartitions"); fflush(NULL);
  For(i,n_trees) Alloc_Bip(list_tree[i]);
  For(i,n_trees) 
    {
      if(!(i%10)) printf("\n. Tree %d",i);
      Get_Bip(list_tree[i]->noeud[0],
	      list_tree[i]->noeud[0]->v[0],
	      list_tree[i]);
    }


  PhyML_Printf("\n. Getting ranks"); fflush(NULL);
  Get_Tips_Y_Rank(ref_tree);
  For(i,n_trees) 
    {
/*       PhyML_Printf("\n. Tree %d\n",i+1); */
      For(j,ref_tree->n_otu) list_tree[i]->noeud[j]->y_rank = list_tree[i]->noeud[j]->ext_node->y_rank;
      Get_All_Y_Rank(list_tree[i]);
/*       Print_Tip_Ordered(list_tree[i]); */
    }

/*   Untangle_Tree_List(n_trees,list_tree,ref_tree); */

  PhyML_Printf("\n. Minimizing"); fflush(NULL);
  Minimize_Tip_Order_Score(n_trees,list_tree,ref_tree);

  For(i,n_trees) Get_All_Y_Rank(list_tree[i]);
/*   For(i,n_trees) Print_Tip_Ordered(list_tree[i]); */

    
  fclose(ps_tree);
  ps_tree = (FILE *)fopen(argv[4],"w");
  DR_Get_Tree_Coord(ref_tree);
  DR_Print_Postscript_Header(1,ps_tree);
  For(i,n_trees)
    {
      DR_Init_Tdraw_Struct(list_tree[i]->ps_tree);
      DR_Get_Tree_Box_Width(list_tree[i]->ps_tree,list_tree[i]);
      if(!list_tree[i]->n_root) Add_Root(list_tree[i]->t_edges[0],list_tree[i]);
      Dist_To_Root(list_tree[i]->n_root,list_tree[i]);
      list_tree[i]->ps_tree->max_dist_to_root = DR_Get_Max_Dist_To_Root(list_tree[i]);
      For(j,ref_tree->n_otu) list_tree[i]->ps_tree->ycoord[j] = ref_tree->ps_tree->ycoord[list_tree[i]->noeud[j]->ext_node->num];
      For(j,ref_tree->n_otu) list_tree[i]->ps_tree->xcoord[j] = ref_tree->ps_tree->xcoord[list_tree[i]->noeud[j]->ext_node->num];
      DR_Get_X_Coord(YES,list_tree[i]->ps_tree,list_tree[i]);
      DR_Get_Y_Coord(YES,list_tree[i]->ps_tree,list_tree[i]);
      if(!i) DR_Print_Tree_Postscript(1,YES,ps_tree,list_tree[i]);
      else   DR_Print_Tree_Postscript(1, NO,ps_tree,list_tree[i]);
    }
  DR_Print_Postscript_EOF(ps_tree);

  fclose(fp_ref_tree);
  fclose(fp_list_tree);
  fclose(ps_tree);
}

/*********************************************************/

void Get_Tips_Y_Rank(t_tree *tree)
{
  phydbl curr_rank;

  curr_rank = .0;
  Get_Tips_Y_Rank_Pre(tree->n_root,tree->n_root->v[0],&curr_rank,tree);
  Get_Tips_Y_Rank_Pre(tree->n_root,tree->n_root->v[1],&curr_rank,tree);
  
  if(curr_rank != tree->n_otu)
    {
      PhyML_Printf("\n. tree->n_otu = %d curr_rank = %d",tree->n_otu,curr_rank);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
}

/*********************************************************/

void Get_Tips_Y_Rank_Pre(t_node *a, t_node *d, phydbl *curr_rank, t_tree *tree)
{
  if(d->tax) 
    {
      d->y_rank = *curr_rank;
      *curr_rank += 1.;
      return;
    }
  else
    {
      int i;
      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      Get_Tips_Y_Rank_Pre(d,d->v[i],curr_rank,tree);
	    }
	}
    }
}

/*********************************************************/

void Get_All_Y_Rank(t_tree *tree)
{
  tree->sum_y_dist_sq = .0;
  tree->sum_y_dist    = .0;
  Get_All_Y_Rank_Pre(tree->n_root,tree->n_root->v[0],tree);
  Get_All_Y_Rank_Pre(tree->n_root,tree->n_root->v[1],tree);
  tree->n_root->y_rank = (tree->n_root->v[0]->y_rank+tree->n_root->v[1]->y_rank)/2.;
  tree->sum_y_dist_sq += POW(tree->n_root->v[0]->y_rank-tree->n_root->v[1]->y_rank,2);
  tree->sum_y_dist    += FABS(tree->n_root->v[0]->y_rank-tree->n_root->v[1]->y_rank);
/*   tree->tip_order_score = tree->sum_y_dist_sq/(tree->n_otu-1)  - POW(tree->sum_y_dist / (tree->n_otu-1),2); */
/*   tree->tip_order_score = tree->sum_y_dist; */
/*   tree->tip_order_score = tree->sum_y_dist_sq; */

  int i;
  tree->tip_order_score = -FLT_MAX;
  For(i,2*tree->n_otu-1)
    {
      if(!tree->noeud[i]->tax && tree->noeud[i]->y_width > tree->tip_order_score)
	tree->tip_order_score = tree->noeud[i]->y_width;
    }
}

/*********************************************************/

void Get_All_Y_Rank_Pre(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;
      phydbl v1,v2;

      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      Get_All_Y_Rank_Pre(d,d->v[i],tree);
	    }
	}

      v1 = v2 = -1.;
      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      if(v1 < .0) v1 = d->v[i]->y_rank;
	      else        v2 = d->v[i]->y_rank;
	    }
	}

      d->y_rank            = (v1+v2)/2.;
      d->y_width           = FABS(v1-v2);
      tree->sum_y_dist_sq += POW(v1-v2,2);
      tree->sum_y_dist    += FABS(v1-v2);
    }
}

/*********************************************************/

phydbl Get_Tip_Order_Score(int n_trees, t_tree **list_tree, t_tree *ref_tree)
{
  int i,j,k;
  phydbl score;

  Get_Tips_Y_Rank(ref_tree);
  
  score = .0;
  For(i,n_trees)
    {
      For(j,ref_tree->n_otu) list_tree[i]->noeud[j]->y_rank = list_tree[i]->noeud[j]->ext_node->y_rank;      
      Get_All_Y_Rank(list_tree[i]);
      score += list_tree[i]->tip_order_score;
    }

  return score;
}

/*********************************************************/

void Swap_One_Node(t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;
      int dir1, dir2;
      t_node *tmp_n;
      t_edge *tmp_e;

      if(d != tree->n_root)
	{
	  dir1 = dir2 = -1;
	  For(i,3)
	    {
	      if((d->v[i] != d->anc) && (d->b[i] != tree->e_root))
		{
		  if(dir1 < 0) dir1 = i;
		  else         dir2 = i;
		}
	    }
	}
      else
	{
	  dir1 = 0;
	  dir2 = 1;
	}

      tmp_n      = d->v[dir2];
      d->v[dir2] = d->v[dir1];
      d->v[dir1] = tmp_n;

      tmp_e      = d->b[dir2];
      d->b[dir2] = d->b[dir1];
      d->b[dir1] = tmp_e;
    }
}

/*********************************************************/

void Minimize_Tip_Order_Score(int n_trees, t_tree **list_tree, t_tree *ref_tree)
{
  int i;
  phydbl score,min_score,old_min_score;
  phydbl diff,eps;
  t_node **node_table;
  int swapped;
  t_node *tmp;

  eps = 1.E-3;
  old_min_score = min_score = score = INT_MAX;
  /*   Print_Tip_Ordered(ref_tree); */

  do
    {
      For(i,2*ref_tree->n_otu-1)
	{	  
	  Swap_One_Node(ref_tree->noeud[i],ref_tree);
	  Get_Tips_Y_Rank(ref_tree);
	  score = (phydbl)Untangle_Tree_List(n_trees,list_tree,ref_tree);
	  if(score == -1) 
	    {
	      return;
	    }

	  if(score < min_score) 
	    {
	      min_score = score;
	      PhyML_Printf("\n- Score = %f",score);
	    }
	  else 
	    {
	      Swap_One_Node(ref_tree->noeud[i],ref_tree);    
	      Get_Tips_Y_Rank(ref_tree);
/* 	      PhyML_Printf("\n+ Score = %f",score); */
	    }
	}
      diff = fabs(old_min_score - min_score);
      old_min_score = min_score;
    }while(diff > eps);

  PhyML_Printf("\n");


  node_table = (t_node **)mCalloc(ref_tree->n_otu,sizeof(t_node *));

  For(i,ref_tree->n_otu) node_table[i] = ref_tree->noeud[i];

/*       bubble sort of conflict nodes according to their y_rank */
  do
    {
      swapped = NO;
      For(i,ref_tree->n_otu-1)
	{
	  if(node_table[i]->y_rank > node_table[i+1]->y_rank)
	    {
	      swapped = YES;
	      tmp             = node_table[i];
	      node_table[i]   = node_table[i+1];
	      node_table[i+1] = tmp;
	    }
	}
    }while(swapped == YES);

  For(i,ref_tree->n_otu)
    {
      PhyML_Printf("\n%s",node_table[i]->name,node_table[i]->y_rank);
    }



  Free(node_table);

}

/*********************************************************/

void Print_Tip_Ordered(t_tree *tree)
{
  Print_Tip_Ordered_Pre(tree->n_root,tree->n_root->v[0],tree);
  Print_Tip_Ordered_Pre(tree->n_root,tree->n_root->v[1],tree);
}

/*********************************************************/

void Print_Tip_Ordered_Pre(t_node *a, t_node *d, t_tree *tree)
{
  
  if(d->tax)
    {
      PhyML_Printf("\n. %f \"%s\"",d->y_rank,d->name);
    }
  else
    {
      int i,dir1,dir2;

      dir1 = dir2 = -1;
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      if(dir1 < 0) dir1 = i;
	      else         dir2 = i;
	    }
	}
      if(d->v[dir1]->y_rank < d->v[dir2]->y_rank)
	{
	  Print_Tip_Ordered_Pre(d,d->v[dir1],tree);
	  Print_Tip_Ordered_Pre(d,d->v[dir2],tree);
	}
      else
	{
	  Print_Tip_Ordered_Pre(d,d->v[dir2],tree);
	  Print_Tip_Ordered_Pre(d,d->v[dir1],tree);
	}
    }
}

/*********************************************************/

int Untangle_Tree_List(int n_trees, t_tree **list_tree, t_tree *ref_tree)
{
  int i,j;
  int tree_score,score;

  score = 0;
  For(i,n_trees) 
    {
/*       PhyML_Printf("\n. Untangling tree %d",i); */
      For(j,ref_tree->n_otu) list_tree[i]->noeud[j]->y_rank = list_tree[i]->noeud[j]->ext_node->y_rank;
      tree_score = Untangle_Tree(list_tree[i]);
      score += tree_score;
      if(tree_score < 0) 
	{
	  return -1;
	}
    }

  return score;
}

/*********************************************************/

int Untangle_Tree(t_tree *tree)
{
  int conflict;
  int n_trials;
  t_node **node_table;
  int i,swapped;
  t_node *tmp;

  node_table = (t_node **)mCalloc(tree->n_otu,sizeof(t_node *));


  For(i,tree->n_otu) node_table[i] = tree->noeud[i];

/*       bubble sort of conflict nodes according to their y_rank */
  do
    {
      swapped = NO;
      For(i,tree->n_otu-1)
	{
	  if(node_table[i]->y_rank > node_table[i+1]->y_rank)
	    {
	      swapped = YES;
	      tmp             = node_table[i];
	      node_table[i]   = node_table[i+1];
	      node_table[i+1] = tmp;
	    }
	}
    }while(swapped == YES);
  

  tree->tip_order_score = 0;

/*   PhyML_Printf("\n\n Enter"); */
  n_trials = 0;
  do
    {
/*       PhyML_Printf("\n >>>>>>>>>> Pass"); */
      conflict= NO;
      Get_All_Y_Rank(tree);
      Untangle_Node(tree->n_root,tree->n_root->v[0],node_table,&conflict,tree);
      Untangle_Node(tree->n_root,tree->n_root->v[1],node_table,&conflict,tree);
      n_trials++;
      if(n_trials > 100) 
	{
	  int i;
	  FILE *ps_tree;
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  ps_tree = (FILE *)fopen("failed_tree.ps","w");
	  DR_Print_Postscript_Header(1,ps_tree);
	  tree->ps_tree = DR_Make_Tdraw_Struct(tree);
	  DR_Init_Tdraw_Struct(tree->ps_tree);
	  DR_Get_Tree_Box_Width(tree->ps_tree,tree);
	  Dist_To_Root(tree->n_root,tree);
	  tree->ps_tree->max_dist_to_root = DR_Get_Max_Dist_To_Root(tree);
	  For(i,tree->n_otu) tree->ps_tree->ycoord[i] = tree->noeud[i]->y_rank * (int)(tree->ps_tree->page_height / (tree->n_otu));
	  DR_Get_X_Coord(NO,tree->ps_tree,tree);
	  DR_Get_Y_Coord(YES,tree->ps_tree,tree);
	  DR_Print_Tree_Postscript(1,YES,ps_tree,tree);
	  DR_Print_Postscript_EOF(ps_tree);
	  fclose(ps_tree);
	  Warn_And_Exit("");	  
	}
    }while(conflict == YES);  

/*   PhyML_Printf("\n <<<<<<<<<< Leave"); */

  Free(node_table);
  return tree->tip_order_score;
}

/*********************************************************/

void Untangle_Node(t_node *a, t_node *d, t_node **node_table, int *conflict, t_tree *tree)
{
  if(*conflict == YES) return;

  if(d->tax) return;
  else
    {
      int i,j;
      int d_a;
      int beg,end;
      phydbl min,max;
      int n_down, n_up;
      t_node *lca;
      t_node **conflict_nodes;
      int n_conflicts;
      int swapped;
      phydbl eps,tmp_rank;
      t_node *tmp_node;
      int n_moved;
      eps = 1.E-8;

      For(i,3)
	{
	  if((d->v[i] != d->anc) && (d->b[i] != tree->e_root))
	    {
	      Untangle_Node(d,d->v[i],node_table,conflict,tree);
	    }
	}
      
      if(*conflict == YES) return;

      lca = NULL;

      For(i,3)
	{
	  if((d->v[i] == d->anc) || (d->b[i] == tree->e_root))
	    {
	      d_a = i;
	      break;
	    }
	}

      min =  FLT_MAX;
      max = -FLT_MAX;
      For(i,d->bip_size[d_a])
	{
	  if(d->bip_node[d_a][i]->y_rank < min)
	    min = d->bip_node[d_a][i]->y_rank;

	  if(d->bip_node[d_a][i]->y_rank > max)
	    max = d->bip_node[d_a][i]->y_rank;	  
	}


      n_conflicts = 0;
      For(i,tree->n_otu)
	{
	  if((node_table[i]->y_rank > min - eps) && (node_table[i]->y_rank < max + eps))
	    {
/* 	      conflict_nodes[n_conflicts] = tree->noeud[i]; */
	      n_conflicts++;	      
	    }
	}

      For(i,tree->n_otu)
	{
	  if(node_table[i]->y_rank > min - eps)
	    {
	      conflict_nodes = node_table+i;
	      break;
	    }
	}



      /* bubble sort of conflict nodes according to their y_rank */
      /*       do */
      /* 	{ */
      /* 	  swapped = NO; */
      /* 	  For(i,n_conflicts-1)  */
      /* 	    { */
      /* 	      if(conflict_nodes[i]->y_rank > conflict_nodes[i+1]->y_rank) */
      /* 		{ */
      /* 		  swapped = YES; */
      /* 		  lca                 = conflict_nodes[i]; */
      /* 		  conflict_nodes[i]   = conflict_nodes[i+1]; */
      /* 		  conflict_nodes[i+1] = lca; */
      /* 		} */
      /* 	    }	     */
      /* 	}while(swapped == YES); */
      

/*       if(n_conflicts) */
/* 	{ */
/* 	  PhyML_Printf("\n"); */
/* 	  For(j,d->bip_size[d_a]) printf("\n. d:%d in bip: %s",d->num,d->bip_node[d_a][j]->name); */
/* 	} */

      beg = 0;
      end = n_conflicts;
      n_moved = 0;
      do
	{
/* 	  printf("\n. v0?%s v1?%s %3d n_moved = %d n_sons=%d n_conflicts=%d min=%f max=%f", */
/* 		 d == tree->n_root->v[0]?"YES":"NO", */
/* 		 d == tree->n_root->v[1]?"YES":"NO", */
/* 		 d->num,n_moved,d->bip_size[d_a],n_conflicts,min,max); */
/* 	  For(i,d->bip_size[d_a]) printf("\n. %d %f",d->bip_node[d_a][i]->num,d->bip_node[d_a][i]->y_rank); */

	  for(i=beg;i<end;i++)
	    {
	      For(j,d->bip_size[d_a]) 
		if(conflict_nodes[i] == d->bip_node[d_a][j]) 
		  break;
	      
	      if(j == d->bip_size[d_a])
		{
		  *conflict = YES;

		  n_moved++;
		  
		  lca = Find_Lca(d,conflict_nodes[i],tree);
		  
/* 		  PhyML_Printf("\n. Detected conflict for ``%s'' (rank:%f min=%f max=%f lca=%f)", */
/* 			       conflict_nodes[i]->name, */
/* 			       conflict_nodes[i]->y_rank, */
/* 			       min,max,lca->y_rank); */
		  

		  if(lca->y_rank > d->y_rank)
		    {
		      end--;
		      for(j=i;j<n_conflicts-1;j++)
			{
/* 			  PhyML_Printf("\n+ Moved (%s,%s) from (%f,%f)", */
/* 				       conflict_nodes[j]->name,conflict_nodes[j+1]->name, */
/* 				       conflict_nodes[j]->y_rank,conflict_nodes[j+1]->y_rank); */

			  tmp_rank                  = conflict_nodes[j]->y_rank;
			  conflict_nodes[j]->y_rank = conflict_nodes[j+1]->y_rank;
			  conflict_nodes[j+1]->y_rank = tmp_rank;

			  tmp_node          = conflict_nodes[j];
			  conflict_nodes[j] = conflict_nodes[j+1];
			  conflict_nodes[j+1] = tmp_node;		     

/* 			  PhyML_Printf(" to (%s,%s) (%f,%f)", */
/* 				       conflict_nodes[j]->name,conflict_nodes[j+1]->name, */
/* 				       conflict_nodes[j]->y_rank,conflict_nodes[j+1]->y_rank); */

			  tree->tip_order_score++;
			}
		    }
		  else
		    {
		      beg++;
		      for(j=i;j>0;j--)
			{
/* 			  PhyML_Printf("\n- Moved (%s,%s) from (%f,%f)", */
/* 				       conflict_nodes[j]->name,conflict_nodes[j-1]->name, */
/* 				       conflict_nodes[j]->y_rank,conflict_nodes[j-1]->y_rank); */

			  tmp_rank                  = conflict_nodes[j]->y_rank;
			  conflict_nodes[j]->y_rank = conflict_nodes[j-1]->y_rank;
			  conflict_nodes[j-1]->y_rank = tmp_rank;
			  

			  tmp_node          = conflict_nodes[j];
			  conflict_nodes[j] = conflict_nodes[j-1];
			  conflict_nodes[j-1] = tmp_node;		     
			  
/* 			  PhyML_Printf(" to (%f,%f)",conflict_nodes[j]->y_rank,conflict_nodes[j-1]->y_rank); */

			  tree->tip_order_score++;
			}
		    }
		  		 
/* 		  if(!Check_Tip_Ranks(tree)) */
/* 		    { */
/* 		      int i; */
/* 		      FILE *ps_tree; */
/* 		      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__); */
/* 		      ps_tree = (FILE *)fopen("failed_tree.ps","w"); */
/* 		      DR_Print_Postscript_Header(1,ps_tree); */
/* 		      tree->ps_tree = DR_Make_Tdraw_Struct(tree); */
/* 		      DR_Init_Tdraw_Struct(tree->ps_tree); */
/* 		      DR_Get_Tree_Box_Width(tree->ps_tree,tree); */
/* 		      Dist_To_Root(tree->n_root,tree); */
/* 		      tree->ps_tree->max_dist_to_root = DR_Get_Max_Dist_To_Root(tree); */
/* 		      For(i,tree->n_otu) tree->ps_tree->ycoord[i] = tree->noeud[i]->y_rank * (int)(tree->ps_tree->page_height / (tree->n_otu)); */
/* 		      DR_Get_X_Coord(NO,tree->ps_tree,tree); */
/* 		      DR_Get_Y_Coord(YES,tree->ps_tree,tree); */
/* 		      DR_Print_Tree_Postscript(1,YES,ps_tree,tree); */
/* 		      DR_Print_Postscript_EOF(ps_tree); */
/* 		      fclose(ps_tree); */
/* 		      Warn_And_Exit(""); */
/* 		    } */

		  break;
		}      
	    }
	}while(n_moved + d->bip_size[d_a] != n_conflicts);

/*       Free(conflict_nodes); */
      return;
    }
}

/*********************************************************/

int Check_Tip_Ranks(t_tree *tree)
{
  int i,j;
  phydbl eps;

  eps = 1.E-6;

  For(i,tree->n_otu-1)
    {
      for(j=i+1;j<tree->n_otu;j++)
	{
	  if(fabs(tree->noeud[i]->y_rank - tree->noeud[j]->y_rank) < eps)
	    {
	      return 0;
	    }
	}
    }
  return 1;
}

/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/

#endif
