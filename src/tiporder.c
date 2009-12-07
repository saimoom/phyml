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
  t_tree **list_tree,*ref_tree,*tree;
  FILE *fp_ref_tree,*fp_list_tree,*fp_coord,*ps_tree;
  int i,j,k;
  int n_trees;
  option *ref_io,*list_io;
  char **name_table;

  ref_io  = (option *)Make_Input();
  list_io = (option *)Make_Input();

  fp_ref_tree  = (FILE *)fopen(argv[1],"r");
  fp_list_tree = (FILE *)fopen(argv[2],"r");
  fp_coord     = (FILE *)fopen(argv[3],"r");

  if(!fp_ref_tree) 
    {
      PhyML_Printf("\n. Can't find %s",argv[1]);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(!fp_list_tree) 
    {
      PhyML_Printf("\n. Can't find %s",argv[2]);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  ref_io->fp_in_tree = fp_ref_tree;
  Read_Tree_File(ref_io);
  fclose(ref_io->fp_in_tree);
  ref_tree = ref_io->tree;
  ref_tree->io = ref_io;

  list_io->fp_in_tree = fp_list_tree;
  Read_Tree_File(list_io);
  fclose(list_io->fp_in_tree);
  list_tree = list_io->treelist->tree;
  n_trees = list_io->treelist->list_size;
  PhyML_Printf("\n. Read %d trees\n",n_trees);
 
  For(i,n_trees) list_tree[i]->io = list_io;

  name_table = (char **)mCalloc(ref_tree->n_otu,sizeof(char **));
  For(i,ref_tree->n_otu) name_table[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));


  /* Sort translation table such that tree->noeud[i]->name == tree->io->short_tax_name[i] for all i */
  Sort_Translation_Table(ref_tree);
  Read_Taxa_Zscores(fp_coord,ref_tree);


  /* Find matching tips */
  For(i,n_trees)
    {
      Sort_Translation_Table(list_tree[i]);

      For(j,ref_tree->n_otu) 
	{
	  For(k,ref_tree->n_otu) 
	    {
	      if(!strcmp(ref_io->long_tax_names[j],list_io->long_tax_names[k]))
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

  PhyML_Printf("\n. Getting ancestors"); fflush(NULL);
  For(i,n_trees)   
    {
      Update_Ancestors(list_tree[i]->n_root,list_tree[i]->n_root->v[0],list_tree[i]);
      Update_Ancestors(list_tree[i]->n_root,list_tree[i]->n_root->v[1],list_tree[i]);
      list_tree[i]->n_root->anc = NULL;
    }

  PhyML_Printf("\n. Getting bipartitions"); fflush(NULL);

  Free_Bip(ref_tree);
  Alloc_Bip(ref_tree);
  Get_Bip(ref_tree->noeud[0],
	  ref_tree->noeud[0]->v[0],
	  ref_tree);

  For(i,n_trees) 
    {
      if(!(i%10))
      PhyML_Printf("\n. Getting bipartition for tree %d",i);
      Free_Bip(list_tree[i]);
      Alloc_Bip(list_tree[i]);
      Get_Bip(list_tree[i]->noeud[0],
	      list_tree[i]->noeud[0]->v[0],
	      list_tree[i]);
    }


  PhyML_Printf("\n. Getting tip ranks"); fflush(NULL);
/*   Get_Tips_Y_Rank(ref_tree); */

  Get_Tips_Y_Rank_From_Zscores(ref_tree);

/*   PhyML_Printf("\n. Minimizing"); fflush(NULL); */
/*   Minimize_Tip_Order_Score(n_trees,list_tree,ref_tree); */



  ps_tree  = (FILE *)fopen("order_tree.ps","w");

  Test_Node_Table_Consistency(ref_tree);

  ref_tree->ps_tree = DR_Make_Tdraw_Struct(ref_tree);
  DR_Get_Tree_Coord(ref_tree);
  For(j,ref_tree->n_otu) 
    {
      ref_tree->ps_tree->ycoord[j] = 
	(ref_tree->noeud[j]->y_rank/ref_tree->n_otu)*
	ref_tree->ps_tree->page_height;
    }

  list_io->z_scores = (phydbl *)mCalloc(ref_tree->n_otu,sizeof(phydbl));

  DR_Print_Postscript_Header(1,ps_tree);
  For(i,n_trees)
    {
      tree = list_tree[i];

      Test_Node_Table_Consistency(tree);
      tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
      RATES_Init_Rate_Struct(tree->rates,tree->n_otu);
      MC_Least_Square_Node_Times(tree->e_root,tree);
      MC_Adjust_Node_Times(tree);
      RATES_Update_Cur_Bl(tree);

      tree->ps_tree = DR_Make_Tdraw_Struct(tree);
      DR_Init_Tdraw_Struct(tree->ps_tree);
      DR_Get_Tree_Box_Width(tree->ps_tree,tree);
      Dist_To_Root(tree->n_root,tree);
      tree->ps_tree->max_dist_to_root = DR_Get_Max_Dist_To_Root(tree);
 
      For(j,ref_tree->n_otu) tree->io->z_scores[j] = ref_tree->io->z_scores[tree->noeud[j]->ext_node->num];
      Get_Tips_Y_Rank_From_Zscores(tree);
      Untangle_Tree(tree);
      For(j,ref_tree->n_otu) tree->ps_tree->ycoord[j] =  (tree->noeud[j]->y_rank/tree->n_otu)*tree->ps_tree->page_height;

/*       For(j,ref_tree->n_otu) tree->ps_tree->ycoord[j] = ref_tree->ps_tree->ycoord[tree->noeud[j]->ext_node->num]; */
      For(j,ref_tree->n_otu) list_io->z_scores[j] = ref_io->z_scores[tree->noeud[j]->ext_node->num];

      DR_Get_Y_Coord(YES,tree->ps_tree,tree);
      DR_Get_X_Coord( NO,tree->ps_tree,tree);

      if(!i) DR_Print_Tree_Postscript(1,YES,ps_tree,tree);
      else   DR_Print_Tree_Postscript(1, NO,ps_tree,tree);
    }
  DR_Print_Postscript_EOF(ps_tree);
  fclose(ps_tree);



  ps_tree  = (FILE *)fopen("ref_tree.ps","w");

  Test_Node_Table_Consistency(ref_tree);


  DR_Print_Postscript_Header(1,ps_tree);
  tree = ref_tree;
  tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
  RATES_Init_Rate_Struct(tree->rates,tree->n_otu);
  MC_Least_Square_Node_Times(tree->e_root,tree);
  MC_Adjust_Node_Times(tree);
  RATES_Update_Cur_Bl(tree);
  tree->ps_tree = DR_Make_Tdraw_Struct(tree);
  DR_Init_Tdraw_Struct(tree->ps_tree);
  DR_Get_Tree_Box_Width(tree->ps_tree,tree);
  Dist_To_Root(tree->n_root,tree);
  tree->ps_tree->max_dist_to_root = DR_Get_Max_Dist_To_Root(tree);
  
  Get_Tips_Y_Rank_From_Zscores(tree);
  Untangle_Tree(tree);
  For(j,tree->n_otu) tree->ps_tree->ycoord[j] =  (tree->noeud[j]->y_rank/tree->n_otu)*tree->ps_tree->page_height;

  DR_Get_Y_Coord(YES,tree->ps_tree,tree);
  DR_Get_X_Coord( NO,tree->ps_tree,tree);
  DR_Print_Tree_Postscript(1,YES,ps_tree,tree);
  DR_Print_Postscript_EOF(ps_tree);
  fclose(ps_tree);

  fclose(fp_ref_tree);
  fclose(fp_list_tree);
  fclose(fp_coord);
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
  int i,j;
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
      for(i=ref_tree->n_otu;i<2*ref_tree->n_otu-1;i++)
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


  For(i,ref_tree->n_otu)
    {
      For(j,ref_tree->n_otu)
	{
	  if(!strcmp(ref_tree->io->short_tax_names[i],ref_tree->noeud[j]->name))
	    {
	      Free(ref_tree->noeud[j]->name);
	      ref_tree->noeud[j]->name = (char *)mCalloc((int)strlen(ref_tree->io->long_tax_names[i])+1,sizeof(char));
	      strcpy(ref_tree->noeud[j]->name,ref_tree->io->long_tax_names[i]);
	      break;
	    }
	}
    }

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
/*       PhyML_Printf("\n. Untangling tree %3d",i); */
      For(j,ref_tree->n_otu) list_tree[i]->noeud[j]->y_rank = list_tree[i]->noeud[j]->ext_node->y_rank;
      tree_score = Untangle_Tree(list_tree[i]);
/*       PhyML_Printf(" score = %3d",tree_score); */
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
      if(n_trials > 1000) 
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
      eps = (phydbl)1./(2.*tree->n_otu);

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


      beg = 0;
      end = n_conflicts;
      n_moved = 0;
      do
	{
	  for(i=beg;i<end;i++)
	    {
	      For(j,d->bip_size[d_a]) if(conflict_nodes[i] == d->bip_node[d_a][j]) break;
	      
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

			  tmp_rank                    = conflict_nodes[j]->y_rank;
			  conflict_nodes[j]->y_rank   = conflict_nodes[j+1]->y_rank;
			  conflict_nodes[j+1]->y_rank = tmp_rank;

			  tmp_node            = conflict_nodes[j];
			  conflict_nodes[j]   = conflict_nodes[j+1];
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

			  tmp_rank                    = conflict_nodes[j]->y_rank;
			  conflict_nodes[j]->y_rank   = conflict_nodes[j-1]->y_rank;
			  conflict_nodes[j-1]->y_rank = tmp_rank;
			  
			  tmp_node            = conflict_nodes[j];
			  conflict_nodes[j]   = conflict_nodes[j-1];
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

void Read_Taxa_Zscores(FILE *fp_coord, t_tree *tree)
{
  char *name,*line;
  phydbl z;
  int i;

  name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  fgets(line,T_MAX_LINE,fp_coord);
  Free(line);

  tree->io->z_scores = (phydbl *)mCalloc(tree->n_otu,sizeof(phydbl));

  do
    {
      if(fscanf(fp_coord,"%s\t%f\n",name,&z) == EOF) break;
      PhyML_Printf("\n. Read %s. Z-score: %f",name,z);

      For(i,tree->n_otu) if(!strcmp(tree->io->long_tax_names[i],name)) break;
      
      if(i == tree->n_otu)
	{
	  PhyML_Printf("\n. Could not find taxon '%s' in coordinate file.",name);
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      tree->io->z_scores[i] = z;
      
    }while(1);
	
  Free(name);
}

/*********************************************************/

void Read_Taxa_Coordinates(FILE *fp_coord, t_tree *tree)
{
  char *name,*line;
  phydbl lon, lat;
  int i;

  name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  fgets(line,T_MAX_LINE,fp_coord);
  Free(line);

  tree->io->lat = (phydbl *)mCalloc(tree->n_otu,sizeof(phydbl));
  tree->io->lon = (phydbl *)mCalloc(tree->n_otu,sizeof(phydbl));

  do
    {
      if(fscanf(fp_coord,"%s\t%f\t%f\n",name,&lat,&lon) == EOF) break;
      PhyML_Printf("\n. Read %s %f %f",name,lat,lon);

      For(i,tree->n_otu) if(!strcmp(tree->io->long_tax_names[i],name)) break;
      
      if(i == tree->n_otu)
	{
	  PhyML_Printf("\n. Could not find taxon '%s' in coordinate file.",name);
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      tree->io->lat[i] = lat;
      tree->io->lon[i] = lon;
      
    }while(1);
	
  Free(name);
}

/*********************************************************/

void Get_Tips_Y_Rank_From_Zscores(t_tree *tree)
{
  int i,j;

  For(i,tree->n_otu) tree->noeud[i]->y_rank = .0;

  For(i,tree->n_otu-1)
    {
      for(j=i+1;j<tree->n_otu;j++)
	{
	  if(tree->io->z_scores[i] > tree->io->z_scores[j])
	    {
	      tree->noeud[i]->y_rank += 1.0;
	    }
	  else
	    {
	      tree->noeud[j]->y_rank += 1.0;
	    }
	}
    }
}

/*********************************************************/
/* Sort translation table such that tree->noeud[i]->name == tree->io->short_tax_name[i] for all i */
void  Sort_Translation_Table(t_tree *tree)
{
  int i,j;
  char *s;


  Test_Node_Table_Consistency(tree);

  For(i,tree->n_otu-1)
    {
      for(j=i+1;j<tree->n_otu;j++)
	{
	  if(!strcmp(tree->noeud[i]->name,tree->io->short_tax_names[j]))
	    {
	      s = tree->io->short_tax_names[i];
	      tree->io->short_tax_names[i] = tree->io->short_tax_names[j];
	      tree->io->short_tax_names[j] = s;

	      s = tree->io->long_tax_names[i];
	      tree->io->long_tax_names[i] = tree->io->long_tax_names[j];
	      tree->io->long_tax_names[j] = s;
	      
	      break;
	    }
	}
    }
}

/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/

#endif
