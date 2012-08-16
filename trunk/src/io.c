/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "io.h"


/* Tree parser function. We need to pass a pointer to the string of characters
   since this string might be freed and then re-allocated by that function (i.e.,
   its address in memory might change)
*/
t_tree *Read_Tree(char **s_tree)
{
  char **subs;
  int i,n_ext,n_int,n_otu;
  t_tree *tree;
  int degree,len;
  t_node *root_node;
  
  n_int = n_ext = 0;
  
  n_otu=0;
  For(i,(int)strlen((*s_tree))) if((*s_tree)[i] == ',') n_otu++;
  n_otu+=1;

  tree = Make_Tree_From_Scratch(n_otu,NULL);
  subs = Sub_Trees((*s_tree),&degree);
  Clean_Multifurcation(subs,degree,3);

  if(degree == 2) 
    {
      /* Unroot_Tree(subs); */
      /* degree = 3; */
      /* root_node = tree->t_nodes[n_otu]; */
      root_node      = tree->t_nodes[2*n_otu-2];
      root_node->num = 2*n_otu-2;
      tree->n_root   = root_node;
      n_int         -= 1;
    }
  else
    {      
      root_node      = tree->t_nodes[n_otu];
      root_node->num = n_otu;
      tree->n_root   = NULL;
   }
  
  if(degree > 3) /* Multifurcation at the root. Need to re-assemble the subtrees
		    since Clean_Multifurcation added sets of parenthesis and
		    the corresponding NULL edges */
    {
      degree = 3;
      Free((*s_tree));
      len = 0;
      For(i,degree) len += (strlen(subs[i])+1);
      len += 5;

      (*s_tree) = (char *)mCalloc(len,sizeof(char));

      (*s_tree)[0] = '('; (*s_tree)[1] = '\0';
      For(i,degree) 
	{
	  strcat((*s_tree),subs[i]);
	  strcat((*s_tree),",\0");
	}

      sprintf((*s_tree)+strlen((*s_tree))-1,"%s",");\0");
      
      For(i,NODE_DEG_MAX) Free(subs[i]);
      Free(subs);

      subs = Sub_Trees((*s_tree),&degree);
    }

  root_node->tax = 0;

  tree->has_branch_lengths = 0;
  tree->num_curr_branch_available = 0;
  For(i,degree) R_rtree((*s_tree),subs[i],root_node,tree,&n_int,&n_ext);

  for(i=degree;i<NODE_DEG_MAX;i++) Free(subs[i]);
  Free(subs);

  if(tree->n_root)
    {
      tree->e_root = tree->t_edges[tree->num_curr_branch_available];
            
      For(i,3) if(tree->n_root->v[0]->v[i] == tree->n_root) { tree->n_root->v[0]->v[i] = tree->n_root->v[1]; break; }
      For(i,3) if(tree->n_root->v[1]->v[i] == tree->n_root) { tree->n_root->v[1]->v[i] = tree->n_root->v[0]; break; }

      Connect_One_Edge_To_Two_Nodes(tree->n_root->v[0],
      				    tree->n_root->v[1],
      				    tree->e_root,
      				    tree);

      tree->e_root->l->v = tree->n_root->l[0] + tree->n_root->l[1];
      if(tree->e_root->l->v > 0.0)
	tree->n_root_pos = tree->n_root->l[0] / tree->e_root->l->v;
      else
	tree->n_root_pos = .5;
    }
  
  return tree;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* 'a' in t_node a stands for ancestor. 'd' stands for descendant */ 
void R_rtree(char *s_tree_a, char *s_tree_d, t_node *a, t_tree *tree, int *n_int, int *n_ext)
{
  int i;
  t_node *d;
  int n_otu = tree->n_otu;

  if(strstr(s_tree_a," ")) 
    {
      PhyML_Printf("\n. [%s]",s_tree_a);
      Warn_And_Exit("\n. Err: the tree must not contain a ' ' character\n");
    }

  if(s_tree_d[0] == '(')
    {
      char **subs;
      int degree;
      t_edge *b;

      (*n_int)+=1;

      if((*n_int + n_otu) == (2*n_otu-1))
	{
	  PhyML_Printf("\n== The number of internal nodes in the tree exceeds the number of taxa minus one.");
	  PhyML_Printf("\n== There probably is a formating problem in the input tree.");
	  PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
	  Exit("\n");	  
	}

      d      = tree->t_nodes[n_otu+*n_int];
      d->num = n_otu+*n_int;
      d->tax = 0;
      b      = tree->t_edges[tree->num_curr_branch_available];

      Read_Branch_Label(s_tree_d,s_tree_a,tree->t_edges[tree->num_curr_branch_available]);
      Read_Branch_Length(s_tree_d,s_tree_a,tree);      

      For(i,3)
	{
	  if(!a->v[i])
	    {
	      a->v[i]=d;
	      d->l[0]=tree->t_edges[tree->num_curr_branch_available]->l->v;
	      a->l[i]=tree->t_edges[tree->num_curr_branch_available]->l->v;
	      break;
	    }
	}
      d->v[0]=a;

      if(a != tree->n_root)
	{
	  Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
	  tree->num_curr_branch_available++;
	}
      
      subs=Sub_Trees(s_tree_d,&degree);

      if(degree > 2)
	{
	  Clean_Multifurcation(subs,degree,2);

	  Free(s_tree_d);

	  s_tree_d = (char *)mCalloc(strlen(subs[0])+strlen(subs[1])+5,sizeof(char));

	  strcat(s_tree_d,"(");
	  strcat(s_tree_d,subs[0]);
	  strcat(s_tree_d,",");
	  strcat(s_tree_d,subs[1]);
	  strcat(s_tree_d,")");
	  For(i,b->n_labels)
	    {
	      strcat(s_tree_d,"#");
	      strcat(s_tree_d,b->labels[i]);
	    }
	  
	  For(i,NODE_DEG_MAX) Free(subs[i]);
	  Free(subs);

	  subs=Sub_Trees(s_tree_d,&degree);
	}

      R_rtree(s_tree_d,subs[0],d,tree,n_int,n_ext);
      R_rtree(s_tree_d,subs[1],d,tree,n_int,n_ext);
  
      for(i=2;i<NODE_DEG_MAX;i++) Free(subs[i]);
      Free(subs);
    }

  else
    {
      int i;

      d      = tree->t_nodes[*n_ext];
      d->tax = 1;

      Read_Node_Name(d,s_tree_d,tree);
      Read_Branch_Label(s_tree_d,s_tree_a,tree->t_edges[tree->num_curr_branch_available]); 
      Read_Branch_Length(s_tree_d,s_tree_a,tree);
      
      For(i,3)
	{
	 if(!a->v[i])
	   {
	     a->v[i]=d;
	     d->l[0]=tree->t_edges[tree->num_curr_branch_available]->l->v;
	     a->l[i]=tree->t_edges[tree->num_curr_branch_available]->l->v;
	     break;
	   }
	}
      d->v[0]=a;

      if(a != tree->n_root)
	{
	  Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
	  tree->num_curr_branch_available++;
	}

      d->num=*n_ext;
      (*n_ext)+=1;
    }
  
  Free(s_tree_d);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Read_Branch_Label(char *s_d, char *s_a, t_edge *b)
{
  char *sub_tp;
  char *p;
  int posp,posl;

  sub_tp = (char *)mCalloc(3+(int)strlen(s_d)+1,sizeof(char));
  /* sub_tp = (char *)mCalloc(T_MAX_LINE,sizeof(char)); */

  sub_tp[0] = '(';
  sub_tp[1] = '\0';
  strcat(sub_tp,s_d);
  strcat(sub_tp,"#");
  p = strstr(s_a,sub_tp);
  
  if(!p)
    {
      sub_tp[0] = ',';
      sub_tp[1] = '\0';
      strcat(sub_tp,s_d);
      strcat(sub_tp,"#");
      p = strstr(s_a,sub_tp);
    }

  b->n_labels = 0;
  if(p)
    {
      if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
      b->n_labels++;

      posp = strlen(s_d);
      while(p[posp] != '#') posp++;
      posp++;

      posl = 0;
      do 
	{
	  b->labels[b->n_labels-1][posl] = p[posp];
	  posl++;
	  posp++;
	  if(p[posp] == '#') 
	    { 
	      b->labels[b->n_labels-1][posl] = '\0';
	      b->n_labels++;
	      if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
	      posp++;
	      posl=0;
	    }
	}
      while((p[posp] != ':') &&
      	    (p[posp] != ',') &&
      	    (p[posp] != '('));

      b->labels[b->n_labels-1][posl] = '\0';
    }

  if(p)
    {
      /* if(b->n_labels == 1) */
      /* 	PhyML_Printf("\n. Read label '%s' on t_edge %3d.",b->labels[0],b->num); */
      /* else */
      /* 	{ */
      /* 	  PhyML_Printf("\n. Read labels "); */
      /* 	  For(i,b->n_labels) PhyML_Printf("'%s' ",b->labels[i]); */
      /* 	  PhyML_Printf("on t_edge %3d.",b->num); */
      /* 	} */

      if(!strcmp(b->labels[0],"NULL"))
	{
	  b->does_exist = NO;
	}
    }
  /* else */
  /*   { */
  /*     PhyML_Printf("\n. No label found on %s",s_d); */
  /*   } */
  Free(sub_tp);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Read_Branch_Length(char *s_d, char *s_a, t_tree *tree)
{
  char *sub_tp;
  char *p;
  t_edge *b;
  int i;

  b = tree->t_edges[tree->num_curr_branch_available];

  /* sub_tp = (char *)mCalloc(T_MAX_LINE,sizeof(char)); */
  sub_tp = (char *)mCalloc(10+strlen(s_d)+1,sizeof(char));

  For(i,b->n_labels)
    {
      strcat(s_d,"#");
      strcat(s_d,b->labels[i]);
    }

  sub_tp[0] = '(';
  sub_tp[1] = '\0';
  strcat(sub_tp,s_d);
  strcat(sub_tp,":");
  p = strstr(s_a,sub_tp);

  if(!p)
    {
      sub_tp[0] = ',';
      sub_tp[1] = '\0';
      strcat(sub_tp,s_d);
      strcat(sub_tp,":");
      p = strstr(s_a,sub_tp);
    }


  if(p)
    {
      b->l->v = atof((char *)p+(int)strlen(sub_tp));
      tree->has_branch_lengths = YES;
      b->does_exist = YES;
    }      
  else
    {
      b->l->v = -1.;
    }


  Free(sub_tp);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Read_Node_Name(t_node *d, char *s_tree_d, t_tree *tree)
{
  int i;
  
  if(!tree->t_edges[tree->num_curr_branch_available]->n_labels)
    {
      d->name = (char *)mCalloc(strlen(s_tree_d)+1,sizeof(char ));
      strcpy(d->name,s_tree_d);
    }
  else
    {
      i = 0;
      do
	{
	  d->name = (char *)realloc(d->name,(i+1)*sizeof(char ));
	  d->name[i] = s_tree_d[i];
	  i++;
	}
      while(s_tree_d[i] != '#');
      d->name[i] = '\0';
    }
  d->ori_name = d->name;

}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Clean_Multifurcation(char **subtrees, int current_deg, int end_deg)
{

  if(current_deg <= end_deg) return;
  else
    {
      char *s_tmp;
      int i;

      /* s_tmp = (char *)mCalloc(T_MAX_LINE,sizeof(char)); */
      s_tmp = (char *)mCalloc(10+
      			      (int)strlen(subtrees[0])+1+
      			      (int)strlen(subtrees[1])+1,
      			      sizeof(char));
      
      strcat(s_tmp,"(\0");
      strcat(s_tmp,subtrees[0]);
      strcat(s_tmp,",\0");
      strcat(s_tmp,subtrees[1]);
      strcat(s_tmp,")#NULL\0"); /* Add the label 'NULL' to identify a non-existing edge */
      Free(subtrees[0]);
      subtrees[0] = s_tmp;

      for(i=1;i<current_deg-1;i++) strcpy(subtrees[i],subtrees[i+1]);

      Clean_Multifurcation(subtrees,current_deg-1,end_deg);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


char **Sub_Trees(char *tree, int *degree)
{
  char **subs;
  int posbeg,posend;
  int i;

  if(tree[0] != '(') {*degree = 1; return NULL;}

  subs=(char **)mCalloc(NODE_DEG_MAX,sizeof(char *));

  For(i,NODE_DEG_MAX) subs[i]=(char *)mCalloc(strlen(tree)+1,sizeof(char));


  posbeg=posend=1;
  (*degree)=0;
  do
    {
      posbeg = posend;
      if(tree[posend] != '(')
	{
	  while((tree[posend] != ',' ) &&
		(tree[posend] != ':' ) &&
		(tree[posend] != '#' ) &&
		(tree[posend] != ')' )) 
	    {
	      posend++ ;
	    }
	  posend -= 1;
	}
      else posend=Next_Par(tree,posend);

      while((tree[posend+1] != ',') &&
	    (tree[posend+1] != ':') &&
	    (tree[posend+1] != '#') &&
	    (tree[posend+1] != ')')) {posend++;}


      strncpy(subs[(*degree)],tree+posbeg,posend-posbeg+1);
/*       strcat(subs[(*degree)],"\0"); */
      subs[(*degree)][posend-posbeg+1]='\0'; /* Thanks to Jean-Baka Domelevo-Entfellner */

      posend += 1;
      while((tree[posend] != ',') &&
	    (tree[posend] != ')')) {posend++;}
      posend+=1;


      (*degree)++;
      if((*degree) == NODE_DEG_MAX)
	{
	  For(i,(*degree))
	    PhyML_Printf("\n. Subtree %d : %s\n",i+1,subs[i]);

	  PhyML_Printf("\n. The degree of a t_node cannot be greater than %d\n",NODE_DEG_MAX);
	  Warn_And_Exit("\n");
	}
    }
  while(tree[posend-1] != ')');

  return subs;
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Next_Par(char *s, int pos)
{
  int curr;

  curr=pos+1;

  while(*(s+curr) != ')')
    {
      if(*(s+curr) == '(') curr=Next_Par(s,curr);
      curr++;
    }

  return curr;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Print_Tree(FILE *fp, t_tree *tree)
{
  char *s_tree;
  int i;

  s_tree = (char *)Write_Tree(tree,NO);

  if(OUTPUT_TREE_FORMAT == NEWICK) PhyML_Fprintf(fp,"%s\n",s_tree);
  else if(OUTPUT_TREE_FORMAT == NEXUS)
    {
      PhyML_Fprintf(fp,"#NEXUS\n");
      PhyML_Fprintf(fp,"BEGIN TREES;\n");
      PhyML_Fprintf(fp,"\tTRANSLATE\n");
      For(i,tree->n_otu) PhyML_Fprintf(fp,"\t%3d\t%s,\n",i+1,tree->t_nodes[i]->name);
      PhyML_Fprintf(fp,"\tUTREE PAUP_1=\n");
      PhyML_Fprintf(fp,"%s\n",s_tree);
      PhyML_Fprintf(fp,"ENDBLOCK;");
    }
  Free(s_tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


char *Write_Tree(t_tree *tree, int custom)
{
  char *s;
  int i,available;
  int init_len;

  init_len = 3*(int)T_MAX_NAME;

#ifndef MPI
  s=(char *)mCalloc(init_len,sizeof(char));
  available = init_len;
#elif defined MPI
  s=(char *)mCalloc(T_MAX_LINE,sizeof(char));
#endif

  
  s[0]='(';

  if(custom == NO)
    {
      if(!tree->n_root)
	{
	  i = 0;
	  while((!tree->t_nodes[tree->n_otu+i]->v[0]) ||
		(!tree->t_nodes[tree->n_otu+i]->v[1]) ||
		(!tree->t_nodes[tree->n_otu+i]->v[2])) i++;
	  
	  R_wtree(tree->t_nodes[tree->n_otu+i],tree->t_nodes[tree->n_otu+i]->v[0],&available,&s,tree);
	  R_wtree(tree->t_nodes[tree->n_otu+i],tree->t_nodes[tree->n_otu+i]->v[1],&available,&s,tree);
	  R_wtree(tree->t_nodes[tree->n_otu+i],tree->t_nodes[tree->n_otu+i]->v[2],&available,&s,tree);
	}
      else
	{
	  R_wtree(tree->n_root,tree->n_root->v[0],&available,&s,tree);
	  R_wtree(tree->n_root,tree->n_root->v[1],&available,&s,tree);
	}
    }
  else
    {
      PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  s[(int)strlen(s)-1]=')';
  s[(int)strlen(s)]=';';

  return s;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void R_wtree(t_node *pere, t_node *fils, int *available, char **s_tree, t_tree *tree)
{
  int i,p;
  char *format;
  int last_len;

  format = (char *)mCalloc(100,sizeof(char));

  sprintf(format,"%%.%df",tree->bl_ndigits);

  p = -1;
  if(fils->tax)
    {
      /* printf("\n- Writing on %p",*s_tree); */
      /* ori_len = *pos; */

      last_len = (int)strlen(*s_tree);

      if(OUTPUT_TREE_FORMAT == NEWICK)
	{
	  if(tree->write_tax_names == YES)
	    {
	      if(tree->io && tree->io->long_tax_names) 
		{
		  strcat(*s_tree,tree->io->long_tax_names[fils->num]);
		}
	      else
		{
		  strcat(*s_tree,fils->name);
		}
	    }
	  else if(tree->write_tax_names == NO)
	    {
	      sprintf(*s_tree+(int)strlen(*s_tree),"%d",fils->num);
	    }
	}
      else if(OUTPUT_TREE_FORMAT == NEXUS)
	{
	  sprintf(*s_tree+(int)strlen(*s_tree),"%d",fils->num+1);
	}
      else
	{
	  PhyML_Printf("\n. Unknown tree format.");
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  PhyML_Printf("\n. s=%s\n",*s_tree);
	}

      if((fils->b) && (fils->b[0]) && (fils->b[0]->l->v > -1.))
	{
	  (*s_tree)[(int)strlen(*s_tree)] = ':';

#ifndef PHYTIME
	  if(!tree->n_root)
	    {
	      sprintf(*s_tree+(int)strlen(*s_tree),format,fils->b[0]->l->v);
	    }
	  else
	    {
	      if(pere == tree->n_root)
		{
		  phydbl root_pos = (fils == tree->n_root->v[0])?(tree->n_root_pos):(1.-tree->n_root_pos);
		  sprintf(*s_tree+(int)strlen(*s_tree),format,tree->e_root->l->v * root_pos);
		}
	      else
		{
		  sprintf(*s_tree+(int)strlen(*s_tree),format,fils->b[0]->l->v);
		}
	    }		
#else
	  if(!tree->n_root)
	    {
	      sprintf(*s_tree+(int)strlen(*s_tree),format,fils->b[0]->l->v);
	    }
	  else
	    {
	      sprintf(*s_tree+(int)strlen(*s_tree),format,tree->rates->cur_l[fils->num]);
	    }
#endif
	}

      /* strcat(*s_tree,","); */
      (*s_tree)[(int)strlen(*s_tree)] = ',';


#ifndef MPI      
      (*available) -= ((int)strlen(*s_tree) - last_len);

      /* printf("\n0 Available = %d [%d %d]",(*available),(int)strlen(*s_tree),last_len); */
      /* printf("\n0 %s [%d,%d]",*s_tree,(int)(int)strlen(*s_tree),*available); */

      if(*available < 0)
	{
	  PhyML_Printf("\n== s=%s\n",*s_tree);
	  PhyML_Printf("\n== len=%d\n",(int)strlen(*s_tree));
	  PhyML_Printf("\n== The sequence names in your input file might be too long.");
	  PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      if(*available < (int)T_MAX_NAME)
      	{
      	  (*s_tree) = (char *)mRealloc(*s_tree,(int)strlen(*s_tree)+3*(int)T_MAX_NAME,sizeof(char));
	  For(i,3*(int)T_MAX_NAME) (*s_tree)[(int)strlen(*s_tree)+i] = '\0';
      	  (*available) = 3*(int)T_MAX_NAME;
	  /* printf("\n. ++ 0 Available = %d",(*available)); */
      	}
#endif

    }
  else
    {
      (*s_tree)[(int)strlen(*s_tree)]='(';

#ifndef MPI
      (*available) -= 1;

      /* printf("\n1 Available = %d [%d %d]",(*available),(int)strlen(*s_tree),last_len); */
      /* printf("\n1 %s [%d,%d]",*s_tree,(int)(int)strlen(*s_tree),*available); */

      if(*available < (int)T_MAX_NAME)
      	{
      	  (*s_tree) = (char *)mRealloc(*s_tree,(int)strlen(*s_tree)+3*(int)T_MAX_NAME,sizeof(char));
	  For(i,3*(int)T_MAX_NAME) (*s_tree)[(int)strlen(*s_tree)+i] = '\0';
      	  (*available) = 3*(int)T_MAX_NAME;
	  /* printf("\n. ++ 1 Available = %d",(*available)); */
      	}
#endif
      /* (*available)--; */

      /* if(*available < (int)T_MAX_NAME/2) */
      /* 	{ */
      /* 	  (*s_tree) = (char *)mRealloc(*s_tree,*pos+(int)T_MAX_NAME,sizeof(char)); */
      /* 	  (*available) = (int)T_MAX_NAME; */
      /* 	} */


      if(tree->n_root)
	{
	  For(i,3)
	    {
	      if((fils->v[i] != pere) && (fils->b[i] != tree->e_root))
		R_wtree(fils,fils->v[i],available,s_tree,tree);
	      else p=i;
	    }
	}
      else
	{
	  For(i,3)
	    {
	      if(fils->v[i] != pere)
		R_wtree(fils,fils->v[i],available,s_tree,tree);
	      else p=i;
	    }
	}
      
      if(p < 0)
	{
	  PhyML_Printf("\n== fils=%p root=%p root->v[0]=%p root->v[1]=%p",fils,tree->n_root,tree->n_root->v[0],tree->n_root->v[1]);
	  PhyML_Printf("\n== tree->e_root=%p fils->b[0]=%p fils->b[1]=%p fils->b[2]=%p",tree->e_root,fils->b[0],fils->b[1],fils->b[2]);		       
	  PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
	  Exit("\n");
	}

      last_len = (int)strlen(*s_tree);

      (*s_tree)[last_len-1] = ')';

      if((fils->b) && (fils->b[p]->l->v > -1.))
	{
	  if(tree->print_boot_val)
	    {
	      sprintf(*s_tree+(int)strlen(*s_tree),"%d",fils->b[p]->bip_score);
	    }
	  else if(tree->print_alrt_val)
	    {
	      sprintf(*s_tree+(int)strlen(*s_tree),"%f",fils->b[p]->ratio_test);
	    }
	  
	  fflush(NULL);

	  (*s_tree)[(int)strlen(*s_tree)] = ':';

#ifndef PHYTIME
	  if(!tree->n_root)
	    {
	      sprintf(*s_tree+(int)strlen(*s_tree),format,fils->b[p]->l->v);
	    }
	  else
	    {
	      if(pere == tree->n_root)
		{
		  phydbl root_pos = (fils == tree->n_root->v[0])?(tree->n_root_pos):(1.-tree->n_root_pos);
		  sprintf(*s_tree+(int)strlen(*s_tree),format,tree->e_root->l->v * root_pos);
		}
	      else
		{
		  sprintf(*s_tree+(int)strlen(*s_tree),format,fils->b[p]->l->v);
		}
	    }
#else
	  if(!tree->n_root)
	    {
	      sprintf(*s_tree+(int)strlen(*s_tree),format,fils->b[p]->l->v);
	    }
	  else
	    {
	      sprintf(*s_tree+(int)strlen(*s_tree),format,tree->rates->cur_l[fils->num]);
	    }
#endif	  
	}

      /* strcat(*s_tree,","); */
      (*s_tree)[(int)strlen(*s_tree)] = ',';
      

#ifndef MPI
      (*available) -= ((int)strlen(*s_tree) - last_len);
      
      /* printf("\n2 Available = %d [%d %d]",(*available),(int)strlen(*s_tree),last_len); */
      /* printf("\n2 %s [%d,%d]",*s_tree,(int)(int)strlen(*s_tree),*available); */

      if(*available < 0)
	{
	  PhyML_Printf("\n. available = %d",*available);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      if(*available < (int)T_MAX_NAME)
      	{
      	  (*s_tree) = (char *)mRealloc(*s_tree,(int)strlen(*s_tree)+3*(int)T_MAX_NAME,sizeof(char));
	  For(i,3*(int)T_MAX_NAME) (*s_tree)[(int)strlen(*s_tree)+i] = '\0';
      	  (*available) = 3*(int)T_MAX_NAME;
	  /* printf("\n. ++ 2 Available = %d",(*available)); */
      	}
#endif

      /* if(*available < (int)T_MAX_NAME/2) */
      /* 	{ */
      /* 	  (*s_tree) = (char *)mRealloc(*s_tree,*pos+(int)T_MAX_NAME,sizeof(char)); */
      /* 	  (*available) = (int)T_MAX_NAME; */
      /* 	} */
    }

  Free(format);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void R_wtree_Custom(t_node *pere, t_node *fils, int *available, char **s_tree, int *pos, t_tree *tree)
{
  int i,p,ori_len;
  char *format;

  format = (char *)mCalloc(100,sizeof(char));

  sprintf(format,"%%.%df",tree->bl_ndigits);
  /* strcpy(format,"%f"); */

  p = -1;
  if(fils->tax)
    {
/*       printf("\n- Writing on %p",*s_tree); */
      ori_len = *pos;

      if(OUTPUT_TREE_FORMAT == NEWICK)
	{
	  if(tree->write_tax_names == YES)
	    {
	      if(tree->io && tree->io->long_tax_names) 
		{
		  strcat(*s_tree,tree->io->long_tax_names[fils->num]);
		  (*pos) += (int)strlen(tree->io->long_tax_names[fils->num]);
		}
	      else
		{
		  strcat(*s_tree,fils->name);
		  (*pos) += (int)strlen(fils->name);
		}	  
	    }
	  else if(tree->write_tax_names == NO)
	    {
	      (*pos) += sprintf(*s_tree+*pos,"%d",fils->num);
	    }
	}
      else if(OUTPUT_TREE_FORMAT == NEXUS)
	{
	  (*pos) += sprintf(*s_tree+*pos,"%d",fils->num+1);
	}
      else
	{
	  PhyML_Printf("\n. Unknown tree format.");
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  PhyML_Printf("\n. s=%s\n",*s_tree);
	}

      if((fils->b) && (fils->b[0]) && (fils->b[0]->l->v > -1.))
	{
	  /* if(tree->print_labels) */
	  /*   { */
	  /*     if(fils->b[0]->n_labels < 10) */
	  /* 	For(i,fils->b[0]->n_labels)  */
	  /* 	  { */
	  /* 	    (*pos) += sprintf(*s_tree+*pos,"#%s",fils->b[0]->labels[i]); */
	  /* 	  } */
	  /*     else */
	  /* 	{ */
	  /* 	  (*pos) += sprintf(*s_tree+*pos,"#%d_labels",fils->b[0]->n_labels); */
	  /* 	} */
	  /*   } */

	  strcat(*s_tree,":");
	  (*pos)++;

#ifndef PHYTIME
	  if(!tree->n_root)
	    {
	      (*pos) += sprintf(*s_tree+*pos,format,fils->b[0]->l->v);
	    }
	  else
	    {
	      if(pere == tree->n_root)
		{
		  phydbl root_pos = (fils == tree->n_root->v[0])?(tree->n_root_pos):(1.-tree->n_root_pos);
		  (*pos) += sprintf(*s_tree+*pos,format,tree->e_root->l->v * root_pos);
		}
	      else
		{
		  (*pos) += sprintf(*s_tree+*pos,format,fils->b[0]->l->v);
		}
	    }		
#else
	  if(!tree->n_root)
	    {
	      (*pos) += sprintf(*s_tree+*pos,format,fils->b[0]->l->v);
	    }
	  else
	    {
	      (*pos) += sprintf(*s_tree+*pos,format,tree->rates->cur_l[fils->num]);
	    }
#endif

	  /* !!!!!!!!!!!!!!!!!!!!1 */
	  if(tree->print_labels)
	    {
	      if(fils->b[0]->n_labels < 10)
		For(i,fils->b[0]->n_labels) 
		  {
		    (*pos) += sprintf(*s_tree+*pos,"::%s",fils->b[0]->labels[i]);
		  }
	      else
		{
		  (*pos) += sprintf(*s_tree+*pos,"::%d_labels",fils->b[0]->n_labels);
		}
	    }

	}

      strcat(*s_tree,",");
      (*pos)++;

      (*available) = (*available) - (*pos - ori_len);

      if(*available < 0)
	{
	  PhyML_Printf("\n. s=%s\n",*s_tree);
	  PhyML_Printf("\n. len=%d\n",strlen(*s_tree));
	  PhyML_Printf("\n. The sequence names in your input file might be too long.");
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      if(*available < (int)T_MAX_NAME/2)
	{
	  (*s_tree) = (char *)mRealloc(*s_tree,*pos+(int)T_MAX_NAME,sizeof(char));
	  (*available) = (int)T_MAX_NAME;
	}
/*       printf(" %s [%d,%d]",*s_tree,(int)strlen(*s_tree),*available); */
    }
  else
    {

      (*s_tree)[(*pos)]='(';
      (*s_tree)[(*pos)+1]='\0';
      (*pos)++;
      (*available)--;

      if(*available < (int)T_MAX_NAME/2)
	{
	  (*s_tree) = (char *)mRealloc(*s_tree,*pos+(int)T_MAX_NAME,sizeof(char));
	  (*available) = (int)T_MAX_NAME;
	}

      if(tree->n_root)
	{
	  For(i,3)
	    {
	      if((fils->v[i] != pere) && (fils->b[i] != tree->e_root))
		R_wtree_Custom(fils,fils->v[i],available,s_tree,pos,tree);
	      else p=i;
	    }
	}
      else
	{
	  For(i,3)
	    {
	      if(fils->v[i] != pere)
		R_wtree_Custom(fils,fils->v[i],available,s_tree,pos,tree);
	      else p=i;
	    }
	}

      ori_len = *pos;
      
      if(p < 0)
	{
	  PhyML_Printf("\n. fils=%p root=%p root->v[0]=%p root->v[1]=%p",fils,tree->n_root,tree->n_root->v[0],tree->n_root->v[1]);
	  PhyML_Printf("\n. tree->e_root=%p fils->b[0]=%p fils->b[1]=%p fils->b[2]=%p",tree->e_root,fils->b[0],fils->b[1],fils->b[2]);		       
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

/*       printf("\n+ Writing on %p",*s_tree); */
      (*s_tree)[(*pos)-1] = ')';
      (*s_tree)[(*pos)]   = '\0';

      if((fils->b) && (fils->b[p]->l->v > -1.))
	{
	  if(tree->print_boot_val)
	    {
	      (*pos) += sprintf(*s_tree+*pos,"%d",fils->b[p]->bip_score);
	    }
	  else if(tree->print_alrt_val)
	    {
	      (*pos) += sprintf(*s_tree+*pos,"%f",fils->b[p]->ratio_test);
	    }
	  
	  fflush(NULL);

	  /* if((tree->print_labels) && (fils->b[p]->labels != NULL)) */
	  /*   { */
	  /*     if(fils->b[p]->n_labels < 10) */
	  /* 	For(i,fils->b[p]->n_labels)  */
	  /* 	  { */
	  /* 	    (*pos) += sprintf(*s_tree+*pos,"#%s",fils->b[p]->labels[i]); */
	  /* 	  } */
	  /*     else */
	  /* 	{ */
	  /* 	  (*pos) += sprintf(*s_tree+*pos,"#%d_labels",fils->b[p]->n_labels); */
	  /* 	} */
	  /*   } */

	  strcat(*s_tree,":");
	  (*pos)++;

#ifndef PHYTIME
	  if(!tree->n_root)
	    {
	      (*pos) += sprintf(*s_tree+*pos,format,fils->b[p]->l->v);
	    }
	  else
	    {
	      if(pere == tree->n_root)
		{
		  phydbl root_pos = (fils == tree->n_root->v[0])?(tree->n_root_pos):(1.-tree->n_root_pos);
		  (*pos) += sprintf(*s_tree+*pos,format,tree->e_root->l->v * root_pos);
		}
	      else
		{
		  (*pos) += sprintf(*s_tree+*pos,format,fils->b[p]->l->v);
		}
	    }
#else
	  if(!tree->n_root)
	    {
	      (*pos) += sprintf(*s_tree+*pos,format,fils->b[p]->l->v);
	    }
	  else
	    {
	      (*pos) += sprintf(*s_tree+*pos,format,tree->rates->cur_l[fils->num]);
	    }
#endif

	  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!1 */
	  if((tree->print_labels) && (fils->b[p]->labels != NULL))
	    {
	      if(fils->b[p]->n_labels < 10)
		For(i,fils->b[p]->n_labels) 
		  {
		    (*pos) += sprintf(*s_tree+*pos,"::%s",fils->b[p]->labels[i]);
		  }
	      else
		{
		  (*pos) += sprintf(*s_tree+*pos,"::%d_labels",fils->b[p]->n_labels);
		}
	    }

	}
      strcat(*s_tree,",");
      (*pos)++;
      (*available) = (*available) - (*pos - ori_len);

      if(*available < 0)
	{
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      if(*available < (int)T_MAX_NAME/2)
	{
	  (*s_tree) = (char *)mRealloc(*s_tree,*pos+(int)T_MAX_NAME,sizeof(char));
	  (*available) = (int)T_MAX_NAME;
	}
/*       printf(" %s [%d,%d]",*s_tree,(int)strlen(*s_tree),*available); */
    }

  Free(format);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Detect_Align_File_Format(option *io)
{
  int c;
  fpos_t curr_pos;
  
  fgetpos(io->fp_in_align,&curr_pos);
  
  errno = 0;

  while((c=fgetc(io->fp_in_align)) != EOF)
    {
      if(errno) io->data_file_format = PHYLIP;
      else if(c == '#')
	{
	  char s[10],t[6]="NEXUS";
	  if(!fgets(s,6,io->fp_in_align))
	    {     
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	  if(!strcmp(t,s)) 
	    {
	      fsetpos(io->fp_in_align,&curr_pos);
	      io->data_file_format = NEXUS;
	      return;
	    }
	}
    }
  
  fsetpos(io->fp_in_align,&curr_pos);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Detect_Tree_File_Format(option *io)
{
  int c;
  fpos_t curr_pos;
  
  fgetpos(io->fp_in_tree,&curr_pos);

  errno = 0;

  while((c=fgetc(io->fp_in_tree)) != EOF)
    {
      if(errno) 
	{
	  io->tree_file_format = PHYLIP;
	  PhyML_Printf("\n. Detected PHYLIP tree file format.");
	}
      else if(c == '#')
	{
	  char s[10],t[6]="NEXUS";
	  if(!fgets(s,6,io->fp_in_tree))
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	  if(!strcmp(t,s))
	    {
	      fsetpos(io->fp_in_tree,&curr_pos);
	      io->tree_file_format = NEXUS;
	      PhyML_Printf("\n. Detected NEXUS tree file format.");
	      return;
	    }
	}
    }
  
  fsetpos(io->fp_in_tree,&curr_pos);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

align **Get_Seq(option *io)
{
  io->data = NULL;


  if(!io->fp_in_align)
    {
      PhyML_Printf("\n== Filehandle to '%s' seems to be closed.",io->in_align_file);
      Exit("\n");
    }

  Detect_Align_File_Format(io);

  switch(io->data_file_format)
    {
    case PHYLIP: 
      {
	io->data = Get_Seq_Phylip(io);
	break;
      }
    case NEXUS:
      {
	io->nex_com_list = Make_Nexus_Com();
	Init_Nexus_Format(io->nex_com_list);
	Get_Nexus_Data(io->fp_in_align,io);
	Free_Nexus(io);
	break;
      }
    default:
      {
	PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	Warn_And_Exit("");
	break;
      }
    }

  if(!io->data)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  else
    {
      int i,j;
      char **buff;
      int *remove;
      int n_unkn,n_removed,pos;

      buff = (char **)mCalloc(io->n_otu,sizeof(char *));
      For(i,io->n_otu) buff[i] = (char *)mCalloc(io->data[0]->len,sizeof(char));
      remove = (int *)mCalloc(io->data[0]->len,sizeof(int));

      n_removed = 0;

      For(i,io->data[0]->len)
	{
	  For(j,io->n_otu)
	    {
	      if((io->data[j]->state[i] == '?') || (io->data[j]->state[i] == '-')) io->data[j]->state[i] = 'X';
	      if((io->datatype == NT) && (io->data[j]->state[i] == 'N')) io->data[j]->state[i] = 'X';
	      if(io->data[j]->state[i] == 'U') io->data[j]->state[i] = 'T';
	    }

	  n_unkn = 0;
	  For(j,io->n_otu) if(io->data[j]->state[i] == 'X') n_unkn++;

	  if(n_unkn == io->n_otu)
	    {
	      remove[i] = 1;
	      n_removed++;
	    }

	  For(j,io->n_otu) buff[j][i] = io->data[j]->state[i];
	}

      pos = 0;
      For(i,io->data[0]->len)
	{
/* 	  if(!remove[i]) */
/* 	    { */
	      For(j,io->n_otu) io->data[j]->state[pos] = buff[j][i];
	      pos++;
/* 	    } */
	}

      For(i,io->n_otu) io->data[i]->len = pos;
      For(i,io->n_otu) Free(buff[i]);
      Free(buff);
      Free(remove);
    }


  return io->data;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* align **Get_Seq_Nexus(option *io) */
/* { */
/*   char *s,*ori_s; */
/*   char *token; */
/*   int in_comment; */
/*   nexcom *curr_com; */
/*   nexparm *curr_parm; */
/*   int nxt_token_t,cur_token_t; */

/*   s = (char *)mCalloc(T_MAX_LINE,sizeof(char)); */
/*   token = (char *)mCalloc(T_MAX_TOKEN,sizeof(char)); */
      	  
/*   ori_s      = s; */
/*   in_comment = NO; */
/*   curr_com   = NULL; */
/*   curr_parm  = NULL; */
/*   nxt_token_t = NEXUS_COM;  */
/*   cur_token_t = -1;  */

/*   while(fgets(s,T_MAX_LINE,io->fp_in_align)) */
/*     {       */
/*       do */
/* 	{	   */
/* 	  Get_Token(&s,token);	   */

/* /\* 	  PhyML_Printf("\n. Token: '%s' next_token=%d cur_token=%d",token,nxt_token_t,cur_token_t); *\/ */

/* 	  if(token[0] == '\0') break; */

/* 	  if(token[0] == ';')  */
/* 	    { */
/* 	      curr_com   = NULL; */
/* 	      curr_parm  = NULL; */
/* 	      nxt_token_t = NEXUS_COM; */
/* 	      cur_token_t = -1; */
/* 	      break; /\* End of command *\/  */
/* 	    } */

/* 	  if(nxt_token_t == NEXUS_EQUAL)  */
/* 	    { */
/* 	      cur_token_t = NEXUS_VALUE; */
/* 	      nxt_token_t = NEXUS_PARM; */
/* 	      continue; */
/* 	    } */

/* 	  if((nxt_token_t == NEXUS_COM) && (cur_token_t != NEXUS_VALUE))  */
/* 	    { */
/* 	      Find_Nexus_Com(token,&curr_com,&curr_parm,io->nex_com_list); */
/* 	      if(curr_com)  */
/* 		{ */
/* 		  nxt_token_t = curr_com->nxt_token_t; */
/* 		  cur_token_t = curr_com->cur_token_t; */
/* 		} */
/* 	      if(cur_token_t != NEXUS_VALUE) continue; */
/* 	    } */

/* 	  if((nxt_token_t == NEXUS_PARM) && (cur_token_t != NEXUS_VALUE))  */
/* 	    { */
/* 	      Find_Nexus_Parm(token,&curr_parm,curr_com); */
/* 	      if(curr_parm)  */
/* 		{ */
/* 		  nxt_token_t = curr_parm->nxt_token_t; */
/* 		  cur_token_t = curr_parm->cur_token_t; */
/* 		} */
/* 	      if(cur_token_t != NEXUS_VALUE) continue; */
/* 	    } */

/* 	  if(cur_token_t == NEXUS_VALUE) */
/* 	    { */
/* 	      if((curr_parm->fp)(token,curr_parm,io))  /\* Read in parameter value *\/ */
/* 		{ */
/* 		  nxt_token_t = NEXUS_PARM; */
/* 		  cur_token_t = -1; */
/* 		} */
/* 	    } */
/* 	} */
/*       while(strlen(token) > 0); */
/*     } */

/*   Free(ori_s); */
/*   Free(token); */

/*   return io->data; */
/* } */

/* /\*********************************************************\/ */

void Get_Nexus_Data(FILE *fp, option *io)
{
  char *token;
  nexcom *curr_com;
  nexparm *curr_parm;
  int nxt_token_t,cur_token_t;

  token = (char *)mCalloc(T_MAX_TOKEN,sizeof(char));
      	  
  curr_com   = NULL;
  curr_parm  = NULL;
  nxt_token_t = NEXUS_COM; 
  cur_token_t = -1; 

  do
    {
      if(!Get_Token(fp,token)) break;

/*       PhyML_Printf("\n+ Token: '%s' next_token=%d cur_token=%d",token,nxt_token_t,cur_token_t); */

      if(token[0] == ';') 
	{
	  curr_com    = NULL;
	  curr_parm   = NULL;
	  nxt_token_t = NEXUS_COM;
	  cur_token_t = -1;
	}
      
      if(nxt_token_t == NEXUS_EQUAL) 
	{
	  cur_token_t = NEXUS_VALUE;
	  nxt_token_t = NEXUS_PARM;
	  continue;
	}
      
      if((nxt_token_t == NEXUS_COM) && (cur_token_t != NEXUS_VALUE)) 
	{
	  Find_Nexus_Com(token,&curr_com,&curr_parm,io->nex_com_list);
	  if(curr_com) 
	    {
	      nxt_token_t = curr_com->nxt_token_t;
	      cur_token_t = curr_com->cur_token_t;
	    }
	  if(cur_token_t != NEXUS_VALUE) continue;
	}
      
      if((nxt_token_t == NEXUS_PARM) && (cur_token_t != NEXUS_VALUE)) 
	{
	  Find_Nexus_Parm(token,&curr_parm,curr_com);
	  if(curr_parm) 
	    {
	      nxt_token_t = curr_parm->nxt_token_t;
	      cur_token_t = curr_parm->cur_token_t;
	    }
	  if(cur_token_t != NEXUS_VALUE) continue;
	}
      
      if(cur_token_t == NEXUS_VALUE)
	{
	  if((curr_parm->fp)(token,curr_parm,io))  /* Read in parameter value */
	    {
	      nxt_token_t = NEXUS_PARM;
	      cur_token_t = -1;
	    }
	}
    }
  while(strlen(token) > 0);
  
  Free(token);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Get_Token(FILE *fp, char *token)
{
  char c;
  
  c = ' ';
  while(c == ' ' || c == '\t' || c == '\n') 
    {
      c = fgetc(fp);
      if(c == EOF) return 0;
    }

  if(c == '"')
    {
      do
	{
	  *token = c;
	  token++;
	  c = fgetc(fp);
	  if(c == EOF) return 0;
	}
      while(c != '"');
      *token = c;
      c = fgetc(fp);
      if(c == EOF) return 0;
      *(token+1) = '\0';
      return 1;
    }

  if(c == '[')
    {
      Skip_Comment(fp);
      c = fgetc(fp);
      if(c == EOF) return 0;
      return 1;
    }

  if(c == '#')      { *token = c; token++; }
  else if(c == ';') { *token = c; token++; }
  else if(c == ',') { *token = c; token++; }
  else if(c == '.') { *token = c; token++; }
  else if(c == '=') { *token = c; token++; }
  else if(c == '(') { *token = c; token++; }
  else if(c == ')') { *token = c; token++; }
  else if(c == '{') { *token = c; token++; }
  else if(c == '}') { *token = c; token++; }
  else if(c == '?') { *token = c; token++; }
  else if(c == '-') { *token = c; token++; }
  else
    {
      while(isgraph(c) && c != ';' && c != '-' && c != ',')
	{
	  *(token++) = c;
	  c = fgetc(fp);
	  if(c == EOF) return 0;
	}

      fseek(fp,-1*sizeof(char),SEEK_CUR);

    }
  *token = '\0';
  return 1;
}


/* void Get_Token(char *line, char *token) */
/* { */
/*   while(**line == ' ' || **line == '\t') (*line)++; */


/*   if(**line == '"')  */
/*     { */
/*       do { *token = **line; (*line)++; token++; } while(**line != '"'); */
/*       *token = **line; */
/*       (*line)++; */
/*       *(token+1) = '\0'; */
/*       return; */
/*     } */

/*   if(**line == '[')  */
/*     { */
/*       int in_comment; */

/*       in_comment = 1; */
/*       do  */
/* 	{  */
/* 	  (*line)++;  */
/* 	  if(**line == '[')  */
/* 	    { */
/* 	      in_comment++; */
/* 	    } */
/* 	  else if(**line == ']') in_comment--;	   */
/* 	} */
/*       while(in_comment); */
/*       (*line)++; */
/*       return; */
/*     } */


/*   if(**line == '#')      {*token = **line; (*line)++; token++; } */
/*   else if(**line == ';') {*token = **line; (*line)++; token++; } */
/*   else if(**line == ',') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '.') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '=') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '(') {*token = **line; (*line)++; token++; } */
/*   else if(**line == ')') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '{') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '}') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '?') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '-') {*token = **line; (*line)++; token++; } */
/*   else */
/*     { */
/*       while(isgraph(**line) && **line != ';' && **line != '=' && **line != ',')  */
/* 	{ */
/* 	  *(token++) = **line; */
/* 	  (*line)++;  */
/* 	} */
/*     } */
/*   *token = '\0'; */
/* } */

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


align **Get_Seq_Phylip(option *io)
{

  Read_Ntax_Len_Phylip(io->fp_in_align,&io->n_otu,&io->init_len);

  if(io->n_otu > N_MAX_OTU)
    {
      PhyML_Printf("\n. The number of taxa should not exceed %d",N_MAX_OTU);
      Exit("\n");
    }

  if(io->interleaved) io->data = Read_Seq_Interleaved(io);
  else                io->data = Read_Seq_Sequential(io);

  return io->data;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Read_Ntax_Len_Phylip(FILE *fp ,int *n_otu, int *n_tax)
{
  char *line;
  int readok;
  
  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));  

  readok = 0;
  do
    {
      if(fscanf(fp,"%s",line) == EOF)
	{
	  Free(line); 
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      else
	{
	  if(strcmp(line,"\n") && strcmp(line,"\r") && strcmp(line,"\t"))
	    {
	      sscanf(line,"%d",n_otu);
	      if(*n_otu <= 0) Warn_And_Exit("\n. The number of taxa cannot be negative.\n");

	      if(!fscanf(fp,"%s",line)) Exit("\n");
	      sscanf(line,"%d",n_tax);
	      if(*n_tax <= 0) Warn_And_Exit("\n. The sequence length cannot be negative.\n");
	      else readok = 1;
	    }
	}
    }while(!readok);
  
  Free(line);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


align **Read_Seq_Sequential(option *io)
{
  int i;
  char *line;
  align **data;
/*   char c; */
  char *format;

  format = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  line   = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  data   = (align **)mCalloc(io->n_otu,sizeof(align *));

/*   while((c=fgetc(in))!='\n'); */
 /*  while(((c=fgetc(io->fp_in_align))!='\n') && (c != ' ') && (c != '\r') && (c != '\t')); */

  For(i,io->n_otu)
    {
      data[i]        = (align *)mCalloc(1,sizeof(align));
      data[i]->name  = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      data[i]->state = (char *)mCalloc(io->init_len*io->state_len+1,sizeof(char));

      data[i]->is_ambigu = NULL;
      data[i]->len = 0;

      sprintf(format, "%%%ds", T_MAX_NAME);

      if(!fscanf(io->fp_in_align,format,data[i]->name)) Exit("\n");

      Check_Sequence_Name(data[i]->name);

      while(data[i]->len < io->init_len * io->state_len) Read_One_Line_Seq(&data,i,io->fp_in_align);

      if(data[i]->len != io->init_len * io->state_len)
	{
	  PhyML_Printf("\n. Err: Problem with species %s's sequence (check the format).\n",data[i]->name);
	  PhyML_Printf("\n. Observed sequence length: %d, expected length: %d\n",data[i]->len, io->init_len * io->state_len);
	  Warn_And_Exit("");
	}
    }

  For(i,io->n_otu) data[i]->state[data[i]->len] = '\0';
  
  Restrict_To_Coding_Position(data,io);

  Free(format);
  Free(line);
  
  return data;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


align **Read_Seq_Interleaved(option *io)
{
  int i,end,num_block;
  char *line;
  align **data;
/*   char c; */
  char *format;



  line   = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  format = (char *)mCalloc(T_MAX_NAME, sizeof(char));
  data   = (align **)mCalloc(io->n_otu,sizeof(align *));


/*   while(((c=fgetc(io->fp_in_align))!='\n') && (c != ' ') && (c != '\r') && (c != '\t')); */


  end = 0;
  For(i,io->n_otu)
    {
      data[i]        = (align *)mCalloc(1,sizeof(align));
      data[i]->name  = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      data[i]->state = (char *)mCalloc(io->init_len*io->state_len+1,sizeof(char));

      data[i]->len       = 0;
      data[i]->is_ambigu = NULL;

      sprintf(format, "%%%ds", T_MAX_NAME);
/*       sprintf(format, "%%%ds", 10); */
      
      if(!fscanf(io->fp_in_align,format,data[i]->name)) Exit("\n");

      Check_Sequence_Name(data[i]->name);

      if(!Read_One_Line_Seq(&data,i,io->fp_in_align))
	{
	  end = 1;
	  if((i != io->n_otu) && (i != io->n_otu-1))
	    {
	      PhyML_Printf("\n. Err: Problem with species %s's sequence.\n",data[i]->name);
	      PhyML_Printf("\n. Observed sequence length: %d, expected length: %d\n",data[i]->len, io->init_len * io->state_len);
	      Warn_And_Exit("");
	    }
	  break;
	}
    }
  
  if(data[0]->len == io->init_len * io->state_len) end = 1;

/*   if(end) printf("\n. finished yet '%c'\n",fgetc(io->fp_in_align)); */
  if(!end)
    {

      end = 0;

      num_block = 1;
      do
	{
	  num_block++;

	  /* interblock */
	  if(!fgets(line,T_MAX_LINE,io->fp_in_align)) break;

	  if(line[0] != 13 && line[0] != 10)
	    {
	      PhyML_Printf("\n. One or more missing sequences in block %d.\n",num_block-1);
	      Warn_And_Exit("");
	    }

	  For(i,io->n_otu) if(data[i]->len != io->init_len * io->state_len) break;

	  if(i == io->n_otu) break;

	  For(i,io->n_otu)
	    {
	      if(data[i]->len > io->init_len * io->state_len)
		{
		  PhyML_Printf("\n. Observed length=%d expected length=%d.\n",data[i]->len,io->init_len * io->state_len);
		  PhyML_Printf("\n. Err: Problem with species %s's sequence.\n",data[i]->name);
		  Warn_And_Exit("");
		}
	      else if(!Read_One_Line_Seq(&data,i,io->fp_in_align))
		{
		  end = 1;
		  if((i != io->n_otu) && (i != io->n_otu-1))
		    {
		      PhyML_Printf("\n. Err: Problem with species %s's sequence.\n",data[i]->name);
		      PhyML_Printf("\n. Observed sequence length: %d, expected length: %d.\n",data[i]->len, io->init_len * io->state_len);
		      Warn_And_Exit("");
		    }
		  break;
		}
	    }
	}while(!end);
    }

  For(i,io->n_otu) data[i]->state[data[i]->len] = '\0';

  For(i,io->n_otu)
    {
      if(data[i]->len != io->init_len * io->state_len)
	{
	  PhyML_Printf("\n. Check sequence '%s' length (expected length: %d, observed length: %d) [OTU %d].\n",data[i]->name,io->init_len,data[i]->len,i+1);
	  Warn_And_Exit("");
	}
    }

  Restrict_To_Coding_Position(data,io);

  Free(format);
  Free(line);

  return data;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int Read_One_Line_Seq(align ***data, int num_otu, FILE *in)
{
  char c = ' ';
  int nchar = 0;
  
  while(1)
    {
/*       if((c == EOF) || (c == '\n') || (c == '\r')) break; */
      
      if((c == 13) || (c == 10)) 
	{
/* 	  PhyML_Printf("[%d %d]\n",c,nchar); fflush(NULL); */
	  if(!nchar)
	    {
	      c=(char)fgetc(in);
	      continue;
	    }
	  else 
	    { 
/* 	      PhyML_Printf("break\n");  */
	      break; 
	    }
	}
      else if(c == EOF)
	{
/* 	  PhyML_Printf("EOL\n"); */
	  break;
	}
      else if((c == ' ') || (c == '\t') || (c == 32)) 
	{
/* 	  PhyML_Printf("[%d]",c); */
	  c=(char)fgetc(in); 
	  continue;
	}

      nchar++;
      Uppercase(&c);
      
      if(c == '.')
	{
	  c = (*data)[0]->state[(*data)[num_otu]->len];
	  if(!num_otu)
	    Warn_And_Exit("\n. Err: Symbol \".\" should not appear in the first sequence\n");
	}
      (*data)[num_otu]->state[(*data)[num_otu]->len]=c;
      (*data)[num_otu]->len++;
/*       PhyML_Printf("%c",c); */
      c = (char)fgetc(in);
      if(c == ';') break;
    }

  if(c == EOF) return 0;
  else return 1;  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_tree *Read_Tree_File(option *io)
{
  t_tree *tree;

  if(!io->fp_in_tree)
    {
      PhyML_Printf("\n== Filehandle to '%s' seems to be closed.",io->in_tree_file);
      Exit("\n");
    }


  Detect_Tree_File_Format(io);

  io->treelist->list_size = 0;

  switch(io->tree_file_format)
    {
    case PHYLIP: 
      {
	do
	  {
	    io->treelist->tree = (t_tree **)realloc(io->treelist->tree,(io->treelist->list_size+1)*sizeof(t_tree *));
	    io->tree = Read_Tree_File_Phylip(io->fp_in_tree);
	    if(!io->tree) break;
	    if(io->treelist->list_size > 1) PhyML_Printf("\n. Reading tree %d",io->treelist->list_size+1);
	    io->treelist->tree[io->treelist->list_size] = io->tree;
	    io->treelist->list_size++;
	  }while(io->tree);
	break;
      }
    case NEXUS:
      {
	io->nex_com_list = Make_Nexus_Com();
	Init_Nexus_Format(io->nex_com_list);
	Get_Nexus_Data(io->fp_in_tree,io);
	Free_Nexus(io);
	break;
      }
    default:
      {
	PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	Warn_And_Exit("");
	break;
      }
    }
  
  if(!io->long_tax_names)
    {
      int i;

      tree = io->treelist->tree[0];

      io->long_tax_names  = (char **)mCalloc(tree->n_otu,sizeof(char *));
      io->short_tax_names = (char **)mCalloc(tree->n_otu,sizeof(char *));

      For(i,tree->n_otu)
	{
	  io->long_tax_names[i] = (char *)mCalloc(strlen(tree->t_nodes[i]->name)+1,sizeof(char));
	  io->short_tax_names[i] = (char *)mCalloc(strlen(tree->t_nodes[i]->name)+1,sizeof(char));
	  strcpy(io->long_tax_names[i],tree->t_nodes[i]->name);
	  strcpy(io->short_tax_names[i],tree->t_nodes[i]->name);
	}
    }
  return NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


char *Return_Tree_String_Phylip(FILE *fp_input_tree)
{
  char *line;
  int i;
  char c;
  int open,maxopen;

  if(fp_input_tree == NULL)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");      
    }
 
  do
    {
      c=fgetc(fp_input_tree);
    }
  while((c != '(') && (c != EOF));
  
  
  if(c==EOF) return NULL;

  line = (char *)mCalloc(1,sizeof(char));
  open = 1;
  maxopen = open;

  i=0;
  for(;;)
    {
      if((c == ' ') || (c == '\n'))
	{
	  c=fgetc(fp_input_tree);
	  if(c == EOF || c == ';') break;
	  else continue;
	}
      
      if(c == '[')
	{
	  Skip_Comment(fp_input_tree);
	  c = fgetc(fp_input_tree);
	  if(c == EOF || c == ';') break;
	}

      line = (char *)mRealloc(line,i+2,sizeof(char));

      line[i]=c;
      i++;
      c=fgetc(fp_input_tree);
      if(c==EOF || c==';') break;
      if(c=='(') open++;
      if(c==')') open--;
      if(open>maxopen) maxopen = open;
    }
  line[i] = '\0';
  

  /* if(maxopen == 1) return NULL; */
  return line;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


t_tree *Read_Tree_File_Phylip(FILE *fp_input_tree)
{
  char *line;
  t_tree *tree;

  line = Return_Tree_String_Phylip(fp_input_tree);
  tree = Read_Tree(&line);
  Free(line);
  
  return tree;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Site(calign *cdata, int num, int n_otu, char *sep, int stepsize)
{
  int i,j;
  For(i,n_otu)
    {
      PhyML_Printf("%s   ",cdata->c_seq[i]->name);
      For(j,stepsize)
	PhyML_Printf("%c",cdata->c_seq[i]->state[num+j]);
      PhyML_Printf("%s",sep);
    }
  PhyML_Fprintf(stderr,"%s",sep);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Site_Lk(t_tree *tree, FILE *fp)
{
  int site;
  int catg;
  char *s;
  phydbl postmean;

  if(!tree->io->print_site_lnl)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(!tree->io->print_trace)
    {
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      
      PhyML_Fprintf(fp,"Note : P(D|M) is the probability of site D given the model M (i.e., the site likelihood)\n");
      if(tree->mod->ras->n_catg > 1 || tree->mod->ras->invar)
	PhyML_Fprintf(fp,"P(D|M,rr[x]) is the probability of site D given the model M and the relative rate\nof evolution rr[x], where x is the class of rate to be considered.\nWe have P(D|M) = \\sum_x P(x) x P(D|M,rr[x]).\n");
      PhyML_Fprintf(fp,"\n\n");
      
      sprintf(s,"Site");
      PhyML_Fprintf(fp, "%-7s",s);
      
      sprintf(s,"P(D|M)");
      PhyML_Fprintf(fp,"%-16s",s);
      
      if(tree->mod->ras->n_catg > 1)
	{
	  For(catg,tree->mod->ras->n_catg)
	    {
	      sprintf(s,"P(D|M,rr[%d]=%5.4f)",catg+1,tree->mod->ras->gamma_rr->v[catg]);
	      PhyML_Fprintf(fp,"%-22s",s);
	    }
	  
	  sprintf(s,"Posterior mean");
	  PhyML_Fprintf(fp,"%-22s",s);
	}
      
      
      if(tree->mod->ras->invar)
	{
	  sprintf(s,"P(D|M,rr[0]=0)");
	  PhyML_Fprintf(fp,"%-16s",s);
	}
      PhyML_Fprintf(fp,"\n");
      
      For(site,tree->data->init_len)
	{
	  PhyML_Fprintf(fp,"%-7d",site+1);
	  PhyML_Fprintf(fp,"%-16g",tree->cur_site_lk[tree->data->sitepatt[site]]);      
	  if(tree->mod->ras->n_catg > 1)
	    {
	      For(catg,tree->mod->ras->n_catg)
		PhyML_Fprintf(fp,"%-22g",(phydbl)EXP(tree->log_site_lk_cat[catg][tree->data->sitepatt[site]]));

	      postmean = .0;
	      For(catg,tree->mod->ras->n_catg) 
		postmean += 
		tree->mod->ras->gamma_rr->v[catg] * 
		EXP(tree->log_site_lk_cat[catg][tree->data->sitepatt[site]]) * 
		tree->mod->ras->gamma_r_proba->v[catg];
	      postmean /= tree->cur_site_lk[tree->data->sitepatt[site]];

	      PhyML_Fprintf(fp,"%-22g",postmean);
	    }
	  if(tree->mod->ras->invar)
	    {
	      if((phydbl)tree->data->invar[tree->data->sitepatt[site]] > -0.5)
		PhyML_Fprintf(fp,"%-16g",tree->mod->e_frq->pi->v[tree->data->invar[tree->data->sitepatt[site]]]);
	      else
		PhyML_Fprintf(fp,"%-16g",0.0);
	    }
	  PhyML_Fprintf(fp,"\n");
	}
      Free(s);
    }
  else
    {
      For(site,tree->data->init_len)
	PhyML_Fprintf(fp,"%.2f\t",LOG(tree->cur_site_lk[tree->data->sitepatt[site]]));
      PhyML_Fprintf(fp,"\n");
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Print_Seq(align **data, int n_otu)
{
  int i,j;

  PhyML_Printf("%d\t%d\n",n_otu,data[0]->len);
  For(i,n_otu)
    {
      For(j,20)
	{
	  if(j<(int)strlen(data[i]->name))
	     putchar(data[i]->name[j]);
	  else putchar(' ');
	}
/*       PhyML_Printf("%10d  ",i); */
      For(j,data[i]->len)
	{
	  PhyML_Printf("%c",data[i]->state[j]);
	}
      PhyML_Printf("\n");
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Print_CSeq(FILE *fp, int compressed, calign *cdata)
{
  int i,j;
  int n_otu;
  
  n_otu = cdata->n_otu;
  if(cdata->format == 0)
    {
      PhyML_Fprintf(fp,"%d\t%d\n",n_otu,cdata->init_len);
    }
  else
    {
      PhyML_Fprintf(fp,"#NEXUS\n");
      PhyML_Fprintf(fp,"begin data\n");
      PhyML_Fprintf(fp,"dimensions ntax=%d nchar=%d;\n",n_otu,cdata->init_len);
      PhyML_Fprintf(fp,"format sequential datatype=dna;\n");
      PhyML_Fprintf(fp,"matrix\n");
    }
  
  For(i,n_otu)
    {
      For(j,50)
	{
	  if(j<(int)strlen(cdata->c_seq[i]->name))
	    fputc(cdata->c_seq[i]->name[j],fp);
	  else fputc(' ',fp);
	}
      
      if(compressed == YES) /* Print out compressed sequences */
	PhyML_Fprintf(fp,"%s",cdata->c_seq[i]->state);
      else /* Print out uncompressed sequences */
	{
	  For(j,cdata->init_len)
	    {
	      PhyML_Fprintf(fp,"%c",cdata->c_seq[i]->state[cdata->sitepatt[j]]);
	    }
	}
      PhyML_Fprintf(fp,"\n");
    }
  PhyML_Fprintf(fp,"\n");

  if(cdata->format == 1)
    {
      PhyML_Fprintf(fp,";\n");
      PhyML_Fprintf(fp,"END;\n");
    }


/*   PhyML_Printf("\t"); */
/*   For(j,cdata->crunch_len) */
/*     PhyML_Printf("%.0f ",cdata->wght[j]); */
/*   PhyML_Printf("\n"); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Print_CSeq_Select(FILE *fp, int compressed, calign *cdata, t_tree *tree)
{
  int i,j;
  int n_otu;
  phydbl eps;

  int slice = 14;

  eps = 1.E-6;
  n_otu = 0;
  For(i,cdata->n_otu)
    if(tree->rates->nd_t[i] < tree->rates->time_slice_lims[slice] + eps)
      n_otu++;
  
  PhyML_Fprintf(fp,"%d\t%d\n",n_otu,cdata->init_len);

  n_otu = cdata->n_otu;

  For(i,n_otu)
    {
      if(tree->rates->nd_t[i] < tree->rates->time_slice_lims[slice] + eps)
	{
	  For(j,50)
	    {
	      if(j<(int)strlen(cdata->c_seq[i]->name))
		fputc(cdata->c_seq[i]->name[j],fp);
	      else fputc(' ',fp);
	    }
	  
	  if(compressed == YES) /* Print out compressed sequences */
	    PhyML_Fprintf(fp,"%s",cdata->c_seq[i]->state);
	  else /* Print out uncompressed sequences */
	    {
	      For(j,cdata->init_len)
		{
		  PhyML_Fprintf(fp,"%c",cdata->c_seq[i]->state[cdata->sitepatt[j]]);
		}
	    }
	  PhyML_Fprintf(fp,"\n");
	}
    }

  if(cdata->format == 1)
    {
      PhyML_Fprintf(fp,";\n");
      PhyML_Fprintf(fp,"END;\n");
    }


/*   PhyML_Printf("\t"); */
/*   For(j,cdata->crunch_len) */
/*     PhyML_Printf("%.0f ",cdata->wght[j]); */
/*   PhyML_Printf("\n"); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Dist(matrix *mat)
{
  int i,j;

  For(i,mat->n_otu)
    {
      PhyML_Printf("%s ",mat->name[i]);

      For(j,mat->n_otu)
	PhyML_Printf("%9.6f ",mat->dist[i][j]);
      PhyML_Printf("\n");
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Print_Node(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  int dir;
  dir = -1;
  For(i,3) if(a->v[i] == d) {dir = i; break;}
  PhyML_Printf("Node nums = %3d %3d  (dir=%d);",a->num,d->num,dir);
  PhyML_Printf("Node names = '%s' '%s' ; ",a->name,d->name);
  For(i,3) if(a->v[i] == d)
    {
      if(a->b[i])
	{
	  PhyML_Printf("Branch num = %3d%c (%d %d) %f",
		       a->b[i]->num,a->b[i]==tree->e_root?'*':' ',a->b[i]->left->num,
		       a->b[i]->rght->num,a->b[i]->l->v);
	  if(a->b[i]->left->tax) PhyML_Printf(" WARNING LEFT->TAX!");
	  break;
	}
    }
  PhyML_Printf("\n");

  if(d->tax) return;
  else
    For(i,3)
      if(d->v[i] != a) Print_Node(d,d->v[i],tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Model(t_mod *mod)
{
  int i,j,k;

  PhyML_Printf("\n. name=%s",mod->modelname);
  PhyML_Printf("\n. string=%s",mod->custom_mod_string);
  PhyML_Printf("\n. mod_num=%d",mod->mod_num);
  PhyML_Printf("\n. ns=%d",mod->ns);
  PhyML_Printf("\n. n_catg=%d",mod->ras->n_catg);
  PhyML_Printf("\n. kappa=%f",mod->kappa->v);
  PhyML_Printf("\n. alpha=%f",mod->ras->alpha->v);
  PhyML_Printf("\n. lambda=%f",mod->lambda->v);
  PhyML_Printf("\n. pinvar=%f",mod->pinvar->v);
  PhyML_Printf("\n. br_len_multiplier=%f",mod->br_len_multiplier->v);
  PhyML_Printf("\n. whichmodel=%d",mod->whichmodel);
  PhyML_Printf("\n. update_eigen=%d",mod->update_eigen);
  PhyML_Printf("\n. bootstrap=%d",mod->bootstrap);
  PhyML_Printf("\n. n_diff_rr=%d",mod->r_mat->n_diff_rr);
  PhyML_Printf("\n. invar=%d",mod->ras->invar);
  PhyML_Printf("\n. use_m4mod=%d",mod->use_m4mod);
  PhyML_Printf("\n. gamma_median=%d",mod->ras->gamma_median);
  PhyML_Printf("\n. state_len=%d",mod->io->state_len);
  PhyML_Printf("\n. log_l=%d",mod->log_l);
  PhyML_Printf("\n. l_min=%f",mod->l_min);
  PhyML_Printf("\n. l_max=%f",mod->l_max);
  PhyML_Printf("\n. free_mixt_rates=%d",mod->ras->free_mixt_rates);
  PhyML_Printf("\n. gamma_mgf_bl=%d",mod->gamma_mgf_bl);
  
  PhyML_Printf("\n. Pi\n");
  For(i,mod->ns) PhyML_Printf(" %f ",mod->e_frq->pi->v[i]);
  PhyML_Printf("\n");
  For(i,mod->ns) PhyML_Printf(" %f ",mod->e_frq->pi_unscaled->v[i]);
  
  PhyML_Printf("\n. Rates\n");
  For(i,mod->ras->n_catg) PhyML_Printf(" %f ",mod->ras->gamma_r_proba->v[i]);
  PhyML_Printf("\n");
  For(i,mod->ras->n_catg) PhyML_Printf(" %f ",mod->ras->gamma_r_proba_unscaled->v[i]);
  PhyML_Printf("\n");
  For(i,mod->ras->n_catg) PhyML_Printf(" %f ",mod->ras->gamma_rr->v[i]);
  PhyML_Printf("\n");
  For(i,mod->ras->n_catg) PhyML_Printf(" %f ",mod->ras->gamma_rr_unscaled->v[i]);
  
  
  
  PhyML_Printf("\n. Qmat \n");
  if(mod->whichmodel == CUSTOM)
    {
      fflush(NULL);
      For(i,6) {PhyML_Printf(" %12f ",mod->r_mat->rr->v[i]); fflush(NULL);}
      For(i,6) {PhyML_Printf(" %12f ",mod->r_mat->rr_val->v[i]); fflush(NULL);}
      For(i,6) {PhyML_Printf(" %12d ",mod->r_mat->rr_num->v[i]); fflush(NULL);}
      For(i,6) {PhyML_Printf(" %12d ",mod->r_mat->n_rr_per_cat->v[i]); fflush(NULL);}
    }
  For(i,mod->ns)
    {
      PhyML_Printf("  ");
      For(j,4)
	PhyML_Printf("%8.5f  ",mod->r_mat->qmat->v[i*4+j]);
      PhyML_Printf("\n");
    }

  PhyML_Printf("\n. Freqs");
  PhyML_Printf("\n");
  For(i,mod->ns) PhyML_Printf(" %12f ",mod->user_b_freq->v[i]);
  PhyML_Printf("\n");
  For(i,mod->ns) PhyML_Printf(" %12f ",mod->e_frq->pi->v[i]);
  PhyML_Printf("\n");
  For(i,mod->ns) PhyML_Printf(" %12f ",mod->e_frq->pi_unscaled->v[i]);

  PhyML_Printf("\n. Eigen\n");
  For(i,2*mod->ns)       PhyML_Printf(" %f ",mod->eigen->space[i]);      
  PhyML_Printf("\n");
  For(i,2*mod->ns)       PhyML_Printf(" %f ",mod->eigen->space_int[i]); 
  PhyML_Printf("\n");
  For(i,mod->ns)         PhyML_Printf(" %f ",mod->eigen->e_val[i]);      
  PhyML_Printf("\n");
  For(i,mod->ns)         PhyML_Printf(" %f ",mod->eigen->e_val_im[i]);   
  PhyML_Printf("\n");
  For(i,mod->ns*mod->ns) PhyML_Printf(" %f ",mod->eigen->r_e_vect[i]);   
  PhyML_Printf("\n");
  For(i,mod->ns*mod->ns) PhyML_Printf(" %f ",mod->eigen->r_e_vect[i]);   
  PhyML_Printf("\n");
  For(i,mod->ns*mod->ns) PhyML_Printf(" %f ",mod->eigen->r_e_vect_im[i]);
  PhyML_Printf("\n");
  For(i,mod->ns*mod->ns) PhyML_Printf(" %f ",mod->eigen->l_e_vect[i]);   
  PhyML_Printf("\n");
  For(i,mod->ns*mod->ns) PhyML_Printf(" %f ",mod->eigen->q[i]);          
  PhyML_Printf("\n");
  
  PhyML_Printf("\n. Pij");
  For(k,mod->ras->n_catg)
    {
      PMat(0.01*mod->ras->gamma_rr->v[k],mod,mod->ns*mod->ns*k,mod->Pij_rr->v);
      PhyML_Printf("\n. l=%f\n",0.01*mod->ras->gamma_rr->v[k]);
      For(i,mod->ns)
	{
	  PhyML_Printf("  ");
	  For(j,mod->ns)
	    PhyML_Printf("%8.5f  ",mod->Pij_rr->v[k*mod->ns*mod->ns+i*mod->ns+j]);
	  PhyML_Printf("\n");
	}
    }
  

  PhyML_Printf("\n");

  fflush(NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Mat(matrix *mat)
{
  int i,j;

  PhyML_Printf("%d",mat->n_otu);
  PhyML_Printf("\n");

  For(i,mat->n_otu)
    {
      For(j,13)
	{
	  if(j>=(int)strlen(mat->name[i])) putchar(' ');
	  else putchar(mat->name[i][j]);
	}

      For(j,mat->n_otu)
	{
	  char s[2]="-";
	  if(mat->dist[i][j] < .0)
	    PhyML_Printf("%12s",s);
	  else
	    PhyML_Printf("%12f",mat->dist[i][j]);
	}
      PhyML_Printf("\n");
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

FILE *Openfile(char *filename, int mode)
{
  /* mode = 0 -> read */
  /* mode = 1 -> write */
  /* mode = 2 -> append */

  FILE *fp;
  char *s;
  int open_test=0;

/*   s = (char *)mCalloc(T_MAX_FILE,sizeof(char)); */

/*   strcpy(s,filename); */

  s = filename;

  fp = NULL;

  switch(mode)
    {
    case 0 :
      {
	while(!(fp = (FILE *)fopen(s,"r")) && ++open_test<10)
	  {
	    PhyML_Printf("\n. Can't open file '%s', enter a new name : ",s);
	    Getstring_Stdin(s);
	  }
	break;
      }
    case 1 :
      {
	fp = (FILE *)fopen(s,"w");
	break;
      }
    case 2 :
      {
	fp = (FILE *)fopen(s,"a");
	break;
      }

    default : break;

    }

/*   Free(s); */

  return fp;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Print_Fp_Out(FILE *fp_out, time_t t_beg, time_t t_end, t_tree *tree, option *io, int n_data_set, int num_tree)
{
  char *s;
  div_t hour,min;
  int i;

/*   For(i,2*tree->n_otu-3) fprintf(fp_out,"\n. * Edge %3d: %f",i,tree->t_edges[i]->l); */
  
  if((!n_data_set) || (!num_tree))
    {
      rewind(fp_out);
      Print_Banner_Small(fp_out);
    }

  PhyML_Fprintf(fp_out,"\n. Sequence filename: \t\t\t%s", Basename(io->in_align_file));
  PhyML_Fprintf(fp_out,"\n. Data set: \t\t\t\t#%d",n_data_set);

  if(io->mod->s_opt->random_input_tree)
    PhyML_Fprintf(fp_out,"\n. Random init tree: \t\t\t#%d",num_tree+1);
  else if(io->n_trees > 1)
    PhyML_Fprintf(fp_out,"\n. Starting tree number: \t\t\t#%d",num_tree+1);
  
  if(io->mod->s_opt->opt_topo)
    {
      if(io->mod->s_opt->topo_search == NNI_MOVE) PhyML_Fprintf(fp_out,"\n. Tree topology search : \t\tNNIs");
      else if(io->mod->s_opt->topo_search == SPR_MOVE) PhyML_Fprintf(fp_out,"\n. Tree topology search : \t\tSPRs");
      else if(io->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR) PhyML_Fprintf(fp_out,"\n. Tree topology search : \t\tBest of NNIs and SPRs");
    }
  else
    {
      PhyML_Fprintf(fp_out,"\n. Tree topology: \t\t\tfixed");
    }


  /* was after Sequence file ; moved here FLT */
  s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  if(io->in_tree == 2)
    {
      strcat(strcat(strcat(s,"user tree ("),io->in_tree_file),")");
    }
  else
    {
      if(!io->mod->s_opt->random_input_tree)
	{
	  if(io->in_tree == 0)
	    strcat(s,"BioNJ");
	  if(io->in_tree == 1)
	    strcat(s,"parsimony");
	}
      else
	{
	  strcat(s,"random tree");
	}
    }

  PhyML_Fprintf(fp_out,"\n. Initial tree: \t\t\t%s",s);
  Free(s);

  if(tree->io->datatype == NT)
    {
      PhyML_Fprintf(fp_out,"\n. Model of nucleotides substitution: \t%s",io->mod->modelname);
      if(io->mod->whichmodel == CUSTOM)
      PhyML_Fprintf(fp_out," (%s)",io->mod->custom_mod_string);
    }
  else if(tree->io->datatype == AA)
    {
      PhyML_Fprintf(fp_out,"\n. Model of amino acids substitution: \t%s",io->mod->modelname);
      if(io->mod->whichmodel == CUSTOMAA) PhyML_Fprintf(fp_out," (%s)",io->aa_rate_mat_file);
    }
  else
    {
      fprintf(fp_out,"\n. Substitution model: \t\t\t%s",io->mod->modelname);
    }


  PhyML_Fprintf(fp_out,"\n. Number of taxa: \t\t\t%d",tree->n_otu);/*added FLT*/

  PhyML_Fprintf(fp_out,"\n. Log-likelihood: \t\t\t%.5f",tree->c_lnL);/*was last ; moved here FLT*/

  Unconstraint_Lk(tree);
  PhyML_Fprintf(fp_out,"\n. Unconstrained likelihood: \t\t%.5f",tree->unconstraint_lk);

  PhyML_Fprintf(fp_out,"\n. Parsimony: \t\t\t\t%d",tree->c_pars);

  PhyML_Fprintf(fp_out,"\n. Tree size: \t\t\t\t%.5f",tree->size);

  if(tree->mod->ras->n_catg > 1 && tree->mod->ras->free_mixt_rates == NO)
    {
      PhyML_Fprintf(fp_out,"\n. Discrete gamma model: \t\t%s","Yes");
      PhyML_Fprintf(fp_out,"\n  - Number of categories: \t\t%d",tree->mod->ras->n_catg);
      PhyML_Fprintf(fp_out,"\n  - Gamma shape parameter: \t\t%.3f",tree->mod->ras->alpha->v);
    }
  else if(tree->mod->ras->free_mixt_rates == YES)
    {
      PhyML_Fprintf(fp_out,"\n. Discrete gamma model: \t\t%s","No");
      PhyML_Fprintf(fp_out,"\n  - Number of categories: \t\t%d",tree->mod->ras->n_catg);
      For(i,tree->mod->ras->n_catg)
	{
	  PhyML_Fprintf(fp_out,"\n  - Relative rate in class %d: \t\t%.5f [prop=%4f] \t\t",i+1,tree->mod->ras->gamma_rr->v[i],tree->mod->ras->gamma_r_proba->v[i]);
	}
    }

  if(tree->mod->ras->invar) PhyML_Fprintf(fp_out,"\n. Proportion of invariant: \t\t%.3f",tree->mod->pinvar->v);

  /*was before Discrete gamma model ; moved here FLT*/
  if((tree->mod->whichmodel == K80)   ||
     (tree->mod->whichmodel == HKY85) ||
     (tree->mod->whichmodel == F84))
    PhyML_Fprintf(fp_out,"\n. Transition/transversion ratio: \t%.3f",tree->mod->kappa->v);
  else if(tree->mod->whichmodel == TN93)
    {
      PhyML_Fprintf(fp_out,"\n. Transition/transversion ratio for purines: \t\t\t%.3f",
		    tree->mod->kappa->v*2.*tree->mod->lambda->v/(1.+tree->mod->lambda->v));
      PhyML_Fprintf(fp_out,"\n. Transition/transversion ratio for pyrimidines: \t\t\t%.3f",
	      tree->mod->kappa->v*2./(1.+tree->mod->lambda->v));
    }

  if(tree->io->datatype == NT)
    {
      PhyML_Fprintf(fp_out,"\n. Nucleotides frequencies:");
      PhyML_Fprintf(fp_out,"\n  - f(A)=%8.5f",tree->mod->e_frq->pi->v[0]);
      PhyML_Fprintf(fp_out,"\n  - f(C)=%8.5f",tree->mod->e_frq->pi->v[1]);
      PhyML_Fprintf(fp_out,"\n  - f(G)=%8.5f",tree->mod->e_frq->pi->v[2]);
      PhyML_Fprintf(fp_out,"\n  - f(T)=%8.5f",tree->mod->e_frq->pi->v[3]);
    }

  /*****************************************/
  if((tree->mod->whichmodel == GTR) ||
     (tree->mod->whichmodel == CUSTOM))
    {
      int i,j;

      Update_Qmat_GTR(tree->mod->r_mat->rr->v,
		      tree->mod->r_mat->rr_val->v,
		      tree->mod->r_mat->rr_num->v,
		      tree->mod->e_frq->pi->v,
		      tree->mod->r_mat->qmat->v);

      PhyML_Fprintf(fp_out,"\n");
      PhyML_Fprintf(fp_out,". GTR relative rate parameters : \n");
      PhyML_Fprintf(fp_out,"  A <-> C   %8.5f\n",  tree->mod->r_mat->rr->v[0]);
      PhyML_Fprintf(fp_out,"  A <-> G   %8.5f\n",  tree->mod->r_mat->rr->v[1]);
      PhyML_Fprintf(fp_out,"  A <-> T   %8.5f\n",  tree->mod->r_mat->rr->v[2]);
      PhyML_Fprintf(fp_out,"  C <-> G   %8.5f\n",  tree->mod->r_mat->rr->v[3]);
      PhyML_Fprintf(fp_out,"  C <-> T   %8.5f\n",  tree->mod->r_mat->rr->v[4]);
      PhyML_Fprintf(fp_out,"  G <-> T   %8.5f\n",tree->mod->r_mat->rr->v[5]);


      PhyML_Fprintf(fp_out,"\n. Instantaneous rate matrix : ");
      PhyML_Fprintf(fp_out,"\n  [A---------C---------G---------T------]\n");
      For(i,4)
	{
	  PhyML_Fprintf(fp_out,"  ");
	  For(j,4)
	    PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->r_mat->qmat->v[i*4+j]);
	  PhyML_Fprintf(fp_out,"\n");
	}
      PhyML_Fprintf(fp_out,"\n");
    }
  /*****************************************/


  if(io->ratio_test == 1)
    {
      PhyML_Fprintf(fp_out,". aLRT statistics to test branches");
    }
  else if(io->ratio_test == 2)
    {
      PhyML_Fprintf(fp_out,". aLRT branch supports (cubic approximation, mixture of Chi2s distribution)");
    }


  PhyML_Fprintf(fp_out,"\n");
  PhyML_Fprintf(fp_out,"\n. Run ID:\t\t\t\t%s", (io->append_run_ID) ? (io->run_id_string): ("none"));
  PhyML_Fprintf(fp_out,"\n. Random seed:\t\t\t\t%d", io->r_seed);
  PhyML_Fprintf(fp_out,"\n. Subtree patterns aliasing:\t\t%s",io->do_alias_subpatt?"yes":"no");
  PhyML_Fprintf(fp_out,"\n. Version:\t\t\t\t%s", VERSION);

  hour = div(t_end-t_beg,3600);
  min  = div(t_end-t_beg,60  );

  min.quot -= hour.quot*60;

  PhyML_Fprintf(fp_out,"\n. Time used:\t\t\t\t%dh%dm%ds (%d seconds)", hour.quot,min.quot,(int)(t_end-t_beg)%60,(int)(t_end-t_beg));

  PhyML_Fprintf(fp_out,"\n\n");
  PhyML_Fprintf(fp_out," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
  PhyML_Fprintf(fp_out," Suggested citations:\n");
  PhyML_Fprintf(fp_out," S. Guindon, JF. Dufayard, V. Lefort, M. Anisimova, W. Hordijk, O. Gascuel\n");
  PhyML_Fprintf(fp_out," \"New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0.\"\n");
  PhyML_Fprintf(fp_out," Systematic Biology. 2010. 59(3):307-321.\n");
  PhyML_Fprintf(fp_out,"\n");
  PhyML_Fprintf(fp_out," S. Guindon & O. Gascuel\n");
  PhyML_Fprintf(fp_out," \"A simple, fast, and accurate algorithm to estimate large phylogenies by maximum likelihood\"\n");
  PhyML_Fprintf(fp_out," Systematic Biology. 2003. 52(5):696-704.\n");
  PhyML_Fprintf(fp_out," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/*FLT wrote this function*/
void Print_Fp_Out_Lines(FILE *fp_out, time_t t_beg, time_t t_end, t_tree *tree, option *io, int n_data_set)
{
  char *s;
  /*div_t hour,min;*/

  if (n_data_set==1)
      {

	PhyML_Fprintf(fp_out,". Sequence file : [%s]\n\n", Basename(io->in_align_file));

	if((tree->io->datatype == NT) || (tree->io->datatype == AA))
	  {
	    (tree->io->datatype == NT)?
	      (PhyML_Fprintf(fp_out,". Model of nucleotides substitution : %s\n\n",io->mod->modelname)):
	      (PhyML_Fprintf(fp_out,". Model of amino acids substitution : %s\n\n",io->mod->modelname));
	  }

	s = (char *)mCalloc(100,sizeof(char));

	switch(io->in_tree)
	  {
	  case 0: { strcpy(s,"BioNJ");     break; }
	  case 1: { strcpy(s,"parsimony"); break; }
	  case 2: { strcpy(s,"user tree ("); 
	            strcat(s,io->in_tree_file); 
	            strcat(s,")");         break; }
	  }

	PhyML_Fprintf(fp_out,". Initial tree : [%s]\n\n",s);

	Free(s);

	PhyML_Fprintf(fp_out,"\n");

	/*headline 1*/
	PhyML_Fprintf(fp_out, ". Data\t");

	PhyML_Fprintf(fp_out,"Nb of \t");

	PhyML_Fprintf(fp_out,"Likelihood\t");

	PhyML_Fprintf(fp_out, "Discrete   \t");

	if(tree->mod->ras->n_catg > 1)
	  PhyML_Fprintf(fp_out, "Number of \tGamma shape\t");

	PhyML_Fprintf(fp_out,"Proportion of\t");

	if(tree->mod->whichmodel <= 6)
	  PhyML_Fprintf(fp_out,"Transition/ \t");

	PhyML_Fprintf(fp_out,"Nucleotides frequencies               \t");

	if((tree->mod->whichmodel == GTR) ||
	   (tree->mod->whichmodel == CUSTOM))
	  PhyML_Fprintf(fp_out,"Instantaneous rate matrix              \t");

	/*    PhyML_Fprintf(fp_out,"Time\t");*/

	PhyML_Fprintf(fp_out, "\n");


	/*headline 2*/
	PhyML_Fprintf(fp_out, "  set\t");

	PhyML_Fprintf(fp_out,"taxa\t");

	PhyML_Fprintf(fp_out,"loglk     \t");

	PhyML_Fprintf(fp_out, "gamma model\t");

	if(tree->mod->ras->n_catg > 1)
	  PhyML_Fprintf(fp_out, "categories\tparameter  \t");

	PhyML_Fprintf(fp_out,"invariant    \t");

	if(tree->mod->whichmodel <= 6)
	  PhyML_Fprintf(fp_out,"transversion\t");

	PhyML_Fprintf(fp_out,"f(A)      f(C)      f(G)      f(T)    \t");

	if((tree->mod->whichmodel == GTR) ||
	   (tree->mod->whichmodel == CUSTOM))
	  PhyML_Fprintf(fp_out,"[A---------C---------G---------T------]\t");

	/*    PhyML_PhyML_Fprintf(fp_out,"used\t");*/

	PhyML_Fprintf(fp_out, "\n");


	/*headline 3*/
	if(tree->mod->whichmodel == TN93)
	  {
	    PhyML_Fprintf(fp_out,"    \t      \t          \t           \t");
	    if(tree->mod->ras->n_catg > 1) PhyML_Fprintf(fp_out,"         \t         \t");
	    PhyML_Fprintf(fp_out,"             \t");
	    PhyML_Fprintf(fp_out,"purines pyrimid.\t");
	    PhyML_Fprintf(fp_out, "\n");
          }

          PhyML_Fprintf(fp_out, "\n");
      }


  /*line items*/

  PhyML_Fprintf(fp_out,"  #%d\t",n_data_set);

  PhyML_Fprintf(fp_out,"%d   \t",tree->n_otu);

  PhyML_Fprintf(fp_out,"%.5f\t",tree->c_lnL);

  PhyML_Fprintf(fp_out,"%s        \t",
	  (tree->mod->ras->n_catg>1)?("Yes"):("No "));
  if(tree->mod->ras->n_catg > 1)
    {
      PhyML_Fprintf(fp_out,"%d        \t",tree->mod->ras->n_catg);
      PhyML_Fprintf(fp_out,"%.3f    \t",tree->mod->ras->alpha->v);
    }

  /*if(tree->mod->ras->invar)*/
    PhyML_Fprintf(fp_out,"%.3f    \t",tree->mod->pinvar->v);

  if(tree->mod->whichmodel <= 5)
    {
      PhyML_Fprintf(fp_out,"%.3f     \t",tree->mod->kappa->v);
    }
  else if(tree->mod->whichmodel == TN93)
    {
      PhyML_Fprintf(fp_out,"%.3f   ",
		    tree->mod->kappa->v*2.*tree->mod->lambda->v/(1.+tree->mod->lambda->v));
      PhyML_Fprintf(fp_out,"%.3f\t",
	      tree->mod->kappa->v*2./(1.+tree->mod->lambda->v));
    }


  if(tree->io->datatype == NT)
    {
      PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->e_frq->pi->v[0]);
      PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->e_frq->pi->v[1]);
      PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->e_frq->pi->v[2]);
      PhyML_Fprintf(fp_out,"%8.5f\t",tree->mod->e_frq->pi->v[3]);
    }
  /*
  hour = div(t_end-t_beg,3600);
  min  = div(t_end-t_beg,60  );

  min.quot -= hour.quot*60;

  PhyML_Fprintf(fp_out,"%dh%dm%ds\t", hour.quot,min.quot,(int)(t_end-t_beg)%60);
  if(t_end-t_beg > 60)
    PhyML_Fprintf(fp_out,". -> %d seconds\t",(int)(t_end-t_beg));
  */

  /*****************************************/
  if((tree->mod->whichmodel == GTR) || (tree->mod->whichmodel == CUSTOM))
    {
      int i,j;

      For(i,4)
	{
	  if (i!=0) {
	    /*format*/
	    PhyML_Fprintf(fp_out,"      \t     \t          \t           \t");
	    if(tree->mod->ras->n_catg > 1) PhyML_Fprintf(fp_out,"          \t           \t");
	    PhyML_Fprintf(fp_out,"             \t                                      \t");
	  }
	  For(j,4)
	    PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->r_mat->qmat->v[i*4+j]);
	  if (i<3) PhyML_Fprintf(fp_out,"\n");
	}
    }
  /*****************************************/

  PhyML_Fprintf(fp_out, "\n\n");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Print_Freq(t_tree *tree)
{

  switch(tree->io->datatype)
    {
    case NT:
      {
	PhyML_Printf("A : %f\n",tree->mod->e_frq->pi->v[0]);
	PhyML_Printf("C : %f\n",tree->mod->e_frq->pi->v[1]);
	PhyML_Printf("G : %f\n",tree->mod->e_frq->pi->v[2]);
	PhyML_Printf("T : %f\n",tree->mod->e_frq->pi->v[3]);

	/* PhyML_Printf("U : %f\n",tree->mod->e_frq->pi->v[4]); */
	/* PhyML_Printf("M : %f\n",tree->mod->e_frq->pi->v[5]); */
	/* PhyML_Printf("R : %f\n",tree->mod->e_frq->pi->v[6]); */
	/* PhyML_Printf("W : %f\n",tree->mod->e_frq->pi->v[7]); */
	/* PhyML_Printf("S : %f\n",tree->mod->e_frq->pi->v[8]); */
	/* PhyML_Printf("Y : %f\n",tree->mod->e_frq->pi->v[9]); */
	/* PhyML_Printf("K : %f\n",tree->mod->e_frq->pi->v[10]); */
	/* PhyML_Printf("B : %f\n",tree->mod->e_frq->pi->v[11]); */
	/* PhyML_Printf("D : %f\n",tree->mod->e_frq->pi->v[12]); */
	/* PhyML_Printf("H : %f\n",tree->mod->e_frq->pi->v[13]); */
	/* PhyML_Printf("V : %f\n",tree->mod->e_frq->pi->v[14]); */
	/* PhyML_Printf("N : %f\n",tree->mod->e_frq->pi->v[15]); */
	break;
      }
    case AA:
      {
	PhyML_Printf("A : %f\n",tree->mod->e_frq->pi->v[0]);
	PhyML_Printf("R : %f\n",tree->mod->e_frq->pi->v[1]);
	PhyML_Printf("N : %f\n",tree->mod->e_frq->pi->v[2]);
	PhyML_Printf("D : %f\n",tree->mod->e_frq->pi->v[3]);
	PhyML_Printf("C : %f\n",tree->mod->e_frq->pi->v[4]);
	PhyML_Printf("Q : %f\n",tree->mod->e_frq->pi->v[5]);
	PhyML_Printf("E : %f\n",tree->mod->e_frq->pi->v[6]);
	PhyML_Printf("G : %f\n",tree->mod->e_frq->pi->v[7]);
	PhyML_Printf("H : %f\n",tree->mod->e_frq->pi->v[8]);
	PhyML_Printf("I : %f\n",tree->mod->e_frq->pi->v[9]);
	PhyML_Printf("L : %f\n",tree->mod->e_frq->pi->v[10]);
	PhyML_Printf("K : %f\n",tree->mod->e_frq->pi->v[11]);
	PhyML_Printf("M : %f\n",tree->mod->e_frq->pi->v[12]);
	PhyML_Printf("F : %f\n",tree->mod->e_frq->pi->v[13]);
	PhyML_Printf("P : %f\n",tree->mod->e_frq->pi->v[14]);
	PhyML_Printf("S : %f\n",tree->mod->e_frq->pi->v[15]);
	PhyML_Printf("T : %f\n",tree->mod->e_frq->pi->v[16]);
	PhyML_Printf("W : %f\n",tree->mod->e_frq->pi->v[17]);
	PhyML_Printf("Y : %f\n",tree->mod->e_frq->pi->v[18]);
	PhyML_Printf("V : %f\n",tree->mod->e_frq->pi->v[19]);

	PhyML_Printf("N : %f\n",tree->mod->e_frq->pi->v[20]);
	break;
      }
    default : {break;}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Print_Settings(option *io)
{
  int answer;
  char *s;

  s = (char *)mCalloc(100,sizeof(char));
  
  PhyML_Printf("\n\n\n");
  PhyML_Printf("\n\n");

  PhyML_Printf("                                 ..........................                                      \n");
  PhyML_Printf(" ooooooooooooooooooooooooooooo        CURRENT SETTINGS        ooooooooooooooooooooooooooooooooooo\n");
  PhyML_Printf("                                 ..........................                                      \n");

  PhyML_Printf("\n                . Sequence filename:\t\t\t\t %s", Basename(io->in_align_file));

  if(io->datatype == NT) strcpy(s,"dna");
  else if(io->datatype == AA) strcpy(s,"aa");
  else strcpy(s,"generic");

  PhyML_Printf("\n                . Data type:\t\t\t\t\t %s",s);
  PhyML_Printf("\n                . Alphabet size:\t\t\t\t %d",io->mod->ns);

  PhyML_Printf("\n                . Sequence format:\t\t\t\t %s", io->interleaved ? "interleaved": "sequential");
  PhyML_Printf("\n                . Number of data sets:\t\t\t\t %d", io->n_data_sets);

  PhyML_Printf("\n                . Nb of bootstrapped data sets:\t\t\t %d", io->mod->bootstrap);

  if (io->mod->bootstrap > 0)
    PhyML_Printf("\n                . Compute approximate likelihood ratio test:\t no");
  else
    {
      if(io->ratio_test == 1)
	PhyML_Printf("\n                . Compute approximate likelihood ratio test:\t yes (aLRT statistics)");
      else if(io->ratio_test == 2)
	PhyML_Printf("\n                . Compute approximate likelihood ratio test:\t yes (Chi2-based parametric branch supports)");
      else if(io->ratio_test == 3)
	PhyML_Printf("\n                . Compute approximate likelihood ratio test:\t yes (Minimum of SH-like and Chi2-based branch supports)");
      else if(io->ratio_test == 4)
	PhyML_Printf("\n                . Compute approximate likelihood ratio test:\t yes (SH-like branch supports)");
      else if(io->ratio_test == 5)
	PhyML_Printf("\n                . Compute approximate likelihood ratio test:\t yes (aBayes branch supports)");
    }

  PhyML_Printf("\n                . Model name:\t\t\t\t\t %s", io->mod->modelname);

  if(io->datatype == AA && io->mod->whichmodel == CUSTOMAA) PhyML_Printf(" (%s)",io->aa_rate_mat_file);

  if (io->datatype == NT)
    {
      if ((io->mod->whichmodel == K80)  ||
	  (io->mod->whichmodel == HKY85)||
	  (io->mod->whichmodel == F84)  ||
	  (io->mod->whichmodel == TN93))
	{
	  if (io->mod->s_opt->opt_kappa)
	    PhyML_Printf("\n                . Ts/tv ratio:\t\t\t\t\t estimated");
	  else
	    PhyML_Printf("\n                . Ts/tv ratio:\t\t\t\t\t %f", io->mod->kappa->v);
	}
    }

  if (io->mod->s_opt->opt_pinvar)
    PhyML_Printf("\n                . Proportion of invariable sites:\t\t estimated");
  else
    PhyML_Printf("\n                . Proportion of invariable sites:\t\t %f", io->mod->pinvar->v);


  PhyML_Printf("\n                . Number of subst. rate categs:\t\t\t %d", io->mod->ras->n_catg);
  if(io->mod->ras->n_catg > 1)
    {
      if(io->mod->ras->free_mixt_rates == NO)
	{
	  if(io->mod->s_opt->opt_alpha)
	    PhyML_Printf("\n                . Gamma distribution parameter:\t\t\t estimated");
	  else
	    PhyML_Printf("\n                . Gamma distribution parameter:\t\t\t %f", io->mod->ras->alpha->v);
	  PhyML_Printf("\n                . 'Middle' of each rate class:\t\t\t %s",(io->mod->ras->gamma_median)?("median"):("mean"));
	}
    }
    
  
  if(io->datatype == AA)
    PhyML_Printf("\n                . Amino acid equilibrium frequencies:\t\t %s", (io->mod->s_opt->opt_state_freq) ? ("empirical"):("model"));
  else if(io->datatype == NT)
    {
      if((io->mod->whichmodel != JC69) &&
	 (io->mod->whichmodel != K80)  &&
	 (io->mod->whichmodel != F81))
	{
	  if(!io->mod->s_opt->user_state_freq)
	    {
	      PhyML_Printf("\n                . Nucleotide equilibrium frequencies:\t\t %s", (io->mod->s_opt->opt_state_freq) ? ("ML"):("empirical"));
	    }
	  else
	    {
	      PhyML_Printf("\n                . Nucleotide equilibrium frequencies:\t\t %s","user-defined");
	    }
	}
    }

  PhyML_Printf("\n                . Optimise tree topology:\t\t\t %s", (io->mod->s_opt->opt_topo) ? "yes": "no");

  switch(io->in_tree)
    {
    case 0: { strcpy(s,"BioNJ");     break; }
    case 1: { strcpy(s,"parsimony"); break; }
    case 2: { strcpy(s,"user tree ("); 
	strcat(s,Basename(io->in_tree_file)); 
	strcat(s,")");         break; }
    }

  if(io->mod->s_opt->opt_topo)
    {
      if(io->mod->s_opt->topo_search == NNI_MOVE) PhyML_Printf("\n                . Tree topology search:\t\t\t\t NNIs");
      else if(io->mod->s_opt->topo_search == SPR_MOVE) PhyML_Printf("\n                . Tree topology search:\t\t\t\t SPRs");
      else if(io->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR) PhyML_Printf("\n                . Tree topology search:\t\t\t\t Best of NNIs and SPRs");



      PhyML_Printf("\n                . Starting tree:\t\t\t\t %s",s);

      PhyML_Printf("\n                . Add random input tree:\t\t\t %s", (io->mod->s_opt->random_input_tree) ? "yes": "no");
      if(io->mod->s_opt->random_input_tree)
	PhyML_Printf("\n                . Number of random starting trees:\t\t %d", io->mod->s_opt->n_rand_starts);	
    }
  else
    if(!io->mod->s_opt->random_input_tree)
      PhyML_Printf("\n                . Evaluated tree:\t\t\t\t \"%s\"",s);

  PhyML_Printf("\n                . Optimise branch lengths:\t\t\t %s", (io->mod->s_opt->opt_bl) ? "yes": "no");

  answer = 0;
  if(io->mod->s_opt->opt_alpha  ||
     io->mod->s_opt->opt_kappa  ||
     io->mod->s_opt->opt_lambda ||
     io->mod->s_opt->opt_pinvar ||
     io->mod->s_opt->opt_rr) answer = 1;
  
  PhyML_Printf("\n                . Optimise substitution model parameters:\t %s", (answer) ? "yes": "no");

  PhyML_Printf("\n                . Run ID:\t\t\t\t\t %s", (io->append_run_ID) ? (io->run_id_string): ("none"));
  PhyML_Printf("\n                . Random seed:\t\t\t\t\t %d", io->r_seed);
  PhyML_Printf("\n                . Subtree patterns aliasing:\t\t\t %s",io->do_alias_subpatt?"yes":"no");
  PhyML_Printf("\n                . Version:\t\t\t\t\t %s", VERSION);


  PhyML_Printf("\n\n oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");

  PhyML_Printf("\n\n");
  fflush(NULL);
  
  Free(s);
}

void Print_Banner(FILE *fp)
{
  PhyML_Fprintf(fp,"\n");
  PhyML_Fprintf(fp," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
  PhyML_Fprintf(fp,"                                                                                                  \n");
  PhyML_Fprintf(fp,"                                 ---  PhyML %s  ---                                             \n",VERSION);
  PhyML_Fprintf(fp,"                                                                                                  \n");
  PhyML_Fprintf(fp,"    A simple, fast, and accurate algorithm to estimate large phylogenies by maximum likelihood    \n");
  PhyML_Fprintf(fp,"                            Stephane Guindon & Olivier Gascuel                                      \n");
  PhyML_Fprintf(fp,"                                                                                                  \n");
  PhyML_Fprintf(fp,"                           http://www.atgc-montpellier.fr/phyml                                          \n");
  PhyML_Fprintf(fp,"                                                                                                  \n");
  PhyML_Fprintf(fp,"                         Copyright CNRS - Universite Montpellier II                                 \n");
  PhyML_Fprintf(fp,"                                                                                                  \n");
  PhyML_Fprintf(fp," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Print_Banner_Small(FILE *fp)
{
  PhyML_Fprintf(fp,"\n");
  PhyML_Fprintf(fp," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
  PhyML_Fprintf(fp,"                                  ---  PhyML %s  ---                                             \n",VERSION);
  PhyML_Fprintf(fp,"                            http://www.atgc-montpellier.fr/phyml                                          \n");
  PhyML_Fprintf(fp,"                         Copyright CNRS - Universite Montpellier II                                 \n");
  PhyML_Fprintf(fp," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



void Print_Data_Set_Number(option *io, FILE *fp)
{
  PhyML_Fprintf(fp,"\n");
  PhyML_Fprintf(fp,"                                                                                                  \n");
  PhyML_Fprintf(fp,"                                 [ Data set number %3d ]                                           \n",io->curr_gt+1);
  PhyML_Fprintf(fp,"                                                                                                  \n");
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Lk(t_tree *tree, char *string)
{
  time(&(tree->t_current));
  PhyML_Printf("\n. (%5d sec) [%15.4f] %s",
	       (int)(tree->t_current-tree->t_beg),tree->c_lnL,
	       string);
#ifndef QUIET 
  fflush(NULL);
#endif
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Print_Pars(t_tree *tree)
{
  time(&(tree->t_current));
  PhyML_Printf("\n. (%5d sec) [%5d]",(int)(tree->t_current-tree->t_beg),tree->c_pars);
#ifndef QUIET
  fflush(NULL);
#endif
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Lk_And_Pars(t_tree *tree)
{	
  time(&(tree->t_current));

  PhyML_Printf("\n. (%5d sec) [%15.4f] [%5d]",
	 (int)(tree->t_current-tree->t_beg),
	 tree->c_lnL,tree->c_pars);

#ifndef QUIET
  fflush(NULL);
#endif
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Read_Qmat(phydbl *daa, phydbl *pi, FILE *fp)
{
  int i,j;
  phydbl sum;
  double val;

  rewind(fp);

  for(i=1;i<20;i++)
    {
      For(j,19)
	{
/* 	  if(!fscanf(fp,"%lf",&(daa[i*20+j]))) Exit("\n"); */
	  if(!fscanf(fp,"%lf",&val)) 
	    {
	      PhyML_Printf("\n== Rate matrix file does not appear to have a proper format. Please refer to the documentation.");
	      Exit("\n");
	    }
	  daa[i*20+j] = (phydbl)val;
	  daa[j*20+i] = daa[i*20+j];
	  if(j == i-1) break; 
	}
    }


  For(i,20) 
    { 
      if(!fscanf(fp,"%lf",&val)) Exit("\n");
      pi[i] = (phydbl)val;
    }
  sum = .0;
  For(i,20) sum += pi[i];
  if(FABS(sum - 1.) > 1.E-06)
    {
      PhyML_Printf("\n. Sum=%f",sum);
      PhyML_Printf("\n. Scaling amino-acid frequencies...\n");
      For(i,20) pi[i] /= sum;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Print_Qmat_AA(phydbl *daa, phydbl *pi)
{
  int i,j,cpt;

  cpt = 0;
  For(i,20)
    {
      for(j=0;j<i;j++)
	{
	  PhyML_Printf("daa[%2d*20+%2d] = %10f;  ",i,j,daa[i*20+j]);
	  cpt++;
	  if(!(cpt%4)) PhyML_Printf("\n");
	}
    }

  PhyML_Printf("\n\n");
  PhyML_Printf("for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j];\n\n");
  For(i,20) PhyML_Printf("pi[%d] = %f; ",i,pi[i]);
  PhyML_Printf("\n");
  PhyML_Printf("Ala\tArg\tAsn\tAsp\tCys\tGln\tGlu\tGly\tHis\tIle\tLeu\tLys\tMet\tPhe\tPro\tSer\tThr\tTrp\tTyr\tVal\n");
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Square_Matrix_Generic(int n, phydbl *mat)
{
  int i,j;

  PhyML_Printf("\n");
  For(i,n)
    {
      PhyML_Printf("[%3d]",i);
      For(j,n)
	{
	  PhyML_Printf("%12.5f ",mat[i*n+j]);
	}
      PhyML_Printf("\n");
    }
  PhyML_Printf("\n");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Print_Diversity(FILE *fp, t_tree *tree)
{
  
  Print_Diversity_Pre(tree->t_nodes[0],
		      tree->t_nodes[0]->v[0],
		      tree->t_nodes[0]->b[0],
		      fp,
		      tree);

/*       mean_div_left = .0; */
/*       For(k,ns)  */
/* 	{ */
/* 	  mean_div_left += (k+1) * tree->t_edges[j]->div_post_pred_left[k]; */
/* 	} */
/*       mean_div_rght = .0; */
/*       For(k,ns) mean_div_rght += (k+1) * tree->t_edges[j]->div_post_pred_rght[k]; */

/*       mean_div_left /= (phydbl)tree->data->init_len; */
/*       mean_div_rght /= (phydbl)tree->data->init_len; */

/*       PhyML_Fprintf(fp,"%4d 0 %f\n",j,mean_div_left); */
/*       PhyML_Fprintf(fp,"%4d 1 %f\n",j,mean_div_rght); */


/*       mean_div_left = .0; */
/*       For(k,ns) mean_div_left += tree->t_edges[j]->div_post_pred_left[k]; */

/*       mean_div_rght = .0; */
/*       For(k,ns)  */
/* 	{ */
/* 	  mean_div_rght += tree->t_edges[j]->div_post_pred_rght[k]; */
/* 	} */

/*       if((mean_div_left != tree->data->init_len) || (mean_div_rght != tree->data->init_len)) */
/* 	{ */
/* 	  PhyML_Printf("\n. mean_div_left = %f mean_div_rght = %f init_len = %d", */
/* 		 mean_div_left,mean_div_rght,tree->data->init_len); */
/* 	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
/* 	  Warn_And_Exit(""); */
/* 	} */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Print_Diversity_Pre(t_node *a, t_node *d, t_edge *b, FILE *fp, t_tree *tree)
{
  int k,ns;

  ns = -1;

  if(d->tax) return;
  else
    {

      if(tree->io->datatype == NT)      ns = 4;
      else if(tree->io->datatype == AA) ns = 20;

      if(d == b->left) For(k,ns) PhyML_Fprintf(fp,"%4d 0 %2d %4d\n",b->num,k,b->div_post_pred_left[k]);
      else             For(k,ns) PhyML_Fprintf(fp,"%4d 1 %2d %4d\n",b->num,k,b->div_post_pred_rght[k]);

      For(k,3) if(d->v[k] != a) Print_Diversity_Pre(d,d->v[k],d->b[k],fp,tree);
    }

}

t_tree *Read_User_Tree(calign *cdata, t_mod *mod, option *io)
{
  t_tree *tree;

  
  PhyML_Printf("\n. Reading tree..."); fflush(NULL);
  if(io->n_trees == 1) rewind(io->fp_in_tree);
  tree = Read_Tree_File_Phylip(io->fp_in_tree);
  if(!tree) Exit("\n. Input tree not found...");
  /* Add branch lengths if necessary */
  if(!tree->has_branch_lengths) Add_BioNJ_Branch_Lengths(tree,cdata,mod);

  return tree;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Print_Time_Info(time_t t_beg, time_t t_end)
{
  div_t hour,min;

  hour = div(t_end-t_beg,3600);
  min  = div(t_end-t_beg,60  );
  min.quot -= hour.quot*60;

  PhyML_Printf("\n. Time used %dh%dm%ds\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
  PhyML_Printf("\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PhyML_Printf(char *format, ...)
{
  va_list ptr;

 #ifdef MPI
  if(Global_myRank == 0)
    {
      va_start (ptr, format);
      vprintf (format, ptr);
      va_end(ptr);
    }
#else
      va_start (ptr, format);
      vprintf (format, ptr);
      va_end(ptr);
#endif
  
  fflush (NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void PhyML_Fprintf(FILE *fp, char *format, ...)
{
  va_list ptr;

#ifdef MPI
  if(Global_myRank == 0)
    {
      va_start (ptr, format);
      vfprintf (fp,format, ptr);
      va_end(ptr);
    }
#else
      va_start (ptr, format);
      vfprintf (fp,format, ptr);
      va_end(ptr);
#endif
  
  fflush (NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Read_Clade_Priors(char *file_name, t_tree *tree)
{
  FILE *fp;
  char *s,*line;
  int n_clade_priors;
  int clade_size;
  char **clade_list;
  int i,pos;
  phydbl prior_low,prior_up;
  int node_num;

  PhyML_Printf("\n");
  PhyML_Printf("\n. Reading prior on node ages.\n");

  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  s    = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  clade_list = (char **)mCalloc(tree->n_otu,sizeof(char *));
  For(i,tree->n_otu) clade_list[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));

  fp = Openfile(file_name,0);
  
  n_clade_priors = 0;  
  do
    {
      if(!fgets(line,T_MAX_LINE,fp)) break;

      clade_size = 0;
      pos = 0;
      do
	{
	  i = 0;

	  while(line[pos] == ' ') pos++;

	  while((line[pos] != ' ') && (line[pos] != '\n') && line[pos] != '#')
	    {
	      s[i] = line[pos];
	      i++;
	      pos++;
	    }
	  s[i] = '\0';

	  /* PhyML_Printf("\n. s = %s\n",s); */
	  
	  if(line[pos] == '\n' || line[pos] == '#') break;
	  pos++;

	  if(strcmp(s,"|"))
	    {
	      strcpy(clade_list[clade_size],s);
	      clade_size++;
	    }
	  else 
	    break;	    
	}
      while(1);

      
      if(line[pos] != '#' && line[pos] != '\n')
	{
	  double val1, val2;
/* 	  sscanf(line+pos,"%lf %lf",&prior_up,&prior_low); */
	  sscanf(line+pos,"%lf %lf",&val1,&val2);
	  prior_up = (phydbl)val1;
	  prior_low = (phydbl)val2;
	  node_num = -1;
	  if(!strcmp("@root@",clade_list[0])) node_num = tree->n_root->num;
	  else node_num = Find_Clade(clade_list, clade_size, tree);

	  n_clade_priors++;  

	  if(node_num < 0)
	    {
	      PhyML_Printf("\n");
	      PhyML_Printf("\n");
	      PhyML_Printf("\n. .................................................................");
	      PhyML_Printf("\n. WARNING: could not find any clade in the tree referred to with the following taxon names:");
	      For(i,clade_size) PhyML_Printf("\n. \"%s\"",clade_list[i]);	      
	      PhyML_Printf("\n. .................................................................");
	      /* sleep(3); */
	    }
	  else
	    {	      
	      tree->rates->t_has_prior[node_num] = YES;
	      tree->rates->t_prior_min[node_num] = MIN(prior_low,prior_up);
	      tree->rates->t_prior_max[node_num] = MAX(prior_low,prior_up);

	      if(FABS(prior_low - prior_up) < 1.E-6 && tree->t_nodes[node_num]->tax == YES)
		tree->rates->nd_t[node_num] = prior_low;

	      PhyML_Printf("\n");
	      PhyML_Printf("\n. [%3d]..................................................................",n_clade_priors);
	      PhyML_Printf("\n. Node %4d matches the clade referred to with the following taxon names:",node_num);
	      For(i,clade_size) PhyML_Printf("\n. - \"%s\"",clade_list[i]);
	      PhyML_Printf("\n. Lower bound set to: %15f time units.",MIN(prior_low,prior_up));
	      PhyML_Printf("\n. Upper bound set to: %15f time units.",MAX(prior_low,prior_up));
	      PhyML_Printf("\n. .......................................................................");
	    }
	}
    }
  while(1);
  
  
  PhyML_Printf("\n. Read prior information on %d %s.",n_clade_priors,n_clade_priors > 1 ? "clades":"clade");

  if(!n_clade_priors)
    {
      PhyML_Printf("\n. PhyTime could not find any prior on node age.");
      PhyML_Printf("\n. This is likely due to a problem in the calibration ");
      PhyML_Printf("\n. file format. Make sure, for instance, that there is ");
      PhyML_Printf("\n. a blank character between the end of the last name");
      PhyML_Printf("\n. of each clade and the character `|'. Otherwise, ");
      PhyML_Printf("\n. please refer to the example file.\n");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  For(i,tree->n_otu) Free(clade_list[i]);
  Free(clade_list);
  Free(line);
  Free(s);
  fclose(fp);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

option *Get_Input(int argc, char **argv)
{

  option *io;
  t_mod *mod;
  t_opt *s_opt;
  m4 *m4mod;

  io    = (option *)Make_Input();
  mod   = (t_mod *)Make_Model_Basic();
  s_opt = (t_opt *)Make_Optimiz();
  m4mod = (m4 *)M4_Make_Light();

  Set_Defaults_Input(io);
  Set_Defaults_Model(mod);
  Set_Defaults_Optimiz(s_opt);

  io->mod        = mod;
  io->mod->m4mod = m4mod;

  mod->io        = io;
  mod->s_opt     = s_opt;


#ifdef MPI
  Read_Command_Line(io,argc,argv);
#elif defined (PHYTIME)
  Read_Command_Line(io,argc,argv);
#else
  putchar('\n');

  switch (argc)
    {
    case 1:
      {
	Launch_Interface(io);
	break;
      }
      /*
	case 2:
	Usage();
	break;
      */
    default:
      Read_Command_Line(io,argc,argv);
    }
#endif
  
  return io;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Set_Whichmodel(int select)
{
  int wm;

  wm = -1;

  switch(select)
    {
    case 1:
      {
	wm = JC69;
	break;
      }
    case 2:
      {
	wm = K80;
	break;
      }
    case 3:
      {
	wm = F81;
	break;
      }
    case 4:
      {
	wm = HKY85;
	break;
      }
    case 5:
      {
	wm = F84;
	break;
      }
    case 6:
      {
	wm = TN93;
	break;
      }
    case 7:
      {
	wm = GTR;
	break;
      }
    case 8:
      {
	wm = CUSTOM;
	break;
      }
    case 11:
      {
	wm = WAG;
	break;
      }
    case 12:
      {
	wm = DAYHOFF;
	break;
      }
    case 13:
      {
	wm = JTT;
	break;
      }
    case 14:
      {
	wm = BLOSUM62;
	break;
      }
    case 15:
      {
	wm = MTREV;
	break;
      }
    case 16:
      {
	wm = RTREV;
	break;
      }
    case 17:
      {
	wm = CPREV;
	break;
      }
    case 18:
      {
	wm = DCMUT;
	break;
      }
    case 19:
      {
	wm = VT;
	break;
      }
    case 20:
      {
	wm = MTMAM;
	break;
      }
    case 21:
      {
	wm = MTART;
	break;
      }
    case 22:
      {
	wm = HIVW;
	break;
      }
    case 23:
      {
	wm = HIVB;
	break;
      }
    case 24:
      {
	wm = CUSTOMAA;
	break;
      }
    case 25:
      {
	wm = LG;
	break;
      }
    default:
      {
	PhyML_Printf("\n== Model number %d is unknown. Please use a valid model name",select);
	Exit("\n");
	break;
      }
    }

  return wm;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Print_Data_Structure(t_tree *root)
{
  int n_partition_elem;
  void *buff;
  char *s;

  PhyML_Printf("\n\n.:: Printing the data structure ::. \n\n");

  PhyML_Printf("\n. Starting tree:\t%12s",
	       root->io->in_tree == 2?root->io->in_tree_file:"BioNJ");

  PhyML_Printf("\n. Tree topology search:\t%12s",
	       root->io->mod->s_opt->opt_topo==YES?
	       root->io->mod->s_opt->topo_search==SPR_MOVE?"spr":
	       root->io->mod->s_opt->topo_search==NNI_MOVE?"nni":
	       "spr+nni":"no");

  n_partition_elem = 1;
  buff = root;
  do
    {
      buff = ((t_tree *)buff)->next;
      if(!buff) break;
      n_partition_elem++;
    }
  while(1);
  
  s = (char *)mCalloc(2,sizeof(char));
  s[0] = ' ';
  s[1] = '\0';
  buff = root;
  do
    {
      s = (char *)mRealloc(s,(int)(strlen(s)+strlen(((t_tree *)buff)->io->in_align_file))+2+2,sizeof(char));
      strcat(s,((t_tree *)buff)->io->in_align_file);
      strcat(s,", ");
      buff = ((t_tree *)buff)->next;
      if(!buff) break;      
    }
  while(1);
  s[(int)strlen(s)-2]=' ';
  s[(int)strlen(s)-1]='\0';

  PhyML_Printf("\n. Processing %d data %s (%s)",n_partition_elem,n_partition_elem>1?"sets":"set",s);

  buff = root;
  
  do
    {
      int class = 0;
      t_tree *tree = ((t_tree *)buff)->child;
      
      PhyML_Printf("\n\n. Mixture model for data set %s",((t_tree *)buff)->io->in_align_file);
      PhyML_Printf("\n.");
      PhyML_Printf("\n. Number of rate classes:\t\t\t%12d",((t_tree *)buff)->mod->ras->n_catg);
      if(((t_tree *)buff)->mod->ras->n_catg > 1)
	{
	  PhyML_Printf("\n. Model of rate variation:\t\t\t%12s",
		       ((t_tree *)buff)->mod->ras->free_mixt_rates?"FreeRates":
		       ((t_tree *)buff)->mod->ras->invar?"Gamma+Inv":"Gamma");
	  if(((t_tree *)buff)->mod->ras->free_mixt_rates == NO)
	    {
	      PhyML_Printf("\n. Gamma shape parameter:\t\t\t%12f",((t_tree *)buff)->mod->ras->alpha->v);
	      PhyML_Printf("\n. Optimize parameter:\t\t\t\t%12s",((t_tree *)buff)->mod->s_opt->opt_alpha==YES?"yes":"no");
	    }
	  if(((t_tree *)buff)->mod->ras->invar == YES)
	    {
	      PhyML_Printf("\n. Proportion of invariable sites:\t\t%12f",((t_tree *)buff)->mod->pinvar->v);
	    }
	}
	  
      do
	{
	  PhyML_Printf("\n. Class %d",class+1);
	  PhyML_Printf("\n.\t\t Substitution model:\t\t%12s [%p]",tree->mod->modelname,(void *)tree->mod);
	  
	  if(tree->mod->whichmodel == K80 ||
	     tree->mod->whichmodel == HKY85 ||
	     tree->mod->whichmodel == TN93)
	    {
	      PhyML_Printf("\n.\t\t Ts/tv ratio:\t\t\t%12p",(void *)tree->mod->kappa);
	      PhyML_Printf("\n.\t\t Optimise ts/tv ratio:\t\t%12s",tree->mod->s_opt->opt_kappa?"yes":"no");
	      PhyML_Printf("\n.\t\t Value of the ts/tv ratio:\t%12f",tree->mod->kappa->v);
	    }
	  else if(tree->mod->whichmodel == GTR ||
		  tree->mod->whichmodel == CUSTOM)	    
	    PhyML_Printf("\n.\t\t Optimise subst. rates:\t\t%12s",tree->mod->s_opt->opt_rr?"yes":"no");

	  class++;
	  
	  tree = tree->next;
	  if(tree && tree->is_mixt_tree == YES) break;
	}
      while(tree);

      buff = ((t_tree *)buff)->next;
      if(!buff) break;
    }
  while(1);  
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
