#include "sergeii.h"

int My_Function(char *s_tree)
{
  t_tree *tree;
  int i;
  
  printf("\n. string: %s\n",s_tree);

  tree = Read_Tree(&s_tree);
  PhyML_Printf("\n. Current tree: %s",Write_Tree(tree,NO));
  
  For(i,2*tree->n_otu-3)
    {
      tree->a_edges[i]->l->v *= 10.;      
    }

  PhyML_Printf("\n. New tree: %s",Write_Tree(tree,NO));
  
  return(0);
}
