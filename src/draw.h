#ifndef DRAW_H
#define DRAW_H

#include "utilities.h"

void DR_Dist_To_Root_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void DR_Dist_To_Root(t_node *n_root, t_tree *tree);
void DR_Get_X_Coord_Pre(t_node *a, t_node *d, t_edge *b, tdraw *w, t_tree *tree);
void DR_Get_X_Coord(tdraw *w, t_tree *tree);
tdraw *DR_Make_Tdraw_Struct(t_tree *tree);
void DR_Init_Tdraw_Struct(tdraw *d);
void DR_Get_Tree_Box_Width(tdraw *w, t_tree *tree);
void DR_Get_Y_Coord_Post(t_node *a, t_node *d, t_edge *b, int *next_y_slot, tdraw *w, t_tree *tree);
void DR_Get_Y_Coord(tdraw *w, t_tree *tree);
void DR_Get_Tree_Coord(t_tree *tree);
phydbl DR_Get_Max_Dist_To_Root(t_tree *tree);
void DR_Print_Tree_Postscript(int tree_num, FILE *fp, t_tree *tree);
void DR_Print_Tree_Postscript_Pre(t_node *a, t_node *d, FILE *fp, tdraw *w, t_tree *tree);
void DR_Print_Postscript_EOF(FILE *fp);
void DR_Print_Postscript_Header(int n_pages, FILE *fp);


#endif
