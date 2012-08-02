/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef MIXT_H
#define MIXT_H

#include "utilities.h"


void MIXT_Connect_Edges_To_Next_Prev_Child_Parent(t_tree *tree);
void MIXT_Connect_Nodes_To_Next_Prev_Child_Parent(t_tree *tree);
void MIXT_Turn_Branches_OnOff(int onoff, t_tree *tree);
void MIXT_Connect_Sprs_To_Next_Prev_Child_Parent(t_tree *tree);
phydbl *MIXT_Get_Lengths_Of_This_Edge(t_edge *mixt_b);
void MIXT_Set_Lengths_Of_This_Edge(phydbl *lens, t_edge *mixt_b);
void MIXT_Post_Order_Lk(t_node *mixt_a, t_node *mixt_d, t_tree *mixt_tree);
void MIXT_Pre_Order_Lk(t_node *mixt_a, t_node *mixt_d, t_tree *mixt_tree);

#endif