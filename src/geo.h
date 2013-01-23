/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef GEO_H
#define GEO_H

#include "utilities.h"

int GEO_Main(int argc, char **argv);
t_geo *GEO_Make_Geo_Basic();
void GEO_Make_Geo_Complete(int ldscape_sz,int n_dim,int n_tax,t_geo *t);
void Free_Geo(t_geo *t);
void GEO_Update_Fmat(t_geo *t);
void GEO_Update_Sorted_Nd(t_geo *t,t_tree *tree);
void GEO_Update_Occup(t_geo *t,t_tree *tree);
void GEO_Update_Rmat(t_node *n,t_geo *t,t_tree *tree);
phydbl GEO_Lk(t_geo *t,t_tree *tree);
void GEO_Init_Tloc_Tips(t_geo *t,t_tree *tree);
phydbl GEO_Total_Migration_Rate(t_node *n,t_geo *t);
int GEO_Get_Arrival_Location(t_node *n,t_geo *t,t_tree *tree);
void GEO_Simulate_Coordinates(int n, t_geo *t);
t_tree *GEO_Simulate(t_geo *t, int n_otu);
void GEO_Optimize_Sigma(t_geo *t, t_tree *tree);
phydbl GEO_Wrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree);
void GEO_Optimize_Lambda(t_geo *t, t_tree *tree);

#endif
