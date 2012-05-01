/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef INIT_H
#define INIT_H

#include "utilities.h"

void Init_Eigen_Struct(eigen *this);
void Init_Scalar_Dbl(scalar_dbl *p);
void Init_Scalar_Int(scalar_int *p);
void Init_Vect_Dbl(int len,vect_dbl *p);
void Init_Vect_Int(int len,vect_int *p);
void Init_Tree(t_tree *tree,int n_otu);
void Init_Edge_Light(t_edge *b,int num);
void Init_Node_Light(t_node *n,int num);
void Init_NNI(nni *a_nni);
void Init_Nexus_Format(nexcom **com);
void Init_Mat(matrix *mat,calign *data);
void Set_Defaults_Input(option *io);
void Set_Defaults_Model(t_mod *mod);
void Set_Defaults_Optimiz(optimiz *s_opt);
void XML_Init_Node(xml_node *parent, xml_node *new, char *name);

#endif
