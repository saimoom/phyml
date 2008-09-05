/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#ifndef MCMC_H
#define MCMC_H

void MCMC_Lexp(arbre *tree);
void MCMC_Print_Param(FILE *fp, arbre *tree);
void MCMC_Rates_Std(arbre *tree);
void MCMC_Rates(arbre *tree);
void MCMC_Rates_Pre(node *a, node *d, arbre *tree);
tmcmc *MCMC_Make_MCMC_Struct();
void MCMC_Init_MCMC_Struct(tmcmc *mcmc);
void MCMC(arbre *tree);
void MCMC_Alpha(arbre *tree);
void MCMC_Times_Pre(node *a, node *d, arbre *tree);
void MCMC_Times(arbre *tree);
void MCMC_Randomize_Branch_Lengths(arbre *tree);
void MCMC_Randomize_Node_Times(arbre *tree);
void MCMC_Randomize_Node_Times_Pre(node *a, node *d, arbre *tree);
void MCMC_Randomize_Lexp(arbre *tree);
void MCMC_Randomize_Alpha(arbre *tree);
void MCMC_One_Node_Time(node *d, arbre *tree);
void MCMC_One_Rate(edge *b, arbre *tree);
void MCMC_One_Rate(edge *b, arbre *tree);
void MCMC_No_Change(arbre *tree);
void MCMC_Nu(arbre *tree);
void MCMC_Randomize_Nu(arbre *tree);
node *MCMC_Select_Random_Node_Pair(phydbl t_sup, arbre *tree);
void MCMC_Modify_Rates(arbre *tree);
void MCMC_Modify_Subtree_Rate(node *a, node *d, phydbl new_rate, arbre *tree);
void MCMC_Randomize_Rates(arbre *tree);
void MCMC_Step_Rate(arbre *tree);
void MCMC_Stick_Rates(arbre *tree);
void MCMC_Stick_Rates_Pre(node *a, node *d, arbre *tree);



#endif
