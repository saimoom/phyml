/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef MODELS_H
#define MODELS_H

#include "utilities.h"
#include "eigen.h"
#include "free.h"
#include "stats.h"

void PMat(phydbl l, t_mod *mod, int pos, phydbl *Pij);
void  PMat_K80(phydbl l,phydbl kappa, int pos, phydbl *Pij);
void  PMat_TN93(phydbl l, t_mod *mod, int pos, phydbl *Pij);
void  PMat_Empirical(phydbl l, t_mod *mod, int pos, phydbl *Pij);
void PMat_Zero_Br_Len(t_mod *mod, int pos, phydbl *Pij);
void PMat_Gamma(phydbl l, t_mod *mod, int pos, phydbl *Pij);

int GetDaa (phydbl *daa, phydbl *pi, char *file_name);
void Init_Model(calign *data, t_mod *mod, option *io);
void Update_Qmat_GTR(phydbl *rr, phydbl *rr_val, int *rr_num, phydbl *pi, phydbl *qmat);
void Update_Qmat_HKY(phydbl kappa, phydbl *pi, phydbl *qmat);
void Update_Qmat_Generic(phydbl *rr, phydbl *pi, int ns, phydbl *qmat);
void Translate_Custom_Mod_String(t_mod *mod);
void Set_Model_Parameters(t_mod *mod);
phydbl GTR_Dist(phydbl *F, phydbl alpha, eigen *eigen_struct);
phydbl General_Dist(phydbl *F, t_mod *mod, eigen *eigen_struct);

int Init_Qmat_WAG(phydbl *daa, phydbl *pi);
int Init_Qmat_Dayhoff(phydbl *daa, phydbl *pi);
int Init_Qmat_JTT(phydbl *daa, phydbl *pi);
int Init_Qmat_RtREV(phydbl *daa, phydbl *pi);
int Init_Qmat_CpREV(phydbl *daa, phydbl *pi);
int Init_Qmat_VT(phydbl *daa, phydbl *pi);
int Init_Qmat_Blosum62(phydbl *daa, phydbl *pi);
int Init_Qmat_MtMam(phydbl *daa, phydbl *pi);
int Init_Qmat_MtArt(phydbl *daa, phydbl *pi); // Added by Federico Abascal
int Init_Qmat_HIVb(phydbl *daa, phydbl *pi);  // Added by Federico Abascal
int Init_Qmat_HIVw(phydbl *daa, phydbl *pi);  // Added by Federico Abascal
void Switch_From_Mod_To_M4mod(t_mod *mod);
void Switch_From_M4mod_To_Mod(t_mod *mod);
void PMat_JC69(phydbl l, int pos, phydbl *Pij, t_mod *mod);
phydbl Get_Lambda_F84(phydbl *pi, phydbl *kappa);
void Update_Eigen(t_mod *mod);

#endif
