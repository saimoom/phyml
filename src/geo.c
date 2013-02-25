/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/


#include "geo.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int GEO_Main(int argc, char **argv)
{
  t_geo *t;
  int n_tax;
  t_tree *tree;
  int seed;
  int n_vars;
  phydbl *res;
  FILE *fp;
  int pid;
  char *s;
  int rand_loc;

  s = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  
  strcpy(s,"geo.out");
  pid = getpid();
  sprintf(s+strlen(s),".%d",pid);

  fp = fopen(s,"w");

  seed = getpid();
  /* seed = 28224; */
  printf("\n. Seed = %d",seed);
  srand(seed);

  t = GEO_Make_Geo_Basic();
  GEO_Init_Geo_Struct(t);


  t->tau        = Uni()*(t->max_tau-t->min_tau)  + t->min_tau;
  t->lbda       = EXP(Uni()*(LOG(t->max_lbda)-LOG(t->min_lbda))  + LOG(t->min_lbda));
  t->sigma      = Uni()*(t->max_sigma-t->min_sigma) + t->min_sigma;

  /* t->tau        = 3.0; */
  /* t->lbda       = 0.0001; */
  /* t->sigma      = 10.; */


  t->ldscape_sz = (int)atoi(argv[1]);
  t->n_dim      = 2;
  n_tax         = (int)atoi(argv[2]);

  PhyML_Fprintf(fp,"\n# SigmaTrue\t SigmaThresh\t LbdaTrue\t TauTrue\txTrue\t yTrue\t xRand\t yRand\t Sigma5\t Sigma50\t Sigma95\t Lbda5\t Lbda50\t Lbda95\t ProbLbdaInf1\t Tau5\t Tau50\t Tau95\t X5\t X50\t X95\t Y5\t Y50\t Y95\t RandX5\t RandX50\t RandX95\t RandY5\t RandY50\t RandY95\t");
  PhyML_Fprintf(fp,"\n");


  GEO_Make_Geo_Complete(t->ldscape_sz,t->n_dim,n_tax,t);

  t->cov[0*t->n_dim+0] = t->sigma;
  t->cov[1*t->n_dim+1] = t->sigma;
  t->cov[0*t->n_dim+1] = 0.0;
  t->cov[1*t->n_dim+0] = 0.0;

  GEO_Simulate_Coordinates(t->ldscape_sz,t);

  tree = GEO_Simulate(t,n_tax);

  GEO_Get_Sigma_Max(t);

  rand_loc = Rand_Int(0,t->ldscape_sz-1);

  PhyML_Printf("\nSigma: %f\t Lbda: %f\t Tau:%f\t x:%f\t y:%f rand.x:%f\t rand.y:%f\t",
               t->sigma,
               t->lbda,
               t->tau,
               t->ldscape[t->loc[tree->n_root->num]*t->n_dim+0],
               t->ldscape[t->loc[tree->n_root->num]*t->n_dim+1],
               t->ldscape[rand_loc*t->n_dim+0],
               t->ldscape[rand_loc*t->n_dim+1]);

  PhyML_Fprintf(fp,"%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t",
                t->sigma,
                t->sigma_thresh,
                t->lbda,
                t->tau,
                t->ldscape[t->loc[tree->n_root->num]*t->n_dim+0],
                t->ldscape[t->loc[tree->n_root->num]*t->n_dim+1],
                t->ldscape[rand_loc*t->n_dim+0],
                t->ldscape[rand_loc*t->n_dim+1]);

  GEO_Get_Locations_Beneath(t,tree);

  GEO_Randomize_Locations(tree->n_root,
                          tree->geo,
                          tree);
  GEO_Update_Occup(t,tree);
  GEO_Lk(t,tree);
  PhyML_Printf("\n. Init loglk: %f",tree->geo->c_lnL);

  /* Exit("\n"); */
  /* int c = 0; */
  /* do */
  /*   { */
  /*     GEO_Optimize_Sigma(t,tree); */
  /*     GEO_Optimize_Lambda(t,tree); */
  /*     GEO_Optimize_Tau(t,tree); */
  /*     printf("\n. Sigma = %f Lambda = %f Tau = %f Loglk = %f",tree->geo->sigma,tree->geo->lbda,tree->geo->tau,GEO_Lk(t,tree)); */
  /*     c++; */
  /*   } */
  /* while(c < 30); */


  tree->mcmc = MCMC_Make_MCMC_Struct();
  MCMC_Complete_MCMC(tree->mcmc,tree);

  tree->mcmc->io               = NULL;
  tree->mcmc->is               = NO;
  tree->mcmc->use_data         = YES;
  tree->mcmc->run              = 0;
  tree->mcmc->sample_interval  = 1E+3;
  tree->mcmc->chain_len        = 1E+6;
  tree->mcmc->chain_len_burnin = 1E+5;
  tree->mcmc->randomize        = YES;
  tree->mcmc->norm_freq        = 1E+3;
  tree->mcmc->max_tune         = 1.E+20;
  tree->mcmc->min_tune         = 1.E-10;
  tree->mcmc->print_every      = 2;
  tree->mcmc->is_burnin        = NO;
  tree->mcmc->nd_t_digits      = 1;


  t->tau   = 1.0;
  t->lbda  = 1.0;
  t->sigma = 1.0;

  n_vars = 10;
  tree->mcmc->chain_len = 1.E+8;
  tree->mcmc->sample_interval = 50;
  
  MCMC_Complete_MCMC(tree->mcmc,tree);

  GEO_Update_Occup(t,tree);
  GEO_Lk(t,tree);
  PhyML_Printf("\n. Init loglk: %f",tree->geo->c_lnL);
  
  res = (phydbl *)mCalloc(tree->mcmc->chain_len / tree->mcmc->sample_interval * n_vars,sizeof(phydbl));

  tree->mcmc->start_ess[tree->mcmc->num_move_geo_sigma]  = YES;
  tree->mcmc->start_ess[tree->mcmc->num_move_geo_lambda] = YES;
  tree->mcmc->start_ess[tree->mcmc->num_move_geo_tau]    = YES;

  tree->mcmc->run = 0;
  do
    {
      MCMC_Geo_Lbda(tree);
      MCMC_Geo_Sigma(tree);
      MCMC_Geo_Tau(tree);
      MCMC_Geo_Loc(tree);
      MCMC_Geo_Updown_Tau_Lbda(tree);

      
      if(tree->mcmc->run%tree->mcmc->sample_interval == 0)
        {
          MCMC_Copy_To_New_Param_Val(tree->mcmc,tree);
          
          MCMC_Update_Effective_Sample_Size(tree->mcmc->num_move_geo_lambda,tree->mcmc,tree);
          MCMC_Update_Effective_Sample_Size(tree->mcmc->num_move_geo_sigma,tree->mcmc,tree);
          MCMC_Update_Effective_Sample_Size(tree->mcmc->num_move_geo_tau,tree->mcmc,tree);

          PhyML_Printf("\n. Run %6d Sigma: %12f [%4.0f] Lambda: %12f [%4.0f] Tau: %12f [%4.0f] LogLk: %12f x: %12f y:%12f",
                       tree->mcmc->run,

                       tree->geo->sigma,
                       tree->mcmc->ess[tree->mcmc->num_move_geo_sigma],

                       tree->geo->lbda,
                       tree->mcmc->ess[tree->mcmc->num_move_geo_lambda],

                       tree->geo->tau,
                       tree->mcmc->ess[tree->mcmc->num_move_geo_tau],

                       tree->geo->c_lnL,
 
                       t->ldscape[t->loc[tree->n_root->num]*t->n_dim+0],
                       t->ldscape[t->loc[tree->n_root->num]*t->n_dim+1]);

          rand_loc = Rand_Int(0,t->ldscape_sz-1);

          res[0 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = tree->geo->sigma; 
          res[1 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = tree->geo->lbda; 
          res[2 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = tree->geo->tau; 
          res[3 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = tree->geo->c_lnL; 
          
                
          res[4 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = t->ldscape[t->loc[tree->n_root->num]*t->n_dim+0];
          res[5 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = t->ldscape[t->loc[tree->n_root->num]*t->n_dim+1];

          res[6 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = t->ldscape[rand_loc*t->n_dim+0];
          res[7 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = t->ldscape[rand_loc*t->n_dim+1];
        }

      tree->mcmc->run++;
      MCMC_Get_Acc_Rates(tree->mcmc);


      if(tree->mcmc->ess[tree->mcmc->num_move_geo_sigma] > 100. &&
         tree->mcmc->ess[tree->mcmc->num_move_geo_tau]   > 100. &&
         tree->mcmc->ess[tree->mcmc->num_move_geo_lambda]> 100.) break;

    }
  while(tree->mcmc->run < tree->mcmc->chain_len);
  
  
  PhyML_Fprintf(fp,"%f\t %f\t %f\t  %f\t %f\t %f\t  %f\t  %f\t %f\t %f\t  %f\t %f\t %f\t  %f\t %f\t %f\t  %f\t %f\t %f\t  %f\t %f\t %f\t  \n",
                Quantile(res+0*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.025),
                Quantile(res+0*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.50),
                Quantile(res+0*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.975),
                
                Quantile(res+1*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.025),
                Quantile(res+1*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.50),
                Quantile(res+1*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.975),
                
                Prob(res+1*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,1.0),
                
                Quantile(res+2*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.025),
                Quantile(res+2*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.50),
                Quantile(res+2*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.975),

                Quantile(res+4*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.025),
                Quantile(res+4*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.50),
                Quantile(res+4*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.975),

                Quantile(res+5*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.025),
                Quantile(res+5*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.50),
                Quantile(res+5*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.975),

                Quantile(res+6*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.025),
                Quantile(res+6*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.50),
                Quantile(res+6*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.975),

                Quantile(res+7*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.025),
                Quantile(res+7*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.50),
                Quantile(res+7*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.975)
                );
  
  Free(s);

  fclose(fp);

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_geo *GEO_Make_Geo_Basic()
{
  t_geo *t;
  t = (t_geo *)mCalloc(1,sizeof(t_geo));
  return(t);  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Make_Geo_Complete(int ldscape_sz, int n_dim, int n_tax, t_geo *t)
{
  // F matrix
  t->f_mat = (phydbl *)mCalloc(ldscape_sz*ldscape_sz,sizeof(phydbl));

  // R matrix
  t->r_mat = (phydbl *)mCalloc(ldscape_sz*ldscape_sz,sizeof(phydbl));

  // Occupation vectors: one vector for each node
  t->occup = (int *)mCalloc((int)(2*n_tax-1)*ldscape_sz,sizeof(int)); 

  // Locations
  t->ldscape = (phydbl *)mCalloc((int)(ldscape_sz*n_dim),sizeof(phydbl));

  // Lineage locations
  t->loc = (int *)mCalloc((int)(2*n_tax-1),sizeof(int));

  // Sorted node heights
  t->sorted_nd = (t_node **)mCalloc((int)(2*n_tax-1),sizeof(t_node *));

  // Covariance matrix
  t->cov = (phydbl *)mCalloc((int)(n_dim*n_dim),sizeof(phydbl));

  // gives the location occupied beneath each node in the tree
  t->loc_beneath = (int *)mCalloc((int)(2*n_tax-1)*ldscape_sz,sizeof(int));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Geo(t_geo *t)
{
  Free(t->f_mat);
  Free(t->r_mat);
  Free(t->occup);
  Free(t->loc);
  Free(t->ldscape);
  Free(t->sorted_nd);
  Free(t->cov);
  Free(t->loc_beneath);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Update F matrix. Assume a diagonal covariance matrix.
void GEO_Update_Fmat(t_geo *t)
{
  phydbl *loc1, *loc2;
  int i,j,k;
  int err;
  phydbl lognloc;
  
  For(i,t->n_dim) t->cov[i*t->n_dim+i] = t->sigma; // Diagonal covariance matrix. Same variance in every direction

  lognloc = LOG(t->ldscape_sz);

  // Fill in F matrix;
  for(i=0;i<t->ldscape_sz;i++)
    {
      loc1 = t->ldscape + i*t->n_dim;

      for(j=i;j<t->ldscape_sz;j++)
        {
          loc2 = t->ldscape + j*t->n_dim;

          t->f_mat[i*t->ldscape_sz+j] = .0;

          // Calculate log(f(l_i,l_j)) - log(f(l_i,l_i) - log(m) 
          For(k,t->n_dim) t->f_mat[i*t->ldscape_sz+j] += Log_Dnorm(loc2[k],loc1[k],SQRT(t->cov[k*t->n_dim+k]),&err);
          t->f_mat[i*t->ldscape_sz+j] -= lognloc;
          For(k,t->n_dim) t->f_mat[i*t->ldscape_sz+j] -= Log_Dnorm(loc1[k],loc1[k],SQRT(t->cov[k*t->n_dim+k]),&err);

          // Take the exponential
          t->f_mat[i*t->ldscape_sz+j] = EXP(t->f_mat[i*t->ldscape_sz+j]);
          
          /* printf("\n. i=%d j=%d f=%f %f %f %f %f",i,j,t->f_mat[i*t->ldscape_sz+j],loc1[0],loc2[0],loc1[1],loc2[1]); */

          // Matrix is symmetric
          t->f_mat[j*t->ldscape_sz+i] = t->f_mat[i*t->ldscape_sz+j];
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Sort node heights from oldest to youngest age.
void GEO_Update_Sorted_Nd(t_geo *t, t_tree *tree)
{
  int i;
  int swap;
  t_node *buff;

  buff = NULL;

  For(i,2*tree->n_otu-1) t->sorted_nd[i] = tree->a_nodes[i];

  // Bubble sort of the node heights
  do
    {
      swap = NO;
      For(i,2*tree->n_otu-2) 
        {
          if(tree->rates->nd_t[t->sorted_nd[i+1]->num] < tree->rates->nd_t[t->sorted_nd[i]->num])
            {
              buff              = t->sorted_nd[i];
              t->sorted_nd[i]   = t->sorted_nd[i+1];
              t->sorted_nd[i+1] = buff;

              swap = YES;
            }
        }
    }
  while(swap == YES);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Update the set of vectors of occupation along the tree
void GEO_Update_Occup(t_geo *t, t_tree *tree)
{
  int i,j;
  t_node *v1, *v2;

  GEO_Update_Sorted_Nd(t,tree);

  For(i,t->ldscape_sz*(2*tree->n_otu-1)) t->occup[i] = 0;

  t->occup[tree->n_root->num*t->ldscape_sz + t->loc[tree->n_root->num]] = 1;
  
  for(i=1;i<2*tree->n_otu-1;i++)
    {
      For(j,t->ldscape_sz) 
        {
          t->occup[t->sorted_nd[i]->num*t->ldscape_sz + j] = 
            t->occup[t->sorted_nd[i-1]->num*t->ldscape_sz + j];
        }

      
      if(t->sorted_nd[i-1]->tax == NO)
        {
          v1 = v2 = NULL;
          if(t->sorted_nd[i-1] == tree->n_root)
            {
              v1 = tree->n_root->v[1];
              v2 = tree->n_root->v[2];
            }
          else
            {
              For(j,3)
                {
                  if(t->sorted_nd[i-1]->v[j] != t->sorted_nd[i-1]->anc &&
                     t->sorted_nd[i-1]->b[j] != tree->e_root)
                    {
                      if(!v1) v1 = t->sorted_nd[i-1]->v[j];
                      else    v2 = t->sorted_nd[i-1]->v[j];
                    }
                }
            }

          
          if(t->loc[v1->num] != t->loc[t->sorted_nd[i-1]->num])
            {
              t->occup[t->sorted_nd[i]->num * t->ldscape_sz + t->loc[v1->num]]++;
            }
          else
            {
              t->occup[t->sorted_nd[i]->num * t->ldscape_sz + t->loc[v2->num]]++;
            }
        }
    }

  /* printf("\n"); */
  /* For(i,2*tree->n_otu-1) */
  /*   { */
  /*     printf("\n. Node %3d: ",t->sorted_nd[i]->num); */
  /*     For(j,t->ldscape_sz) */
  /*       { */
  /*         printf("%3d [%12f;%12f]   ", */
  /*                t->occup[t->sorted_nd[i]->num*t->ldscape_sz + j], */
  /*                t->ldscape[j*t->n_dim+0],t->ldscape[j*t->n_dim+1]); */
  /*       } */
  /*   } */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Calculate R mat at node n
void GEO_Update_Rmat(t_node *n, t_geo *t, t_tree *tree)
{
  int i,j;
  phydbl lbda_j;

  GEO_Update_Fmat(t);

  For(i,t->ldscape_sz)
    {
      For(j,t->ldscape_sz)
        {
          lbda_j = ((t->occup[n->num*t->ldscape_sz + j]==0) ? (1.0) : (t->lbda));
          t->r_mat[i*t->ldscape_sz+j] = t->f_mat[i*t->ldscape_sz+j] * lbda_j * t->tau;          
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl GEO_Lk(t_geo *t, t_tree *tree)
{
  int i;
  phydbl loglk;
  phydbl R;
  int dep,arr; // departure and arrival location indices;
  t_node *curr_n,*prev_n;
  phydbl sum;

  GEO_Update_Occup(t,tree);     // Same here.

  prev_n = NULL;
  curr_n = NULL;
  loglk = .0;
  for(i=1;i<tree->n_otu-1;i++) // Consider all the time slices, from oldest to youngest. 
                             // Start at first node below root
    {
      prev_n = t->sorted_nd[i-1]; // node just above
      curr_n = t->sorted_nd[i];   // current node
      
      GEO_Update_Rmat(curr_n,t,tree); // NOTE: don't need to do that every time. Add check later.

      R = GEO_Total_Migration_Rate(curr_n,t); // Total migration rate calculated at node n

      dep = t->loc[curr_n->num]; // departure location
      arr =                      // arrival location
        (t->loc[curr_n->v[1]->num] == t->loc[curr_n->num] ? 
         t->loc[curr_n->v[2]->num] : 
         t->loc[curr_n->v[1]->num]);
      
      /* printf("\n%f\t%f", */
      /*        t->ldscape[arr*t->n_dim+0]-t->ldscape[dep*t->n_dim+0], */
      /*        t->ldscape[arr*t->n_dim+1]-t->ldscape[dep*t->n_dim+1]); */

      loglk -= R * FABS(tree->rates->nd_t[curr_n->num] - tree->rates->nd_t[prev_n->num]);
      loglk += LOG(t->r_mat[dep * t->ldscape_sz + arr]);

      /* printf("\n. R = %f r_mat = %f dt = %f loglk = %f", */
      /*        R, */
      /*        t->r_mat[dep * t->ldscape_sz + arr], */
      /*        FABS(t->sorted_nd[i] - t->sorted_nd[i-1]),loglk); */

    }


  // Likelihood for the first 'slice' (i.e., the part just below the root down to
  // the next node)
  GEO_Update_Rmat(tree->n_root,t,tree);

  loglk -= LOG(t->ldscape_sz); 
  dep = t->loc[tree->n_root->num];
  arr = 
    (t->loc[tree->n_root->num] != t->loc[tree->n_root->v[1]->num] ? 
     t->loc[tree->n_root->v[1]->num] :
     t->loc[tree->n_root->v[2]->num]);
  
  /* printf("\n %f %f",t->ldscape[dep],t->ldscape[arr]); */

  loglk += LOG(t->r_mat[dep * t->ldscape_sz + arr]);
    
  sum = .0;
  For(i,t->ldscape_sz) sum += t->r_mat[dep * t->ldscape_sz + i];
  
  loglk -= LOG(sum);

  tree->geo->c_lnL = loglk;

  return loglk;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Init_Tloc_Tips(t_geo *t, t_tree *tree)
{
  int i;

  // TO DO
  For(i,tree->n_otu)
    {
      t->loc[i] = i;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Do not forget to call GEO_Update_Rmat (with node n) before calling this function
phydbl GEO_Total_Migration_Rate(t_node *n, t_geo *t)
{
  phydbl R;
  int i,j;

  R = .0;

  For(i,t->ldscape_sz)
    {
      For(j,t->ldscape_sz)
        {
          R += 
            t->r_mat[i * t->ldscape_sz + j] * 
            t->occup[n->num * t->ldscape_sz + i];
        }
    }

  return R;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Find the arrival location for the migration leaving from n
int GEO_Get_Arrival_Location(t_node *n, t_geo *t, t_tree *tree)
{
  int i;
  t_node *v1, *v2; // the two daughters of n

  v1 = v2 = NULL;

  For(i,3)
    {
      if(n->v[i] && n->v[i] != n->anc)
        {
          if(!v1) v1 = n->v[i];
          else    v2 = n->v[i];
        }
    }

  if(t->loc[v1->num] == t->loc[v2->num]) // Migrated to the same location as that of n
    {
      if(t->loc[n->num] != t->loc[v1->num])
        {
          PhyML_Printf("\n== Error detected in location labeling.");
          PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
          Exit("\n");    
        }
      else
        return t->loc[n->num];
    }
  else // Migrated to a different spot
    {
      if((t->loc[v1->num] != t->loc[n->num]) && (t->loc[v2->num] != t->loc[n->num]))
        {
          PhyML_Printf("\n== Error detected in location labeling.");
          PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
          Exit("\n");    
        }
      else
        {
          if(t->loc[v1->num] == t->loc[n->num]) return t->loc[v2->num]; // v2 gets the new location, v1 inheritates from n
          else                                  return t->loc[v1->num]; // v1 gets the new location, v2 inheritates from n
        }
    }
  return -1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_tree *GEO_Simulate(t_geo *t, int n_otu)
{
  t_tree *tree;
  int n_branching_nodes;
  t_node **branching_nodes; // vector of nodes available for branching out
  phydbl *p_branch; // p_branch[i]: proba that the i-th node in branching_nodes branches out (p_x vector in the article)
  phydbl *p_mig; // p_branch[i]: proba of migrating to location i from the location of the edge branching out (q_i vector in the article)
  int hit;
  phydbl time;
  int dep, arr;
  int i,j;
  phydbl sum;
  phydbl R;
  int *occup; // occupation vector. Updated as we move down the tree
  int nd_idx;
  t_node *buff_nd;
  phydbl buff_t;
  int buff_l;
  int swap;


  tree = Make_Tree_From_Scratch(n_otu,NULL);
  tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
  RATES_Init_Rate_Struct(tree->rates,NULL,tree->n_otu);
  tree->n_root = tree->a_nodes[2*tree->n_otu-2]; // Set the root node to the last element in the list of nodes
  tree->geo = t;

  For(i,2*tree->n_otu-2) tree->rates->nd_t[i] = -1.;

  occup = (int *)mCalloc(t->ldscape_sz,sizeof(int));
  
  GEO_Update_Fmat(t);

  branching_nodes = (t_node **)mCalloc(tree->n_otu,sizeof(t_node *));
  branching_nodes[0] = tree->n_root;
  n_branching_nodes  = 1;
  nd_idx = 0;

  
  p_branch = (phydbl *)mCalloc(tree->n_otu,sizeof(phydbl ));
  p_branch[0] = 1.0;

  p_mig = (phydbl *)mCalloc(t->ldscape_sz,sizeof(phydbl ));

  time = 0.0; // Time at the root node (this will be changed afterward)

  // Sample a location uniformly for the root
  t->loc[tree->n_root->num] = Rand_Int(0,t->ldscape_sz-1);
  
  // Update the occupancy vector
  occup[t->loc[tree->n_root->num]] = 1;

  dep = arr = -1;

 // total migration rate
  R = 0.0;
  For(i,t->ldscape_sz) 
    {
      R += 
        t->f_mat[t->loc[tree->n_root->num]*t->ldscape_sz+i] * 
        ((occup[i] == 0) ? (1.0) : (t->lbda)) * 
        t->tau;
    }

  do
    {      
      // Select the node that branches out
      hit = Sample_i_With_Proba_pi(p_branch,n_branching_nodes);
      
      /* printf("\n. [%d] Select node %d (location %d)",n_branching_nodes,branching_nodes[hit]->num,t->loc[branching_nodes[hit]->num]); */

      // Set the time for the branching node
      tree->rates->nd_t[branching_nodes[hit]->num] = time;


      /* printf("\n. Set its time to %f",time); */

      // Select the destination location
      dep = t->loc[branching_nodes[hit]->num]; // Departure point
           
      sum = .0;
      For(i,t->ldscape_sz) // Total rate of migration out of departure point
        {
          p_mig[i] = 
            t->f_mat[dep*t->ldscape_sz+i] * 
            ((occup[i] == 0) ? (1.0) : (t->lbda)) * 
            t->tau;

          sum += p_mig[i];
        }      
      For(i,t->ldscape_sz) p_mig[i] /= sum;

      arr = Sample_i_With_Proba_pi(p_mig,t->ldscape_sz);

      printf("\n. Migrate from %d [%5.2f,%5.2f] to %d [%5.2f,%5.2f]",
             dep,
             t->ldscape[dep*t->n_dim+0],
             t->ldscape[dep*t->n_dim+1],
             arr,
             t->ldscape[arr*t->n_dim+0],
             t->ldscape[arr*t->n_dim+1]);

      /* printf("\n%f\t%f", */
      /*        t->ldscape[arr*t->n_dim+0]-t->ldscape[dep*t->n_dim+0], */
      /*        t->ldscape[arr*t->n_dim+1]-t->ldscape[dep*t->n_dim+1]); */
             


      // Update vector of occupation
      occup[arr]++;
      
      /* printf("\n. Remove %d. Add %d and %d",branching_nodes[hit]->num,tree->a_nodes[nd_idx]->num,tree->a_nodes[nd_idx+1]->num); */
      // Connect two new nodes to the node undergoing a branching event
      tree->a_nodes[nd_idx]->v[0]   = branching_nodes[hit];
      tree->a_nodes[nd_idx+1]->v[0] = branching_nodes[hit];
      branching_nodes[hit]->v[1] = tree->a_nodes[nd_idx];
      branching_nodes[hit]->v[2] = tree->a_nodes[nd_idx+1];

      // update branching_nodes vector. Element 'hit' is being replaced so that the corresponding node can no longer branch out
      branching_nodes[hit]                 = tree->a_nodes[nd_idx];
      branching_nodes[n_branching_nodes]   = tree->a_nodes[nd_idx+1];

      // Update t_loc vector.
      t->loc[tree->a_nodes[nd_idx]->num]   = dep;
      t->loc[tree->a_nodes[nd_idx+1]->num] = arr;

      // Update total migration rate 
      R = .0;
      For(i,t->ldscape_sz)
        {
          if(occup[i] > 0)
            {
              For(j,t->ldscape_sz)
                {
                  R += 
                    occup[i] *
                    t->f_mat[i*t->ldscape_sz+j] * 
                    ((occup[j] == 0) ? (1.0) : (t->lbda)) *
                    t->tau;
                }
            }
        }

      // Set the time until next branching event
      time = time + Rexp(R);
    
      // Update p_branch vector
      For(i,n_branching_nodes+1)
        {
          dep = t->loc[branching_nodes[i]->num];
          p_branch[i] = 0.0;
          For(j,t->ldscape_sz)
            {
              p_branch[i] +=
                t->f_mat[dep*t->ldscape_sz+j] * 
                ((occup[j] == 0) ? (1.0) : (t->lbda)) * 
                t->tau / R;

              /* printf("\n. %f %f %f %f", */
              /*        R, */
              /*        t->f_mat[dep*t->ldscape_sz+j], */
              /*        ((occup[j]>0) ? (t->lbda) : (1.0)), */
              /*        t->tau); */
            }
          /* printf("\n. %f ",p_branch[i]); */
        }

              
      // Increase the number of branching nodes by one (you just added 2 new and removed 1 old)
      n_branching_nodes++;
      nd_idx += 2;

    }
  while(n_branching_nodes < n_otu);

  printf("\n");
  For(i,t->ldscape_sz)
    printf("%d,",occup[i]);
  printf("\n");


  // Set the times at the tips
  For(i,2*tree->n_otu-1) if(tree->rates->nd_t[i] < 0.0) tree->rates->nd_t[i] = time;
  
  // Reverse time scale
  For(i,2*tree->n_otu-1) tree->rates->nd_t[i] -= time;
  /* For(i,2*tree->n_otu-1) tree->rates->nd_t[i] = FABS(tree->rates->nd_t[i]); */

  //  Bubble sort to put all the tips at the top of the tree->a_nodes array
  do
    {
      swap = NO;
      For(i,2*tree->n_otu-2)
        {
          if(!tree->a_nodes[i+1]->v[1] && tree->a_nodes[i]->v[1])
            {
              buff_nd            = tree->a_nodes[i+1];
              tree->a_nodes[i+1] = tree->a_nodes[i];
              tree->a_nodes[i]   = buff_nd;

              buff_t                 = tree->rates->nd_t[i+1];
              tree->rates->nd_t[i+1] = tree->rates->nd_t[i];
              tree->rates->nd_t[i]   = buff_t;

              buff_l      = t->loc[i+1];
              t->loc[i+1] = t->loc[i];
              t->loc[i]   = buff_l;

              swap = YES;
            }
        }
    }
  while(swap == YES);


  // The rest below is just bookeeping...


  For(i,2*tree->n_otu-1) tree->a_nodes[i]->num = i;
  For(i,2*tree->n_otu-1) 
    {
      if(i < tree->n_otu) tree->a_nodes[i]->tax = YES;
      else                tree->a_nodes[i]->tax = NO;
    }

  /* printf("\n++++++++++++++++++\n"); */
  /* For(i,2*tree->n_otu-1) */
  /*   { */
  /*     printf("\n. Node %3d [%p] anc:%3d v1:%3d v2:%3d time: %f", */
  /*            tree->a_nodes[i]->num, */
  /*            (void *)tree->a_nodes[i], */
  /*            tree->a_nodes[i]->v[0] ? tree->a_nodes[i]->v[0]->num : -1, */
  /*            tree->a_nodes[i]->v[1] ? tree->a_nodes[i]->v[1]->num : -1, */
  /*            tree->a_nodes[i]->v[2] ? tree->a_nodes[i]->v[2]->num : -1, */
  /*            tree->rates->nd_t[i]);              */
    /* } */


  For(i,tree->n_otu) 
    {
      if(!tree->a_nodes[i]->name) tree->a_nodes[i]->name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      strcpy(tree->a_nodes[i]->name,"x");
      sprintf(tree->a_nodes[i]->name+1,"%d",i);      
    }

  tree->n_root->v[1]->v[0] = tree->n_root->v[2];
  tree->n_root->v[2]->v[0] = tree->n_root->v[1];
  
  tree->num_curr_branch_available = 0;
  Connect_Edges_To_Nodes_Recur(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);

  tree->e_root = tree->n_root->v[1]->b[0];

  For(i,2*tree->n_otu-3)
    {
      tree->a_edges[i]->l->v = FABS(tree->rates->nd_t[tree->a_edges[i]->left->num] - 
                                    tree->rates->nd_t[tree->a_edges[i]->rght->num]);
    }

  tree->e_root->l->v =
    FABS(tree->rates->nd_t[tree->n_root->v[1]->num] -
         tree->rates->nd_t[tree->n_root->num]) +
    FABS(tree->rates->nd_t[tree->n_root->v[2]->num] -
         tree->rates->nd_t[tree->n_root->num]);

  tree->n_root_pos = 
    FABS(tree->rates->nd_t[tree->n_root->v[2]->num] -
         tree->rates->nd_t[tree->n_root->num]) / tree->e_root->l->v;

  /* printf("\n. %s ",Write_Tree(tree,NO)); */

  DR_Draw_Tree("essai.ps",tree);

  /* For(i,tree->n_otu) */
  /*   printf("\n. %4s %4d [%5.2f %5.2f]",tree->a_nodes[i]->name, */
  /*          t->loc[i], */
  /*          t->ldscape[t->loc[i]*t->n_dim+0], */
  /*          t->ldscape[t->loc[i]*t->n_dim+1]); */

  
  Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);		

  Free(branching_nodes);
  Free(p_branch);
  Free(p_mig);

  return(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Simualte n coordinates (in 2D)
void GEO_Simulate_Coordinates(int n, t_geo *t)
{
  int i;
  phydbl width;

  width = 5.;

  For(i,n)
    {
      t->ldscape[i*t->n_dim+0] = -width/2. + Uni()*width;
      t->ldscape[i*t->n_dim+1] = -width/2. + Uni()*width;
    }


  /* t->ldscape[0*t->n_dim+0] = 0.0; */
  /* t->ldscape[0*t->n_dim+1] = 0.0; */

  /* t->ldscape[1*t->n_dim+0] = 0.1; */
  /* t->ldscape[1*t->n_dim+1] = 0.1; */

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Optimize_Sigma(t_geo *t, t_tree *tree)
{
  Generic_Brent_Lk(&(t->sigma),
                   t->min_sigma,
                   t->max_sigma,
                   1.E-5,
                   1000,
                   NO,
                   GEO_Wrap_Lk,NULL,tree,NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Optimize_Lambda(t_geo *t, t_tree *tree)
{
  Generic_Brent_Lk(&(t->lbda),
                   t->min_lbda,
                   t->max_lbda,
                   1.E-5,
                   1000,
                   NO,
                   GEO_Wrap_Lk,NULL,tree,NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Optimize_Tau(t_geo *t, t_tree *tree)
{
  Generic_Brent_Lk(&(t->tau),
                   t->min_tau,
                   t->max_tau,
                   1.E-5,
                   1000,
                   NO,
                   GEO_Wrap_Lk,NULL,tree,NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl GEO_Wrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree)
{
  return GEO_Lk(tree->geo,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Init_Geo_Struct(t_geo *t)
{
  t->c_lnL        = UNLIKELY;

  t->sigma        = 1.0;
  t->min_sigma    = 1.E-3;
  t->max_sigma    = 10.;
  t->sigma_thresh = t->max_sigma;

  t->lbda         = 1.0;
  t->min_lbda     = 1.E-2;
  t->max_lbda     = 1.E+1;
  
  t->tau          = 1.0;
  t->min_tau      = 1.E-3;
  t->max_tau      = 1.E+1;

  t->tau          = 1.0;

  t->n_dim        = -1;
  t->ldscape_sz   = -1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Randomize_Locations_Pre(t_node *n, t_geo *t, t_tree *tree)
{

  if(n->tax == YES) return;
  else
    {
      int i;
      
      GEO_Randomize_Locations(n,t,tree);

      For(i,3)
        if(n->v[i] != n->anc && n->b[i] != tree->e_root) 
          GEO_Randomize_Locations_Pre(n->v[i],t,tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Randomize_Locations(t_node *n, t_geo *t, t_tree *tree)
{
  t_node *v1, *v2;
  int i;
  phydbl *probs; // vector of probability of picking each location
  phydbl sum;
  
  probs = (phydbl *)mCalloc(t->ldscape_sz,sizeof(phydbl));

  v1 = v2 = NULL;
  For(i,3)
    {
      if(n->v[i] != n->anc && n->b[i] != tree->e_root)
        {
          if(!v1) v1 = n->v[i];
          else    v2 = n->v[i];
        }
    }

  if(t->loc[v1->num] != t->loc[n->num] && v1->tax == NO)
    {
      sum = 0.0;
      For(i,t->ldscape_sz) sum += t->loc_beneath[v1->num * t->ldscape_sz + i];
      For(i,t->ldscape_sz) probs[i] = t->loc_beneath[v1->num * t->ldscape_sz + i]/sum;
      
      t->loc[v1->num] = Sample_i_With_Proba_pi(probs,t->ldscape_sz);      
    }
  else if(t->loc[v2->num] != t->loc[n->num] && v2->tax == NO)
    {
      sum = 0.0;
      For(i,t->ldscape_sz) sum += t->loc_beneath[v2->num * t->ldscape_sz + i];
      For(i,t->ldscape_sz) probs[i] = t->loc_beneath[v2->num * t->ldscape_sz + i]/sum;
      
      Sample_i_With_Proba_pi(probs,t->ldscape_sz);

      t->loc[v2->num] = Sample_i_With_Proba_pi(probs,t->ldscape_sz);      
    }

  Free(probs);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Get_Locations_Beneath(t_geo *t, t_tree *tree)
{
  int i;

  GEO_Get_Locations_Beneath_Post(tree->n_root,tree->n_root->v[1],t,tree);
  GEO_Get_Locations_Beneath_Post(tree->n_root,tree->n_root->v[2],t,tree);

  For(i,t->ldscape_sz)
    t->loc_beneath[tree->n_root->num*t->ldscape_sz+i] =
    t->loc_beneath[tree->n_root->v[1]->num*t->ldscape_sz+i] +
    t->loc_beneath[tree->n_root->v[2]->num*t->ldscape_sz+i] ;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Get_Locations_Beneath_Post(t_node *a, t_node *d, t_geo *t, t_tree *tree)
{

  if(d->tax) 
    {
      t->loc_beneath[d->num*t->ldscape_sz+t->loc[d->num]] = 1;
      return;
    }
  else
    {
      int i;
      t_node *v1, *v2;

      For(i,3) 
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            {
              GEO_Get_Locations_Beneath_Post(d,d->v[i],t,tree);
            }
        }
          

      v1 = v2 = NULL;
      For(i,3) 
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            {
              if(!v1) v1 = d->v[i];
              else    v2 = d->v[i];
            }
        }
          
      For(i,t->ldscape_sz) 
        t->loc_beneath[ d->num*t->ldscape_sz+i] = 
        t->loc_beneath[v1->num*t->ldscape_sz+i] + 
        t->loc_beneath[v2->num*t->ldscape_sz+i] ;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Get_Sigma_Max(t_geo *t)
{
  int i,j;
  phydbl max_dist,dist,inv_max_dist;
  phydbl sigma_a, sigma_b, sigma_c;
  phydbl overlap_a, overlap_b, overlap_c;
  phydbl d_intersect;
  phydbl overlap_target;
  phydbl eps;
  int n_iter,n_iter_max;

  eps = 1.E-6;
  overlap_target = 0.95;
  n_iter_max = 100;

  dist = .0;
  max_dist = -1.;
  inv_max_dist = -1.;
  For(i,t->ldscape_sz-1)
    {
      for(j=i+1;j<t->ldscape_sz;j++)
        {
          dist = POW(t->ldscape[i*t->n_dim+0] - t->ldscape[j*t->n_dim+0],1);
          if(dist > max_dist) max_dist = dist;
          dist = POW(t->ldscape[i*t->n_dim+1] - t->ldscape[j*t->n_dim+1],1);
          if(dist > max_dist) max_dist = dist;
        }
    }
  
  inv_max_dist = 1./max_dist;
  


  sigma_a = t->min_sigma; sigma_b = 1.0; sigma_c = t->max_sigma;
  /* sigma_a = t->min_sigma; sigma_b = 1.0; sigma_c = 10.; */
  n_iter = 0;
  do
    {
      d_intersect = Inverse_Truncated_Normal(inv_max_dist,0.0,sigma_a,0.0,max_dist);
      overlap_a = 
        (Pnorm(max_dist,0.0,sigma_a) - Pnorm(d_intersect,0.0,sigma_a))/
        (Pnorm(max_dist,0.0,sigma_a) - Pnorm(0.0,0.0,sigma_a)) + 
        d_intersect / max_dist;
      /* printf("\n. inter: %f %f [%f]",d_intersect,max_dist,d_intersect / max_dist); */

      d_intersect = Inverse_Truncated_Normal(inv_max_dist,0.0,sigma_b,0.0,max_dist);
      overlap_b = 
        (Pnorm(max_dist,0.0,sigma_b) - Pnorm(d_intersect,0.0,sigma_b))/
        (Pnorm(max_dist,0.0,sigma_b) - Pnorm(0.0,0.0,sigma_b)) + 
        d_intersect / max_dist;
      /* printf("\n. inter: %f %f [%f]",d_intersect,max_dist,d_intersect / max_dist); */
      
      d_intersect = Inverse_Truncated_Normal(inv_max_dist,0.0,sigma_c,0.0,max_dist);
      overlap_c = 
        (Pnorm(max_dist,0.0,sigma_c) - Pnorm(d_intersect,0.0,sigma_c))/
        (Pnorm(max_dist,0.0,sigma_c) - Pnorm(0.0,0.0,sigma_c)) + 
        d_intersect / max_dist;

      /* printf("\n. inter: %f %f [%f]",d_intersect,max_dist,d_intersect / max_dist); */
      
      /* printf("\n. sigma_a:%f overlap_a:%f sigma_b:%f overlap_b:%f sigma_c:%f overlap_c:%f", */
      /*        sigma_a,overlap_a, */
      /*        sigma_b,overlap_b, */
      /*        sigma_c,overlap_c); */

      if(overlap_target > overlap_a && overlap_target < overlap_b)
        {
          sigma_c = sigma_b;
          sigma_b = sigma_a + (sigma_c - sigma_a)/2.;
        }
      else if(overlap_target > overlap_b && overlap_target < overlap_c)
        {
          sigma_a = sigma_b;
          sigma_b = sigma_a + (sigma_c - sigma_a)/2.;
        }
      else if(overlap_target < overlap_a)
        {
          sigma_a /= 2.;
        }
      else if(overlap_target > overlap_c)
        {
          sigma_c *= 2.;
        }

      n_iter++;

    }
  while(sigma_c - sigma_a > eps && n_iter < n_iter_max);

  /* if(sigma_c - sigma_a > eps) */
  /*   { */
  /*     PhyML_Printf("\n== Error detected in getting maximum value of sigma."); */
  /*     PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__); */
  /*     Exit("\n");     */
  /*   } */
  /* else */
  /*   { */
  /*     PhyML_Printf("\n== Threshold for sigma: %f",sigma_b); */
  /*   } */

  t->sigma_thresh = sigma_b;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Geo_Updown_Tau_Lbda(t_tree *tree)
{
  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL,new_lnL;
  phydbl cur_tau,new_tau;
  phydbl cur_lbda,new_lbda;
  
  K        = tree->mcmc->tune_move[tree->mcmc->num_move_geo_updown_tau_lbda];
  cur_lnL  = tree->geo->c_lnL;
  new_lnL  = tree->geo->c_lnL;
  cur_tau  = tree->geo->tau;
  new_tau  = tree->geo->tau;
  cur_lbda = tree->geo->lbda;
  new_lbda = tree->geo->lbda;

  u = Uni();
  mult = EXP(K*(u-0.5));

  /* Multiply tau by K */
  new_tau = cur_tau * K;
  
  /* Divide lbda by same amount */
  new_lbda = cur_lbda / K;


  if(
     new_lbda < tree->geo->min_lbda || new_lbda > tree->geo->max_lbda ||
     new_tau  < tree->geo->min_tau  || new_tau  > tree->geo->max_tau
     )
    {
      tree->mcmc->run_move[tree->mcmc->num_move_geo_updown_tau_lbda]++;
      return;
    }
  
  tree->geo->tau  = new_tau;
  tree->geo->lbda = new_lbda;

  if(tree->mcmc->use_data) new_lnL = GEO_Lk(tree->geo,tree);

  ratio = 0.0;
  /* Proposal ratio: 2n-2=> number of multiplications, 1=>number of divisions */
  ratio += 0.0*LOG(mult); /* (1-1)*LOG(mult); */
  /* Likelihood density ratio */
  ratio += (new_lnL - cur_lnL);

  /* printf("\n. new_tau: %f new_lbda:%f cur_lnL:%f new_lnL:%f",new_tau,new_lbda,cur_lnL,new_lnL); */


  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  if(u > alpha) /* Reject */
    {
      tree->geo->tau   = cur_tau;
      tree->geo->lbda  = cur_lbda;
      tree->geo->c_lnL = cur_lnL;
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_geo_updown_tau_lbda]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_geo_updown_tau_lbda]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


