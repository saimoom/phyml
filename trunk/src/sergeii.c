#include "spr.h"
#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "bionj.h"
#include "models.h"
#include "free.h"
#include "help.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"
#include "mixt.h"
#include "sergeii.h"
#ifdef MPI
#include "mpi_boot.h"
#endif
#define MAX 500




int My_main(int argc, char **argv){

  FILE *f;
  int i, j, s, n_taxa, n_clade, *clade_size, *node_num, *monitor, n_monitor;
  phydbl *low, *up, *t_prior_min, *t_prior_max;
  xml_node *n;
  t_tree *tree;
  t_rate *rates;
  t_mcmc *mcmc;
  option *io;
  align **data;
  char *topology;
  char **seq_name_list, **seq_list, ***clade, **bound_l, **bound_u, **clade_name, **monitor_list;
  i = 0;
  j = 0;
  s = 0; 
 
  //memory for list of clades to be monitored:
  monitor_list = (char **)mCalloc(MAX,sizeof(char *));
  For(i, MAX)
    {
      monitor_list[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    }

  //memory for options:
  io = (option *)Make_Input();

  //memory for data:
  data = (align **)mCalloc(MAX,sizeof(align *));
  For(i, MAX)
    {
      data[i] = (align *)mCalloc(MAX,sizeof(align));
    }

  //memory for tree:
  tree = (t_tree *)mCalloc(MAX,sizeof(t_tree));
  
  //memory for mcmc:
  mcmc = (t_mcmc *)mCalloc(MAX,sizeof(t_mcmc));

  //memory for rates:
  rates = (t_rate *)mCalloc(MAX,sizeof(t_rate));

  //memory for clade size:
  clade_size = (int *)mCalloc(MAX,sizeof(int));

  //memory for nodes numbers:
  node_num = (int *)mCalloc(MAX,sizeof(int));
  
  //memory for topology:
  topology = (char *)mCalloc(T_MAX_SEQ,sizeof(char));

  //memory for bounds:
  t_prior_min = (phydbl *)mCalloc(MAX,sizeof(phydbl));
  t_prior_max = (phydbl *)mCalloc(MAX,sizeof(phydbl));

  //memory for sequences:
  seq_name_list = (char **)mCalloc(MAX,sizeof(char *));
  For(i, MAX)
    {
      seq_name_list[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    }
 
  seq_list = (char **)mCalloc(MAX,sizeof(char *));
  For(i, MAX)
    {
      seq_list[i] = (char *)mCalloc(T_MAX_SEQ,sizeof(char));
    }
  
  //memory for monitor flag:
  monitor = (int *)mCalloc(MAX,sizeof(int));
  
  //memory for clade name:
  clade_name = (char **)mCalloc(MAX,sizeof(char *));
  For(j, MAX)
    {
      clade_name[j] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    }
 
  //memory for calibration:
  clade = (char ***)mCalloc(MAX,sizeof(char **));
  For(j, MAX)
    {
      clade[j] = (char **)mCalloc(MAX,sizeof(char *));
      For(i, MAX)
	{
	  clade[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
	}
    }

  bound_l = (char **)mCalloc(MAX,sizeof(char *));
  For(i, MAX)
    {
      bound_l[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    }

  bound_u = (char **)mCalloc(MAX,sizeof(char *));
  For(i, MAX)
    {
      bound_u[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    }

  //memory for bounds (read from file):
  low = (phydbl *)mCalloc(MAX,sizeof(phydbl));
  up = (phydbl *)mCalloc(MAX,sizeof(phydbl));
  
  //file can/cannot be open:
  if ((f =(FILE *)fopen(argv[1], "r")) == NULL)
    {
      PhyML_Printf("\n==File can not be open...\n");
      Exit("\n");
      return 0;
    }

  io -> data = data;
  n = XML_Load_File(f);
  n = n -> child;
  i = 0;
  j = 0;
  n_clade = 0;

  do
    {
      if(!strcmp(n -> name, "topology"))//looking for a node <topology>.
	{
	  strcpy(topology, n -> value);
      	  n = n -> next;
	}
     
      else if(!strcmp(n -> name, "alignment"))//looking for a node <alignment>.
	{
	  i = 0;
	  do
	    {
	      strcpy(seq_name_list[i], n -> child -> attr -> value);
	      strcpy(seq_list[i], n -> child -> value);
	      i++;
	      if(n -> child -> next) n -> child = n -> child -> next;
	      else break;
	    }
	  while(n -> child);
	  n_taxa = i;
	  i = 1;
	  For(i, n_taxa) if(strlen(seq_list[0]) != strlen(seq_list[i]))
		{
		  printf("\n. Sequences are of different length. Please check your data...\n");
		  break;
		}
     
	  i = 0;
	  n = n -> next;
	}
      else if(!strcmp(n -> name, "clade"))//looking for a node <clade>.
	{
	  strcpy(clade_name[j], n -> attr -> value);
	  i = 0;
	  do
	    {
	      clade[j][i] = n -> child -> attr -> value;
	      i++;
	      if(n -> child -> next) n -> child = n -> child -> next;
	      else break;
	    }
	  while(n -> child);
	  clade_size[j] = i;
	  strcpy(bound_l[j], n -> attr -> next -> value);
	  strcpy(bound_u[j],  n -> attr -> next -> next -> value);
	  j++;
	  i = 0;
	  n_clade++;
	  n = n -> next;	  
	}
      else if(!strcmp(n -> name, "monitor"))//looking for a node <monitor>.
	{
	  do
	    {
	      strcpy(monitor_list[i], n -> child -> attr -> value);
	      i++;	      
	      if(n -> child -> next) n -> child = n -> child -> next;
	      else break;
	    }
	  while(n -> child);
	  n_monitor = i;
	  n = n -> next;
	}
      else
	n = n -> next;
    }
  while(n); 
  
  i = 0;
  j = 0;  

  /*
  PhyML_Printf("\n. Topology:'%s'\n", topology);
  For(i, n_taxa)
    {
      PhyML_Printf("\n. Seq_name:'%s' Sequence:'%s'\n", seq_name_list[i], seq_list[i]);
    }

  For(j, n_clade)
    {
      For(i, clade_size[j])
	{
	  PhyML_Printf("\n. Clade name:'%s' Clade:'%s'\n", clade_name[j], clade[j][i]);
	}
      PhyML_Printf("\n. Lower:'%s' Upper:'%s'\n", bound_l[j], bound_u[j]);
    }
 
  */
  //monitor - setting flags to clades:
  For(j, n_clade){
    For(i, n_monitor){
      if(!strcmp(clade_name[j], monitor_list[i])) monitor[j] = 1;
	 }
    }

  //check if a number of flags is equal to a number of clades to be monitored:
  For(j, n_clade)
    {   
      s = s + monitor[j];      
    }
  if(s != n_monitor)
    {
      PhyML_Printf("\n==Check your data on clades ... \n");
      Exit("\n");
    }
  
  //checking sequence names:
  For(i, n_taxa) Check_Sequence_Name(seq_name_list[i]);
  
  //setting tree:
  tree = Read_Tree(&topology);
  io -> tree = tree;

  //setting alignments to data:
  For(i, n_taxa)
    {
      io -> data[i] -> name = seq_name_list[i];
      io -> data[i] -> state = seq_list[i];
    } 
  
  //setting node numbers:
  For(j, n_clade)
    {
      low[j] = String_To_Dbl(bound_l[j]);
      up[j] = String_To_Dbl(bound_u[j]);
    }

  tree -> mcmc = mcmc;
  mcmc -> monitor = monitor;

  For(i, n_monitor)
    {
      PhyML_Printf("\n. Clades to be monitored:'%s'\n", monitor_list[i]);
    }
 
  //setting upper and lower bounds for clades:
  For(j, n_clade)
    {
      node_num[j] = Find_Clade(clade[j], clade_size[j], tree);
    }
  /*
  //setting priors on bounds:
  For(j, n_clade)
    {
      t_prior_min[node_num[j]] = low[j];
      t_prior_max[node_num[j]] = up[j];
    }

  tree -> rates = rates;
  rates -> t_prior_min = t_prior_min;
  rates -> t_prior_max = t_prior_max;
  */

  For(i, n_monitor)
    {
      PhyML_Printf("\n. Clades to be monitored:'%s'\n", monitor_list[i]);
    }

  For(i, MAX)
    {
      free(monitor_list[i]);
    }
  free(monitor_list);  

  For(i, MAX)
    {
      free(bound_l[i]);
    }
  free(bound_l);

  For(i, MAX)
    {
      free(bound_u[i]);
    }
  free(bound_u);

  return 1;   
}

	





