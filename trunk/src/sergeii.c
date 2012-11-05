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



int PhyTime_XML(int argc, char **argv){

  FILE *f;
  char **clade, *clade_name, **mon_list;
  phydbl low, up;
  int i, j, n_taxa, clade_size, node_num, n_mon;
  xml_node *n_r, *n_t, *n_m, *n_cur;
  t_cal *last_calib, *cur;
  align **data;
  option *io;

  i = 0;
  j = 0;

  last_calib = NULL;

  //file can/cannot be open:
  if ((f =(FILE *)fopen(argv[1], "r")) == NULL)
    {
      PhyML_Printf("\n==File can not be open...\n");
      Exit("\n");
      return 0;
    }

  n_r = XML_Load_File(f);


  ////////////////////////////////////////////////////////////////////////////
  //////////////////////reading tree topology:////////////////////////////////

  //memory for options:
  io = (option *)Make_Input();

  //looking for a node <topology>
  n_t = XML_Search_Node_Name("topology", YES, n_r);

  //setting tree:
  io -> tree = Read_Tree(&n_t -> value);
 
  //memory for rates:
  io -> tree -> rates = RATES_Make_Rate_Struct(io -> tree -> n_otu);
 
  // TO DO: check that the tree has a root...
  Update_Ancestors(io -> tree -> n_root, io -> tree -> n_root -> v[0], io -> tree);
  Update_Ancestors(io -> tree -> n_root, io -> tree -> n_root -> v[1], io -> tree);
		  

  ////////////////////////////////////////////////////////////////////////////
  //////////////////////reading sequences and clades://///////////////////////

  //memory for monitor flag:
  io -> mcmc -> monitor = (int *)mCalloc(2 * io -> tree -> n_otu - 1,sizeof(int));
 
  //memory for sequences:
  n_cur = XML_Search_Node_Name("alignment", YES, n_r);

  data   = (align **)mCalloc(io -> tree -> n_otu,sizeof(align *));
  For(i, io -> tree -> n_otu)
    {
      data[i]          = (align *)mCalloc(1,sizeof(align));
      data[i] -> name  = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      data[i] -> state = (char *)mCalloc(strlen(n_cur -> child -> value) + 1,sizeof(char));
    }
  io -> data = data;
 
  //memory for clade:
  clade_name = (char *)mCalloc(T_MAX_NAME,sizeof(char));

  clade = (char **)mCalloc(T_MAX_FILE,sizeof(char *));
  For(i, T_MAX_FILE)
    {
      clade[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    }

  //memory for list of clades to be monitored
  mon_list = (char **)mCalloc(T_MAX_FILE,sizeof(char *));
  For(i, T_MAX_FILE)
    {
      mon_list[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    }

 
  //reading monitor node:
  i = 0;
  n_m = XML_Search_Node_Name("monitor", YES, n_r);
  if(n_m != NULL)
    {
      do
	{
	  strcpy(mon_list[i], n_m -> child -> attr -> value);
	  i++;
	  if(n_m -> child) n_m -> child = n_m -> child -> next;
	  else break;
	}
      while(n_m -> child);
      n_mon = i;
    }
  else
    {
      n_mon = 0;
      PhyML_Printf("\n. There is no clade to be monitored. \n");
    }

  //chekcing for calibration node (upper or lower bound) to exist:
  n_cur = XML_Search_Node_Name("calibration", YES, n_r);

  if(n_cur == NULL)
    {
      PhyML_Printf("\n==There is no calibration information provided. \n");
      PhyML_Printf("\n==Please check your data. \n");
      Exit("\n");
    }
  else
    {
      if(XML_Search_Node_Name("upper", NO, n_cur -> child) == NULL && XML_Search_Node_Name("lower", NO, n_cur -> child) == NULL)
	{
	  PhyML_Printf("\n==There is no calibration information provided. \n");
	  PhyML_Printf("\n==Please check your data. \n");
	  Exit("\n");
	}
    }   

  n_r = n_r -> child; 

  do
    {
      if(!strcmp(n_r -> name, "alignment"))//looking for a node <alignment>.
	{
	  i = 0;
	  do
	    {
	      strcpy(io -> data[i] -> name, n_r -> child -> attr -> value);
	      strcpy(io -> data[i] -> state, n_r -> child -> value);
	      i++;
	      if(n_r -> child -> next) n_r -> child = n_r -> child -> next;
	      else break;
	    }
	  while(n_r -> child);
	  n_taxa = i;

	  //checking if a sequences of the same lengths:
	  i = 1;
	  For(i, n_taxa) if(strlen(io -> data[0] -> state) != strlen(io -> data[i] -> state))
	    {
	      printf("\n. Sequences are of different length. Please check your data...\n");
	      Exit("\n");
	      break;
	    }  

	  //checking sequence names:
	  i = 0;
	  For(i, n_taxa) Check_Sequence_Name(io -> data[i] -> name);

	  //check if a number of tips is equal to a number of taxa:
	  if(n_taxa != io -> tree -> n_otu)
	    {
	      PhyML_Printf("\n==Number of taxa is not the same as a number of tips. Check your data...\n");
	      Exit("\n");
	    }   

	  n_r = n_r -> next;
	}
      else if(!strcmp(n_r -> name, "calibration"))//looking for a node <calibration>.
	{
	  low = BIG;
	  up  = BIG;
	  n_cur = XML_Search_Node_Name("lower", NO, n_r -> child);
	  if(n_cur != NULL) low = String_To_Dbl( n_cur -> value);
	  n_cur = XML_Search_Node_Name("upper", NO, n_r -> child);
	  if(n_cur != NULL) up = String_To_Dbl(n_cur -> value);
	  do
	    {
	      if(!strcmp("appliesto", n_r -> child -> name)) 
		{
		  strcpy(clade_name, n_r -> child -> child -> attr -> value);//reached clade names
		  if (io -> tree -> rates -> calib == NULL) io -> tree -> rates -> calib = Make_Calib();
		  if(last_calib)
		    {
		      last_calib -> next = io -> tree -> rates -> calib;
		      io -> tree -> rates -> calib -> prev = last_calib;
		    }
		  last_calib = io -> tree -> rates -> calib;
		  io -> tree -> rates -> calib -> lower = low;
		  //PhyML_Printf("\n. '%f'\n", io -> tree -> rates -> calib -> lower);
		  io -> tree -> rates -> calib -> upper = up;
		  //PhyML_Printf("\n. '%f'\n", io -> tree -> rates -> calib -> upper);
		  if(!strcmp("@root@", clade_name))
		    {
		      node_num = io -> tree -> n_root -> num;
		      For(j, n_mon)
			{
			  if(!strcmp("@root@", mon_list[j])) io -> mcmc -> monitor[node_num] = YES;
			}	
		    }
		  else		  
		    {
		      i = 0; 
		      do
			{
			  strcpy(clade[i], n_r -> child -> child -> child -> attr -> value); 
			  i++;
			  if(n_r -> child -> child -> child -> next) n_r -> child -> child -> child = n_r -> child -> child -> child -> next;
			  else break;
			}
		      while(n_r -> child -> child -> child);
		      clade_size = i;
		      node_num = Find_Clade(clade, clade_size, io -> tree);
			{
			  if(!strcmp(clade_name, mon_list[j])) io -> mcmc -> monitor[node_num] = YES;
			}	
		    }

		      io -> tree -> rates -> calib -> node_num = node_num;
		      //PhyML_Printf("\n. '%d'\n", io -> tree -> rates -> calib -> node_num);
		      io -> tree -> rates -> t_prior_min[node_num] = low;
		      io -> tree -> rates -> t_prior_max[node_num] = up;
		      io -> tree -> rates -> t_has_prior[node_num] = YES;
		      io -> tree -> rates -> calib -> proba = String_To_Dbl(n_r -> child -> attr -> value);
		      //PhyML_Printf("\n. '%f'\n", io -> tree -> rates -> calib -> proba);
		      io -> tree -> rates -> calib = io -> tree -> rates -> calib -> next;
 		      if(n_r -> child -> next) n_r -> child = n_r -> child -> next;
		      else break;
		}
	      else if(n_r -> child -> next) n_r -> child = n_r -> child -> next;
	      else break;	      
	    }
	  while(n_r -> child);	   
	  n_r = n_r -> next;

	}
      else if (n_r -> next) n_r = n_r -> next;
      else break;
    }
  while(n_r);

  io -> tree -> rates -> calib = last_calib;
  while(io -> tree -> rates -> calib -> prev) io -> tree -> rates -> calib = io -> tree -> rates -> calib -> prev;
  cur = io -> tree -> rates -> calib;
  /*
  do
    {
      PhyML_Printf("\n '%d' \n", cur -> node_num);
      if(cur -> next) cur = cur -> next;
      else break;
    }
  while(cur);  
  */

  free(clade_name);
  For(i, T_MAX_FILE)
    {
      free(clade[i]);
    }
  free(clade);

  For(i, T_MAX_FILE)
    {
      free(mon_list[i]);
    }
  free(mon_list);

  //PhyML_Printf("\n. '%f' '%f'\n", io -> tree -> rates -> t_prior_min[11], io -> tree -> rates -> t_prior_max[11]);
  //PhyML_Printf("\n. '%f' '%f'\n", io -> tree -> rates -> t_prior_min[18], io -> tree -> rates -> t_prior_max[18]);


  fclose(f);

  return 1;   
}

