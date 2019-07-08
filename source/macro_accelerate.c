#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_blas.h>
#include "my_linalg.h"


/*******/

void
calc_matom_matrix(xplasma, matom_matix)
     PlasmaPtr xplasma;
     double *matom_matrox;
{
  MacroPtr mplasma;
  double t_e, ne;
  int nbbd, nbbu, nbfd, nbfu;

  mplasma = &macromain[xplasma->nplasma]; //telling us where in the matom structure we are

  t_e = xplasma->t_e;           //electron temperature 
  ne = xplasma->ne;             //electron number density


  for (uplvl = 0; uplvl < nlevels_macro; uplvl++)
    {
      nbbd = config[uplvl].n_bbd_jump;    //store these for easy access -- number of bb downward jumps
      nbbu = config[uplvl].n_bbu_jump;    // number of bb upward jump from this configuration
      nbfd = config[uplvl].n_bfd_jump;    // number of bf downward jumps from this transition
      nbfu = config[uplvl].n_bfu_jump;    // number of bf upward jumps from this transiion


      /* bb*/
      for (n = 0; n < nbbd; n++)
	{
   
	  line_ptr = &line[config[uplvl].bbd_jump[n]];

	  rad_rate = (a21 (line_ptr) * p_escape (line_ptr, xplasma));
	  coll_rate = q21 (line_ptr, t_e);        // this is multiplied by ne below
	  
	  if (coll_rate < 0)
	    {
	      coll_rate = 0;
	    }

	  bb_cont = rad_rate + (coll_rate * ne);
      
	  target_level = line_ptr->nconfigl;
      
	  Q_matrix[uplvl][target_level] = Qcont = bb_cont * config[target_level].ex;     //energy of lower state    
	  R_matrix[uplvl][uplvl] = Rcont = bb_cont * (config[uplvl].ex - config[target_level].ex);  //energy difference
	  
	  Q_norm[uplvl] += Qcont+Rcont;
	}

      /* bf */
      for (n = 0; n < nbfd; n++)
	{
	  
	  cont_ptr = &phot_top[config[uplvl].bfd_jump[n]];        //pointer to continuum
	  if (n < 25) //?
	    {
	      sp_rec_rate = mplasma->recomb_sp[config[uplvl].bfd_indx_first + n];   //need this twice so store it
	      bf_cont = (sp_rec_rate + q_recomb (cont_ptr, t_e) * ne) * ne;
	    }
	  else
	    {
	      bf_cont = 0.0;
	    }


	target_level = config[phot_top[config[uplvl].bfd_jump[n]].nlev;
			       
        Q_matrix[uplvl][target_level] = Qcont = bf_cont * config[target_level].ex;       //energy of lower state
        R_matrix[uplvl][uplvl] = Rcont = bf_cont * (config[uplvl].ex - config[target_lvl].ex);  //energy difference

	Q_norm[uplvl] += Qcont+Rcont;	
      }


      
    }

  
}
  
