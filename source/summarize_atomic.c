/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 
Arguments:		

Returns:
 
Description:	
		
Notes:

History:

**************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"


#include "python.h"
#define NSPEC	20

int
main (argc, argv)
     int argc;
     char *argv[];
{

  char atomic_filename[LINELENGTH];
  int n, m,i;
  int xlines;
  FILE *fptr, *fopen ();
  

  write_atomicdata = 1;
  strcat (atomic_filename, argv[1]);

  get_atomic_data (atomic_filename);

  Log ("The number of elements is %d\n", nelements);
  Log ("The number of ions is %d\n", nions);

  fptr = fopen ("ele.adat", "w");
  fprintf(fptr,"name\tz\tfirstion\tnions\tabun\t\n");
  for (i=0;i<nelements;i++)
  {
	  fprintf(fptr,"%s\t%i\t%i\t%i\t%e\n",ele[i].name,ele[i].z,ele[i].firstion,ele[i].nions,ele[i].abun);
  }
  fclose(fptr);
  
  
  fptr = fopen ("ion.adat", "w");
  fprintf(fptr,"z\tistate\tip\tg\tnmax\tn_lte_max\tfirstlevel\tnlevels\tfirst_nlte_level\tfirst_levden\t\
	  nlte\tphot_info\tmacro_info\tntop_first\tntop_ground\tntop\tnxphot\tlev_type\tdrflag\tnxdrecomb\tcpartflag\t\
	  nxcpart\ttotal_rrflag\tnxtotalrr\tbad_gs_rr_t_flag\tbad_gs_rr_r_flag\tnxbadgsrr\tdere_di_flag\tnxderedi\n");
  for (i=0;i<nions;i++)
  {
	  fprintf(fptr,"%i\t%i\t%e\t%e\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n",
	  ion[i].z,ion[i].istate,ion[i].ip,ion[i].g,ion[i].nmax,ion[i].n_lte_max,ion[i].firstlevel,ion[i].nlevels,ion[i].first_nlte_level,
	  ion[i].first_levden,ion[i].nlte,ion[i].phot_info,ion[i].macro_info,ion[i].ntop_first,ion[i].ntop_ground,ion[i].ntop,ion[i].nxphot,
	  ion[i].lev_type,ion[i].drflag,ion[i].nxdrecomb,ion[i].cpartflag,ion[i].nxcpart,ion[i].total_rrflag,ion[i].nxtotalrr,ion[i].bad_gs_rr_t_flag,
	  ion[i].bad_gs_rr_r_flag,ion[i].nxbadgsrr,ion[i].dere_di_flag,ion[i].nxderedi);
  }
  fclose(fptr);
  
  
  fptr = fopen ("config.adat","w");
	  fprintf(fptr,"z\tistate\tnion\tnden\tisp\tilv\tmacro_info\tn_bbu_jump\tn_bbd_jump\tn_bfu_jump\tn_bfd_jump\tg\tq_num\tex\trad_rate\t\
		  bbu_indx_first\tbfu_indx_first\tbfd_indx_first\n");
  for (i=0;i<nlevels;i++)
  {
	  fprintf(fptr,"%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%e\t%e\t%e\t%e\t%i\t%i\t%i\n",config[i].z,config[i].istate,config[i].nion,
  	  config[i].nden,config[i].isp,config[i].ilv,config[i].macro_info,config[i].n_bbu_jump,config[i].n_bbd_jump,config[i].n_bfu_jump,
	  config[i].n_bfd_jump,config[i].g,config[i].q_num,config[i].ex,config[i].rad_rate,config[i].bbu_indx_first,config[i].bfu_indx_first,
	  config[i].bfd_indx_first);
  }
  fclose(fptr);
  
  fptr=fopen("lines.adat","w");
  fprintf(fptr,"nion\tz\tistate\tgl\tgu\tnconfigl\tnconfigu\tlevl\tlevu\tmacro_info\tfreq\tf\tel\teu\tpow\twhere_in_list\tdown_index\tup_index\n");
  for (i=0;i<nlines;i++)
  {
	  fprintf(fptr,"%i\t%i\t%i\t%e\t%e\t%i\t%i\t%i\t%i\t%i\t%e\t%e\t%e\t%e\t%e\t%i\t%i\t%i\n",
      line[i].nion,line[i].z,line[i].istate,line[i].gl,line[i].gu,line[i].nconfigl,line[i].nconfigu,line[i].levl,line[i].levu,
	  line[i].macro_info,line[i].freq,line[i].f,line[i].el,line[i].eu,line[i].pow,line[i].where_in_list,line[i].down_index,line[i].up_index);
  }
  fclose(fptr);
  
  fptr=fopen("topbase_phot.adat","w");
  fprintf(fptr,"nlev\tuplev\tnion\tz\tistate\tnp\tnlast\tmacro_info\tdown_index\tup_index\tf0\tsigma0\tf\tsigma\n");
  for (i=0;i<nlevels;i++)
  {
	  fprintf(fptr,"%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%e\t%e\t%e\t%e\n",phot_top[i].nlev,phot_top[i].uplev,
	  phot_top[i].nion,phot_top[i].z,phot_top[i].istate,phot_top[i].np,phot_top[i].nlast,phot_top[i].macro_info,
	  phot_top[i].down_index,phot_top[i].up_index,phot_top[i].freq[0],phot_top[i].x[0],phot_top[i].f,phot_top[i].sigma);
  }
  fclose(fptr);

  fptr=fopen("xphot_tab.adat","w");
  fprintf(fptr,"nlev\tuplev\tnion\tz\tistate\tnp\tnlast\tmacro_info\tdown_index\tup_index\tf0\tsigma0\tf\tsigma\n");
  for (i=0;i<nions;i++)
  {
	  fprintf(fptr,"%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%e\t%e\t%e\t%e\n",xphot_tab[i].nlev,xphot_tab[i].uplev,
	  xphot_tab[i].nion,xphot_tab[i].z,xphot_tab[i].istate,xphot_tab[i].np,xphot_tab[i].nlast,xphot_tab[i].macro_info,
	  xphot_tab[i].down_index,xphot_tab[i].up_index,xphot_tab[i].freq[0],xphot_tab[i].x[0],xphot_tab[i].f,xphot_tab[i].sigma);
  }
  fclose(fptr);





  n = 0;
  while (n < nions)
    {
      xlines = 0;
      m=0;
      while (m < nlines)
	{
		
		if (line[m].nion == n)
	    {
	      xlines++;
	    }
		if (line[m].nion ==176)
		{
			printf ("Line data %i %i %e\n",line[m].nconfigl,line[m].nconfigu,line[m].freq);
		}
	  m++;
	}
	  Log ("ion %3d %2d %2d %4d  %4d\n", n, ion[n].z, ion[n].istate,
	       ion[n].nlevels, xlines);
	  n = n + 1;
	}
    }


/* 
   a21 alculates and returns the Einstein A coefficient 
   History:
   98aug        ksl     Coded and debugged
   99jan        ksl Modified so would shortcircuit calculation if 
   called multiple times for same a
 */
#define A21_CONSTANT 7.429297e-22	// 8 * PI * PI * E * E / (MELEC * C * C * C)

struct lines *a21_line_ptr;
double a21_a;

double
a21 (line_ptr)
     struct lines *line_ptr;
{
  double freq;

  if (a21_line_ptr != line_ptr)
    {
      freq = line_ptr->freq;
      a21_a =
	A21_CONSTANT * line_ptr->gl / line_ptr->gu * freq * freq *
	line_ptr->f;
      a21_line_ptr = line_ptr;
    }

  return (a21_a);
}
