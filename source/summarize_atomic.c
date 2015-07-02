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
  int n, m;
  int xlines;

  write_atomicdata = 1;
  strcat (atomic_filename, argv[1]);

  get_atomic_data (atomic_filename);

  Log ("The number of elements is %d\n", nelements);
  Log ("The number of ions is %d", nions);

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
