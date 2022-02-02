
/***********************************************************/
/** @file  unit_test.c
 * @author ksl
 * @date   January, 2021
 *
 * @brief  This file contains main and various related routines
 * to test one or more low level routines in Python
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/stat.h>
#include <string.h>

#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"



int
main (argc, argv)
     int argc;
     char *argv[];
{
  int n;
  double velocity[3];
  double t;
  double v;
  double speed = 0;
  double speed2 = 0;

//  FILE *fptr, *fopen ();

  int my_rank;                  // these two variables are used regardless of parallel mode
  int np_mpi;                   // rank and number of processes, 0 and 1 in non-parallel



#ifdef MPI_ON
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &np_mpi);
#else
  my_rank = 0;
  np_mpi = 1;
#endif

  np_mpi_global = np_mpi;       // Global variable which holds the number of MPI processes
  rank_global = my_rank;        // Global variable which holds the rank of the active MPI process
  Log_set_mpi_rank (my_rank, np_mpi);   // communicates my_rank to kpar

  init_rand (1084515760 + (13 * rank_global));



  t = 1.0e4;

  for (n = 0; n < 1000000; n++)
  {
    compton_get_thermal_velocity (t, velocity);
    printf ("%10.3e %10.3e %10.3e\n", velocity[0], velocity[1], velocity[2]);
    v = sqrt (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]);
    speed += v;
    speed2 += v * v;

  }


  printf (" Average speed %e  average v**2  %e\n", speed / n, speed2 / n);

  printf (" The expected average is %e\n", 2. / sqrt (3.141592) * sqrt (2. * BOLTZMANN * t / MELEC));


  exit (0);

}
