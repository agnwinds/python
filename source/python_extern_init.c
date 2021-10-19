#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "atomic.h"
#include "python.h" 

 
int np_mpi_global;              // Global variable which holds the number of MPI processes

int rank_global;

int verbosity;                  /* verbosity level. 0 low, 10 is high */

int cell_phot_stats;            //1=do  it, 0=dont do it
int ncstat;                     // the actual number we are going to log
int ncell_stats[NCSTAT];        //the numbers of the cells we are going to log

int nerr_no_Jmodel;
int nerr_Jmodel_wrong_freq;


int xxxbound;

FILE *pstatptr=NULL;

