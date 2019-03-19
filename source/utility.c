/***********************************************************/
/** @file  utility.c
 * @author ejp
 * @date   March 2019
 *
 * @brief
 * Functions which are used for general utility purposes in Python.
 *
 *
 * ### Notes ###
 *
 * Created originally to store the Exit wrapper function as it felt more natural
 * placing it here than trying to kludge some global variables into xlog.c.
 *
 ***********************************************************/


#include <stdio.h>
#include <stdlib.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/**
 *  @brief Wrapper function to exit MPI/Python.
 *
 *  @param[in] int error_code. An integer classifying the error. Note that this
 *  number should be a positive non-zero integer.
 *
 *  @details
 *  When MPI is in use, using the standard C library exit function can be somewhat
 *  dangerous, as if a process exits with return code 0, MPI assumes that the
 *  process has exited correctly and the program will continue. This results in
 *  a deadlock and without any safety mechanism can result in an MPI run never
 *  finishing. Hence, in multiprocessor mode, one should use use MPI_Abort
 *  (which isn't a very graceful) to ensure that if a process does need to exit,
 *  then all processes will exit with it.
 *
 *  Before exiting, an error summary is printed to screen and and log is flushed
 *  to ensure everything is up to date to aid with error diagnosis.
 *
 *  ### Programming Notes ###
 *
 *  Maybe this function belongs elsewhere.
 *
 *  The error summary is only printed to the master processes' stdout, but the
 *  error summary can also be found in each MPI processes' diag file.
 *
 **********************************************************/

void
Exit (int error_code)
{
#ifdef MPI_ON

  Log_parallel ("Rank %i: Python exiting all processes with error %i\n", rank_global, error_code);
  Log_parallel ("Check diag/%s_%i.diag for more information\n", files.root, rank_global);

  if (np_mpi_global > 1)
  {
    error_summary_parallel ("summary prior to MPI abort");
    MPI_Abort (MPI_COMM_WORLD, error_code);
  }
  else
  {
    error_summary ("summary prior to abort");
    exit (error_code);
  }

#else
  Log ("\n--------------------------------------------------------------------------\n");
  Log ("Python exiting with error %i\n", error_code);
  Log ("Check diag/%s.diag for more information\n", files.root);
  error_summary ("summary prior to abort");
  exit (error_code);
#endif
}
