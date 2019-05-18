
/***********************************************************/
/** @file  signal.c
 * @author ksl
 * @date   March, 2018
 *
 * @brief  Routines indicating  how  many ionization or spectal
 * generation cycles have been completed in 
 * in a file named root.sig where 
 * root is the rootname of the .pf file and which stop the
 * program if one has set a maximum execution time and this
 * has been exceeded. 
 *
 * The purpose of these routines are to give a very high level
 * view of where the program is at any time, and in particular
 * whether the program was complete at the time the program exited.
 *
 * The were written for running models on systems whre there 
 * are time limits on how long a program can run, such as 
 * beowulf clusters.  
 *
 * They permit you to stop (checkpoint) 
 * the program after a certain time, and then using scripts
 * to easily determine whether the program ran to completion
 * and if not
 * to restart the program to finish a run
 *
 *
 ***********************************************************/




#include <stdio.h>
#include <strings.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "log.h"

#include "atomic.h"
#include "python.h"




/**********************************************************/
/** 
 * @brief      xsignal generates a single line message to a file names root.sig
 *
 * @param [in] char *  root   root name of the file to wirte to
 * @param [in] char *  format   A format string
 * @param [in]   ...   The remaining inputs for the fprintf statement
 * @return     Always  returns 0
 *
 * If one cannot write to the .sig file, Python will exit
 *
 * @details
 * 
 * ### Notes ###
 * In principle almost anything can be written to the .sig file, as this 
 * uses vfprintf to write to the file.  In pracitce all of the messages
 * have the format 
 *
 * 
 * Mon Nov 10 09:05:34 2008     10.0  message 
 * 
 * where the message is determined by the format and the extra variables 
 * that are passed to the program.  
 *
 * The messages are all written from the rank 0 thread.
 * 
 *
 **********************************************************/

int
xsignal (char *root, char *format, ...)
{

  va_list ap, ap2;

  char curtime[LINELENGTH];
  char message[LINELENGTH];
  FILE *fopen (), *sptr;
  char filename[LINELENGTH];
  double elapsed_time;


  /* if we are in MPI mode we only want the master process to write to sig file */
#ifdef MPI_ON
  if (rank_global == 0)
  {
#endif


    /* Make the filemne */
    strcpy (filename, "");
    strcpy (filename, root);
    strcat (filename, ".sig");

    /* Open the file so that it will append if the file exists */


    if ((sptr = fopen (filename, "a")) == NULL)
    {
      Error ("xsignal: Could not even open signal file %s\n", filename);
      Exit (0);
    }

    /* Now generate the message */

    /* Get the current time */
    get_time (curtime);


    elapsed_time = timer ();

    /* Get the time since the time was initiated */


    fprintf (sptr, "%s %8.1f ", curtime, elapsed_time);


    va_start (ap, format);
    va_copy (ap2, ap);          /* Added because vfprintf can change ap */
    vfprintf (sptr, format, ap);
    va_end (ap);


    vsprintf (message, format, ap2);
    Log ("xxx %s %8.1f %s", curtime, elapsed_time, message);


    fclose (sptr);

#ifdef MPI_ON
  }
#endif

  return (0);
}



/**********************************************************/
/** 
 * @brief      Remove the old signal file so that one can begin again
 *
 * @param [in, out] char *  root   Root name of the .sig file
 * @return     Always returns 0
 *
 * @details
 * The routine checks whether a .sig file exists, and if so 
 * removes it
 *
 * ### Notes ###
 *
 * Thread 0 is responsible for executing the removal.
 *
 **********************************************************/

int
xsignal_rm (char *root)
{

#ifdef MPI_ON
  if (rank_global == 0)         // only remove sig file if 0th thread
  {
#endif

    char filename[LINELENGTH];
    char command[LINELENGTH];
    FILE *tmp_ptr;
    /* Make the filemne */
    strcpy (filename, "");
    strcpy (filename, root);
    strcat (filename, ".sig");

    /* first check if the file exists */

    if ((tmp_ptr = fopen (filename, "r")) == NULL)
    {
      return (0);
    }


    strcpy (command, "rm ");
    strcat (command, filename);
    system (command);

#ifdef MPI_ON
  }
#endif

  return (0);

}




double max_time = -1.0;


/**********************************************************/
/** 
 * @brief      Set the maximum time the program should run
 *
 * @param [in, out] char *  root   Root name of the .sig file where the max time is recorded
 * @param [in, out] double  t   The maximum time one wants to proposal to run
 * @return     Always returns 0
 *
 * @details
 * The routine simple sets the max_time (an external variable) and
 * writs information about this to .sig and log files.
 *
 * ### Notes ###
 *
 **********************************************************/

int
set_max_time (char *root, double t)
{
  Log ("Setting maximum time for program to run to be %d s\n", t);
  xsignal (root, "MAXTIME set to %.0f\n", t);
  max_time = t;
  return 0;
}




/**********************************************************/
/** 
 * @brief      check whether the elapsed time is greater than the max_time
 *   and terminate the program if that is the case.
 *
 * @param [in] char *  root   root name of the files where information is logged
 * @return     Returns 0, unless the time has exeeded the maximum time in which case
 * the program writes a comment to .sig file and exits
 *
 * @details
 *
 *  If the current time is greater than the allowed time. Then the program is
 *  stopped, and a signal is sent to the .sig file that the program can be 
 *  restarted
 *
 *  If the maximum time has not been set then this routine is a NOP
 *
 *
 * ### Notes ###
 *
 * Generally speaking one should send a message to xsignal before
 * invoking check_time to record the status of the signal file
 *
 *
 **********************************************************/

int
check_time (char *root)
{
  double t;
  if (max_time > 0.0 && (t = timer () > max_time))
  {
    error_summary ("Time allowed has expired expired\n");
    xsignal (root, "COMMENT max_time %.1f exceeded\n", max_time);
    Exit (0);  // At present this allows the code to exit compltely, but produces MPI warninngs see #518
    //OLD Exit (1);
  };

  return (0);
}
