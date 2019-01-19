


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

  Description:

  Arguments:  (Input via .pf file)


  Returns:

  Notes:

  This is roughly what the Makefile should look like. Uncomment the mv
  statement to store binary in ~/bin

CC = gcc
CFLAGS = -c -g -I$$HOME/include -Wall
LDFLAGS= -L$$HOME/lib -lm -lkpar
BIN = $$HOME/bin

simple: simple.o
	gcc simple.o  $(LDFLAGS) -o simple
#	mv $@ $(BIN)



  History:
2004	ksl	Coded as better.c

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#include "log.h"
#include "atomic.h"
#include "python.h"




/* the form of the structure which contains
 * the map from words to values
    struct  xchoices {
        char choices[10][LINELENGTH];
        int  vals[10];
        int  n;
    } zz_spec; 
*/



/* this is a routine which is intended to initialize a structure 
 * used in the special instance when we are going to use the rdchoice
 * command, but when we want different options to be present in different
 * calls.  init_choices only needs to be called once
 */

int xinit_choices = 0;

int
init_choices ()
{
  /* Initialize the structure that contains all of the types of possilbe radiations */

  char *xchoices[] = { "bb", "uniform", "power", "cloudy", "brems", "none", "models" };
  int xvals[] = { SPECTYPE_BB, SPECTYPE_UNIFORM, SPECTYPE_POW, SPECTYPE_CL_TAB, SPECTYPE_BREM, SPECTYPE_NONE, SPECTYPE_MODEL };
  int num_choices = 7;

  if (xinit_choices)
    return (0);

  int i = 0;
  for (i = 0; i < num_choices; i++)
  {
    strcpy (zz_spec.choices[i], xchoices[i]);
    zz_spec.vals[i] = xvals[i];
  }
  zz_spec.n = num_choices;

  /* To initialize more than one sets of choices one needs essentially to repeat the lines aboe
   * making sure that the names of the input arrays changed.  Alternatively, one could create
   * a routine for each structure one wanted to initialize
   *
   * If one had large numbers of variables to intialize by this method, one could put all oft he
   * relevant data into a file
   */

  xinit_choices = 1;


  //OLD printf ("Verify %s\n", zz_spec.choices[1]);


  return (0);
}


int
get_choices (question, choices, qstruct)
     char *question;
     char *choices;
     struct rdpar_choices *qstruct;
{
  //char *xchoices[] = {"bb", "uniform", "power", "cloudy", "brems", "none", "model"};
  // int xvals      [] = {SPECTYPE_BB, SPECTYPE_UNIFORM, SPECTYPE_POW, SPECTYPE_CL_TAB, SPECTYPE_BREM, SPECTYPE_NONE, SPECTYPE_MODEL};
  int num_choices = 7;


  char cur_choices[10][LINELENGTH];
  int cur_values[10] = { -999 };
  int cur_num;
  char cur_string[LINELENGTH];


  if (xinit_choices == 0)
  {
    Error ("init_choices needs to have been called before get_choices.  Programming error\n");
    exit (0);
  }

  /*
   * A typical question will be disk.spec_type(bb,power), that is to
   * say * it will not include all of the possibileities
   */


  /* Now match the current possibilities to the integer choices */

  int n, nstart, nstop;
  char dummy[LINELENGTH];
  nstart = nstop = 0;
  for (n = 0; n < strlen (question); n++)
  {
    if (question[n] == '(')
    {
      nstart = n;
    }
    if (question[n] == ')')
    {
      nstop = n;
    }
  }


  snprintf (dummy, nstop - nstart, "%s", &question[nstart + 1]);
  //OLDprintf ("dummy2 %s\n", dummy);

  //OLD strcpy (dummy, " ");
  //OLD  strncpy (dummy, &question[nstart + 1], nstop - nstart - 1);

  //OLD printf ("question %s\n", question);
  //OLD printf ("dummy %d %d %lu %s\n", nstart, nstop, strlen (question), dummy);


  strcpy (dummy, " ");
  strcpy (dummy, question);
  //OLD printf ("question %s\n", dummy);
  int ncommas = 0;
  int nparen = 0;
  int i = 0;
  int j = 0;
  for (i = 0; i < strlen (dummy); i++)
  {
    if (dummy[i] == '(')
    {
      nparen += 1;
      dummy[i] = ' ';
    }
    if (dummy[i] == ')')
    {
      nparen += 1;
      dummy[i] = ' ';
    }
    if (dummy[i] == ',')
    {
      ncommas += 1;
      dummy[i] = ' ';
    }
  }
  //OLD printf ("question %s\n", dummy);

  cur_num = sscanf (dummy, "%*s %s %s %s %s %s %s %s %s %s %s",
                    cur_choices[0], cur_choices[1],
                    cur_choices[2], cur_choices[3],
                    cur_choices[4], cur_choices[5], cur_choices[6], cur_choices[7], cur_choices[8], cur_choices[9]);

  //OLD printf ("The current number of choices is %d\n", cur_num);

  /* Now we have to find out the values associated with these choices */
  //OLD printf ("Test1 %s\n", qstruct->choices[1]);
  //OLD printf ("Test2 %d\n", qstruct->vals[1]);
  //OLD printf ("Test3 %s\n", zz_spec.choices[1]);
  //OLD printf ("Test4 %d\n", zz_spec.vals[1]);


  for (i = 0; i < cur_num; i++)
  {
    for (j = 0; j < num_choices; j++)
    {
      if (strncmp (cur_choices[i], qstruct->choices[j], strlen (qstruct->choices[j])) == 0)
      {
        //OLD printf ("gotcha  %s  s %d\n", qstruct->choices[j], qstruct->vals[j]);
        cur_values[i] = qstruct->vals[j];
        break;
      }
      if (j == num_choices)
      {
        Error ("Did not match %s\n", cur_choices[i]);
        Error ("This is a programming error\n");
        exit (0);
      }
    }
  }
  /*
   * If we have gotten to this point, we can construct the integer
   * string that represents the choices that are currently available
   */

  sprintf (cur_string, "%d", cur_values[0]);
  for (i = 1; i < cur_num; i++)
  {
    sprintf (cur_string, "%s,%d", cur_string, cur_values[i]);
  }

  //OLD printf ("The new string is %s\n", cur_string);
  strcpy (choices, cur_string);






  return (0);
}
