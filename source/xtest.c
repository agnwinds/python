

/***********************************************************/
/** @file  xtest.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  A test routine to check a part of the program 
 * for any special test.  This is not intended to be 
 * used for anything permanent.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/** 
 * @brief      Conduct a diagnostic test for some purpose
 *
 *
 * @details
 *  This is a totally uncommented routine for use 
 *  with Python to perform some specific diagnostic 
 *  test
 *
 * ### Notes ###
 *
 * xtest is called based on a commmand line switch
 * -xtest.  Depending on whether an exit command is
 * included here (the normal situation), the code
 * will continue to run after that.
 *
 * The routine is called after all of the inputs
 * of a run have been read in.  In many cases
 * (namely when the -r (restart option is used,
 * a previous wind file.
 *
 * xtest is intended to be temporary.  Developers should
 * not expect to be the same the next time they want to
 * debug something.  As a result, developers may want
 * to save their version of xtest.c outside of the
 * git repository for later use.
 *
 *
 *
 **********************************************************/

int
xtest ()
{
  double x[] = { 4.4300e+11, 0.0000e+00, 6.0384e+16 };
  double v_grad[3][3];
  int i;
  int j;
  double ds_fraction;



  get_derivative ();

  rel_mode = REL_MODE_LINEAR;
  model_vgrad (1, x, v_grad);
  Log ("\n OLD Linear\n");
  for (i = 0; i < 3; i++)
  {
    Log ("Results %11.4e %11.4e %11.4e \n", v_grad[i][0], v_grad[i][1], v_grad[i][2]);
  }

  rel_mode = REL_MODE_FULL;
  Log ("\n OLD Full  \n");
  model_vgrad (1, x, v_grad);
  for (i = 0; i < 3; i++)
  {
    Log ("Results %11.4e %11.4e %11.4e \n", v_grad[i][0], v_grad[i][1], v_grad[i][2]);
  }


  ds_fraction = 0.01;

  for (j = 0; j <= 10; j++)
  {

    ds_fraction /= 10;

    Log ("ds_fraction: %10.3e\n", ds_fraction);

    rel_mode = REL_MODE_LINEAR;
    Log ("\n NEW Linear\n");
    xmodel_vgrad (ds_fraction, 1, x, v_grad);
    for (i = 0; i < 3; i++)
    {
      Log ("Results %11.4e %11.4e %11.4e \n", v_grad[i][0], v_grad[i][1], v_grad[i][2]);
    }


    rel_mode = REL_MODE_FULL;
    Log ("\n NEW Full\n");
    xmodel_vgrad (ds_fraction, 1, x, v_grad);
    for (i = 0; i < 3; i++)
    {
      Log ("Results %11.4e %11.4e %11.4e \n", v_grad[i][0], v_grad[i][1], v_grad[i][2]);
    }

  }



  Exit (0);

  return (0);

}
