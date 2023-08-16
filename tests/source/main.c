/* ****************************************************************************************************************** */
/**
 *  @file main.c
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date August 2023
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

#include <stdlib.h>

#include "tests/tests.h"

int cuda_init (void);
int cuda_finish (void);

/** *******************************************************************************************************************

 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

int
main (void)
{
#ifdef CUDA_ON
  cuda_init ();
#endif

  test_solve_matrix ();
  test_invert_matrix ();

#ifdef CUDA_ON
  cuda_finish ();
#endif

  return EXIT_SUCCESS;
}
