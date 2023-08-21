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

#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

#include "tests/tests.h"

/** *******************************************************************************************************************
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

int
main (void)
{
  // Initialize the CUnit test registry
  if (CU_initialize_registry () != CUE_SUCCESS)
  {
    return CU_get_error ();
  }

#ifdef CUDA_ON
  cuda_init ();
#endif

  // solve_matrix_small ();
  // solve_matrix_matrix_ion ();
  // invert_matrix_small ();

  // Add a suite to the registry
  CU_pSuite suite = CU_add_suite ("Matrix Functions Suite", NULL, NULL);
  if (suite == NULL)
  {
    CU_cleanup_registry ();
    return CU_get_error ();
  }

  // Add the test cases to the suite
  if ((CU_add_test (suite, "Solve Matrix: small", solve_matrix_small) == NULL) ||
      (CU_add_test (suite, "Solve Matrix: matrix ion", solve_matrix_matrix_ion) == NULL) ||
      (CU_add_test (suite, "Invert Matrix: small", invert_matrix_small) == NULL))
  {
    CU_cleanup_registry ();
    return CU_get_error ();
  }

  // Run the tests
  CU_basic_set_mode (CU_BRM_VERBOSE);
  CU_basic_run_tests ();

  // Clean up the CUnit registry
  CU_cleanup_registry ();

#ifdef CUDA_ON
  cuda_finish ();
#endif

  return CU_get_error ();

  return EXIT_SUCCESS;
}
