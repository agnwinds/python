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
#include "../source/log.h"

/** *******************************************************************************************************************
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

int
main (void)
{
  Log_set_verbosity (0);

  // Initialize the CUnit test registry
  if (CU_initialize_registry () != CUE_SUCCESS)
  {
    return CU_get_error ();
  }

#ifdef CUDA_ON
  cusolver_create ();
#endif

  // Add a suite to the registry
  CU_pSuite suite = CU_add_suite ("Matrix Functions Suite", NULL, NULL);
  if (suite == NULL)
  {
    CU_cleanup_registry ();
    return CU_get_error ();
  }

  // Add the test cases to the suite
  if ((CU_add_test (suite, "Solve Matrix", test_solve_matrix) == NULL) ||
      (CU_add_test (suite, "Invert Matrix", test_invert_matrix) == NULL))
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
  cusolver_destroy ();
#endif

  return CU_get_error ();
}
