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

  // Create test suites
  create_matrix_test_suite ();

  // Run the tests
  CU_basic_set_mode (CU_BRM_VERBOSE);
  CU_basic_run_tests ();

  // Clean up the CUnit registry
  CU_cleanup_registry ();

  return CU_get_error ();
}
