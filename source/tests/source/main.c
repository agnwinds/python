/* ****************************************************************************************************************** */
/**
 *  @file main.c
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date August 2023
 *
 *  @brief Main function for unit tests
 *
 *  ***************************************************************************************************************** */

#include <stdlib.h>
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

/* Test suite prototypes */
#include "tests/tests.h"

/* Python source code */
#include "../source/log.h"


/** *******************************************************************************************************************
 *
 *  @brief Entry point for unit tests
 *
 *  ***************************************************************************************************************** */

int
main (void)
{
  /* I don't want any output from Python, thank you very much */
  Log_set_verbosity (0);

  /* Initialize the CUnit test registry */
  if (CU_initialize_registry () != CUE_SUCCESS)
  {
    return CU_get_error ();
  }

  /* Create test suites */
  create_matrix_test_suite ();

  /* Run the test suites */
  CU_basic_set_mode (CU_BRM_VERBOSE);
  CU_basic_run_tests ();

  /* Clean up the CUnit registry */
  CU_cleanup_registry ();

  return CU_get_error ();
}
