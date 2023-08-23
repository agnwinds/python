/** ********************************************************************************************************************
 *
 * @file assert.h
 * @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 * @date August 2023
 *
 * @brief Additional custom assertions for testing
 *
 * ****************************************************************************************************************** */

#ifndef ASSERT_H
#define ASSERT_H

#define TRUE 1
#define FALSE 0
#define CU_BUFFER_SIZE 1024
#define EPSILON 1.0e-6

#define CU_CHECK_DOUBLE_ARRAY_EQ_FATAL(actual, expected, size, granularity)                                            \
{                                                                                                                      \
    int i;                                                                                                             \
    int not_equal_count = 0;                                                                                           \
    for (i = 0; i < size; ++i)                                                                                         \
    {                                                                                                                  \
        if (fabs((double)(actual[i]) - (double)(expected[i])) >= fabs((double)(granularity)))                          \
        {                                                                                                              \
            not_equal_count += 1;                                                                                      \
        }                                                                                                              \
    }                                                                                                                  \
  CU_assertImplementation(not_equal_count == 0, __LINE__, "CU_CHECK_DOUBLE_ARRAY_EQ_FATAL", __FILE__, "", TRUE);       \
}

#define CU_CHECK_DOUBLE_ARRAY_EQ(actual, expected, size, granularity)                                                  \
{                                                                                                                      \
    int i;                                                                                                             \
    int not_equal_count = 0;                                                                                           \
    for (i = 0; i < size; ++i)                                                                                         \
    {                                                                                                                  \
        if (fabs((double)(actual[i]) - (double)(expected[i])) >= fabs((double)(granularity)))                          \
        {                                                                                                              \
            not_equal_count += 1;                                                                                      \
        }                                                                                                              \
    }                                                                                                                  \
  CU_assertImplementation(not_equal_count == 0, __LINE__, "CU_CHECK_DOUBLE_ARRAY_EQ_FATAL", __FILE__, "", FALSE);      \
}

#define CU_FAIL_MSG_FATAL(...)                                                                                         \
{                                                                                                                      \
    char buffer[CU_BUFFER_SIZE];                                                                                       \
    snprintf(buffer, CU_BUFFER_SIZE, __VA_ARGS__);                                                                     \
    CU_FAIL_FATAL(buffer);                                                                                             \
}

#endif
