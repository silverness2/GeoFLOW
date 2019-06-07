/*
 * Filename:	assert.hpp
 * Author:      bflynt
 * Created:		Jun 11, 2016
 * Copyright:   2016, Bryan Flynt, All rights reserved.
 */
#ifndef STDX_ASSERT_H_
#define STDX_ASSERT_H_

/*
 * Macros for:
 * ASSERT - Assert condition evaluates to true at
 * runtime, otherwise exit program
 *
 * STATIC_ASSERT - Assert condition evaluates to true at
 * compile time.
 */


#if !defined(NDEBUG)
#include <cstdlib>   // std::exit(EXIT_FAILURE)
#include <iostream>  // std::cerr
#define ASSERT(EXP)                                                          \
		do {                                                                 \
			if (!(EXP)) {                                                    \
				std::cerr << std::endl;                                      \
				std::cerr << "***** Failed Assertion *****" << std::endl;    \
				std::cerr << "Failed expression: " << # EXP << std::endl;    \
				std::cerr << "File: " << __FILE__ << std::endl;              \
				std::cerr << "Line: " << __LINE__ << std::endl;              \
				std::cerr << std::endl;                                      \
				std::exit(EXIT_FAILURE);                                     \
			}                                                                \
		} while (0)

#else
#define ASSERT(EXP)
#endif

// ------------------------------------------------------------------- //
//                           STATIC_ASSERT                             //
// ------------------------------------------------------------------- //
//
// Simply forwards to the C++ static_assert function
//
#define STATIC_ASSERT(...) static_assert(__VA_ARGS__)



#endif /* STDX_ASSERT_H_ */
