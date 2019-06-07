/*
 * macros.hpp
 *
 *  Created on: Apr 25, 2019
 *      Author: bflynt
 */

#ifndef NULL_MACROS_HPP_
#define NULL_MACROS_HPP_




// ------------------------------------------------------------------- //
//                           NULL_STATEMENT                            //
// ------------------------------------------------------------------- //

/// Force the execution of a null statement
/**
 * Many compilers will not do a loop if there isn't a
 * statement inside so make a do nothing statement to insure it
 * does the loop.  Otherwise it will optimize the loop away.
 */
#define NULL_STATEMENT do{if(0){int nullstatement = 0;}}while(0)



// ------------------------------------------------------------------- //
//                           NULL_USE Definition                       //
// ------------------------------------------------------------------- //

/// Do nothing with a variable to make warning go away
/**
 * A null use of a variable, use to avoid compiler
 * warnings about unused variables.
 */
#define NULL_USE(VAR) do{(void)VAR;}while(0)



#endif /* NULL_MACROS_HPP_ */
