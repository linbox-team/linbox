/* linbox/util/debug.h
 *
 * Copyright (C) 2001,2010 LinBox
 * Copyright (C) 2001 Bradford Hovinen
 * Copyright (C) 2010 LinBox
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 * Modified by BB.
 *
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *
 */

/*! @file util/debug.h
 * @ingroup util
 * Various utilities for debugging.
 * @todo we should put vector printing elsewhere.
 */

#ifndef __LINBOX_util_debug_H
#define __LINBOX_util_debug_H

#include "linbox/util/error.h"
#include <iostream>
#include <list>
#include <sstream>
#include <vector>

/**
 * Check an assertion (à la \c std::assert).
 *
 * If in DEBUG mode, throws a \ref PreconditionFailed exception.
 * In RELEASE mode, nothing is checked.
 * @param check assertion to be checked.
 */
#ifdef NDEBUG // à la assert.
#define linbox_check(check)
#else
#ifdef __GNUC__
#define linbox_check(check)                                                                                                      \
    if (!(check)) throw ::LinBox::PreconditionFailed(__func__, __FILE__, __LINE__, #check);
#else
#define linbox_check(check)                                                                                                      \
    if (!(check)) throw ::LinBox::PreconditionFailed(__FILE__, __LINE__, #check);
#endif
#endif

// @fixme THIS IS UGLY
#define THIS_CODE_COMPILES_BUT_IS_NOT_TESTED                                                                                     \
    std::cout << "*** Warning *** " << std::endl                                                                                 \
              << __func__ << " in " << __FILE__ << ':' << __LINE__ << " is not tested" << std::endl;

#define THIS_CODE_MAY_NOT_COMPILE_AND_IS_NOT_TESTED                                                                              \
    throw(" *** Warning ***  this piece of code is not compiled by default and may not work")

#define LB_FILE_LOC __func__, __FILE__, __LINE__

#if defined(LinBoxSrcOnly) or defined(LinBoxTestOnly)
// for all-source compilation
#include "linbox/util/debug.C"
#endif

#include <givaro/givprint.h>

// @note Taken from contracts.h
#if _DEBUG == 2

// The debug mode also defines the following macros. Failure of any of these macros leads to
// program termination. The user is notified of the error condition with the right file name
// and line number. The actual failing operation is also printed using the stringizing operator #

#define ASSERT(bool_expression) assert(bool_expression)
#define REQUIRE(bool_expression) ASSERT(bool_expression)
#define STATE(expression) expression

#elif _DEBUG == 1

// @fixme InvalidOperationException does not exists
#define ASSERT(bool_expression) assert(bool_expression)
#define REQUIRE(bool_expression)                                                                                                 \
    if (!bool_expression) throw InvalidOperationException(__FILE__, __LINE__);
#define STATE(expression) expression

#elif _DEBUG == 0

// When built in release mode, the _DEBUG flag would not be defined, thus there will be no overhead
// in the final release from these checks.

#define ASSERT(ignore) ((void)0)
#define REQUIRE(ignore) ((void)0)
#define STATE(ignore) ((void)0)

#endif

#endif // __LINBOX_util_debug_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
