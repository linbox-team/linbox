/* Copyright (C) LinBox
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_util_contracts_H
#define __LINBOX_util_contracts_H

// Object class is the base class for all
// objects in the system. All classes inheriting from this class need 
// to define a method IsValid. This method should perform a
// consistency check on the state of the object. Note that 
// this method needs to be defined only when a debug build is made
#include <cassert>
#include <cstddef>
#if 0
class Object
{
public:
#ifdef _DEBUG
    virtual bool IsValid() const = 0;
#endif
    
};
#endif
#if _DEBUG == 2

// The debug mode also defines the following macros. Failure of any of these macros leads to
// program termination. The user is notified of the error condition with the right file name
// and line number. The actual failing operation is also printed using the stringizing operator #

#define ASSERT(bool_expression) assert(bool_expression)
#define IS_VALID(obj) ASSERT((obj) != NULL && (obj)->IsValid())
#define REQUIRE(bool_expression) ASSERT(bool_expression)
#define ENSURE(bool_expression) ASSERT(bool_expression)
#define STATE(expression) expression

#elif _DEBUG == 1

#define ASSERT(bool_expression) assert(bool_expression)
#define IS_VALID(obj) ASSERT((obj) != NULL && (obj)->IsValid())
#define REQUIRE(bool_expression) if( ! bool_expression ) throw InvalidOperationException(__FILE__,__LINE__);
#define ENSURE(bool_expression) ASSERT(bool_expression)
#define STATE(expression) expression

#elif _DEBUG == 0

// When built in release mode, the _DEBUG flag would not be defined, thus there will be no overhead
// in the final release from these checks.

#define ASSERT(ignore) ((void) 0)
#define IS_VALID(ignore) ((void) 0) 
#define REQUIRE(ignore) ((void) 0)
#define ENSURE(ignore) ((void) 0)
#define STATE(ignore) ((void) 0)

#endif

#endif //__LINBOX_util_contracts_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
