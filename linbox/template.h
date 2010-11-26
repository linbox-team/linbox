/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) <+year+> LinBox
 * Written by <+someone+> <<+her.mail@somewhere.net+>>
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

#ifndef __LINBOX_<+directory_file_name+>_H
#define __LINBOX_<+directory_file_name+>_H

/** @file <+directory/file-name.h+>
  * @brief desc
  * long doxy desc
  */

#include "<++>"

#if 0 // better than commenting out code with big /*   */ 
      // (it shows this code is not really mature yet...)

#define LB_VAR  10  /* local var */
#define LB_DEF      /* local define */
#define LB_MACRO () /* local macro */

#ifdef _LINBOX_DEF   // linbox global define
#ifdef _LINBOX_MACRO // linbox global macro
#if _LINBOX_VAR // linbox global variable
#endif
#endif
#endif

#endif


/*! @file template.h
 * @brief desc
 * what is this (important) file about ?
 */

namespace LinBox 
{
	/**  @brief this function is about...
	 * this important function has comments  !!!
	 */
	template<++>
	void my_func(T & toto) ;
} //LinBox

#include "<+file-name.inl+>" // implementation here

#undef LB_VAR    // environmentalists love us
#undef LB_DEF    // really !
#undef LB_MACRO 

#endif //__LINBOX_<+directory_file_name+>_H

