/* Copyright (C) <+year+> LinBox
 * Written by <+someone+> <<+so.me@on.ne+>>
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

#ifndef __LINBOX_<+the_name+>_H
#define __LINBOX_<+the_name+>_H

/** @file <+the/file.h+>
  * @brief desc
  * long doxy desc
  */

#include <<++>>

#if 0 // better than commenting out code with big /*   */ 
      // (it shows this code is not really mature yet...)

#define LB_VAR  /* local var */
#define LB_DEF  /* local define */

#ifdef _LINBOX_DEF // linbox global define
#endif

#endif

//#if LINBOX_VAR // linbox global variable

/*! @file template.h
 * @brief desc
 * what is this (important) file about ?
 */

namespace LinBox 
{
	/**
	 * this important code has comments so that other people 
	 * will understand it 10 yrs afterwards !!!!!!
	 */
	template<++>
}

#undef LB_VAR // environmentalists love us
#undef LB_DEF // really !

#endif //__LINBOX_<+the_name+>_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
