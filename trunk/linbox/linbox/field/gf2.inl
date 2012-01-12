/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/field/gf2.inl
 * Copyright (C) 2003 Bradford Hovinen
 *
 * Written by Bradford Hovinen, Dumas, bds
 *
 * ------------------------------------
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_field_gf2_INL
#define __LINBOX_field_gf2_INL

#include <iostream>
#include <time.h>

#include <cctype> //isdigit

template<typename Vector>
std::ostream& afficheVector (std::ostream& o, const Vector& C)
{
	for(typename Vector::const_iterator refs =  C.begin();
	    refs != C.end() ;
	    ++refs )
		o << (*refs) << " " ;
	return o;
}



#endif // __LINBOX_field_gf2_INL

