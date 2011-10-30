/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/field/gf2.inl
 * Copyright (C) 2003 Bradford Hovinen
 *
 * Written by Bradford Hovinen, Dumas, bds
 *
 * ------------------------------------
 *
 * See COPYING for license information.
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

