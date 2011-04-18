/* Copyright (C) 2005 LinBox
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr> 
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

#ifndef __LINBOX_rebind_H
#define __LINBOX_rebind_H

namespace LinBox
{

/** \brief used in support of Hom, MatrixHom 

Helps define rebind for vector types.  See blackbox/sparse.h for example of use.
*/
template<class XXX, class U>
struct Rebind 
{
    typedef typename XXX::template rebind<U>::other other;
};

}
#endif //__LINBOX_rebind_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
