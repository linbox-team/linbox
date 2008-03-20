/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/field/rebind.h
 * Copyright (C) 2005 Jean-Guillaume Dumas
 * Time-stamp: <16 Jun 05 09:00:15 Jean-Guillaume.Dumas@imag.fr> 
 */
#ifndef __LINBOX_REBIND_H
#define __LINBOX_REBIND_H

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
#endif
