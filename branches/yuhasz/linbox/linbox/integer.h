/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/** @name linbox/integer.h
 * @memo integer types 
 *
 * Copyright(c)'94-97 by Givaro Team
 * Copyright(c)'2000-2002 by LinBox Team 
 * see the copyright file.
 * Created by M. Samama, T. Gautier
 *
 * Modified Jean-Guillaume.Dumas <Jean-Guillaume.Dumas@imag.fr>
 *          B. David Saunders <saunders@cis.udel.edu>,
 *          Bradford Hovinen <hovinen@cis.udel.edu>
 *          Gilles Villard <Gilles.Villard@ens-lyon.fr>
 *                        JGD Random functions back.                          
 *                        (2002/02/12 16:05:24) 
 *
 */

#ifndef __INTEGER_H
#define __INTEGER_H

#include "linbox-config.h"

#include "gmp++/gmp++.h"

namespace LinBox
{
	///wrapper of GMP
	typedef Integer integer;

	typedef signed __LINBOX_INT8 int8;
	typedef signed __LINBOX_INT16 int16;
	typedef signed __LINBOX_INT32 int32;
	typedef signed __LINBOX_INT64 int64;

	typedef unsigned __LINBOX_INT8 uint8;
	typedef unsigned __LINBOX_INT16 uint16;
	typedef unsigned __LINBOX_INT32 uint32;
	typedef unsigned __LINBOX_INT64 uint64;
}

#endif // __INTEGER_H
