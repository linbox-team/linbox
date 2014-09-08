/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* Copyright (C) 2014  Pascal Giorgi
 *
 * Written by Pascal Giorgi   <pascal.giorgi@lirmm.fr>
 *
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
 */

/*! @file field/Modular/modular-integer.h
 * @ingroup field
 * @brief  representation of <code>Z/mZ</code> over \c integer .
 */
#ifndef __LINBOX_modular_integer_H
#define __LINBOX_modular_integer_H
#include "linbox/integer.h"
#include "linbox/field/modular.h"

#include "fflas-ffpack/field/modular-integer.h"

namespace LinBox {
	

	template<>
	class Modular<integer>	 : public FFPACK::Modular<Givaro::Integer>,
				   public FieldInterface {
	public:
		using FFPACK::Modular<FFPACK::integer>::Modular;
		typedef FFPACK::Modular<FFPACK::integer> Father_t;
			
		inline integer getMaxModulus() const
		{
			return -1 ;
		}
		
		inline integer& convert(integer& x, const integer &y) const {
			return x=y; 
		}
};
} // end of namespace LinBox
#endif
