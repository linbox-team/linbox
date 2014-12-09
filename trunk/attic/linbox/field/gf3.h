/* linbox/field/PID-double.h
 * Copyright (C) 2013 LinBox team
 *
 * Written by : bds
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


/*!  @file field/gf3.h
 * @ingroup field
 * @brief NO DOC
 */

#ifndef __LINBOX_gf3_H
#define __LINBOX_gf3_H

#include "linbox/field/modular.h"
#include "linbox/field/field-traits.h"

namespace LinBox
{


	/*! \ingroup integers
	 * @brief NO DOC
	 */
	class GF3 : public Modular<float> { // could be byte size.

	public:
		GF3(integer p = 3, size_t e = 1) : Modular<float>(3) {}

		typedef float Element;

#if 0
		/// test if unit (1 or -1)
		inline static bool isUnit (const Element& x)
		{ return x != 0; }

		// some specializations and conversions
		double& convert(double& x, const Element& y) const
		{ return x=y; }

		Element& init(Element& x, const float& y) const
		{
			return x=y;
		}

		integer& convert(integer& x, const Element& y) const
		{
			return x=y;
		}

		Element& init(Element& x, const integer& y) const
		{
			return x=y;
		}

#endif
		std::ostream &write (std::ostream &os) const
		{ return os << "GF3"; }
		std::ostream &write (std::ostream &os, const Element &x) const
		{ return Modular<float>::write(os, x); }

	}; //end of class GF3

} //end of namespace LinBox
#endif //__LINBOX_gf3_H



// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
