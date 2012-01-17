/* linbox/switch/cekstv.h
 * Copyright (C) 1999-2001 William J Turner
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>
 *
 * -----------------------------------------------------------
 * 2002-09-26  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Refactoring: The switch object now contains only the information for one 2x2
 * block. A vector of switches is maintained by the butterfly preconditioner
 * instead of the switch object. Since there will be many switch objects, they
 * should be kept very lightweight, so the field is not maintained in the object
 * itself, but instead passed to apply and applyTranspose. Since those methods
 * are inline, this does not create overhead. apply and applyTranspose now take
 * four field elements: the source elements and destination elements. This
 * eliminates the need to keep an additional temporary in the class, and
 * eliminates the need for copying in the butterfly.
 *
 * -----------------------------------------------------------
 * 2002-08-20  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Brought this file into the current Linbox framework:
 *   - Renamed file as cekstv.h
 *   - Renamed class cekstv_switch as CekstvSwitch
 *   - Reindent
 * -----------------------------------------------------------
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
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_cekstv_gf2_H
#define __LINBOX_cekstv_gf2_H

#include <vector>

#include "linbox/field/gf2.h"
#include "linbox/switch/cekstv.h"


namespace LinBox
{
	// Specialization of Butterfly switch object
	template <>
	class CekstvSwitch<GF2>
	{
	public:
		typedef GF2 Field;
		/// Typedef
		typedef Field::Element Element;
		typedef CekstvSwitch<Field> Self_t;
		typedef CekstvSwitchFactory<Field> Factory;

		/** Constructor from a field and a field element.
		 * @param F field in which arithmetic is done
		 * @param switches vector of switches
		 */
		CekstvSwitch (const Field::Element &a) :
			_a (a)
		{}

		~CekstvSwitch () {}

		bool apply (const Field &F, Element &x, Element &y) const
		{
			F.axpyin (x, _a, y);
			F.addin (y, x);
			return true;
		}

		bool applyTranspose (const Field &F, Element &x, Element &y) const
		{
			F.addin (x, y);
			F.axpyin (y, _a, x);
			return true;
		}

		bool apply (const Field &F, std::_Bit_reference x, std::_Bit_reference y) const
		{
			F.axpyin (x, _a, y);
			F.addin (y, x);
			return true;
		}

		bool applyTranspose (const Field &F, std::_Bit_reference x, std::_Bit_reference y) const
		{
			F.addin (x, y);
			F.axpyin (y, _a, x);
			return true;
		}

		template<typename _Tp1>
		struct rebind
		{
			typedef CekstvSwitch<_Tp1> other;

			// special rebind operator() with two fields,
			// indeed local field is not stored in the switch
			void operator() (other *& Ap, const Self_t& A, const _Tp1& T, const Field& F) {
				typename _Tp1::Element u;
				Hom<Field, _Tp1>(F,T).image(u, A._a);
				Ap = new other(u);
			}
		};


	private:

		// Parameter of this 2x2 block
		Field::Element _a;
	};


} // namespace LinBox

#endif // __LINBOX_cekstv_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

