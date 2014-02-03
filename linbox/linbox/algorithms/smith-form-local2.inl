/* linbox/algorithms/localsmith.h
 * Copyright (C) LinBox
 *
 * Written by David Saunders
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
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_smith_form_local2_H
#define __LINBOX_smith_form_local2_H


#include <vector>
#include <list>
//#include <algorithm>

#include "linbox/field/local2_32.h"

namespace LinBox
{

	/**
	  \brief Smith normal form (invariant factors) of a matrix over a local ring.
	  */
	template<class LocalRing>
	class SmithFormLocal;

	template <>
	class SmithFormLocal<Local2_32>
	{
	public:
		typedef Local2_32 LocalPIR;
		typedef LocalPIR::Element Elt;

		template<class Matrix>
		std::list<Elt>& operator()(std::list<Elt>& L, Matrix& A, const LocalPIR& R)
		{   Elt d; R.init(d, 1);
			Elt *p = &(A[0][0]);
			return smithStep(L, d, p, A.rowdim(), A.coldim(), A.getStride(), R);
		}

		std::list<Elt>&
		smithStep(std::list<Elt>& L, Elt& d, Elt* Ap, size_t m, size_t n, size_t stride, const LocalPIR& R)
		{
			if ( m == 0 || n == 0 )
				return L;

			LocalPIR::Exponent g = LocalPIR::Exponent(32); //R.init(g); // must change to 2^31 maybe.
			size_t i, j;
			/* Arguably this search order should be reversed to increase the likelyhood of no col swap,
			   assuming row swaps cheaper.  Not so, however on my example. -bds 11Nov */
			for ( i = 0; i != m; ++i)
			{
				for (j = 0; j != n; ++j)
				{
					R.gcdin(g, Ap[i*stride + j]);
					if ( R.isUnit(g) ) break;
				}
				if ( R.isUnit(g) ) break;
			}
			if ( R.isZero(g) )
			{
				L.insert(L.end(), (m < n) ? m : n, 0);
				return L;
			}
			if ( i != m ) // g is a unit and, because this is a local ring,
			// value at which this first happened also is a unit.
			{ // put pivot in 0,0 position
				if ( i != 0 ) // swap rows
					std::swap_ranges(Ap, Ap+n, Ap + i*stride);
				if ( j != 0 ) // swap cols
					for(size_t k = 0; k != m; ++k)
						std::swap(Ap[k*stride + 0], Ap[k*stride + j]);

				// elimination step - crude and for dense only - fix later
				// Want to use a block method or "left looking" elimination.
				Elt f; R.inv(f, Ap[0*stride + 0] );
				R.negin(f);

				// normalize first row to -1, ...
				for ( j = 0; j != n; ++j)
					R.mulin(Ap[0*stride + j], f);

				// eliminate in subsequent rows
				for ( i = 1; i != m; ++i)
				{
					f = Ap[i*stride + 0];
					for ( j = 0; j != n; ++j)
						R.axpyin( Ap[i*stride +j], f, Ap[0*stride +j] );
				}
				L.push_back(d);
				return smithStep(L, d, Ap + stride+1,m-1, n-1, stride, R);
			}
			else
			{
				for ( i = 0; i != m; ++i)
					for ( j = 0; j != n; ++j)
					{
						R.divin(Ap[i*stride + j], g);
					}
				return smithStep(L, R.mulin(d, g), Ap, m, n, stride, R);
			}
		}

	}; // end SmithFormLocal

} // end LinBox

#endif // __LINBOX_smith_form_local2_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

