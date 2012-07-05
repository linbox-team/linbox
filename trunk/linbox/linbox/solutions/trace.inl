/* linbox/solutions/trace.inl
 * Copyright(C) LinBox
 *  Evolved from an earlier one by Bradford Hovinen <hovinen@cis.udel.edu>
 *  -bds
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

#ifndef __LINBOX_trace_INL
#define __LINBOX_trace_INL

#include <vector>

#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/toeplitz.h"
#include "linbox/solutions/getentry.h"

namespace LinBox
{

	// Trait to show whether or not the BB class has a local trace function.
	template<class BB> struct TraceCategory;

	/// undocumented
	namespace TraceTags
	{
		struct Generic{};
		struct Local{};
	}

	template<class BB> struct TraceCategory		{ typedef TraceTags::Generic Tag; };

	template<class Field>
	struct TraceCategory<ScalarMatrix<Field> > 	{ typedef TraceTags::Local Tag; };

	template<class Field, class PD>
	struct TraceCategory<Toeplitz<Field,PD> >	{ typedef typename TraceTags::Local Tag; };

	template <class BB>
	typename BB::Field::Element& trace(typename BB::Field::Element& t, const BB& A)
	{
		typename TraceCategory<BB>::Tag tt;
		return trace(t, A, tt);
	}

	/* Generic approach.  It will be efficient for BBs with efficient getEntry.
	   If getEntry is constant time on n by n BB's, trace will be in O(n).
	   */
	template <class BB>
	typename BB::Field::Element& trace(typename BB::Field::Element& t, const BB& A, TraceTags::Generic tt)
	{
		typename BB::Field::Element x;
		A.field().init(x, 0);
		A.field().init(t, 0);
		for (size_t i = 0; i < A.coldim(); ++i)
			A.field().addin(t, getEntry(x,A,i,i));
		return t;
	}

	/* Specialization for BB's with local trace function.
	   Allows constant time trace for, eg., ScalarMatrix, Toeplitz.
	   */
	template <class BB>
	typename BB::Field::Element & trace(typename BB::Field::Element & t, const BB& A, TraceTags::Local tt)
	{ return A.trace(t); }

} // namespace LinBox

#endif // __LINBOX_trace_INL


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

