/* linbox/solutions/trace.h
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

#ifndef __LINBOX_trace_H
#define __LINBOX_trace_H
#include "linbox/solutions/solution-tags.h"

namespace LinBox
{

	/** \brief Sum of the eigenvalues.

	  Also it is the sum of the diagonal entries.

	  Runtime on n by n matrix is n times the cost of getEntry().
	  This is linear in n for those classes where getEntry is constant time
	  (eg DenseMatrix and SparseMatrix).
	  Trace is constant time when the diagonal is necessarily constant, 
	  eg. for ScalarMatrix and Toeplitz.
	  Worst case time is cost of n blackbox applies (matrix vector products), 
	  and apply cost typically ranges between O(n) and O(n^2).
	*/
	template <class BB>
	typename BB::Field::Element & trace(typename BB::Field::Element & t, const BB& A);


/************* internal forms *******************/

	template <class BB>
	typename BB::Field::Element& trace(typename BB::Field::Element& t, const BB& A, SolutionTags::Generic tt);
	/* Specialization for BB's with local trace function.
	   Allows constant time trace for, eg., ScalarMatrix, Toeplitz.
	*/
	template <class BB>
	typename BB::Field::Element & trace(typename BB::Field::Element & t, const BB& A, SolutionTags::Local tt);

	/** Trait to show whether or not the BB class has a local trace function.

	 ScalarMatrix<Field> and Toeplitz<Field,PD> have local trace.
	 */
	template<class BB> struct TraceCategory;

} // namespace LinBox

#include "linbox/solutions/trace.inl"

#endif // __LINBOX_trace_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

