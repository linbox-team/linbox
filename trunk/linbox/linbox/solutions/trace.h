/* linbox/solutions/trace.h
 * Copyright(C) LinBox
 *  Evolved from an earlier one by Bradford Hovinen <hovinen@cis.udel.edu>
 *  -bds
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_trace_H
#define __LINBOX_trace_H

#include <vector>

#include "linbox/blackbox/scalar-matrix.h"
//#include "linbox/blackbox/toeplitz.h"
#include "linbox/solutions/getentry.h"

namespace LinBox 
{

// Trait to show whether or not the BB class has a local trace function.
template<class BB> struct TraceCategory;

namespace TraceTags { struct Generic{}; struct Local{}; };

template<class BB> struct TraceCategory		{ typedef TraceTags::Generic Tag; };

template<class Field>
struct TraceCategory<ScalarMatrix<Field> > 	{ typedef TraceTags::Local Tag; };

//template<class Field, class PD>
//struct TraceCategory<Toeplitz<Field,PD> >	{ typedef typename TraceTags::Local Tag; };

// trace
/** \brief Sum of the eigenvalues.
 
Also it is the sum of the diagonal entries.

Runtime on n by n matrix is n times the cost of getEntry().  
This is linear in n for those classes where getEntry is constant time
(eg DenseMatrix and SparseMatrix).  
Trace is constant time when the diagonal is necessarily constant, eg. for ScalarMatrix and Toeplitz.
Worst case time is cost of n blackbox applies (matrix vector products), and apply cost typically ranges between O(n) and O(n^2).

*/
template <class BB> 
typename BB::Field::Element& trace(typename BB::Field::Element& t, const BB& A)
{	typename TraceCategory<BB>::Tag tt; 
	return trace(t, A, tt);	
}

/* Generic approach.  It will be efficient for BBs with efficient getEntry.
If getEntry is constant time on n by n BB's, trace will be in O(n).
*/
template <class BB> 
typename BB::Field::Element& trace(typename BB::Field::Element& t, const BB& A, TraceTags::Generic tt)
{	typename BB::Field::Element x;
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

}; // namespace LinBox

#endif // __LINBOX_trace_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
