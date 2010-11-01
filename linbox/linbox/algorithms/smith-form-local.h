/* linbox/algorithms/localsmith.h
 * Copyright(C) LinBox
 *
 * Written by David Saunders
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_smith_form_local_H
#define __LINBOX_smith_form_local_H

#include <vector>
#include <list>
//#include <algorithm>

#include "linbox/matrix/dense-submatrix.h"

namespace LinBox 
{

/** 
 \brief Smith normal form (invariant factors) of a matrix over a local ring.

  The matrix must be a DenseMatrix over a LocalPID.
  A localPID has the standard ring/field arithmetic functions plus gcdin().

 */
template <class LocalPID> 
class SmithFormLocal{

    public:
	typedef typename LocalPID::Element Elt;

	template<class Matrix>
	std::list<Elt>& operator()(std::list<Elt>& L, Matrix& A, const LocalPID& R) { 
		Elt d; R.init(d, 1);
	    return smithStep(L, d, A, R);
	}

	template<class Matrix>
	std::list<Elt>& smithStep(std::list<Elt>& L, Elt& d, Matrix& A, const LocalPID& R) {

//std::cout << "Dimension: " << A.rowdim() << " " << A.coldim() <<"\n";
	    if ( A.rowdim() == 0 || A.coldim() == 0 ) return L;

	    Elt g; R.init(g, 0);
	    typename Matrix::RowIterator p;
	    typename Matrix::Row::iterator q, r;
    	for ( p = A.rowBegin(); p != A.rowEnd(); ++p) {

        	for (q = p->begin(); q != p->end(); ++q) {
       	       	R.gcdin(g, *q);
		    	if ( R.isUnit(g) ) {R.divin(g, g); break; }
			}

			if ( R.isUnit(g) ) break;
	    }
                
	    if ( R.isZero(g) ) {
			L.insert(L.end(), (A.rowdim() < A.coldim()) ? A.rowdim() : A.coldim(), g);
			return L;
	    }

    	if ( p != A.rowEnd() ) // g is a unit and, 
	    // because this is a local ring, value at which this first happened 
	    // also is a unit.
    	{
	        if ( p != A.rowBegin() ) 
		    swap_ranges(A.rowBegin()->begin(), A.rowBegin()->end(), p->begin());
	        if ( q != p->begin() ) 
		    swap_ranges(A.colBegin()->begin(), A.colBegin()->end(), (A.colBegin() + (q - p->begin()))->begin());

	        // eliminate step - crude and for dense only - fix later
			// Want to use a block method or "left looking" elimination.
			Elt f; R.inv(f, *(A.rowBegin()->begin() ) );
			R.negin(f);
			// normalize first row to -1, ...
	        for ( q = A.rowBegin()->begin() /*+ 1*/; q != A.rowBegin()->end(); ++q)
	            R.mulin(*q, f);
			//
			// eliminate in subsequent rows
	        for ( p = A.rowBegin() + 1; p != A.rowEnd(); ++p)
	            for ( q = p->begin() + 1, r = A.rowBegin()->begin() + 1, f = *(p -> begin()); q != p->end(); ++q, ++r )
		        	R.axpyin( *q, f, *r );

	        DenseSubmatrix<Elt> Ap(A, 1, 1, A.rowdim() - 1, A.coldim() - 1);
			L.push_back(d);
	        return smithStep(L, d, Ap, R); 
   		}
    	else  { 
			typename Matrix::RawIterator p;
	        for (p = A.rawBegin(); p != A.rawEnd(); ++p) R.divin(*p, g); 
	        	return smithStep(L, R.mulin(d, g), A, R);
    	    }
    	}

}; // end SmithFormLocal

} // end LinBox

#include <linbox/algorithms/smith-form-local2.h>
#endif // __LINBOX_smith_form_local_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
