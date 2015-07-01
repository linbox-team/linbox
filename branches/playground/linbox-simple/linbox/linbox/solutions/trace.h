/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/solutions/trace.h
 *  Evolved from an earlier one by Bradford Hovinen <hovinen@cis.udel.edu>
 *  -bds
 *
 * See COPYING for license information.
 */

#ifndef __TRACE_H
#define __TRACE_H

#include <vector>
//#include <algorithm>

// must fix this list...
//#include "linbox/algorithms/wiedemann.h"
//#include "linbox/algorithms/lanczos.h"
//#include "linbox/algorithms/block-lanczos.h"
#include "linbox/util/debug.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/blackbox/dense.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/solutions/methods.h"
#include "linbox/solutions/getentry.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/compose.h"

namespace LinBox 
{


/* for trace we actually use only the blackbox method and local defs for sparse,
 dense, and a few other BB types.
*/
/** \brief sum of eigenvalues
 
Also sum of diagonal entries.
This is the generic one.

testing multiple doc comments.
*/
template <class BB> 
typename BB::Field::Element& trace(typename BB::Field::Element& t, const BB& A)
{ return trace(t, A, Method::Hybrid()); }

// Any BBs that offer a local trace can specialize the BB class for the Hybrid method.
/** \brief our best guess

Hybrid method will choose based on matrix size and type
*/

template <class BB> 
typename BB::Field::Element& trace(typename BB::Field::Element& t, const BB& A, 
		const Method::Hybrid& m)
{ return trace(t, A, Method::Blackbox(m)); }

// DenseMatrix specialization
template <class Field> 
typename Field::Element& trace(typename Field::Element& t, const DenseMatrix<Field>& A, 
		const Method::Hybrid& m)
{	typename Field::Element x;
	A.field().init(t, 0);
	for (size_t i = 0; i < A.coldim(); ++i)  {
		A.getEntry(x,i,i);
		A.field().addin(t, x);
	}
	return t;
}

// SparseMatrix specialization
template <class Field, class Row> 
typename Field::Element& trace(typename Field::Element& t, const SparseMatrix<Field, Row>& A, 
		const Method::Hybrid& m)
{	typename Field::Element x;
	A.field().init(t, 0);
	for (size_t i = 0; i < A.coldim(); ++i) { 
		A.getEntry(x,i,i);
		A.field().addin(t, x);
	}
	return t;
}

// Diagonal specialization
template <class Field, class Trait> 
typename Field::Element& trace(typename Field::Element& t, const Diagonal<Field, Trait>& A, 
		const Method::Hybrid& m)
{	typename Field::Element x;
	A.field().init(t, 0);
	for (size_t i = 0; i < A.coldim(); ++i) { 
		A.getEntry(x,i,i);
		A.field().addin(t, x);
	}
	return t;
}

// scalar matrix specialization 
template <class Field>
typename Field::Element & trace(typename Field::Element & t, const ScalarMatrix<Field>& A, const Method::Hybrid& m) 
{ return A.trace(t); }

/** \brief our elimination (a fake in this case)

Elimination method will go to blackbox.
*/
template <class BB> 
typename BB::Field::Element& trace(typename BB::Field::Element& t, const BB& A, const Method::Elimination& m)
{ return trace(t, A, Method::Blackbox(m)); 
}


    /** Compute the trace of a linear operator A, represented as a black
     * box. This class is parameterized by the black box type so that it can
     * be specialized for different black boxes.
     */

template <class Blackbox>
typename Blackbox::Field::Element &trace (typename Blackbox::Field::Element &res,
                                          const Blackbox          &A, 
                                          const Method::Blackbox& m)
{
    
    typedef typename Blackbox::Field Field;
    typedef std::vector<typename Field::Element> Vector;
    Vector v, w;
    Field F = A.field();
    StandardBasisStream<Field, Vector> stream (F, A.coldim ());
    
    VectorWrapper::ensureDim (v, A.coldim ());
    VectorWrapper::ensureDim (w, A.rowdim ());
    
    F.init (res, 0);
    
    while (stream) {
        stream >> v;
        A.apply (w, v);
        F.addin (res, VectorWrapper::constRef<Field, Vector> (w, stream.pos () - 1));
    }
    
    return res;
}


// Compose< Diagonal, BB > specialization
template <class Field, class Trait, class BlackBox> 
typename Field::Element& trace(typename Field::Element& t, const Compose<Diagonal<Field, Trait>, BlackBox>& A, const Method::Hybrid& m)
{
    typename Field::Element x, y;
    A.field().init(t, 0);
    size_t n = (A.coldim()<A.rowdim()?A.coldim():A.rowdim());
    for (size_t i = 0; i < n; ++i) { 
        getEntry(x, *(A.getRightPtr()), i, i);
        getEntry(y, *(A.getLeftPtr()), i, i);
        A.field().axpyin(t, x, y);
    }
    return t;
}

// Compose< BB, Diagonal > specialization
template <class BlackBox, class Field, class Trait> 
typename Field::Element& trace(typename Field::Element& t, const Compose<BlackBox, Diagonal<Field, Trait> >& A, const Method::Hybrid& m)
{
    typename Field::Element x, y;
    A.field().init(t, 0);
    size_t n = (A.coldim()<A.rowdim()?A.coldim():A.rowdim());
    for (size_t i = 0; i < n; ++i) { 
        getEntry(x, *(A.getRightPtr()), i, i);
        getEntry(y, *(A.getLeftPtr()), i, i);
        A.field().axpyin(t, x, y);
    }
    return t;
}


// Compose< Diagonal, Diagonal > specialization
template <class Field, class T1, class T2> 
typename Field::Element& trace(typename Field::Element& t, const Compose<Diagonal<Field,T1>, Diagonal<Field, T2> >& A, const Method::Hybrid& m)
{
    typename Field::Element x, y;
    A.field().init(t, 0);
    size_t n = (A.coldim()<A.rowdim()?A.coldim():A.rowdim());
    for (size_t i = 0; i < n; ++i) { 
        getEntry(x, *(A.getRightPtr()), i, i);
        getEntry(y, *(A.getLeftPtr()), i, i);
        A.field().axpyin(t, x, y);
    }
    return t;
}


}



#endif // __TRACE_H
