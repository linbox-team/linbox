/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/solutions/getentry.h
 *  Evolved from an earlier one by Bradford Hovinen <hovinen@cis.udel.edu>
 *  -bds
 *
 * See COPYING for license information.
 */

#ifndef __GETENTRY_H
#define __GETENTRY_H

#include <vector>

#include "linbox/util/debug.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/blackbox/dense.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/solutions/methods.h"

namespace LinBox 
{

/** \brief
 * Getting the i,j entry of the blackbox.
 */
template <class BB> 
typename BB::Field::Element& getEntry(typename BB::Field::Element& x, const BB& A, const size_t i, const size_t j)
{ return getEntry(x, A, i, j, Method::Hybrid()); }

// Any BBs that offer a local getEntry can specialize the BB class for the Hybrid method.
/** \brief our best guess

Hybrid method will choose based on matrix size and type
*/

template <class BB> 
typename BB::Field::Element& getEntry(typename BB::Field::Element& x, const BB& A, const size_t i, const size_t j, 
		const Method::Hybrid& m)
{ return getEntry(x, A, i, j, Method::Blackbox(m)); }

// DenseMatrix specialization
template <class Field> 
typename Field::Element& getEntry(typename Field::Element& x, const DenseMatrix<Field>& A, const size_t i, const size_t j, 
		const Method::Hybrid& m)
{	
	return A.getEntry(x,i,j);
}

// SparseMatrix specialization
template <class Field> 
typename Field::Element& getEntry(typename Field::Element& x, const SparseMatrix<Field>& A, const size_t i, const size_t j, 
				  const Method::Hybrid& m)
{
	return A.getEntry(x,i,j);
}

// scalar matrix specialization 
template <class Field>
typename Field::Element & getEntry(typename Field::Element & x, const ScalarMatrix<Field>& A, const size_t i, const size_t j, const Method::Hybrid& m) 
{ return A.getEntry(x, i, j); }


// diagonal specialization 
template <class Field, class Trait>
typename Field::Element & getEntry(typename Field::Element & x, const Diagonal<Field, Trait>& A, const size_t i, const size_t j, const Method::Hybrid& m) 
{ return A.getEntry(x, i, j); }

/** \brief our elimination (a fake in this case)

Elimination method will go to blackbox.
*/
template <class BB> 
typename BB::Field::Element& getEntry(typename BB::Field::Element& x, const BB& A, const size_t i, const size_t j, const Method::Elimination& m)
{ return getEntry(x, A, i, j, Method::Blackbox(m)); 
}


	/** Compute the getEntry of a linear operator A, represented as a black
	 * box. This class is parameterized by the black box type so that it can
	 * be specialized for different black boxes.
	 */

	template <class Blackbox>
	typename Blackbox::Field::Element &getEntry (typename Blackbox::Field::Element &res,
					const Blackbox          &A, const size_t i, const size_t j, 
					const Method::Blackbox& m)
	{

		typedef typename Blackbox::Field Field;
		typedef std::vector<typename Field::Element> Vector;
		Vector v, w;
		VectorWrapper::ensureDim (v, A.coldim ());
		VectorWrapper::ensureDim (w, A.rowdim ());
		const Field& F = A.field();
		typename Field::Element zero; F.init(zero, 0UL);
		typename Vector::iterator it;
		for (it = v.begin (); it != v.end (); ++it)
			F.assign (*it, zero);
		F.init(v[j],1UL);
		F.init (res, 0);
		A.apply (w, v);
		F.assign (res, VectorWrapper::constRef<Field, Vector> (w, i));
		return res;
	}





// Compose< Diagonal, BB > specialization
template <class Field, class Trait, class BlackBox> 
typename Field::Element& getEntry(typename Field::Element& t, const Compose<Diagonal<Field, Trait>, BlackBox>& A, const size_t i, const size_t j, const Method::Hybrid& m)
{
    typename Field::Element y;
    getEntry(y, *(A.getLeftPtr()), i, i);
    getEntry(t, *(A.getRightPtr()), i, j);
    return A.field().mulin(t, y);
}

// Compose< BB, Diagonal > specialization
template <class BlackBox, class Field, class Trait> 
typename Field::Element& getEntry(typename Field::Element& t, const Compose<BlackBox, Diagonal<Field, Trait> >& A, const size_t i, const size_t j, const Method::Hybrid& m)
{
    typename Field::Element y;
    getEntry(y, *(A.getLeftPtr()), i, j);
    getEntry(t, *(A.getRightPtr()), j, j);
    return A.field().mulin(t, y);
}

// Compose< Diagonal, Diagonal > specialization
template <class Field, class T1, class T2> 
typename Field::Element& getEntry(typename Field::Element& t, const Compose<Diagonal<Field,T1>, Diagonal<Field, T2> >& A, const size_t i, const size_t j, const Method::Hybrid& m)
{
    if (i != j) 
        return A.field().init(t, 0UL);
    else {
        typename Field::Element y;
        getEntry(y, *(A.getLeftPtr()), i, i);
        getEntry(t, *(A.getRightPtr()), j, j);
        return A.field().mulin(t, y);
    }
}


}

#endif // __GETENTRY_H
